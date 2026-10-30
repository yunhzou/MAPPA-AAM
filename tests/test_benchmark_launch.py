"""Submission lifecycle contracts; scientific checks stay in existing modules."""
from dataclasses import asdict
import hashlib
import json
from types import SimpleNamespace

import pytest

from rebenchmark import launch


def test_submission_records_dependencies_and_refuses_duplicate_run(tmp_path, monkeypatch):
    repo = tmp_path/'repo'; (repo/'src/mappa').mkdir(parents=True)
    monkeypatch.setattr(launch, 'REPO', repo)
    monkeypatch.setattr(launch.subprocess, 'check_output',
        lambda args, **kw: 'abc123\n' if 'rev-parse' in args else '')
    commands = []
    def submitted(args, **kw):
        commands.append(args)
        return SimpleNamespace(args=args,stdout=str(100+len(commands))+'\n',stderr='',returncode=0,
                               check_returncode=lambda:None)
    monkeypatch.setattr(launch.subprocess, 'run', submitted)
    configuration = tmp_path/'config.json'
    configuration.write_text(json.dumps(dict(slurm=dict(account='test',partition='cpu',long_qos='long',
        **{k:dict(cpus=2,memory_gib=8,time='01:00:00') for k in ('prepare','pilot','coordinator')}))))
    result = launch.submit(tmp_path/'run',configuration)
    assert result['production']=='103'
    assert '--dependency=afterok:101' in commands[1]
    assert '--dependency=afterok:102' in commands[2]
    assert commands[0][-1]=='rebenchmark.shared'
    assert 'bash' not in commands[0]
    assert any(a.endswith('/run_conversion.sbatch') for a in commands[0])
    assert commands[1][-1]=='rebenchmark.launch'
    assert all((tmp_path/'run'/f'{stage}-launch.json').exists() for stage in ('prepare','pilot','production'))
    with pytest.raises(FileExistsError):launch.submit(tmp_path/'run',configuration)
    assert len(commands)==3


def test_capacity_respects_memory_cpu_and_excludes_unavailable_nodes():
    nodes = '\n'.join([
        'NodeName=cpu-0001 State=MIXED CPUTot=64 CPUAlloc=48 RealMemory=65536 AllocMem=20480',
        'NodeName=cpu-0002 State=MIXED CPUTot=64 CPUAlloc=16 RealMemory=65536 AllocMem=32768',
        'NodeName=cpu-0003 State=IDLE+RESERVED CPUTot=64 CPUAlloc=0 RealMemory=65536 AllocMem=0'])
    result = launch.capacity(nodes,2*2**30,2)
    assert result['compute_workers']==29
    assert len(result['groups'])==2
    assert launch.capacity('',2*2**30,4)['groups']==[]


def test_worker_renewal_requires_clean_exit_and_every_node_receipt(tmp_path, monkeypatch):
    job=dict(job='100',group=dict(nodes=1))
    def output(args,**kw):
        if args[0]=='squeue':
            assert '-j' not in args  # Old job IDs disappear from the live scheduler.
            return '100\n'
        pytest.fail('Active job must not be renewed')
    monkeypatch.setattr(launch.subprocess,'check_output',output)
    assert launch.retired_workers(tmp_path,job) is None
    monkeypatch.setattr(launch.subprocess,'check_output',lambda args,**kw:'' if args[0]=='squeue' else 'FAILED|1:0\n')
    with pytest.raises(RuntimeError,match='failed'):launch.retired_workers(tmp_path,job)
    monkeypatch.setattr(launch.subprocess,'check_output',lambda args,**kw:'' if args[0]=='squeue' else 'COMPLETED|0:0\n')
    with pytest.raises(RuntimeError,match='receipts'):launch.retired_workers(tmp_path,job)
    (tmp_path/'retirements').mkdir()
    (tmp_path/'retirements/100-node-42.json').write_text(json.dumps(dict(states=dict(pending=1))))
    assert launch.retired_workers(tmp_path,job)==[dict(states=dict(pending=1))]


def test_small_enzyme_pipeline_prepares_checks_and_reuses_pilot(tmp_path,monkeypatch):
    import csv,gzip
    from mappa import AAMSearchConfig
    from rebenchmark.shared import prepare_campaign,prepare_queue
    monkeypatch.setenv('TMPDIR',str(tmp_path))
    dataset=tmp_path/'dataset.csv.gz'
    with gzip.open(dataset,'wt') as stream:
        writer=csv.DictWriter(stream,fieldnames=['r_id','rxn','ground_truth','original_id'])
        writer.writeheader()
        writer.writerow(dict(r_id=0,rxn='CO>>CO',ground_truth='[CH3:1][OH:2]>>[CH3:1][OH:2]',original_id='tiny'))
    manifest=tmp_path/'manifest.json'
    manifest.write_text(json.dumps(dict(cases=1,dataset_sha256=hashlib.sha256(dataset.read_bytes()).hexdigest())))
    plan=dict(dataset='EnzymeMap',dataset_path=str(dataset),input_manifest=str(manifest),
        input_manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest(),directions=[False,True],
        config=asdict(AAMSearchConfig(seed_count=1,iso_tolerance=1,sweep_cuts=False)),
        seed_batch=1,case_seconds=60,memory_bytes=2**30,slack=5,reference_checks=True,preparation_workers=1,
        identity='test-launch',commit='source',pilot_commits=['source'],pilots=[],expected_cases=1,
        pilot_cases=[0],pilot_workers=2)
    (tmp_path/'plan.json').write_text(json.dumps(plan))
    prepare_campaign(tmp_path)
    launch.pilot_check(tmp_path)
    check=launch.read(tmp_path/'pilot/check.json')
    assert check['verified'] and check['event_outputs']==check['reference_outputs']==2
    scratch=tmp_path/'production-scratch';scratch.mkdir()
    queue=prepare_queue(launch.read(tmp_path/'full/plan.json'),tmp_path/'full',scratch)
    assert set(queue.snapshot())=={'complete'}
    assert queue.claim('no-repeat',2**30,queue.identity)['state']=='finished'
    queue.close()


def test_prepared_recovery_reuses_only_verified_unstarted_work(tmp_path,monkeypatch):
    previous=tmp_path/'previous';(previous/'full').mkdir(parents=True)
    pilot=previous/'pilot';pilot.mkdir()
    (pilot/'check.json').write_text(json.dumps(dict(verified=True)))
    (pilot/'progress.json').write_text(json.dumps(dict(finished=True)))
    index=previous/'index.json';index.write_text('{}')
    plan={k:'same' for k in ('dataset','dataset_path','input_manifest_sha256','directions','config',
        'seed_batch','case_seconds','memory_bytes','slack','reference_checks','expected_cases')}
    plan.update(commit='new',identity='recovery')
    prior=dict(plan,commit='old',index=str(index),index_sha256=hashlib.sha256(index.read_bytes()).hexdigest(),
               pilots=[str(pilot)],pilot_commits=['old'])
    (previous/'full/plan.json').write_text(json.dumps(prior))
    (previous/'production-launch.json').write_text(json.dumps(dict(job='123')))
    monkeypatch.setattr(launch.subprocess,'check_output',lambda args,**kw:'OUT_OF_MEMORY\n' if args[0]=='sacct' else '')
    result=launch.prepared_plan(previous,plan)
    assert result['pilots']==[str(pilot)] and result['pilot_commits']==['old']
    assert result['commit']=='new' and result['index']==str(index)
    changed=dict(plan,config='different')
    with pytest.raises(ValueError,match='settings differ'):launch.prepared_plan(previous,changed)
    (previous/'full/endpoint.json').write_text('{}')
    with pytest.raises(ValueError,match='already started'):launch.prepared_plan(previous,plan)
