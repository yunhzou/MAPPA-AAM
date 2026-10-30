"""Submission lifecycle contracts; scientific checks stay in existing modules."""
from dataclasses import asdict
import hashlib
import json
from types import SimpleNamespace

import pytest

from rebenchmark import launch, execution


def test_scheduler_outage_retries_same_pending_only_request(monkeypatch):
    import subprocess
    calls, pauses = [], []
    def output(arguments, **options):
        calls.append(arguments)
        if len(calls) == 1:
            raise subprocess.CalledProcessError(1, arguments,
                stderr='slurm_load_jobs error: Socket timed out on send/recv operation')
        if len(calls) == 2:
            raise subprocess.TimeoutExpired(arguments, 60)
        return ''
    monkeypatch.setattr(launch.subprocess, 'check_output', output)
    monkeypatch.setattr(launch.time, 'sleep', pauses.append)
    command = ['scancel', '--state=PENDING', '100']
    assert launch.scheduler_output(command) == ''
    assert calls == [command]*3 and pauses == [5, 5]
    def rejected(*args, **kwargs):
        raise subprocess.CalledProcessError(1, args[0], stderr='Access denied')
    monkeypatch.setattr(launch.subprocess, 'check_output', rejected)
    with pytest.raises(subprocess.CalledProcessError):
        launch.scheduler_output(command)


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
    configuration.write_text(json.dumps(dict(slurm=dict(account='test',partition='cpu',
        **{k:dict(qos='short',cpus=2,memory_gib=8,time='01:00:00') for k in ('prepare','pilot','queue','coordinator')}))))
    from rebenchmark import inputs
    monkeypatch.setattr(inputs, 'validate', lambda plan: None)
    result = launch.submit(tmp_path/'run',configuration)
    assert result['production']=='103'
    assert '--dependency=afterok:101' in commands[1]
    assert '--dependency=afterok:102' in commands[2]
    assert 'pilot-check' in commands[0]
    assert not any(a.startswith('--dependency') for a in commands[0])
    assert 'bash' not in commands[0]
    assert any(a.endswith('/run_conversion.sbatch') for a in commands[0])
    assert 'prepare-queue' in commands[1]
    assert all((tmp_path/'run'/f'{stage}-launch.json').exists() for stage in ('pilot','queue','production'))
    with pytest.raises(FileExistsError):launch.submit(tmp_path/'run',configuration)
    assert len(commands)==3


def test_capacity_uses_free_cpu_and_allocatable_memory_with_strict_ceiling():
    nodes = '\n'.join(
        f'NodeName=cpu-{i:04d} State=MIXED CPUTot=96 CPUAlloc=64 RealMemory=207698 AllocMem=131072'
        for i in range(1, 12))
    result = launch.capacity(nodes, 2*2**30, 3, 800)
    assert result['compute_workers'] == 11*31
    assert sum(g['nodes']*g['workers_per_node'] for g in result['groups']) == 341
    idle = nodes.replace('CPUAlloc=64', 'CPUAlloc=0').replace('AllocMem=131072', 'AllocMem=0')
    assert launch.capacity(idle, 2*2**30, 3, 800)['compute_workers'] == 800
    assert launch.capacity(nodes.replace('AllocMem=131072', 'AllocMem=204800'), 2*2**30, 3, 800)['groups'] == []
    assert launch.capacity(nodes+' ReservationName=reserved', 2*2**30, 3, 800)['compute_workers'] == 310
    assert launch.capacity(nodes.replace('State=MIXED', 'State=IDLE+RESERVED'), 2*2**30, 3, 800)['groups'] == []
    assert launch.capacity(idle, 2*2**30, 3, 1)['groups'] == [dict(nodes=1, workers_per_node=1)]


def test_scaling_preserves_live_workers_and_refits_retired_capacity(tmp_path, monkeypatch):
    plan = dict(memory_bytes=2*2**30, slurm=dict(max_worker_jobs=3, target_workers=800))
    live = {'100': 'RUNNING'}
    nodes = 'NodeName=cpu-0001 State=MIXED CPUTot=96 CPUAlloc=64 RealMemory=207698 AllocMem=65536'
    monkeypatch.setattr(launch, 'scheduler_output', lambda args: '\n'.join(k+'|'+v for k,v in live.items()) if args[0]=='squeue' else nodes)
    monkeypatch.setattr(launch, 'retired_workers', lambda full, job: None if job['job'] in live else [])
    calls = []
    def start(full, plan, slot, generation, group):
        calls.append((slot, generation, group)); return dict(job=str(101+len(calls)))
    monkeypatch.setattr(launch, 'start_workers', start)
    jobs = [dict(job='100', generation=0, group=dict(nodes=1, workers_per_node=20))]
    jobs = launch.scale_workers(tmp_path, plan, jobs)
    assert calls == [(1, 0, dict(nodes=1, workers_per_node=31))]
    assert jobs[0]['job'] == '100'
    # Retired allocations are re-sized, not renewed with the old 20-worker shape.
    live.clear()
    jobs = launch.scale_workers(tmp_path, plan, jobs)
    assert calls[-1] == (0, 1, dict(nodes=1, workers_per_node=31))


def test_scaling_cancels_only_pending_and_waits_for_retirement(tmp_path, monkeypatch):
    plan = dict(memory_bytes=2*2**30, slurm=dict(max_worker_jobs=3, target_workers=800))
    calls = []
    def output(args):
        calls.append(args)
        if args[0] == 'squeue': return '100|PENDING'
        if args[0] == 'scontrol': return 'NodeName=cpu-0001 State=MIXED CPUTot=96 CPUAlloc=64 RealMemory=207698 AllocMem=65536'
        return ''
    monkeypatch.setattr(launch, 'scheduler_output', output)
    monkeypatch.setattr(launch, 'retired_workers', lambda *a: None)
    monkeypatch.setattr(launch, 'start_workers', lambda *a: pytest.fail('Must confirm retirement before replacement'))
    jobs = [dict(job='100', generation=0, group=dict(nodes=3, workers_per_node=89))]
    launch.scale_workers(tmp_path, plan, jobs)
    assert ['scancel', '--state=PENDING', '100'] in calls
    assert jobs[0]['cancel_pending']


def test_scaling_counts_existing_workers_against_ceiling(tmp_path, monkeypatch):
    plan = dict(memory_bytes=2*2**30, slurm=dict(max_worker_jobs=3, target_workers=800))
    monkeypatch.setattr(launch, 'retired_workers', lambda *a: None)
    monkeypatch.setattr(launch, 'scheduler_output', lambda args: '100|RUNNING' if args[0]=='squeue' else 'NodeName=cpu-0001 State=IDLE CPUTot=96 CPUAlloc=0 RealMemory=207698 AllocMem=0')
    seen=[]
    def start(full, plan, slot, generation, group):
        seen.append(group);return dict(job='101')
    monkeypatch.setattr(launch, 'start_workers', start)
    launch.scale_workers(tmp_path, plan, [dict(job='100', generation=0, group=dict(nodes=8, workers_per_node=99))])
    assert seen == [dict(nodes=1, workers_per_node=8)]


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
    from rebenchmark.inputs import prepare_campaign
    from rebenchmark.execution import prepare_queue
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
    plan.update(index=str(tmp_path/'index.json'),index_sha256=launch.sha(tmp_path/'index.json'))
    launch.save(tmp_path/'plan.json',plan)
    from rebenchmark.inputs import validate
    validate(plan)
    with pytest.raises(ValueError,match='settings'):
        validate(dict(plan,slack=6))
    launch.pilot_check(tmp_path)
    check=launch.read(tmp_path/'pilot/check.json')
    assert check['verified'] and check['event_outputs']==check['reference_outputs']==2
    scratch=tmp_path/'production-scratch';scratch.mkdir()
    queue=prepare_queue(launch.read(tmp_path/'full/plan.json'),tmp_path/'full',scratch)
    assert set(queue.snapshot())=={'complete'}
    assert queue.claim('no-repeat',2**30,queue.identity)['state']=='finished'
    queue.backup(tmp_path/'full/queue.sqlite')
    launch.save(tmp_path/'full/progress.json',dict(finished=True,states=queue.snapshot()))
    from rebenchmark.report import summarize_original_events
    report=summarize_original_events(tmp_path/'full',tmp_path/'complete-report',1)
    assert report['unresolved_tasks']==0
    assert report['modes']['sweep/bidirectional/seeds1']['recovered']==1
    # Failed reference/decode rows must report unknowns without hiding their keys.
    queue.db.execute("UPDATE tasks SET state='timeout',result=NULL WHERE key LIKE '%/verify/%'")
    queue.db.execute("UPDATE tasks SET state='blocked',result=NULL WHERE key LIKE '%/seed/%'")
    queue.db.commit();queue.backup(tmp_path/'full/queue.sqlite')
    launch.save(tmp_path/'full/progress.json',dict(finished=True,states=queue.snapshot()))
    queue.close()
    from rebenchmark.report import summarize_original_events
    report=summarize_original_events(tmp_path/'full',tmp_path/'report',1)
    assert report['unresolved_tasks']==4
    assert report['modes']['sweep/bidirectional/seeds1']['unknown']==1
    assert report['modes']['sweep/bidirectional/seeds1']['mapping_space_preserved']==0
    assert report['counts']['verify_seed']==2 and report['counts']['convert_seed']==2


def test_retired_checkpoint_handoff_preserves_results_attempts_and_budget(tmp_path,monkeypatch):
    from mappa.mappa_core.core.scheduler import Task,TaskQueue
    from rebenchmark import execution
    import mappa.mappa_core.core.scheduler.queue as queue_module
    clock=[100.]
    monkeypatch.setattr(queue_module.time,'time',lambda:clock[0])
    monkeypatch.setattr(execution.time,'sleep',lambda seconds:None)
    monkeypatch.setattr(execution.socket,'gethostname',lambda:'127.0.0.1')
    monkeypatch.setenv('SLURM_JOB_ID','100')
    full=tmp_path/'full';full.mkdir()
    queue=TaskQueue(tmp_path/'live.sqlite','handoff')
    queue.add([Task('done','execution_helpers:value',{},case='same',case_seconds=60),
               Task('next','execution_helpers:value',{},case='same',case_seconds=60,dependencies=('done',))])
    claim=queue.claim('old',2**30,'handoff');clock[0]=107.
    queue.finish('done',claim['token'],'complete',{'value':42})
    before={table:list(queue.db.execute('SELECT * FROM '+table)) for table in ('tasks','attempts','cases')}
    launch.save(full/'retire.json',{'workers':[]})
    execution.serve_queue(queue,full,{},retire_file=full/'retire.json')
    queue.close()
    assert launch.read(full/'progress.json')['finished'] is False
    plan=dict(identity='handoff',index_sha256='index',campaign_root=str(tmp_path),generation=0,
              claim_until=100,commit='source',python='python',slurm=dict(coordinator=dict(qos='short',
              cpus=2,memory_gib=8,time='01:00:00')))
    submissions=[]
    def submit(*args,**kw):
        submissions.append((args,kw));return {'job':'101'}
    monkeypatch.setattr(launch,'submit_job',submit)
    launch.next_allocation(full,plan,{'complete':1,'pending':1})
    assert submissions[0][1]['dependency']=='afterok:100'
    following=tmp_path/'continuations/0001'
    continued=launch.read(following/'plan.json')
    assert continued['identity']==plan['identity'] and 'claim_until' not in continued
    scratch=tmp_path/'scratch';scratch.mkdir()
    restored=execution.restore_queue(continued,following,scratch)
    assert before=={table:list(restored.db.execute('SELECT * FROM '+table)) for table in before}
    clock[0]=1000.
    claim=restored.claim('new',2**30,'handoff')
    assert claim['task'].key=='next' and claim['deadline']==1053.
    restored.finish('next',claim['token'],'complete',{'value':43})
    assert restored.snapshot()=={'complete':2}
    assert list(restored.db.execute('SELECT task,count(*) FROM attempts GROUP BY task'))==[('done',1),('next',1)]
    restored.close()
    with pytest.raises(FileExistsError):launch.next_allocation(full,plan,{'complete':1,'pending':1})
    assert len(submissions)==1


def test_checkpoint_rejects_active_owners_and_changed_bytes(tmp_path):
    from mappa.mappa_core.core.scheduler import Task,TaskQueue
    queue=TaskQueue(tmp_path/'live.sqlite','handoff')
    queue.add([Task('a','execution_helpers:value',{})])
    queue.claim('active',2**30,'handoff')
    queue.backup(tmp_path/'snapshot.sqlite');queue.close()
    receipt=dict(identity='handoff',index_sha256='index',queue=str(tmp_path/'snapshot.sqlite'),
                 queue_sha256=launch.sha(tmp_path/'snapshot.sqlite'),states={'running':1})
    launch.save(tmp_path/'handoff.json',receipt)
    plan=dict(identity='handoff',index_sha256='index',resume_receipt=str(tmp_path/'handoff.json'))
    scratch=tmp_path/'scratch';scratch.mkdir()
    with pytest.raises(ValueError,match='active owners'):execution.restore_queue(plan,tmp_path,scratch)
    receipt['states']={'pending':1};launch.save(tmp_path/'handoff.json',receipt)
    with pytest.raises(ValueError,match='not retired'):execution.restore_queue(plan,tmp_path,scratch)
    receipt['queue_sha256']='changed';launch.save(tmp_path/'handoff.json',receipt)
    with pytest.raises(ValueError,match='hash differs'):execution.restore_queue(plan,tmp_path,scratch)


def test_pending_cancellation_needs_no_worker_receipt_but_started_cancellation_fails(tmp_path,monkeypatch):
    job=dict(job='100',group=dict(nodes=1),cancel_pending=True)
    started=['None']
    def output(args,**kw):
        if args[0]=='squeue':return ''
        return started[0] if args[-1]=='--format=Start' else 'CANCELLED by 10|0:0\n'
    monkeypatch.setattr(launch.subprocess,'check_output',output)
    assert launch.retired_workers(tmp_path,job)==[]
    started[0]='2026-10-30T00:00:00'
    with pytest.raises(RuntimeError,match='failed'):launch.retired_workers(tmp_path,job)


def test_queue_preparation_publishes_restorable_snapshot_before_coordinator(tmp_path,monkeypatch):
    from mappa.mappa_core.core.scheduler import Task,TaskQueue
    from rebenchmark import execution
    monkeypatch.setenv('TMPDIR',str(tmp_path))
    plan=dict(identity='fresh',index_sha256='index')
    launch.save(tmp_path/'plan.json',plan)
    def prepared(plan,out,scratch):
        queue=TaskQueue(Path(scratch)/'queue.sqlite',plan['identity'])
        queue.add([Task('a','execution_helpers:value',{})])
        return queue
    from pathlib import Path
    monkeypatch.setattr(execution,'prepare_queue',prepared)
    execution.prepare_full_queue(tmp_path)
    receipt=launch.read(tmp_path/'prepared-queue.json')
    assert receipt['states']=={'pending':1}
    assert receipt['queue_sha256']==launch.sha(tmp_path/'queue.sqlite')
    scratch=tmp_path/'restore';scratch.mkdir()
    restored=execution.restore_queue(launch.read(tmp_path/'plan.json'),tmp_path,scratch)
    assert restored.snapshot()=={'pending':1}
    assert restored.db.execute('SELECT count(*) FROM attempts').fetchone()[0]==0
    restored.close()


def test_status_separates_requested_from_running_workers_and_reports_failures(tmp_path,monkeypatch):
    full=tmp_path/'full';full.mkdir()
    launch.save(tmp_path/'submission.json',dict(source='source'))
    launch.save(tmp_path/'plan.json',dict(slurm=dict(target_workers=800)))
    launch.save(tmp_path/'production-launch.json',dict(job='100'))
    launch.save(full/'workers-00-00-launch.json',dict(job='101'))
    launch.save(full/'workers-01-00-launch.json',dict(job='102'))
    launch.save(full/'worker-jobs.json',[dict(job=j,group=dict(nodes=3,workers_per_node=89)) for j in ('101','102')])
    launch.save(full/'progress.json',dict(states=dict(complete=12,pending=20),utc=0))
    def output(args):
        if args[0]=='squeue':return '100|RUNNING|2|1:00|cpu-1\n101|RUNNING|270|1:00|cpu-[2-4]\n102|PENDING|270|0:00|Resources\n999|RUNNING|64|1:00|other\n'
        if args[0]=='sacct':return '100|RUNNING|0:0\n'
        return '100.batch|00:00:04|200M|0|0\n'
    monkeypatch.setattr(launch,'scheduler_output',output)
    result=launch.status(tmp_path)
    assert result['requested_compute_workers']==534 and result['running_compute_workers']==267
    assert result['target_compute_workers']==800 and result['checkpoint'] is None
    assert len(result['scheduler'])==3 and result['progress']['states']['complete']==12
    monkeypatch.setattr(launch,'scheduler_output',lambda args:'' if args[0]=='squeue' else '100|FAILED|1:0')
    assert launch.status(tmp_path)['running_compute_workers']==0


def test_failed_service_retains_local_queue_and_reports_location(tmp_path,monkeypatch):
    from rebenchmark import execution
    from pathlib import Path
    monkeypatch.setenv('TMPDIR',str(tmp_path))
    launch.save(tmp_path/'plan.json',dict(identity='fresh'))
    closed=[]
    def restored(plan,full,scratch):
        (scratch/'queue.sqlite').write_bytes(b'preserve failed queue')
        return SimpleNamespace(close=lambda:closed.append(True))
    monkeypatch.setattr(execution,'restore_queue',restored)
    def failed(*args,**kwargs):raise RuntimeError('transport failed')
    monkeypatch.setattr(execution,'serve_queue',failed)
    monkeypatch.setattr(launch.signal,'signal',lambda *args:None)
    with pytest.raises(RuntimeError,match='transport failed'):execution.service(tmp_path)
    receipt=launch.read(tmp_path/'retained-scratch.json')
    assert (Path(receipt['path'])/'queue.sqlite').read_bytes()==b'preserve failed queue'
    assert closed==[True]


@pytest.mark.parametrize('states,kind', [({'complete':2},'report'),
    ({'complete':1,'timeout':1,'error':1,'blocked':1},'report'),
    ({'complete':1,'memory_limit':1},'memory')])
def test_terminal_routing_and_memory_stage_reports_once(tmp_path,monkeypatch,states,kind):
    monkeypatch.setenv('SLURM_JOB_ID','100')
    full=tmp_path/'full';full.mkdir()
    resource=dict(qos='short',cpus=2,memory_gib=8,time='01:00:00')
    plan=dict(campaign_root=str(tmp_path),identity='test',generation=0,commit='source',python='python',
        memory_bytes=2**30,retry_memory_bytes=8*2**30,retry_workers=1,
        slurm=dict(memory_stage=resource,report=resource))
    submitted=[]
    monkeypatch.setattr(launch,'submit_job',lambda *a,**k:submitted.append((a,k)) or dict(job='101'))
    launch.complete_stage(full,plan,dict(states=states,finished=True))
    assert submitted[0][0][1]==('memory-stage' if kind=='memory' else 'report')
    assert submitted[0][1]['dependency']=='afterok:100'
    if kind=='memory':
        stage=tmp_path/'memory-stage';retry=launch.read(stage/'plan.json')
        assert retry['previous_stage']==str(full) and retry['memory_stage']
        launch.complete_stage(stage,retry,dict(states={'complete':1,'memory_limit':1},finished=True))
        assert submitted[-1][0][1]=='report'  # No repeated escalation loop.
    assert launch.read(tmp_path/'terminal.json')['unresolved']==sum(n for k,n in states.items() if k!='complete')


def test_stage_completion_rejects_live_owners(tmp_path):
    with pytest.raises(ValueError,match='owners'):
        launch.complete_stage(tmp_path,{},dict(states={'running':1},finished=False))


def test_input_submission_is_separate_from_run(tmp_path,monkeypatch):
    plan=dict(commit='source',python='python',slurm=dict(prepare=dict(qos='short',cpus=2,memory_gib=8,time='01:00:00')))
    monkeypatch.setattr(launch,'initialize',lambda *args:plan)
    submitted=[]
    monkeypatch.setattr(launch,'submit_job',lambda *a,**k:submitted.append((a,k)) or {'job':'100'})
    launch.prepare_inputs(tmp_path,tmp_path/'config.json')
    assert len(submitted)==1 and submitted[0][0][1]=='prepare'
    assert launch.read(tmp_path/'submission.json')['kind']=='inputs'


def test_disk_timing_distributions_match_numpy(tmp_path):
    from rebenchmark.timings import DiskSamples,distribution
    samples=DiskSamples(tmp_path/'samples')
    values=[0.,2.,7.,4.,18.,3.]
    for value in values:samples.add('decode','cpu_seconds',value)
    assert samples.summarize()['decode']['cpu_seconds']==distribution(values)


def test_saved_report_rejects_missing_keys_and_changed_specifications(tmp_path):
    from dataclasses import asdict,replace
    from mappa.mappa_core.core.scheduler import Task,TaskQueue,BundleWriter
    from rebenchmark.report import _event_direction_report
    queue=TaskQueue(tmp_path/'queue.sqlite','test')
    task=Task('00000/cut/00000/verify/00','rebenchmark.operations:verify_seed',{'seed':0},case='00000')
    writer=BundleWriter(tmp_path/'bundles',tmp_path)
    direction=dict(direction=0,case=0,blob=asdict(writer.put([task])),reversed=False)
    queue.backup(tmp_path/'snapshot.sqlite')
    with pytest.raises(ValueError,match='grid'):
        _event_direction_report((str(tmp_path/'snapshot.sqlite'),direction,str(tmp_path)))
    queue.add([replace(task,payload={'seed':1})]);queue.backup(tmp_path/'snapshot.sqlite')
    with pytest.raises(ValueError,match='specification'):
        _event_direction_report((str(tmp_path/'snapshot.sqlite'),direction,str(tmp_path)))
    queue.close()


def test_local_resource_stage_publishes_final_checkpoint(tmp_path,monkeypatch):
    from mappa.mappa_core.core.scheduler import Task,TaskQueue
    import time
    queue=TaskQueue(tmp_path/'live.sqlite','local')
    queue.add([Task('value','execution_helpers:value',{})])
    plan=dict(identity='local',memory_bytes=2**30,claim_until=time.time()+60)
    progress=execution.run_local(queue,tmp_path,plan,1,tmp_path)
    assert progress['finished'] and progress['states']=={'complete':1}
    assert launch.read(tmp_path/'runtime.json')['queue']==str(queue.path)
    assert (tmp_path/'queue.sqlite').exists()
    queue.close()


def test_controller_routes_terminal_queue_to_reporting(tmp_path,monkeypatch):
    from datetime import datetime,timedelta
    monkeypatch.setenv('SLURM_JOB_ID','100')
    full=tmp_path/'full';full.mkdir()
    resource=dict(qos='short',cpus=2,memory_gib=8,time='01:00:00')
    plan=dict(campaign_root=str(tmp_path),identity='test',generation=0,commit='source',python='python',
              case_seconds=60,slurm=dict(checkpoint_seconds=30,report=resource))
    launch.save(full/'plan.json',plan)
    launch.save(full/'endpoint.json',{})
    launch.save(full/'progress.json',dict(finished=True,states={'complete':2,'timeout':1}))
    class Process:
        pid=123
        code=None
        def poll(self):return self.code
        def wait(self):self.code=0;return 0
    monkeypatch.setattr(launch.subprocess,'Popen',lambda *a,**k:Process())
    monkeypatch.setattr(launch.signal,'signal',lambda *a:None)
    monkeypatch.setattr(launch,'scheduler_output',lambda args:'EndTime='+(datetime.now()+timedelta(hours=2)).isoformat())
    submitted=[]
    monkeypatch.setattr(launch,'submit_job',lambda *a,**k:submitted.append((a,k)) or {'job':'101'})
    launch.control(full)
    assert len(submitted)==1 and submitted[0][0][1]=='report'
    assert launch.read(tmp_path/'terminal.json')['unresolved']==1
