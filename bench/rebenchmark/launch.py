"""Submit and supervise the existing benchmark queue from a saved configuration.

Scientific work stays in shared.py and the public MAPPA APIs. Slurm handles
stage dependencies; this module renews retired worker capacity, never tasks.
"""
import argparse
import getpass
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import tempfile
import time

from .decoder_campaign import save

REPO = Path(__file__).resolve().parents[2]


def read(path):
    return json.loads(Path(path).read_text())


def active_jobs():
    # Retired jobs disappear from squeue; their terminal state lives in sacct.
    return set(subprocess.check_output(['squeue','-h','-u',getpass.getuser(),
                                       '-o','%i'],text=True).split())


def wrapper(plan, mode, output, workers=1, module='rebenchmark.shared'):
    return ['bash', str(REPO/'bench/rebenchmark/run_conversion.sbatch'), mode,
            str(output), str(REPO), plan['commit'], plan['python'], str(workers), module]


def submit_job(root, name, plan, resources, command, *, dependency=None):
    """An exclusive receipt prevents duplicate stage or worker submissions."""
    receipt = root/f'{name}-launch.json'
    with receipt.open('x') as stream:
        json.dump(dict(state='submitting', utc=time.time()), stream)
    slurm = plan['slurm']
    args = ['sbatch', '--parsable', '--account='+slurm['account'],
            '--partition='+slurm['partition'], '--job-name=mappa-'+name,
            '--output='+str(root/(name+'-%j.log')), *resources]
    if dependency:
        args.append('--dependency='+dependency)
    result = subprocess.run([*args, *command[1:]], text=True, capture_output=True)
    row = dict(command=result.args, stdout=result.stdout, stderr=result.stderr,
               returncode=result.returncode, utc=time.time())
    if result.returncode == 0:
        row['job'] = result.stdout.strip().split(';')[0]
    save(receipt, row)
    result.check_returncode()
    return row


def stage_resources(plan, stage):
    spec = plan['slurm'][stage]
    return ['--qos='+plan['slurm']['long_qos'], '--nodes=1',
            f'--cpus-per-task={spec["cpus"]}', f'--mem={spec["memory_gib"]}G',
            '--time='+spec['time']]


def prepared_plan(previous, plan):
    """Reuse preparation/pilot only when production never exposed a queue."""
    full = previous/'full'
    if (full/'endpoint.json').exists() or (full/'worker-jobs.json').exists():
        raise ValueError('Production already started; reconcile its queue instead')
    job = read(previous/'production-launch.json')['job']
    if job in active_jobs():
        raise ValueError('Previous production allocation is still active')
    state = subprocess.check_output(['sacct','-X','-n','-P','-j',job,
                                    '--format=State'],text=True).strip().strip('|')
    if state not in ('FAILED','OUT_OF_MEMORY','TIMEOUT','CANCELLED'):
        raise ValueError('Need confirmed failed production retirement: '+state)
    prior = read(full/'plan.json')
    for key in ('dataset','dataset_path','input_manifest_sha256','directions','config',
                'seed_batch','case_seconds','memory_bytes','slack','reference_checks','expected_cases'):
        if prior[key] != plan[key]:
            raise ValueError('Prepared scientific settings differ: '+key)
    changed = subprocess.check_output(['git','-C',str(REPO),'diff','--name-only',
        prior['commit'],plan['commit'],'--','src','bench/rebenchmark/shared.py',
        'bench/rebenchmark/decoding.py','bench/rebenchmark/reference.py'],text=True)
    if changed.strip():
        raise ValueError('Prepared scientific implementation differs: '+changed)
    with Path(prior['index']).open('rb') as stream:
        if hashlib.file_digest(stream,'sha256').hexdigest() != prior['index_sha256']:
            raise ValueError('Prepared index changed')
    for pilot in map(Path,prior['pilots']):
        if not read(pilot/'check.json')['verified'] or not read(pilot/'progress.json')['finished']:
            raise ValueError('Pilot has not passed')
    return dict(plan,index=prior['index'],index_sha256=prior['index_sha256'],
                pilots=prior['pilots'],pilot_commits=prior['pilot_commits'])


def submit(output, configuration, prepared=None):
    plan = read(configuration)
    commit = subprocess.check_output(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], text=True).strip()
    if subprocess.check_output(['git', '-C', str(REPO), 'status', '--porcelain', '--untracked-files=no'], text=True).strip():
        raise ValueError('Commit the source before submission')
    if plan.get('commit', commit) != commit:
        raise ValueError('Requested source differs from checkout')
    plan.update(commit=commit, python=sys.executable, identity=output.name+'-'+commit[:7],
                pilots=[], pilot_commits=[commit])
    if prepared is not None:
        plan = prepared_plan(prepared, plan)
    output.mkdir(parents=True, exist_ok=False)
    save(output/'plan.json', plan)
    binaries = {str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest()
                for pattern in ('_engine*.so', '_group_ops*.so') for p in (REPO/'src/mappa').glob(pattern)}
    save(output/'source.json', dict(commit=commit, native_binaries=binaries,
        plan_sha256=hashlib.sha256((output/'plan.json').read_bytes()).hexdigest()))
    preparation = pilot = None
    if prepared is None:
        preparation = submit_job(output, 'prepare', plan, stage_resources(plan, 'prepare'),
            wrapper(plan, 'prepare', output))
        pilot = submit_job(output, 'pilot', plan, stage_resources(plan, 'pilot'),
            wrapper(plan, 'pilot-check', output, module='rebenchmark.launch'),
            dependency='afterok:'+preparation['job'])
    else:
        full = output/'full'; full.mkdir()
        save(full/'plan.json',dict(plan,identity=plan['identity']+'-full'))
        save(output/'prepared-reuse.json',dict(previous=str(prepared),index=plan['index'],
             index_sha256=plan['index_sha256'],pilots=plan['pilots']))
    production = submit_job(output, 'production', plan, stage_resources(plan, 'coordinator'),
        wrapper(plan, 'control', output, module='rebenchmark.launch'),
        dependency='afterok:'+pilot['job'] if pilot else None)
    result = dict(preparation=preparation['job'] if preparation else None,
                  pilot=pilot['job'] if pilot else None, production=production['job'],
                  run=str(output), source=commit)
    save(output/'submission.json', result)
    print(json.dumps(result), flush=True)
    return result


def pilot_check(root):
    from mappa.mappa_core.core.scheduler import read_blob, run
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint
    from .shared import prepare_queue
    plan = read(root/'plan.json')
    index = read(root/'index.json')
    if index['cases'] != plan['expected_cases'] or len(index['directions']) != plan['expected_cases']*len(plan['directions']):
        raise ValueError('Prepared dataset denominator differs')
    selected = [(r['direction'], cut) for r in index['directions']
                if r['case'] in plan['pilot_cases'] for cut in sorted({0, r['cuts']-1})]
    if {r['case'] for r in index['directions'] if r['case'] in plan['pilot_cases']} != set(plan['pilot_cases']) or not selected:
        raise ValueError('Pilot case selection is incomplete')
    pilot = root/'pilot'; pilot.mkdir(exist_ok=False)
    plan.update(identity=plan['identity']+'-pilot', index=str(root/'index.json'),
        index_sha256=hashlib.sha256((root/'index.json').read_bytes()).hexdigest(), selected_cuts=selected)
    save(pilot/'plan.json', plan)
    with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
        queue = prepare_queue(plan, pilot, scratch)
        try:
            run(queue, workers=plan['pilot_workers'], identity=plan['identity'],
                scratch=scratch, bundles=pilot/'bundles', memory_bytes=plan['memory_bytes'])
            states = queue.snapshot()
            queue.backup(pilot/'queue.sqlite')
            save(pilot/'progress.json', dict(states=states,finished=not any(states.get(s,0) for s in ('pending','running'))))
            if set(states) != {'complete'}:
                raise ValueError('Pilot execution requires diagnosis; production stays gated')
            events = references = representatives = 0
            for specification, raw in queue.db.execute('SELECT specification,result FROM tasks'):
                task, resource = _loads_checkpoint(specification), _loads_checkpoint(raw)
                value = read_blob(resource['blob'], scratch)
                if task.operation.endswith(':convert_seed'):
                    d = value['decoded']
                    if d['schema'] != 'mappa.decoded_events/v3' or not d['mapping_space_preserved']:
                        raise ValueError('Pilot decoder policy differs')
                    if value['finalized_sha256'] != value['finalized'].sha256:
                        raise ValueError('Pilot parent identity differs')
                    events += 1; representatives += len(d['representatives'])
                elif task.operation.endswith(':verify_seed'):
                    if value['scope'] != 'original_iso1_aam':
                        raise ValueError('Pilot reference policy differs')
                    references += 1
            if events != len(selected)*plan['config']['seed_count'] or references != events:
                raise ValueError('Pilot output denominator differs')
            save(pilot/'check.json', dict(verified=True,selected_cuts=selected,states=states,
                event_outputs=events,reference_outputs=references,representatives=representatives,
                event_validation='Independent dense recounts in decode_graph; unknown references remain unknown'))
        finally:
            queue.close()
    full = root/'full'; full.mkdir(exist_ok=False)
    plan.pop('selected_cuts')
    plan.update(identity=read(root/'plan.json')['identity']+'-full', pilots=[str(pilot)])
    save(full/'plan.json', plan)


def capacity(node_output, memory_bytes, max_jobs):
    """Group available nodes into a few memory-fitting sizes without pinning."""
    observed = []
    memory_mib = (memory_bytes+2**20-1)//2**20
    for line in node_output.splitlines():
        fields = dict(re.findall(r'(\w+)=([^ ]+)', line))
        if not fields['NodeName'].startswith('cpu-') or any(s in fields['State'] for s in ('DRAIN','DOWN','FAIL','RESERVED')):
            continue
        count = min(int(fields['CPUTot'])-int(fields['CPUAlloc'])-1,
                    (int(fields['RealMemory'])-int(fields['AllocMem'])-4096)//memory_mib)
        if count >= 1:
            observed.append(dict(node=fields['NodeName'],workers=count))
    if not observed:
        return dict(observed=[],groups=[],compute_workers=0)
    levels = sorted({r['workers'] for r in observed})
    sizes = max(itertools.combinations(levels,min(max_jobs,len(levels))),
                key=lambda group:sum(max((n for n in group if n<=r['workers']),default=0) for r in observed))
    groups = []
    for workers in sizes:
        nodes = sum(max((n for n in sizes if n<=r['workers']),default=0)==workers for r in observed)
        if nodes: groups.append(dict(nodes=nodes,workers_per_node=workers))
    return dict(observed=observed,groups=groups,compute_workers=sum(g['nodes']*g['workers_per_node'] for g in groups))


def start_workers(full, plan, slot, generation, group):
    count, workers = group['nodes'], group['workers_per_node']
    memory_mib = (plan['memory_bytes']+2**20-1)//2**20
    return submit_job(full, f'workers-{slot:02d}-{generation:02d}', plan,
        ['--qos='+plan['slurm']['short_qos'], f'--nodes={count}', f'--ntasks={count}',
         '--ntasks-per-node=1', f'--cpus-per-task={workers+1}',
         f'--mem={workers*memory_mib+4096}M', '--time='+plan['slurm']['worker_time']],
        wrapper(plan, 'workers', full, workers))


def retired_workers(full, job):
    """Never replace an active allocation or conceal a failed worker exit."""
    if job['job'] in active_jobs():
        return None
    accounting = subprocess.check_output(['sacct','-X','-n','-P','-j',job['job'],
                                          '--format=State,ExitCode'],text=True).strip().splitlines()
    if not accounting: return None  # Accounting has not caught up with retirement.
    if accounting != ['COMPLETED|0:0']:
        raise RuntimeError(f'Worker {job["job"]} failed: {accounting}')
    paths = list((full/'retirements').glob(job['job']+'-*.json'))
    if len(paths) != job['group']['nodes']:
        raise RuntimeError('Worker retirement receipts incomplete: '+job['job'])
    return [read(path) for path in paths]


def control(root):
    full = root/'full'; plan = read(full/'plan.json')
    # Set the claim horizon when this allocation starts, not while queued.
    plan['claim_until'] = time.time()+plan['slurm']['claim_hours']*3600
    save(full/'plan.json', plan)
    process = subprocess.Popen(wrapper(plan, 'coordinator', full), start_new_session=True)
    def terminate(signum, frame):
        os.killpg(process.pid, signal.SIGTERM)
        raise SystemExit(128+signum)
    signal.signal(signal.SIGTERM, terminate)
    signal.signal(signal.SIGINT, terminate)
    try:
        while not ((full/'endpoint.json').exists() and (full/'progress.json').exists()):
            if process.poll() is not None: raise RuntimeError('Coordinator exited during preparation')
            time.sleep(5)
        while time.time() < plan['claim_until']:
            available = capacity(subprocess.check_output(['scontrol','show','nodes','-o'],text=True),
                                 plan['memory_bytes'],plan['slurm']['max_worker_jobs'])
            if available['groups']: break
            time.sleep(30)
        else: raise RuntimeError('No worker capacity before claim horizon')
        save(full/'capacity.json',available)
        jobs = [dict(group=g,generation=0,**start_workers(full,plan,i,0,g)) for i,g in enumerate(available['groups'])]
        save(full/'worker-jobs.json',jobs)
        while process.poll() is None:
            progress = read(full/'progress.json')
            if progress.get('finished'): break
            for slot, job in enumerate(jobs):
                receipts = retired_workers(full, job)
                if receipts is None or time.time() >= plan['claim_until']: continue
                if any(r['states'].get('pending',0) or r['states'].get('running',0) for r in receipts):
                    generation = job['generation']+1
                    jobs[slot] = dict(group=job['group'],generation=generation,
                        **start_workers(full,plan,slot,generation,job['group']))
                    save(full/'worker-jobs.json',jobs)
            if time.time() >= plan['claim_until'] and not ({j['job'] for j in jobs} & active_jobs()):
                raise RuntimeError('Claim horizon reached; checkpoint retained for explicit continuation')
            time.sleep(30)
        if process.wait(): raise RuntimeError('Coordinator failed')
        progress = read(full/'progress.json'); save(full/'completion.json',progress)
        if progress['states'].get('memory_limit',0):
            stage = root/'memory-stage'; stage.mkdir(exist_ok=False)
            retry = dict(plan, previous_stage=str(full), previous_memory_bytes=plan['memory_bytes'],
                         memory_bytes=plan['retry_memory_bytes'], claim_until=time.time()+6*3600)
            save(stage/'plan.json',retry)
            submit_job(stage,'memory-stage',retry,stage_resources(plan,'memory_stage'),
                       wrapper(retry,'memory-stage',stage,plan['retry_workers']))
        elif set(progress['states']) == {'complete'}:
            submit_job(root,'report',plan,stage_resources(plan,'report'),
                       wrapper(plan,'report',full,plan['slurm']['report']['cpus'],module='rebenchmark.launch'))
    finally:
        if process.poll() is None:
            os.killpg(process.pid, signal.SIGTERM)
            process.wait()


def report(full, workers):
    from .report import summarize_original_events
    summarize_original_events(full,full.parent/'report',workers)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('submit','pilot-check','control','report'))
    parser.add_argument('output',type=Path)
    parser.add_argument('--plan',type=Path)
    parser.add_argument('--prepared',type=Path,help='Reuse a passed pilot and index after pre-worker startup failure')
    parser.add_argument('--workers',type=int,default=1)
    args = parser.parse_args()
    if args.mode == 'submit':
        if args.plan is None: parser.error('submit needs --plan')
        submit(args.output.resolve(),args.plan,args.prepared)
    elif args.mode == 'pilot-check': pilot_check(args.output)
    elif args.mode == 'control': control(args.output)
    else: report(args.output,args.workers)
