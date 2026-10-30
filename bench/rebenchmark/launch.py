"""Submit and supervise the existing benchmark queue from a saved configuration.

Scientific work stays in shared.py and the public MAPPA APIs. Slurm handles
stage dependencies; this module renews retired worker capacity, never tasks.
"""
import argparse
import getpass
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import socket
import shutil
import subprocess
import sys
import tempfile
import time

from .decoder_campaign import save

REPO = Path(__file__).resolve().parents[2]


def read(path):
    return json.loads(Path(path).read_text())


def scheduler_output(arguments):
    """Retry an idempotent scheduler request when the controller is unreachable.

    A scheduler outage must not tear down the independently running task queue.
    Submissions use submit_job instead: retrying sbatch could duplicate a job.
    """
    while True:
        try:
            return subprocess.check_output(arguments, text=True, stderr=subprocess.PIPE, timeout=60)
        except subprocess.CalledProcessError as error:
            if 'Socket timed out on send/recv operation' not in (error.stderr or ''):
                raise
            message = error.stderr.strip()
        except subprocess.TimeoutExpired:
            message = 'Scheduler request timed out'
        print(f'{message}; retaining the queue and retrying {arguments[0]} in 5 seconds', flush=True)
        time.sleep(5)


def active_jobs():
    # Retired jobs disappear from squeue; their terminal state lives in sacct.
    return set(scheduler_output(['squeue','-h','-u',getpass.getuser(), '-o','%i']).split())


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
    return ['--qos='+spec['qos'], '--nodes=1',
            f'--cpus-per-task={spec["cpus"]}', f'--mem={spec["memory_gib"]}G',
            '--time='+spec['time']]


def submit(output, configuration):
    plan = read(configuration)
    commit = subprocess.check_output(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], text=True).strip()
    if subprocess.check_output(['git', '-C', str(REPO), 'status', '--porcelain', '--untracked-files=no'], text=True).strip():
        raise ValueError('Commit the source before submission')
    if plan.get('commit', commit) != commit:
        raise ValueError('Requested source differs from checkout')
    plan.update(commit=commit, python=sys.executable, identity=output.name+'-'+commit[:7],
                campaign_root=str(output), generation=0,
                pilots=[], pilot_commits=[commit])
    output.mkdir(parents=True, exist_ok=False)
    save(output/'plan.json', plan)
    binaries = {str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest()
                for pattern in ('_engine*.so', '_group_ops*.so') for p in (REPO/'src/mappa').glob(pattern)}
    save(output/'source.json', dict(commit=commit, native_binaries=binaries,
        plan_sha256=hashlib.sha256((output/'plan.json').read_bytes()).hexdigest()))
    preparation = submit_job(output, 'prepare', plan, stage_resources(plan, 'prepare'),
        wrapper(plan, 'prepare', output))
    pilot = submit_job(output, 'pilot', plan, stage_resources(plan, 'pilot'),
        wrapper(plan, 'pilot-check', output, module='rebenchmark.launch'),
        dependency='afterok:'+preparation['job'])
    prepared_queue = submit_job(output, 'queue', plan, stage_resources(plan, 'queue'),
        wrapper(plan, 'prepare-queue', output/'full', module='rebenchmark.launch'),
        dependency='afterok:'+pilot['job'])
    production = submit_job(output, 'production', plan, stage_resources(plan, 'coordinator'),
        wrapper(plan, 'control', output/'full', module='rebenchmark.launch'),
        dependency='afterok:'+prepared_queue['job'])
    result = dict(preparation=preparation['job'], pilot=pilot['job'], queue=prepared_queue['job'],
                  production=production['job'], run=str(output), source=commit)
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


def capacity(node_output, memory_bytes, max_jobs, target_workers):
    """Request the configured throughput; current occupancy is not a ceiling."""
    import math
    if min(max_jobs,target_workers) < 1:
        raise ValueError('Worker target and job count must be positive')
    memory_mib = (memory_bytes+2**20-1)//2**20
    observed = []
    for line in node_output.splitlines():
        fields = dict(re.findall(r'(\w+)=([^ ]+)',line))
        if not fields['NodeName'].startswith('cpu-') or any(
                s in fields['State'] for s in ('DRAIN','DOWN','FAIL','RESERVED')):
            continue
        count = min(int(fields['CPUTot'])-1,(int(fields['RealMemory'])-4096)//memory_mib)
        if count > 0:
            observed.append(dict(node=fields['NodeName'],workers=count))
    # Choose a uniform per-node request that the eligible physical nodes can fit.
    options = []
    for limit in {row['workers'] for row in observed}:
        nodes = math.ceil(target_workers/limit)
        if nodes <= sum(row['workers']>=limit for row in observed):
            workers = math.ceil(target_workers/nodes)
            options.append((nodes*workers-target_workers,nodes,workers))
    if not options:
        raise ValueError('Requested worker target exceeds eligible CPU-node capacity')
    _,nodes,workers = min(options)
    jobs = min(max_jobs,nodes)
    groups = [dict(nodes=nodes//jobs+(slot<nodes%jobs),workers_per_node=workers)
              for slot in range(jobs)]
    return dict(observed=observed,groups=groups,target_workers=target_workers,
                compute_workers=nodes*workers)


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
    accounting = scheduler_output(['sacct','-X','-n','-P','-j',job['job'],
                                  '--format=State,ExitCode']).strip().splitlines()
    if not accounting: return None  # Accounting has not caught up with retirement.
    if job.get('cancel_pending') and len(accounting)==1 and accounting[0].startswith('CANCELLED'):
        started = scheduler_output(['sacct','-X','-n','-P','-j',job['job'],
                                    '--format=Start']).strip().strip('|')
        if started in ('Unknown','None'):
            return []
    if accounting != ['COMPLETED|0:0']:
        raise RuntimeError(f'Worker {job["job"]} failed: {accounting}')
    paths = list((full/'retirements').glob(job['job']+'-*.json'))
    if len(paths) != job['group']['nodes']:
        raise RuntimeError('Worker retirement receipts incomplete: '+job['job'])
    return [read(path) for path in paths]


def file_hash(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def restore_queue(plan, full, scratch):
    """Restore a retired snapshot byte-for-byte; do not replan or retry tasks."""
    from mappa.mappa_core.core.scheduler import TaskQueue
    receipt = read(plan['resume_receipt'])
    if receipt['identity'] != plan['identity'] or receipt['index_sha256'] != plan['index_sha256']:
        raise ValueError('Continuation identity differs')
    if receipt['states'].get('running',0):
        raise ValueError('Continuation has active owners')
    path = Path(scratch)/'queue.sqlite'
    shutil.copyfile(receipt['queue'],path)
    if file_hash(path) != receipt['queue_sha256']:
        raise ValueError('Continuation checkpoint hash differs')
    queue = TaskQueue(path,plan['identity'])
    try:
        states = dict(queue.db.execute('SELECT state,count(*) FROM tasks GROUP BY state'))
        if states != receipt['states'] or queue.db.execute('SELECT count(*) FROM active_work').fetchone()[0]:
            raise ValueError('Continuation queue is not retired')
        save(full/'restored.json',dict(receipt=plan['resume_receipt'],states=states,
            queue_sha256=receipt['queue_sha256'],attempts_and_budgets_preserved=True))
        return queue
    except BaseException:
        queue.close()
        raise


def prepare_full_queue(full):
    """Build and publish the task grid before starting the coordinator clock."""
    from .shared import prepare_queue
    plan = read(full/'plan.json')
    with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
        queue = prepare_queue(plan,full,scratch)
        try:
            states = queue.snapshot()
            queue.backup(full/'queue.sqlite')
        finally:
            queue.close()
    receipt = dict(identity=plan['identity'],index_sha256=plan['index_sha256'],
        queue=str(full/'queue.sqlite'),queue_sha256=file_hash(full/'queue.sqlite'),
        states=states,owner=os.environ.get('SLURM_JOB_ID'),utc=time.time())
    save(full/'prepared-queue.json',receipt)
    plan['resume_receipt'] = str(full/'prepared-queue.json')
    save(full/'plan.json',plan)


def status(root):
    """One compact, read-only check; never inspect result archives or live SQLite."""
    now = time.time()
    submission = read(root/'submission.json')
    allocation = read(root/'current-allocation.json') if (root/'current-allocation.json').exists() else None
    stage = Path(allocation['stage']) if allocation else root/'full'
    receipts = list(root.glob('*-launch.json')) + list(stage.glob('*-launch.json'))
    jobs = sorted({str(read(p)['job']) for p in receipts if read(p).get('job')})
    live_rows = scheduler_output(['squeue','-h','-u',getpass.getuser(),'-o','%i|%T|%C|%M|%R']).splitlines()
    live = {r.split('|')[0]:r for r in live_rows if r.split('|')[0] in jobs}
    def observed(name):
        path = stage/name
        if not path.exists(): return None
        return dict(age_seconds=now-path.stat().st_mtime,**read(path))
    checkpoint = stage/'queue.sqlite'
    result = dict(utc=now,run=str(root),stage=str(stage),source=submission['source'],
        target_compute_workers=read(root/'plan.json')['slurm']['target_workers'],
        scheduler=list(live.values()),progress=observed('progress.json'),
        preparation=observed('preparation.json'),runtime=observed('runtime.json'),
        checkpoint=dict(bytes=checkpoint.stat().st_size,age_seconds=now-checkpoint.stat().st_mtime)
                   if checkpoint.exists() else None)
    worker_jobs = read(stage/'worker-jobs.json') if (stage/'worker-jobs.json').exists() else []
    result['requested_compute_workers'] = sum(j['group']['nodes']*j['group']['workers_per_node'] for j in worker_jobs)
    result['running_compute_workers'] = sum(j['group']['nodes']*j['group']['workers_per_node'] for j in worker_jobs
        if live.get(str(j['job']),'').split('|')[1:2]==['RUNNING'])
    result['accounting'] = scheduler_output(['sacct','-X','-n','-P','-j',','.join(jobs),
        '--format=JobID,State,ExitCode,AllocCPUS,Elapsed,Start,End']).splitlines() if jobs else []
    running = [j for j,r in live.items() if r.split('|')[1]=='RUNNING']
    result['resources'] = scheduler_output(['sstat','-a','-n','-P','-j',','.join(running),
        '--format=JobID,AveCPU,MaxRSS,MaxDiskRead,MaxDiskWrite']).splitlines() if running else []
    save(root/'status.json',result)
    return result


def service(full):
    from .decoder_campaign import serve_queue
    plan = read(full/'plan.json')
    def terminate(signum, frame):
        raise SystemExit(128+signum)
    signal.signal(signal.SIGTERM,terminate)
    scratch = Path(tempfile.mkdtemp(dir=os.environ['TMPDIR']))
    queue = restore_queue(plan,full,scratch)
    try:
        serve_queue(queue,full,dict(campaign='matching and original-space events'),
                    retire_file=full/'retire.json')
    except BaseException:
        save(full/'retained-scratch.json',dict(path=str(scratch),host=socket.gethostname(),
            utc=time.time(),reason='Service failed; retain local database and WAL for diagnosis'))
        raise
    finally:
        queue.close()
    shutil.rmtree(scratch)


def next_allocation(full, plan, states):
    if states.get('running',0) or not states.get('pending',0):
        raise ValueError('Continuation requires pending work with no active tasks')
    receipt = dict(identity=plan['identity'],index_sha256=plan['index_sha256'],
        queue=str(full/'queue.sqlite'),queue_sha256=file_hash(full/'queue.sqlite'),
        states=states,owner=os.environ['SLURM_JOB_ID'],utc=time.time())
    save(full/'handoff.json',receipt)
    root = Path(plan['campaign_root'])
    following = root/'continuations'/f'{plan["generation"]+1:04d}'
    following.mkdir(parents=True,exist_ok=False)
    continued = dict(plan,generation=plan['generation']+1,resume_receipt=str(full/'handoff.json'))
    continued.pop('claim_until')
    save(following/'plan.json',continued)
    submitted = submit_job(following,'production',continued,stage_resources(continued,'coordinator'),
        wrapper(continued,'control',following,module='rebenchmark.launch'),
        dependency='afterok:'+os.environ['SLURM_JOB_ID'])
    save(root/'next-allocation.json',dict(stage=str(following),job=submitted['job'],
        generation=continued['generation'],handoff=str(full/'handoff.json')))
    return submitted


def control(full):
    from datetime import datetime
    plan = read(full/'plan.json');root = Path(plan['campaign_root'])
    allocation = scheduler_output(['scontrol','show','job',os.environ['SLURM_JOB_ID'],'-o'])
    end = datetime.fromisoformat(re.search(r'EndTime=(\S+)',allocation)[1]).timestamp()
    # Reserve the longest remaining task, publication, and a disk checkpoint.
    plan['claim_until'] = end-plan['case_seconds']-180-plan['slurm']['checkpoint_seconds']
    save(full/'plan.json',plan)
    save(root/'current-allocation.json',dict(stage=str(full),job=os.environ['SLURM_JOB_ID'],
        generation=plan['generation'],claim_until=plan['claim_until'],allocation_end=end))
    process = subprocess.Popen(wrapper(plan,'service',full,module='rebenchmark.launch'),start_new_session=True)
    def terminate(signum, frame):
        os.killpg(process.pid,signal.SIGTERM)
        raise SystemExit(128+signum)
    signal.signal(signal.SIGTERM,terminate)
    signal.signal(signal.SIGINT,terminate)
    try:
        while not ((full/'endpoint.json').exists() and (full/'progress.json').exists()):
            if process.poll() is not None: raise RuntimeError('Coordinator exited during preparation')
            time.sleep(5)
        jobs=[]
        while process.poll() is None:
            progress=read(full/'progress.json')
            if progress.get('finished'): break
            if time.time() < plan['claim_until']:
                if not jobs:
                    available=capacity(scheduler_output(['scontrol','show','nodes','-o']),
                        plan['memory_bytes'],plan['slurm']['max_worker_jobs'],plan['slurm']['target_workers'])
                    save(full/'capacity.json',available)
                    jobs=[dict(group=g,generation=0,**start_workers(full,plan,i,0,g))
                          for i,g in enumerate(available['groups'])]
                else:
                    for slot,job in enumerate(jobs):
                        if retired_workers(full,job) is not None:
                            generation=job['generation']+1
                            jobs[slot]=dict(group=job['group'],generation=generation,
                                **start_workers(full,plan,slot,generation,job['group']))
                save(full/'worker-jobs.json',jobs)
            else:
                # Atomic Slurm filtering cancels only jobs which have not started.
                live=active_jobs()
                for job in jobs:
                    if job['job'] in live and not job.get('cancel_pending'):
                        scheduler_output(['scancel','--state=PENDING',job['job']])
                        job['cancel_pending']=True
                save(full/'worker-jobs.json',jobs)
                if all(retired_workers(full,job) is not None for job in jobs):
                    from mappa.mappa_core.core.scheduler import connect
                    endpoint=read(full/'endpoint.json')
                    queue=connect((endpoint['host'],endpoint['port']),(full/'worker.key').read_bytes())
                    if queue.snapshot().get('running',0):
                        raise RuntimeError('Workers retired with unresolved running task owners')
                    save(full/'retire.json',dict(utc=time.time(),workers=jobs))
                    break
            time.sleep(30)
        if process.wait(): raise RuntimeError('Coordinator failed')
        progress=read(full/'progress.json');save(full/'completion.json',progress)
        if progress['states'].get('pending',0):
            next_allocation(full,plan,progress['states'])
        elif progress['states'].get('memory_limit',0):
            stage=root/'memory-stage';stage.mkdir(exist_ok=False)
            retry=dict(plan,previous_stage=str(full),previous_memory_bytes=plan['memory_bytes'],
                       memory_bytes=plan['retry_memory_bytes'],claim_until=time.time()+6*3600)
            save(stage/'plan.json',retry)
            submit_job(stage,'memory-stage',retry,stage_resources(plan,'memory_stage'),
                       wrapper(retry,'memory-stage',stage,plan['retry_workers']))
        elif set(progress['states']) == {'complete'}:
            submit_job(root,'report',plan,stage_resources(plan,'report'),
                       wrapper(plan,'report',full,plan['slurm']['report']['cpus'],module='rebenchmark.launch'))
    finally:
        if process.poll() is None:
            os.killpg(process.pid,signal.SIGTERM)
            process.wait()


def report(full, workers):
    from .report import summarize_original_events
    summarize_original_events(full,Path(read(full/'plan.json')['campaign_root'])/'report',workers)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('submit','pilot-check','prepare-queue','control','service','report','status'))
    parser.add_argument('output',type=Path)
    parser.add_argument('--plan',type=Path)
    parser.add_argument('--workers',type=int,default=1)
    args = parser.parse_args()
    if args.mode == 'submit':
        if args.plan is None: parser.error('submit needs --plan')
        submit(args.output.resolve(),args.plan)
    elif args.mode == 'pilot-check': pilot_check(args.output)
    elif args.mode == 'prepare-queue': prepare_full_queue(args.output)
    elif args.mode == 'status': print(json.dumps(status(args.output),indent=2))
    elif args.mode == 'control': control(args.output)
    elif args.mode == 'service': service(args.output)
    else: report(args.output,args.workers)
