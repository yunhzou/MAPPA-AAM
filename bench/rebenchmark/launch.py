"""Submit and supervise the existing benchmark queue from a saved configuration.

Scientific work stays in operations.py and the public MAPPA APIs. Slurm handles
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
import subprocess
import sys
import time

from .evidence import read, save, sha
from . import execution

REPO = Path(__file__).resolve().parents[2]


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


def wrapper(plan, mode, output, workers=1):
    return ['bash', str(REPO/'bench/rebenchmark/run_conversion.sbatch'), mode,
            str(output), str(REPO), plan['commit'], plan['python'], str(workers)]


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


def initialize(output, configuration):
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
    return plan


def prepare_inputs(output, configuration):
    plan = initialize(output, configuration)
    result = submit_job(output, 'prepare', plan, stage_resources(plan, 'prepare'),
        wrapper(plan, 'prepare', output))
    save(output/'submission.json',dict(preparation=result['job'],source=plan['commit'],kind='inputs'))
    return result


def submit(output, configuration):
    from .inputs import validate
    # Validate before allocating a new run or submitting any jobs.
    validate(read(configuration))
    plan = initialize(output, configuration)
    pilot = submit_job(output, 'pilot', plan, stage_resources(plan, 'pilot'),
        wrapper(plan, 'pilot-check', output))
    prepared_queue = submit_job(output, 'queue', plan, stage_resources(plan, 'queue'),
        wrapper(plan, 'prepare-queue', output/'full'),
        dependency='afterok:'+pilot['job'])
    production = submit_job(output, 'production', plan, stage_resources(plan, 'coordinator'),
        wrapper(plan, 'control', output/'full'),
        dependency='afterok:'+prepared_queue['job'])
    result = dict(pilot=pilot['job'], queue=prepared_queue['job'],
                  production=production['job'], run=str(output), source=plan['commit'])
    save(output/'submission.json', result)
    print(json.dumps(result), flush=True)
    return result


def pilot_check(root):
    from mappa.mappa_core.core.scheduler import read_blob, run
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint
    from .execution import prepare_queue
    plan = read(root/'plan.json')
    from .inputs import validate
    index = validate(plan)
    if index['cases'] != plan['expected_cases'] or len(index['directions']) != plan['expected_cases']*len(plan['directions']):
        raise ValueError('Prepared dataset denominator differs')
    selected = [(r['direction'], cut) for r in index['directions']
                if r['case'] in plan['pilot_cases'] for cut in sorted({0, r['cuts']-1})]
    if {r['case'] for r in index['directions'] if r['case'] in plan['pilot_cases']} != set(plan['pilot_cases']) or not selected:
        raise ValueError('Pilot case selection is incomplete')
    pilot = root/'pilot'; pilot.mkdir(exist_ok=False)
    plan.update(identity=plan['identity']+'-pilot', selected_cuts=selected)
    save(pilot/'plan.json', plan)
    with execution.owned_queue(pilot, plan, prepare_queue) as (queue, scratch):
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
    full = root/'full'; full.mkdir(exist_ok=False)
    plan.pop('selected_cuts')
    plan.update(identity=read(root/'plan.json')['identity']+'-full', pilots=[str(pilot)])
    save(full/'plan.json', plan)


def capacity(node_output, memory_bytes, max_jobs, target_workers):
    """Fit currently unallocated CPU/RAM; the worker target is a ceiling."""
    from functools import lru_cache
    if min(max_jobs, target_workers) < 1:
        raise ValueError('Worker ceiling and job count must be positive')
    memory_mib = (memory_bytes+2**20-1)//2**20
    observed = []
    for line in node_output.splitlines():
        fields = dict(re.findall(r'(\w+)=([^ ]+)', line))
        if not fields['NodeName'].startswith('cpu-') or any(
                state in fields['State'] for state in ('DRAIN', 'DOWN', 'FAIL', 'RESERVED', 'MAINT')):
            continue
        if fields.get('ReservationName', '(null)') not in ('(null)', 'N/A'):
            continue
        cpus = int(fields.get('CPUEfctv', fields['CPUTot']))-int(fields['CPUAlloc'])-1
        memory = int(fields['RealMemory'])-int(fields['AllocMem'])-4096
        count = min(cpus, memory//memory_mib)
        if count > 0:
            observed.append(dict(node=fields['NodeName'], workers=count))
    limits = sorted((row['workers'] for row in observed), reverse=True)

    @lru_cache(None)
    def fit(start, slots, remaining):
        if not slots or not remaining:
            return 0, ()
        best = (0, ())
        for count in range(1, len(limits)-start+1):
            workers = min(limits[start+count-1], remaining//count)
            if not workers:
                break
            used = count*workers
            tail, groups = fit(start+count, slots-1, remaining-used)
            if used+tail > best[0]:
                best = used+tail, ((count, workers),)+groups
        return best

    total, groups = fit(0, max_jobs, target_workers)
    return dict(observed=observed,
                groups=[dict(nodes=n, workers_per_node=w) for n, w in groups],
                target_workers=target_workers, compute_workers=total)


def scale_workers(full, plan, jobs):
    """Reconsider free slots and pending requests; never resize running workers."""
    for job in jobs:
        if not job.get('retired') and retired_workers(full, job) is not None:
            job['retired'] = True
    live = dict(line.split('|', 1) for line in scheduler_output(
        ['squeue', '-h', '-u', getpass.getuser(), '-o', '%i|%T']).splitlines())
    active = [job for job in jobs if not job.get('retired')]
    pending = [job for job in active if live.get(job['job']) == 'PENDING']
    if pending and time.time()-max(job['utc'] for job in pending) < 120:
        return jobs  # Give newly submitted requests time to be scheduled.
    retained = [job for job in active if job not in pending]
    occupied = sum(j['group']['nodes']*j['group']['workers_per_node'] for j in retained)
    slots = plan['slurm']['max_worker_jobs']-len(retained)
    remaining = plan['slurm']['target_workers']-occupied
    available = capacity(scheduler_output(['scontrol', 'show', 'nodes', '-o']),
        plan['memory_bytes'], slots, remaining) if slots > 0 and remaining > 0 else dict(groups=[], compute_workers=0)
    save(full/'capacity.json', dict(available, retained_compute_workers=occupied, utc=time.time()))
    signature = lambda groups: sorted((g['nodes'], g['workers_per_node']) for g in groups)
    if pending:
        if signature(j['group'] for j in pending) != signature(available['groups']):
            for job in pending:
                scheduler_output(['scancel', '--state=PENDING', job['job']])
                job['cancel_pending'] = True
            # Cancellation can race with admission. Recheck ownership and resources
            # on the next cycle before filling the slots.
        return jobs
    vacant = [slot for slot in range(plan['slurm']['max_worker_jobs'])
              if slot >= len(jobs) or jobs[slot].get('retired')]
    for slot, group in zip(vacant, available['groups']):
        generation = jobs[slot]['generation']+1 if slot < len(jobs) else 0
        job = dict(group=group, generation=generation,
                   **start_workers(full, plan, slot, generation, group))
        if slot < len(jobs):
            jobs[slot] = job
        else:
            jobs.append(job)
    return jobs


def start_workers(full, plan, slot, generation, group):
    count, workers = group['nodes'], group['workers_per_node']
    memory_mib = (plan['memory_bytes']+2**20-1)//2**20
    return submit_job(full, f'workers-{slot:02d}-{generation:02d}', plan,
        ['--qos='+plan['slurm']['short_qos'], f'--nodes={count}', f'--ntasks={count}',
         '--ntasks-per-node=1', f'--cpus-per-task={workers+1}',
         f'--mem={workers*memory_mib+4096}M', '--time='+plan['slurm']['worker_time']],
        wrapper(plan, 'workers', full, workers))


class WorkerFailure(RuntimeError):
    """An allocation retired unsuccessfully; preserve the independent queue."""


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
        raise WorkerFailure(f'Worker {job["job"]} failed: {accounting}')
    paths = list((full/'retirements').glob(job['job']+'-*.json'))
    if len(paths) != job['group']['nodes']:
        raise WorkerFailure('Worker retirement receipts incomplete: '+job['job'])
    return [read(path) for path in paths]


def status(root):
    """One compact, read-only check; never inspect result archives or live SQLite."""
    now = time.time()
    submission = read(root/'submission.json')
    allocation = read(root/'current-allocation.json') if (root/'current-allocation.json').exists() else None
    stage = Path(allocation['stage']) if allocation else (root if submission.get('kind') == 'inputs' else root/'full')
    receipts = list(root.glob('*-launch.json')) + list(stage.glob('*-launch.json'))
    next_path = root/'next-allocation.json'
    following = read(next_path) if next_path.exists() else None
    if following:
        receipts += list(Path(following['stage']).glob('*-launch.json'))
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
        scheduler=list(live.values()),next_allocation=following,progress=observed('progress.json'),
        preparation=observed('preparation.json'),runtime=observed('runtime.json'),
        report=read(root/'report/progress.json') if (root/'report/progress.json').exists() else None,
        terminal=read(root/'terminal.json') if (root/'terminal.json').exists() else None,
        checkpoint=dict(bytes=checkpoint.stat().st_size,age_seconds=now-checkpoint.stat().st_mtime)
                   if checkpoint.exists() else None)
    worker_jobs = read(stage/'worker-jobs.json') if (stage/'worker-jobs.json').exists() else []
    result['requested_compute_workers'] = sum(j['group']['nodes']*j['group']['workers_per_node'] for j in worker_jobs if not j.get('retired'))
    result['running_compute_workers'] = sum(j['group']['nodes']*j['group']['workers_per_node'] for j in worker_jobs
        if live.get(str(j['job']),'').split('|')[1:2]==['RUNNING'])
    if allocation and allocation.get('kind') == 'memory':
        result['requested_compute_workers'] = allocation['compute_workers']
        result['running_compute_workers'] = (allocation['compute_workers']
            if live.get(str(allocation['job']),'').split('|')[1:2]==['RUNNING'] else 0)
    result['accounting'] = scheduler_output(['sacct','-X','-n','-P','-j',','.join(jobs),
        '--format=JobID,State,ExitCode,AllocCPUS,Elapsed,Start,End']).splitlines() if jobs else []
    running = [j for j,r in live.items() if r.split('|')[1]=='RUNNING']
    result['resources'] = scheduler_output(['sstat','-a','-n','-P','-j',','.join(running),
        '--format=JobID,AveCPU,MaxRSS,MaxDiskRead,MaxDiskWrite']).splitlines() if running else []
    save(root/'status.json',result)
    return result


def next_allocation(full, plan, states):
    if states.get('running',0) or not states.get('pending',0):
        raise ValueError('Continuation requires pending work with no active tasks')
    receipt = dict(identity=plan['identity'],index_sha256=plan['index_sha256'],
        queue=str(full/'queue.sqlite'),queue_sha256=sha(full/'queue.sqlite'),
        states=states,owner=os.environ['SLURM_JOB_ID'],utc=time.time())
    save(full/'handoff.json',receipt)
    root = Path(plan['campaign_root'])
    following = root/'continuations'/f'{plan["generation"]+1:04d}'
    following.mkdir(parents=True,exist_ok=False)
    continued = dict(plan,generation=plan['generation']+1,resume_receipt=str(full/'handoff.json'))
    continued.pop('claim_until', None)
    save(following/'plan.json',continued)
    submitted = submit_job(following,'production',continued,stage_resources(continued,'memory_stage' if continued.get('memory_stage') else 'coordinator'),
        wrapper(continued,'memory-stage' if continued.get('memory_stage') else 'control',following,
                continued['retry_workers'] if continued.get('memory_stage') else 1),
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
    process = subprocess.Popen(wrapper(plan,'service',full),start_new_session=True)
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
        worker_failure = None
        while process.poll() is None:
            progress=read(full/'progress.json')
            if progress.get('finished'): break
            if worker_failure is not None:
                # Healthy workers retain their original claim deadline. Never let
                # one allocation failure tear down their queue or retry its tasks.
                live = active_jobs()
                for job in jobs:
                    if job['job'] in live and not job.get('cancel_pending'):
                        scheduler_output(['scancel', '--state=PENDING', job['job']])
                        job['cancel_pending'] = True
                if not any(job['job'] in live for job in jobs):
                    from mappa.mappa_core.core.scheduler import connect
                    endpoint = read(full/'endpoint.json')
                    queue = connect((endpoint['host'], endpoint['port']), (full/'worker.key').read_bytes())
                    if not queue.snapshot().get('running', 0):
                        save(full/'retire.json', dict(utc=time.time(), workers=jobs, failure=worker_failure))
                        break
            elif time.time() < plan['claim_until']:
                try:
                    jobs = scale_workers(full, plan, jobs)
                except WorkerFailure as error:
                    worker_failure = str(error)
                    save(full/'worker-failure.json', dict(utc=time.time(), error=worker_failure,
                         workers=jobs, action='Drain existing workers; preserve queue and unresolved attempts'))
                save(full/'worker-jobs.json',jobs)
            else:
                # Atomic Slurm filtering cancels only jobs which have not started.
                live=active_jobs()
                for job in jobs:
                    if job['job'] in live and not job.get('cancel_pending'):
                        scheduler_output(['scancel','--state=PENDING',job['job']])
                        job['cancel_pending']=True
                save(full/'worker-jobs.json',jobs)
                try:
                    retired = all(retired_workers(full,job) is not None for job in jobs)
                except WorkerFailure as error:
                    worker_failure = str(error)
                    save(full/'worker-failure.json', dict(utc=time.time(), error=worker_failure,
                         workers=jobs, action='Drain existing workers; preserve queue and unresolved attempts'))
                    retired = False
                if retired:
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
        complete_stage(full, plan, progress)
    finally:
        if process.poll() is None:
            os.killpg(process.pid,signal.SIGTERM)
            process.wait()


def complete_stage(stage, plan, progress):
    """Every retired queue has one explicit next stage, including partial results."""
    states = progress['states']
    if states.get('running', 0):
        raise ValueError('Stage completion requires retired task owners')
    if states.get('pending', 0):
        return next_allocation(stage, plan, states)
    if not progress['finished']:
        raise ValueError('Terminal stage lacks its final checkpoint receipt')
    root = Path(plan['campaign_root'])
    if states.get('memory_limit', 0) and not plan.get('memory_stage'):
        following = root/'memory-stage'
        following.mkdir(exist_ok=False)
        retry = dict(plan, previous_stage=str(stage), previous_memory_bytes=plan['memory_bytes'],
                     memory_bytes=plan['retry_memory_bytes'], memory_stage=True)
        retry.pop('resume_receipt', None)
        retry.pop('claim_until', None)
        save(following/'plan.json', retry)
        job = submit_job(following, 'memory-stage', retry, stage_resources(plan, 'memory_stage'),
            wrapper(retry, 'memory-stage', following, plan['retry_workers']),
            dependency='afterok:'+os.environ['SLURM_JOB_ID'])
        save(root/'next-allocation.json', dict(stage=str(following),job=job['job'],kind='memory'))
        return job
    # Failed/limited tasks remain in the denominator; reporting must not wait forever.
    save(root/'terminal.json', dict(stage=str(stage), states=states,
        unresolved=sum(n for state,n in states.items() if state != 'complete'), utc=time.time()))
    return submit_job(root, 'report', plan, stage_resources(plan, 'report'),
        wrapper(plan, 'report', stage, plan['slurm']['report']['cpus']),
        dependency='afterok:'+os.environ['SLURM_JOB_ID'])


def memory_stage(stage, workers):
    """Run the authorized resource stage, then use the same completion routing."""
    from datetime import datetime
    from mappa.mappa_core.core.scheduler import run
    plan = read(stage/'plan.json')
    allocation = scheduler_output(['scontrol','show','job',os.environ['SLURM_JOB_ID'],'-o'])
    end = datetime.fromisoformat(re.search(r'EndTime=(\S+)', allocation)[1]).timestamp()
    plan['claim_until'] = end-plan['case_seconds']-180-plan['slurm']['checkpoint_seconds']
    save(stage/'plan.json', plan)
    root = Path(plan['campaign_root'])
    save(root/'current-allocation.json',dict(stage=str(stage),job=os.environ['SLURM_JOB_ID'],kind='memory',compute_workers=workers))
    prepare = execution.restore_queue if 'resume_receipt' in plan else execution.prepare_memory_stage
    with execution.owned_queue(stage, plan, prepare) as (queue, scratch):
        progress = execution.run_local(queue, stage, plan, workers, scratch)
    complete_stage(stage,plan,progress)


def report(full, workers):
    from .report import summarize_original_events
    summarize_original_events(full,Path(read(full/'plan.json')['campaign_root'])/'report',workers)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('prepare-inputs','submit','pilot-check','prepare-queue','control','service','memory-stage','report','status'))
    parser.add_argument('output',type=Path)
    parser.add_argument('--plan',type=Path)
    parser.add_argument('--workers',type=int,default=1)
    args = parser.parse_args()
    if args.mode in ('prepare-inputs','submit'):
        if args.plan is None: parser.error(args.mode+' needs --plan')
        action = prepare_inputs if args.mode == 'prepare-inputs' else submit
        action(args.output.resolve(),args.plan)
    elif args.mode == 'pilot-check': pilot_check(args.output)
    elif args.mode == 'prepare-queue': execution.prepare_full_queue(args.output)
    elif args.mode == 'status': print(json.dumps(status(args.output),indent=2))
    elif args.mode == 'control': control(args.output)
    elif args.mode == 'service': execution.service(args.output)
    elif args.mode == 'memory-stage': memory_stage(args.output,args.workers)
    else: report(args.output,args.workers)
