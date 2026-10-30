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
                campaign_root=str(output), generation=0,
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
        wrapper(plan, 'control', output/'full', module='rebenchmark.launch'),
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


def recover_linked_outputs(queue_path, members, scratch):
    """Recover task identities carried by saved decodes and their parent hashes.

    This creates an evidence overlay, not a runnable queue. Unlabelled search
    and reference records remain unresolved; graph similarity is not identity.
    """
    import sqlite3
    from dataclasses import asdict
    from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob, read_blobs
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint
    db = sqlite3.connect(f'file:{Path(queue_path).resolve()}?mode=ro&immutable=1',uri=True)
    blobs = {Blob(**{k:m[k] for k in ('path','offset','size','sha256')}):m for m in members}
    recovered = {}
    def specification(key, operation):
        row = db.execute('SELECT specification,state,result FROM tasks WHERE key=?',(key,)).fetchone()
        if row is None:
            raise ValueError('Saved output has no expected task: '+key)
        task = _loads_checkpoint(row[0])
        if task.operation != operation:
            raise ValueError('Saved output operation differs: '+key)
        return task, row
    def record(key, blob, row, provenance):
        if row[1] == 'complete':
            if _loads_checkpoint(row[2])['blob'] != blob:
                raise ValueError('Saved output conflicts with accepted result: '+key)
            return
        if row[1] not in ('pending','running'):
            raise ValueError('Unexpected recorded task state: '+key)
        value = dict(key=key,blob=asdict(blob),specification_sha256=hashlib.sha256(row[0]).hexdigest(),
                     snapshot_state=row[1],provenance=provenance,
                     attempt_resource_accounting='not present in the surviving checkpoint')
        if key in recovered and recovered[key]['blob'] != value['blob']:
            raise ValueError('Multiple saved outputs for one task: '+key)
        recovered[key] = value
    try:
        decodes = sorted((b for b,m in blobs.items() if m['kind']=='decode'),key=lambda b:(b.path,b.offset))
        seen = set()
        for blob,value in read_blobs(decodes,scratch):
            key = f'{value["task_id"]:05d}/cut/{value["cut"]:05d}/seed/{value["seed"]-1:02d}'
            if key in seen:
                raise ValueError('Duplicate decoder task identity: '+key)
            seen.add(key)
            task,row = specification(key,'rebenchmark.shared:convert_seed')
            p = task.payload
            if (p['input'] != value['input'] or p['task_id'] != value['task_id']
                    or p['number'] != value['cut'] or p['seed'] != value['seed']-1):
                raise ValueError('Saved decoder input identity differs: '+key)
            parent = value['finalized']
            if parent.sha256 != value['finalized_sha256'] or task.dependencies != (p['cut'],):
                raise ValueError('Saved parent identity differs: '+key)
            final,final_row = specification(p['cut'],'mappa.mappa_core.core.scheduler.search:finalize_cut')
            if final.payload['input'] != p['input'] or final.payload['number'] != p['number']:
                raise ValueError('Finalized task input differs: '+key)
            if final_row[1] != 'complete' and (parent not in blobs or blobs[parent]['kind']!='finalize'):
                raise ValueError('Referenced parent is absent from recovery inventory: '+key)
            if p['cut'] not in recovered:
                saved = read_blob(parent,scratch)
                if (saved['cut'] != final.payload['number'] or saved['seeds'] != final.payload['seeds']
                        or any(tuple(c.cuts)!=tuple(final.payload['cut']) for c in saved['graph'].contexts)):
                    raise ValueError('Saved finalized specification differs: '+key)
            record(p['cut'],parent,final_row,dict(decoder_key=key,relation='exact finalized blob reference'))
            record(key,blob,row,dict(relation='saved task/seed/cut/input identity',finalized_key=p['cut']))
        return list(recovered.values())
    finally:
        db.close()


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


def service(full):
    from .shared import prepare_queue
    from .decoder_campaign import serve_queue
    plan = read(full/'plan.json')
    def terminate(signum, frame):
        raise SystemExit(128+signum)
    signal.signal(signal.SIGTERM,terminate)
    with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
        queue = (restore_queue(plan,full,scratch) if 'resume_receipt' in plan
                 else prepare_queue(plan,full,scratch))
        try:
            serve_queue(queue,full,dict(campaign='matching and original-space events'),
                        retire_file=full/'retire.json')
        finally:
            queue.close()


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
    parser.add_argument('mode', choices=('submit','pilot-check','control','service','report'))
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
    elif args.mode == 'service': service(args.output)
    else: report(args.output,args.workers)
