"""One owner for task queues, worker execution and consistent checkpoints."""
from contextlib import contextmanager
import os
from pathlib import Path
import resource
import signal
import secrets
import shutil
import socket
import tempfile
import time
from mappa.mappa_core.core.scheduler import coordinator
from .evidence import read, save, sha


def prepare_queue(plan, out, scratch):
    """Reuse only verified pilot tasks from the identical frozen task grid."""
    from dataclasses import asdict
    import sqlite3
    import time
    from mappa.mappa_core.core.scheduler import Blob, TaskQueue, read_blob
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint
    path = Path(plan['index'])
    if sha(path) != plan['index_sha256']:
        raise ValueError('Frozen task index changed')
    index = read(path)
    queue = TaskQueue(Path(scratch)/'queue.sqlite', plan['identity'])
    # Bulk planning has no live attempts; the published backup is its commit point.
    queue.db.execute('PRAGMA synchronous=NORMAL')
    selected = {tuple(x) for x in plan['selected_cuts']} if 'selected_cuts' in plan else None
    expected = 0
    for row in index['directions']:
        if selected is not None and not any(d == row['direction'] for d, _ in selected):
            continue
        tasks = read_blob(Blob(**row['blob']), scratch)
        if selected is not None:
            tasks = [t for t in tasks if (row['direction'],int(t.key.split('/')[2])) in selected]
        queue.add(tasks)
        expected += len(tasks)
        if row['direction'] % 100 == 0:
            save(out/'preparation.json',dict(direction=row['direction'],tasks=expected,utc=time.time()))
    queue.db.execute('PRAGMA synchronous=FULL')
    reused, charged = [], {}
    for root in map(Path, plan['pilots']):
        prior_plan = read(root/'plan.json')
        if prior_plan['commit'] not in plan['pilot_commits']:
            raise ValueError('Pilot commit is not explicitly accepted')
        for key in ('index_sha256', 'memory_bytes'):
            if prior_plan[key] != plan[key]:
                raise ValueError(f'Pilot {key} differs')
        if not read(root/'progress.json').get('finished'):
            raise ValueError('Pilot still owns live work')
        prior = sqlite3.connect('file:'+str(root/'queue.sqlite')+'?mode=ro', uri=True)
        try:
            if prior.execute("SELECT count(*) FROM tasks WHERE state!='complete'").fetchone()[0]:
                raise ValueError('Diagnose failed pilot tasks before production')
            for key, specification, result in prior.execute('SELECT key,specification,result FROM tasks'):
                expected_row = queue.db.execute('SELECT specification,state FROM tasks WHERE key=?',(key,)).fetchone()
                if expected_row != (specification, 'pending'):
                    raise ValueError(f'Pilot task differs or is duplicated: {key}')
                resource = _loads_checkpoint(result)
                # Verify bytes and successful deserialization, including graph objects.
                read_blob(resource['blob'], scratch)
                queue.db.execute("UPDATE tasks SET state='complete',result=? WHERE key=?",(result,key))
                reused.append(dict(key=key, pilot=str(root), blob=asdict(resource['blob'])))
            for name, spent in prior.execute('SELECT name,spent FROM cases'):
                charged[name] = charged.get(name,0.) + spent
        finally:
            prior.close()
    for name, seconds in charged.items():
        # An exhausted pilot direction remains exhausted; never reset its allowance.
        queue.db.execute('UPDATE cases SET seconds=max(0,seconds-?) WHERE name=?',(seconds,name))
    queue.db.commit()
    queue.optimize()
    save(out/'reused.json', reused)
    save(out/'inventory.json',dict(expected=expected, reused=len(reused), charged_direction_seconds=charged,
        index_sha256=plan['index_sha256'], counts=dict(queue.db.execute('SELECT state,count(*) FROM tasks GROUP BY state'))))
    return queue


def prepare_repeat_queue(plan, out, scratch):
    """A new timed experiment reuses completed stages and repeats explicit ones.

    The original run and its attempts remain immutable. This is an authorized
    benchmark repetition with fresh budgets, not recovery of failed work.
    """
    import sqlite3
    from collections import Counter
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint
    source = Path(plan['repeat_from'])
    prior_plan = read(source/'plan.json')
    progress = read(source/'progress.json')
    if not progress.get('finished') or progress['states'].keys() != {'complete'}:
        raise ValueError('Repeat requires a completed source experiment')
    if prior_plan['index_sha256'] != plan['index_sha256'] or plan['pilots']:
        raise ValueError('Repeat requires the same complete task grid without pilot reuse')
    digest = sha(source/'queue.sqlite')
    if digest != plan['source_queue_sha256']:
        raise ValueError('Source queue hash differs')
    prior = sqlite3.connect('file:'+str(source/'queue.sqlite')+'?mode=ro', uri=True)
    local = sqlite3.connect(Path(scratch)/'source.sqlite')
    try:
        prior.backup(local)
    finally:
        prior.close()
    queue = prepare_queue(plan, out, scratch)
    repeated, reused = Counter(), Counter()
    operations = set(plan['repeat_operations'])
    try:
        for key, specification, state, result in local.execute(
                'SELECT key,specification,state,result FROM tasks'):
            expected = queue.db.execute('SELECT specification FROM tasks WHERE key=?',(key,)).fetchone()
            if expected != (specification,) or state != 'complete' or result is None:
                raise ValueError('Source task differs or is unfinished: '+key)
            operation = _loads_checkpoint(specification).operation
            if operation in operations:
                repeated[operation] += 1
            else:
                queue.db.execute("UPDATE tasks SET state='complete',result=? WHERE key=?",(result,key))
                reused[operation] += 1
        count = queue.db.execute('SELECT count(*) FROM tasks').fetchone()[0]
        if sum(repeated.values())+sum(reused.values()) != count or set(repeated) != operations or not operations:
            raise ValueError('Repeat operation selection or source denominator differs')
        queue.db.commit()
        queue.optimize()
        save(out/'inventory.json',dict(expected=count, repeated=dict(repeated), reused=dict(reused),
            source_queue=str(source/'queue.sqlite'), source_queue_sha256=digest,
            index_sha256=plan['index_sha256'], fresh_experiment_budgets=True,
            counts=queue.snapshot()))
        return queue
    except BaseException:
        queue.close()
        raise
    finally:
        local.close()


def prepare_memory_stage(plan, out, scratch):
    """Continue only authorized memory failures from a retired first stage."""
    import sqlite3
    from mappa.mappa_core.core.scheduler import TaskQueue
    source = Path(plan['previous_stage'])
    progress = read(source/'progress.json')
    if not progress.get('finished') or any(progress['states'].get(k, 0) for k in ('pending', 'running')):
        raise ValueError('Previous memory stage still owns unfinished work')
    path = source/'queue.sqlite'
    digest = sha(path)
    prior = sqlite3.connect('file:'+str(path)+'?mode=ro', uri=True)
    copied = sqlite3.connect(Path(scratch)/'queue.sqlite')
    try:
        prior.backup(copied)
    finally:
        copied.close()
        prior.close()
    queue = TaskQueue(Path(scratch)/'queue.sqlite', plan['identity'])
    keys = queue.retry_memory_limits(plan['previous_memory_bytes'], plan['memory_bytes'],
        reason='User-authorized larger-memory second stage; preserve prior attempts and deadlines')
    save(out/'selection.json', dict(keys=keys, source_queue=str(path), source_sha256=digest,
        previous_memory_bytes=plan['previous_memory_bytes'], memory_bytes=plan['memory_bytes'],
        states=queue.snapshot(), deadlines_preserved=True))
    return queue


def execute(args):
    import os
    import signal
    from mappa.mappa_core.core.scheduler import connect, run
    out = args.output
    plan = read(out/'plan.json')
    if args.mode == 'worker':
        endpoint = read(out/'endpoint.json')
        queue = connect((endpoint['host'],endpoint['port']),(out/'worker.key').read_bytes())
        claim_until = None
        if 'SLURM_JOB_ID' in os.environ:
            import datetime, re, time
            from .launch import scheduler_output
            job = scheduler_output(['scontrol','show','job',os.environ['SLURM_JOB_ID'],'-o'])
            end = datetime.datetime.fromisoformat(re.search(r'EndTime=(\S+)',job)[1]).timestamp()
            claim_until = min(end - plan['case_seconds'] - 180, plan['claim_until'])
        states = run(queue, workers=args.workers, identity=plan['identity'], scratch=os.environ['TMPDIR'],
            bundles=out/'bundles', memory_bytes=plan['memory_bytes'], claim_until=claim_until)
        if claim_until is not None:
            folder = out/'retirements'; folder.mkdir(exist_ok=True)
            save(folder/f'{os.environ["SLURM_JOB_ID"]}-{socket.gethostname()}-{os.getpid()}.json',
                dict(states=states, claim_until=claim_until, utc=time.time(), allocation_end=end,
                     memory_bytes=plan['memory_bytes'], workers=args.workers))
        return


def serve_queue(queue,out,metadata,*,retire_file=None):
    """Publish one endpoint and consistent snapshots for a prepared queue."""
    secret=secrets.token_bytes(32)
    with (out/'worker.key').open('xb') as stream:
        os.chmod(out/'worker.key',0o600);stream.write(secret)
    queue.backup(out/'queue.sqlite')
    with coordinator(queue,(socket.gethostname(),0),secret) as address:
        save(out/'endpoint.json',dict(host=address[0],port=address[1],identity=queue.identity))
        publish_progress(queue,out,metadata,retired=lambda: retire_file is not None and retire_file.exists())
        time.sleep(5)


def publish_progress(queue, output, metadata, *, retired):
    last=time.monotonic()
    try:
        while True:
            states=queue.snapshot();stamp=time.monotonic()
            save(output/'runtime.json',dict(utc=time.time(),host=socket.gethostname(),pid=os.getpid(),
                queue=str(queue.path),disk_free_bytes=shutil.disk_usage(queue.path.parent).free,
                process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
            if stamp-last>=300:
                queue.backup(output/'queue.sqlite');last=stamp
            save(output/'progress.json',dict(metadata,states=states,utc=time.time()))
            if not any(states.get(k,0) for k in ('pending','running')):break
            if retired():
                if states.get('running',0):raise ValueError('Cannot retire a queue with active tasks')
                break
            time.sleep(10)
    finally:
        queue.backup(output/'queue.sqlite')
        states=queue.snapshot()
        save(output/'progress.json',dict(metadata,states=states,utc=time.time(),
            finished=not any(states.get(k,0) for k in ('pending','running'))))


def run_local(queue, output, plan, workers, scratch):
    """Use the same publisher for an allocation containing its own workers."""
    from concurrent.futures import ThreadPoolExecutor
    from mappa.mappa_core.core.scheduler import run
    queue.backup(output/'queue.sqlite')
    with ThreadPoolExecutor(max_workers=1) as pool:
        work=pool.submit(run,queue,workers=workers,identity=plan['identity'],scratch=scratch,
            bundles=output/'bundles',memory_bytes=plan['memory_bytes'],claim_until=plan['claim_until'])
        publish_progress(queue,output,{},retired=work.done)
        work.result()
    return read(output/'progress.json')


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
    if sha(path) != receipt['queue_sha256']:
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
    plan = read(full/'plan.json')
    with owned_queue(full, plan, prepare_queue) as (queue, scratch):
        states = queue.snapshot()
        queue.backup(full/'queue.sqlite')
    receipt = dict(identity=plan['identity'],index_sha256=plan['index_sha256'],
        queue=str(full/'queue.sqlite'),queue_sha256=sha(full/'queue.sqlite'),
        states=states,owner=os.environ.get('SLURM_JOB_ID'),utc=time.time())
    save(full/'prepared-queue.json',receipt)
    plan['resume_receipt'] = str(full/'prepared-queue.json')
    save(full/'plan.json',plan)




@contextmanager
def owned_queue(output, plan, prepare):
    scratch = Path(tempfile.mkdtemp(dir=os.environ['TMPDIR']))
    queue = None
    try:
        queue = prepare(plan, output, scratch)
        yield queue, scratch
    except BaseException:
        save(output/'retained-scratch.json',dict(path=str(scratch),host=socket.gethostname(),
            utc=time.time(),reason='Queue owner failed; retain database and WAL'))
        raise
    else:
        queue.close()
        queue = None
        shutil.rmtree(scratch)
    finally:
        if queue is not None:
            queue.close()


def service(full):
    plan = read(full/'plan.json')
    def terminate(signum, frame):
        raise SystemExit(128+signum)
    signal.signal(signal.SIGTERM,terminate)
    with owned_queue(full, plan, restore_queue) as (queue, scratch):
        serve_queue(queue,full,dict(campaign='matching and original-space events'),
                    retire_file=full/'retire.json')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('worker',))
    parser.add_argument('output', type=Path)
    parser.add_argument('--workers', type=int, default=1)
    execute(parser.parse_args())
