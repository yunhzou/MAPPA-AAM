"""Coordinate saved-graph event decoding through the common CPU task queue."""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
import os
from pathlib import Path
import pickle
import secrets
import shutil
import resource
import signal
import socket
import sqlite3
import tempfile
import time
from mappa.mappa_core.core.scheduler import Blob, TaskQueue, read_blob, connect, coordinator, run
from rebenchmark.decoding import seed_tasks
from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint


def save(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n')
    temporary.replace(path)


def execute(args):
    out=args.output
    plan=json.loads((out/'plan.json').read_text())
    if plan.get('decoder') != 'anchored_fragment_conversion' or 'family_seconds' in plan:
        raise ValueError('A new anchored-fragment conversion plan is required')
    if args.mode=='worker':
        endpoint=json.loads((out/'endpoint.json').read_text())
        queue=connect((endpoint['host'],endpoint['port']),(out/'worker.key').read_bytes())
        run(queue,workers=args.workers,identity=plan['identity'],scratch=os.environ['TMPDIR'],
            bundles=out/'bundles',memory_bytes=plan['memory_bytes'])
        return
    def terminate(signum, frame):
        raise SystemExit(128 + signum)
    signal.signal(signal.SIGTERM, terminate)
    with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
        reused={};charged={}
        for pilot in plan['pilots']:
            root=Path(pilot)
            progress=json.loads((root/'progress.json').read_text())
            assert progress.get('finished') and set(progress['states'])=={'complete'},progress
            db=sqlite3.connect('file:'+str(root/'queue.sqlite')+'?mode=ro',uri=True)
            for case,spent in db.execute('SELECT name,spent FROM cases'):
                charged[case]=charged.get(case,0.)+spent
            for key,specification,state,result in db.execute('SELECT key,specification,state,result FROM tasks'):
                assert key not in reused and state=='complete'
                resource=_loads_checkpoint(result)
                value=read_blob(resource['blob'],scratch)
                assert value['search_rerun'] is False
                reused[key]=(_loads_checkpoint(specification),dict(state=state,resource=dict(resource,blob=asdict(resource['blob'])),
                    pilot=str(root),complete=value['decoded']['complete'],minimum=value['decoded']['observed_minimum']))
            db.close()
        index_path=Path(plan['index'])
        assert hashlib.sha256(index_path.read_bytes()).hexdigest()==plan['index_sha256']
        index=json.loads(index_path.read_text())
        queue=TaskQueue(Path(scratch)/'queue.sqlite',plan['identity'])
        total=0;verified=set();exhausted=[];directions=[]
        selected = {tuple(pair) for pair in plan['selected_cuts']} if 'selected_cuts' in plan else None
        for row in index['tasks']:
            if selected is not None and not any(task == row['task'] for task, _ in selected):
                continue
            cuts=read_blob(Blob(**row['blob']),scratch)
            if selected is not None:
                cuts=[cut for cut in cuts if (row['task'],cut.payload['cut']) in selected]
            tasks=list(seed_tasks(cuts,slack=plan['slack'],branch_limit=plan['branch_limit'],random_seed=plan['random_seed'],
                case_seconds=plan['case_seconds'],memory_bytes=plan['memory_bytes']))
            pending=[]
            for task in tasks:
                total+=1
                if task.key in reused:
                    assert task==reused[task.key][0],task.key
                    verified.add(task.key)
                else:
                    remaining=plan['case_seconds']-charged.get(task.case,0.)
                    if remaining<=0:exhausted.append(task.key)
                    else:pending.append(replace(task,case_seconds=remaining))
            queue.add(pending)
            directions.append(dict(task=row['task'],cuts=len(cuts),expected=len(tasks),pending=len(pending),
                prior_charged_wall_seconds=charged.get(f'{row["task"]:05d}',0.)))
            if len(directions)%100==0:
                save(out/'preparation.json',dict(directions=len(directions),tasks=total,reused=len(verified)))
        assert verified==set(reused)
        save(out/'reused.json',{k:v[1] for k,v in reused.items()})
        save(out/'inventory.json',dict(total=total,reused=len(reused),exhausted=exhausted,directions=directions))
        queue.optimize()
        serve_queue(queue,out,dict(reused=len(reused),exhausted=len(exhausted)))
        queue.close()


def serve_queue(queue,out,metadata,*,retire_file=None):
    """Publish one endpoint and consistent snapshots for a prepared queue."""
    secret=secrets.token_bytes(32)
    with (out/'worker.key').open('xb') as stream:
        os.chmod(out/'worker.key',0o600);stream.write(secret)
    queue.backup(out/'queue.sqlite')
    with coordinator(queue,(socket.gethostname(),0),secret) as address:
        save(out/'endpoint.json',dict(host=address[0],port=address[1],identity=queue.identity))
        last=time.monotonic()
        try:
            while True:
                states=queue.snapshot();stamp=time.monotonic()
                save(out/'runtime.json',dict(utc=time.time(),host=socket.gethostname(),pid=os.getpid(),
                    queue=str(queue.path),disk_free_bytes=shutil.disk_usage(queue.path.parent).free,
                    process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
                if stamp-last>=300:
                    queue.backup(out/'queue.sqlite');last=stamp
                save(out/'progress.json',dict(metadata,states=states,utc=time.time()))
                if not any(states.get(k,0) for k in ('pending','running')):break
                if retire_file is not None and retire_file.exists():
                    if states.get('running',0):
                        raise ValueError('Cannot retire a queue with active tasks')
                    break
                time.sleep(10)
        finally:
            queue.backup(out/'queue.sqlite')
            states=queue.snapshot()
            save(out/'progress.json',dict(metadata,states=states,utc=time.time(),
                finished=not any(states.get(k,0) for k in ('pending','running'))))
        time.sleep(5)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['coordinator','worker'])
    parser.add_argument('output',type=Path)
    parser.add_argument('--workers',type=int,default=1)
    execute(parser.parse_args())
