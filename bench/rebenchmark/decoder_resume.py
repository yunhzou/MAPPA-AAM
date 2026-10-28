"""Continue only never-started decoder tasks after a verified queue retirement."""
import argparse
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import pickle
import shutil
import sqlite3
import subprocess
import tempfile
import time

from mappa.execution import TaskQueue
from .decoder_campaign import save, serve_queue


TERMINAL = {'COMPLETED', 'TIMEOUT', 'FAILED', 'CANCELLED', 'NODE_FAIL', 'OUT_OF_MEMORY', 'PREEMPTED', 'BOOT_FAIL'}


def require_retired(coordinator_job, worker_job):
    """Queued old workers cannot reach the retired endpoint; active owners block."""
    jobs=f'{coordinator_job},{worker_job}'
    accounting=subprocess.run(['sacct','-X','-j',jobs,'--format=JobID,State,End','-n','-P'],
                              check=True,capture_output=True,text=True).stdout
    rows=[line.split('|') for line in accounting.splitlines() if line.strip()]
    parent=[row for row in rows if row[0]==str(coordinator_job)]
    if len(parent)!=1 or parent[0][1].split()[0] not in TERMINAL:
        raise ValueError('Original coordinator has not retired')
    if any(row[1].split()[0] not in TERMINAL|{'PENDING'} for row in rows):
        raise ValueError('Original worker owners remain active')
    live=subprocess.run(['squeue','-h','-r','-j',jobs,'-o','%T'],check=True,capture_output=True,text=True).stdout
    if any(state.strip()!='PENDING' for state in live.splitlines()):
        raise ValueError('Original worker owners remain active')
    return dict(accounting=accounting,queue=live,checked_at=time.time())


def prepare_pending(source, queue, charged_until):
    """Keep every prior outcome outside this queue and carry remaining budgets.

    An attempted task is never eligible, even if its state was reset to pending.
    Interrupted/unknown results require separate diagnosis, not this continuation.
    """
    attempted={row[0] for row in source.execute('SELECT task FROM attempts')}
    cases={name:(seconds,deadline) for name,seconds,deadline in source.execute('SELECT name,seconds,deadline FROM cases')}
    states=dict(source.execute('SELECT state,count(*) FROM tasks GROUP BY state'))
    batch=[];scheduled=0;expired={};directions={}
    for key,data,case,worker,token in source.execute(
            "SELECT key,specification,case_name,worker,token FROM tasks WHERE state='pending' ORDER BY rowid"):
        if key in attempted or worker is not None or token is not None:
            raise ValueError('Pending task has prior ownership: '+key)
        task=pickle.loads(data)
        if task.operation!='rebenchmark.decoding:decode_seed' or task.dependencies or task.key!=key:
            raise ValueError('Unexpected decoder task specification: '+key)
        seconds,deadline=cases[case]
        remaining=seconds if deadline is None else min(seconds,deadline-charged_until)
        if remaining<=0:
            expired[case]=expired.get(case,0)+1
            continue
        batch.append(replace(task,case_seconds=remaining))
        scheduled+=1
        directions.setdefault(case,dict(tasks=0,remaining_seconds=remaining))['tasks']+=1
        if len(batch)==1000:queue.add(batch);batch=[]
    if batch:queue.add(batch)
    assert scheduled+sum(expired.values())==states.get('pending',0)
    return dict(source_states=states,scheduled=scheduled,expired_pending=expired,
                directions=directions,charged_until=charged_until,
                policy='Only never-attempted pending tasks; original remaining directional budgets; all previous outcomes retained in source run')


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(2**20),b''):h.update(block)
    return h.hexdigest()


def execute(out):
    plan=json.loads((out/'plan.json').read_text())
    previous=Path(plan['previous'])
    receipt=json.loads((previous/'retirement-observer/receipt.json').read_text())
    retired=require_retired(plan['previous_coordinator'],plan['previous_workers'])
    source_path=previous/'retirement-observer/queue.sqlite'
    with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
        local=Path(scratch)/'previous.sqlite';shutil.copyfile(source_path,local)
        if digest(local)!=receipt['sha256']:
            raise ValueError('Final queue checksum differs from retirement receipt')
        source=sqlite3.connect('file:'+str(local)+'?mode=ro&immutable=1',uri=True)
        if source.execute('SELECT value FROM identity').fetchone()[0]!=receipt['identity']:
            raise ValueError('Retirement queue identity differs')
        if dict(source.execute('SELECT state,count(*) FROM tasks GROUP BY state'))!=receipt['states']:
            raise ValueError('Retirement states differ')
        queue=TaskQueue(Path(scratch)/'queue.sqlite',plan['identity'])
        summary=prepare_pending(source,queue,retired['checked_at']);source.close()
        save(out/'continuation.json',dict(summary,previous=str(previous),receipt=receipt,ownership=retired))
        queue.optimize()
        serve_queue(queue,out,dict(previous=str(previous)))
        queue.close()



if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('output',type=Path)
    execute(parser.parse_args().output)
