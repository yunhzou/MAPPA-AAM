"""Frozen bidirectional benchmark: prepare, pilot, search, analyze, aggregate.

Search uses the existing public checkpoint pipeline. Every original row stays
in the denominator. Analysis and decoding consume saved evidence; they never
alter search or replace incomplete work with a negative.
"""
import argparse
from dataclasses import asdict
import fcntl
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import traceback
import zipfile

from .evidence import read, save, sha, publish, restore

REPO=Path(__file__).resolve().parents[2]
DIRECTIONS=('R_to_P','P_to_R')
CASE_SECONDS=3600
RDF_SHA='4e287d49b70791ec30acee6d12babc26a8b7d15bd5730a3a4464e78e43eb74ce'


def environment():
    return dict(host=platform.node(),python=sys.version,platform=platform.platform(),
        affinity=sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else None,
        slurm={k:v for k,v in os.environ.items() if k.startswith('SLURM_')},
        dependencies={n:importlib.metadata.version(n) for n in ['numpy','rdkit','pynauty','sympy','z3-solver']})


def load_plan(folder):
    from mappa import AAMProblem,MolecularEndpoint,AAMSearchConfig,AAMSearchPlan
    raw=read(Path(folder)/'input.json');task=read(Path(folder)/'task.json')
    original=AAMProblem(*(MolecularEndpoint(**raw[s]) for s in ('reactant','product')),raw['name'])
    reverse=task['direction']=='P_to_R'
    directed=AAMProblem(original.product,original.reactant,original.name) if reverse else original
    return AAMSearchPlan(original,directed,AAMSearchConfig(**task['config']),reverse)


def prepare(args):
    from mappa import AAMSearchConfig
    from golden_evaluation import prepare as parse
    from prepare_golden_benchmark import audit_block
    assert not args.run.exists(),'Use a new run directory'
    with zipfile.ZipFile(args.dataset) as archive:rdf=archive.read('golden_dataset.rdf')
    assert hashlib.sha256(rdf).hexdigest()==RDF_SHA
    blocks=rdf.decode().split('$RFMT')[1:];assert len(blocks)==1851
    config=asdict(AAMSearchConfig(seed_count=10,branch_limit=100,iso_tolerance=1,random_seed=42))
    args.run.mkdir(parents=True)
    shutil.copy2(args.dataset,args.run/'dataset.zip')
    inputs=[];tasks=[]
    def endpoint(e):return dict(elements=e.elements,coordinates=e.coordinates.tolist(),wbo=e.wbo.tolist(),label=e.label,metadata=dict(e.metadata))
    for index,block in enumerate(blocks):
        folder=args.run/'inputs'/str(index)
        row=audit_block(block,index)
        save(folder/'source.json',row)
        try:
            assert row['status']=='parsed'
            problem,features,mapping=parse(row['mapped_reaction'])
            save(folder/'input.json',dict(name=f'golden_{index}',reactant=endpoint(problem.reactant),product=endpoint(problem.product)))
            save(folder/'reference.json',dict(features=features,mapping=sorted(mapping.items()),mapped_reaction=row['mapped_reaction']))
            inputs.append(dict(index=index,status='prepared',balanced=problem.balanced,
                reactant_atoms=problem.source_atom_count,product_atoms=problem.target_atom_count,
                input_sha256=sha(folder/'input.json'),reference_sha256=sha(folder/'reference.json')))
            tasks.extend(dict(index=index,direction=d,config=config) for d in DIRECTIONS)
        except Exception:
            inputs.append(dict(index=index,status='preparation_error',traceback=traceback.format_exc()))
    source={str(p.relative_to(REPO)):sha(p) for folder in ['src/mappa','native/src','bench/rebenchmark']
            for p in (REPO/folder).rglob('*') if p.is_file() and p.suffix in ('.py','.h','.cpp','.so')}
    for name in ['golden_evaluation.py','publication_analysis.py','publication_timing.py','golden_topk.py','view_golden_remaining.py']:
        source['bench/'+name]=sha(REPO/'bench'/name)
    save(args.run/'manifest.json',dict(schema='mappa.rebenchmark/v1',dataset='Golden',total=1851,
        rdf_sha256=RDF_SHA,dataset_sha256=sha(args.dataset),config=config,inputs=inputs,tasks=tasks,
        git_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
        source_sha256=source,environment=environment(),workers=8,concurrency=100,cpu_target=800,
        partition='cpu_short',memory_mib=22528,stage_seconds=dict(search=7200,analysis=5400),
        direction_policy='R_to_P and P_to_R, plus smaller-first and unconditional bidirectional union',
        seed_policy='One ten-seed search; retain independent contexts and derive every prefix 1..10',
        search_execution='reused_native',preparation='Original Golden adapter; explicit H; no composition exclusions',
        timing='Per-seed shared-cache CPU; full parent+child search CPU, IO and phase CPU separately. Parallel wall time is measured, never summed worker durations.',
        analysis='Original Golden heavy-relation evaluator; all seed/cut outcomes retained; exact query text and SAT/UNSAT/unknown recorded.',
        scope='Search and raw-family/representative metrics first; full decoded event metrics require a separately completed decoder stage on these archives. No search is repeated for downstream metrics.',
        created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
    print(json.dumps(dict(total=1851,prepared=len(tasks)//2,tasks=len(tasks))),flush=True)


def verify_sources(run):
    manifest=read(Path(run)/'manifest.json')
    for name,value in manifest['source_sha256'].items():
        assert sha(REPO/name)==value,name
    return manifest


def child(args):
    plan=load_plan(args.folder)
    if args.stage=='search':
        from .search import search
        search(plan.problem,plan.config,args.folder,args.workers)
    else:
        from .analysis import analyze
        analyze(args.folder,args.workers)


def guarded(command,folder,stage,seconds,memory_mib):
    import psutil
    start=time.perf_counter();usage=resource.getrusage(resource.RUSAGE_CHILDREN);peak=0;status=None
    with (folder/(stage+'.log')).open('w') as log:
        proc=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        while proc.poll() is None:
            try:
                process=psutil.Process(proc.pid)
                processes=[process,*process.children(recursive=True)]
                rss=sum(p.memory_info().rss for p in processes if p.is_running())
                peak=max(peak,rss)
            except (psutil.NoSuchProcess,psutil.AccessDenied):pass
            if time.perf_counter()-start>seconds:status='timeout'
            elif peak>memory_mib*2**20:status='memory_limit'
            if status:
                os.killpg(proc.pid,signal.SIGTERM)
                try:proc.wait(timeout=10)
                except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL)
                break
            time.sleep(.2)
        code=proc.wait()
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    return dict(status=status or ('complete' if code==0 else 'error'),exit_code=code,
        cpu_seconds=after.ru_utime+after.ru_stime-usage.ru_utime-usage.ru_stime,
        wall_seconds=time.perf_counter()-start,peak_tree_rss_bytes=peak,
        hard_seconds=seconds,memory_mib=memory_mib)


def worker(args,manifest=None):
    manifest=manifest or verify_sources(args.run);task=manifest['tasks'][args.task]
    progress=args.run/'progress'/f'{args.task:05d}.json'
    lock=args.run/'locks'/f'{args.task:05d}.lock';lock.parent.mkdir(exist_ok=True)
    with lock.open('a') as stream:
        fcntl.flock(stream,fcntl.LOCK_EX|fcntl.LOCK_NB)
        assert not progress.exists(),'An existing attempt requires explicit diagnosis'
        with tempfile.TemporaryDirectory(prefix='mappa-benchmark-',dir=os.environ.get('TMPDIR','/tmp')) as td:
            folder=Path(td);source=args.run/'inputs'/str(task['index'])
            for name in ('input.json','reference.json','source.json'):shutil.copy2(source/name,folder/name)
            save(folder/'task.json',task);save(folder/'environment.json',environment())
            started=time.perf_counter()
            state=dict(task=args.task,index=task['index'],direction=task['direction'],stages={},complete=False,
                started_utc=time.time(),case_seconds=CASE_SECONDS)
            for stage in ('search','analysis'):
                case_remaining=CASE_SECONDS-(time.perf_counter()-started)
                if case_remaining<=0:
                    state['stop_reason']='case_wall_limit'
                    break
                state['active_stage']=stage;save(progress,state)
                command=[sys.executable,'-m','rebenchmark.campaign','child','--folder',str(folder),
                    '--stage',stage,'--workers',str(manifest['workers'])]
                remaining=float(os.environ.get('MAPPA_ALLOCATION_DEADLINE',str(time.time()+14000)))-time.time()-180
                budget=min(manifest['stage_seconds'][stage],max(1,remaining),case_remaining)
                record=guarded(command,folder,stage,budget,manifest['memory_mib'])
                if record['status']=='timeout':
                    record['limit_scope']='case' if budget==case_remaining else 'stage_or_allocation'
                state['stages'][stage]=record
                save(folder/(stage+'-resources.json'),record)
                save(folder/'status.json',state)
                receipt=publish(folder,args.run/'archives'/f'{args.task:05d}.tar')
                state['archive']=receipt;save(progress,state)
                if record['status']!='complete':break
            state.update(active_stage=None,wall_seconds=time.perf_counter()-started,
                complete=all(state['stages'].get(s,{}).get('status')=='complete' for s in ('search','analysis')))
            save(progress,state)
            print(json.dumps(state),flush=True)


def drain(args):
    """Workers claim independent directions; public search distributes their cuts."""
    from types import SimpleNamespace
    manifest=verify_sources(args.run)
    sizes={row['index']:max(row.get('reactant_atoms',0),row.get('product_atoms',0)) for row in manifest['inputs']}
    order=sorted(range(len(manifest['tasks'])),key=lambda n:(-sizes[manifest['tasks'][n]['index']],n))
    claims=args.run/'claims';claims.mkdir(exist_ok=True)
    deadline=time.time()+3*3600+40*60
    os.environ['MAPPA_ALLOCATION_DEADLINE']=str(deadline)
    for task in order:
        if time.time()>deadline-300:break
        try:(claims/str(task)).mkdir()
        except FileExistsError:continue
        save(claims/str(task)/'owner.json',dict(job=os.environ.get('SLURM_JOB_ID'),pid=os.getpid(),utc=time.time()))
        worker(SimpleNamespace(run=args.run,task=task),manifest)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['prepare','worker','child','drain'])
    p.add_argument('--run',type=Path);p.add_argument('--dataset',type=Path)
    p.add_argument('--task',type=int);p.add_argument('--folder',type=Path)
    p.add_argument('--stage',choices=['search','analysis']);p.add_argument('--workers',type=int,default=8)
    a=p.parse_args();globals()[a.mode](a)


if __name__=='__main__':main()
