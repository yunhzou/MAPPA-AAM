"""Compare reaction snapshot exchange on identical full sweeps in one Slurm pool."""
import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import time

from mappa import AAMSearchConfig
from mappa.execution import BundleStore, StageResults, WorkPool
from mappa.execution.identity import code_identity
from mappa.mappa_core.core.scheduler import Task
from mappa.mappa_core.core.scheduler.bundles import Blob
from enzymemap import prepared_problem


def wait(stage, label):
    next_report = 0
    while True:
        rows = stage._rows()
        if all(row['state'] not in ('pending','running') for row in rows.values()):
            return stage.wait_and_validate()
        stage.pool.check_workers()
        if time.monotonic() >= next_report:
            states = {}
            for row in rows.values():
                states[row['state']] = states.get(row['state'],0)+1
            print(json.dumps(dict(stage=label, tasks=states)), flush=True)
            next_report = time.monotonic()+30
        time.sleep(.2)


def metrics(pool, stage, elapsed):
    rows = stage._rows()
    values = [row['result'] for row in rows.values()]
    cache = [row.get('cache', {}) for row in values]
    summary = dict(wall_seconds=elapsed, tasks=len(values),
        compute_cpu_seconds=sum(row.get('cpu_seconds',0) for row in values),
        total_cpu_seconds=sum(row.get('total_cpu_seconds',0) for row in values),
        publication_seconds=sum(row.get('publish_seconds',0) for row in values),
        peak_worker_rss_bytes=max(row.get('peak_tree_rss_bytes',0) for row in values),
        workers_used=len({row['worker'] for row in values}))
    summary['cache'] = {key: sum(row.get(key,0) for row in cache) for key in
        ('load_seconds','save_seconds','loaded_bytes','saved_bytes','baseline_loads',
         'fragment_calls','fragment_hits','dependency_checks','published')}
    summary['cache']['errors'] = [error for row in cache for error in row.get('errors',[])]
    keys = set(rows)
    attempts = [row for row in pool.queue.db.execute(
        'SELECT task,started,finished FROM attempts WHERE state="complete"') if row[0] in keys]
    active = peak = 0
    for _, change in sorted([(start,1) for _,start,_ in attempts]+[(end,-1) for _,_,end in attempts]):
        active += change; peak = max(peak,active)
    summary['peak_overlapping_tasks'] = peak
    return summary


def run(indices, output, *, count=10, workers=8, local=False):
    output = Path(output); output.mkdir(parents=True, exist_ok=False)
    manifests = [json.loads(Path(path).read_text()) for path in indices]
    if any(item['config'] != manifests[0]['config'] for item in manifests):
        raise ValueError('Pilot input configurations differ')
    all_rows = [row for item in manifests for row in item['directions']]
    cases = sorted({row['case'] for row in all_rows})[:count]
    rows = sorted((row for row in all_rows if row['case'] in cases),key=lambda row:(row['case'],row['reversed']))
    if len(rows) != 2*count or len({(r['case'],r['reversed']) for r in rows}) != len(rows):
        raise ValueError('Pilot requires both directions for each distinct reaction')
    named = {f"{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}": row for row in rows}
    problems = {name: prepared_problem(row) for name,row in named.items()}
    config = AAMSearchConfig(**{**manifests[0]['config'], 'growth_policy':'existence'})
    report = dict(job=os.environ.get('SLURM_JOB_ID'), code=code_identity(), config=asdict(config),
        reaction_ids=cases, directions=len(rows), inputs=rows,
        atoms={name:[len(p.reactant.elements),len(p.product.elements)] for name,p in problems.items()},
        stages={})
    def save():
        (output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    save()
    store = BundleStore(output/'bundles')
    factory = WorkPool.local if local else WorkPool.slurm
    with factory(workers=workers, store=store, memory_bytes=4*2**30, case_seconds=3600.) as pool:
        report['workers'] = pool.queue.ready_workers()
        report['startup_seconds'] = pool.startup_elapsed_seconds
        matched = {}
        for sharing in (False,True):
            label = 'shared_cache' if sharing else 'worker_cache_only'
            started = time.perf_counter()
            matched[label] = wait(pool.search(problems, config, seed_batch=config.seed_count,
                shared_cache=sharing), label)
            report['stages'][label] = metrics(pool,matched[label],time.perf_counter()-started)
            matched[label].save_manifest(output/f'{label}.json')
            save()
            print(json.dumps(dict(stage=label, summary=report['stages'][label])),flush=True)
        equal = 0
        for name in problems:
            previous, current = (matched[label].load(name) for label in ('worker_cache_only','shared_cache'))
            if previous.scope != current.scope:
                raise AssertionError(f'Cut/seed scope differs: {name}')
            for a,b in zip(previous.cuts,current.cuts,strict=True):
                if store.load(a).graph != store.load(b).graph:
                    raise AssertionError(f'Full finalized search graph differs: {name}')
                equal += 1
        report['exact_graph_comparisons'] = equal
        save()
        started = time.perf_counter()
        catalogues = wait(pool.deduplicate(matched['shared_cache']),'deduplication')
        report['stages']['deduplication'] = metrics(pool,catalogues,time.perf_counter()-started)
        catalogues.save_manifest(output/'catalogues.json'); save()
        started = time.perf_counter(); decoded = wait(pool.decode(catalogues),'decoding')
        report['stages']['decoding'] = metrics(pool,decoded,time.perf_counter()-started)
        decoded.save_manifest(output/'decoded.json'); save()
        tasks, keys = [], {}
        for name,row in named.items():
            key = f'verify/{name}'
            tasks.append(Task(key,'enzymemap:verify',dict(catalogue=catalogues.references[name],
                reference=Blob(**row['reference']),reversed=row['reversed']),
                case=key,memory_bytes=4*2**30))
            keys[name] = (key,)
        started=time.perf_counter();pool.queue.add(tasks)
        verified=wait(StageResults('verification',store,pool=pool,tasks=keys),'verification')
        report['verification']={name:verified.load(name) for name in problems}
        report['stages']['verification']=metrics(pool,verified,time.perf_counter()-started)
        verified.save_manifest(output/'verification.json')
        report['cache_files']={path.name:path.stat().st_size for path in (store.directory/'cache').iterdir()}
        report['complete']=True;save()
    print(json.dumps(dict(complete=True,output=str(output))),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index',action='append',required=True)
    parser.add_argument('--output',required=True)
    parser.add_argument('--count',type=int,default=10)
    parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--local',action='store_true')
    args=parser.parse_args()
    run(args.index,args.output,count=args.count,workers=args.workers,local=args.local)
