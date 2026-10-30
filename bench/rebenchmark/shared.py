"""Dataset planning for the portable execution layer; scoring remains in analysis."""
from pathlib import Path
import tarfile
import tempfile

from mappa.mappa_core.core.scheduler import Task
from mappa.mappa_core.core.scheduler.search import plan_search


def plan_direction(case, plan, reference, writer, *, seed_batch=10, case_seconds=3600., memory_bytes=2**30):
    tasks, cuts = plan_search(case, plan.problem, plan.config, writer,
        seed_batch=seed_batch, case_seconds=case_seconds, memory_bytes=memory_bytes)
    evidence = writer.put(dict(reference=reference, reversed=plan.reversed))
    finalization = {t.key: t for t in tasks if t.key in cuts}
    evaluations = []
    for number, key in enumerate(cuts):
        final = finalization[key]
        for seed in range(final.payload['seeds']):
            evaluation = f'{case}/cut/{number:05d}/evaluate/{seed:05d}'
            tasks.append(Task(evaluation, 'rebenchmark.shared:evaluate_seed',
                dict(input=final.payload['input'], reference=evidence, cut=key, seed=seed, number=number),
                dependencies=(key,), case=case, case_seconds=case_seconds, memory_bytes=memory_bytes))
            evaluations.append(evaluation)
    return tasks, cuts, tuple(evaluations)


def evaluate_seed(payload, workspace):
    from mappa import AAMProblem, AAMSearchPlan
    from mappa.mappa_core.core.scheduler.search import _session
    from .analysis import analyze_seed
    from .evidence import read
    session = _session(payload, workspace)
    problem, config, _target, _orbits = session.context
    reference = workspace.load(payload['reference'])
    original = AAMProblem(problem.product, problem.reactant, problem.name) if reference['reversed'] else problem
    plan = AAMSearchPlan(original, problem, config, reference['reversed'])
    finalized = workspace.dependency(payload['cut'])
    with tempfile.TemporaryDirectory(dir=workspace.scratch) as temporary:
        folder = Path(temporary)/'evidence'
        summary = analyze_seed(plan, reference['reference'], finalized['graph'], payload['seed'], folder,
            payload['number'], finalized_bundle_sha256=workspace.dependencies[payload['cut']]['blob'].sha256)
        archive = Path(temporary)/'evidence.tar.gz'
        with tarfile.open(archive, 'w:gz') as tar:
            tar.add(folder, arcname='evidence')
        proof = workspace.writer.put_file(archive)
        return dict(summary=summary, evaluation=read(folder/'evaluation.json'),
                    classes=read(folder/'classes.json'), evidence=proof)


def conversion_tasks(case, tasks, cuts, *, slack=5):
    """Attach the public converter to every seed of each finalized cut."""
    from dataclasses import replace
    finals = {t.key: t for t in tasks if t.key in cuts}
    for key in cuts:
        final = finals[key]
        for seed in range(final.payload['seeds']):
            payload = dict(input=final.payload['input'], cut=key, number=final.payload['number'],
                task_id=int(case), seed=seed, slack=slack, branch_limit=100, random_seed=42,
                decoder='compressed_event_interpretation')
            yield replace(final, key=f'{case}/cut/{payload["number"]:05d}/seed/{seed:02d}',
                operation='rebenchmark.shared:convert_seed', payload=payload, dependencies=(key,))


def convert_seed(payload, workspace):
    from mappa.mappa_core.core.scheduler.search import _session
    from .decoding import decode_graph
    import time
    if payload['decoder'] != 'compressed_event_interpretation':
        raise ValueError('Create a fresh plan for the current decoder')
    start, cpu = time.perf_counter(), time.process_time()
    session = _session(payload, workspace)
    finalized = workspace.dependency(payload['cut'])
    loading = dict(wall_seconds=time.perf_counter()-start, cpu_seconds=time.process_time()-cpu)
    value = decode_graph(session.context[0], finalized['graph'], payload)
    value['timing']['loading'] = loading
    blob = workspace.dependencies[payload['cut']]['blob']
    return dict(value, task_id=payload['task_id'], cut=payload['number'], input=payload['input'],
        finalized=blob, finalized_sha256=blob.sha256, search_rerun=True)


def _initialize_planner(output, scratch):
    from mappa.mappa_core.core.scheduler import BundleWriter
    global _planning_writer
    _planning_writer = BundleWriter(Path(output)/'inputs', scratch)


def _prepare_record(payload):
    """Plan one reaction, publishing only bounded bundles from local scratch."""
    import hashlib
    import json
    import numpy as np
    from collections import Counter
    from dataclasses import asdict
    from mappa import AAMSearchConfig, AAMProblem
    from mappa.mappa_core.utils.io.records import endpoint_record, problem_from_record
    plan, item, output = payload
    case = item['index']
    if plan['dataset'] == 'EnzymeMap':
        from golden_evaluation import prepare
        problem, input_features, _ = prepare(item['rxn'])
        problem = AAMProblem(problem.reactant, problem.product, f'enzymemap_{case}')
        reference_problem, features, mapping = prepare(item['ground_truth'])
        same_index = all(a.elements == b.elements and np.array_equal(a.wbo, b.wbo)
            for a, b in zip((problem.reactant, problem.product),
                            (reference_problem.reactant, reference_problem.product)))
        raw = dict(input=dict(name=problem.name,
            reactant=endpoint_record(problem.reactant, array_lists=True),
            product=endpoint_record(problem.product, array_lists=True)),
            reference=dict(features=features, mapping=sorted(mapping.items()),
                mapped_reaction=item['ground_truth'], input_reaction=item['rxn'],
                input_features=input_features, input_reference_index_agrees=same_index,
                original_id=item['original_id']))
    else:
        folder = Path(plan['input_manifest']).parent/'inputs'/str(case)
        raw = {}
        for name in ('input', 'reference'):
            data = (folder/f'{name}.json').read_bytes()
            if hashlib.sha256(data).hexdigest() != item[f'{name}_sha256']:
                raise ValueError(f'Input {case} {name} hash differs')
            raw[name] = json.loads(data)
        problem = problem_from_record(raw['input'])
    config = AAMSearchConfig(**plan['config'])
    writer = _planning_writer
    reference = asdict(writer.put(raw['reference']))
    rows = []
    for reverse in plan['directions']:
        direction = 2*case+int(reverse)
        directed = AAMProblem(problem.product,problem.reactant,problem.name) if reverse else problem
        tasks, cuts = plan_search(f'{direction:05d}', directed, config, writer,
            seed_batch=plan['seed_batch'], case_seconds=plan['case_seconds'], memory_bytes=plan['memory_bytes'])
        tasks += list(conversion_tasks(f'{direction:05d}',tasks,cuts,slack=plan['slack']))
        if plan.get('reference_checks', False):
            tasks += reference_tasks(f'{direction:05d}', tasks, cuts, reference, reverse)
        order = {'search_batch':0,'finalize_cut':1,'verify_seed':2,'convert_seed':3}
        tasks.sort(key=lambda t: (t.key.split('/')[2], order[t.operation.rsplit(':',1)[1]], t.key))
        rows.append(dict(direction=direction,case=case,reversed=reverse,cuts=len(cuts),
            tasks=len(tasks),counts=dict(Counter(t.operation.rsplit(':',1)[1] for t in tasks)),
            blob=asdict(writer.put(tasks)),reference=reference,
            input_sha256=hashlib.sha256(json.dumps(raw['input'],sort_keys=True).encode()).hexdigest(),
            reference_sha256=hashlib.sha256(json.dumps(raw['reference'],sort_keys=True).encode()).hexdigest()))
    return rows


def reference_tasks(case, tasks, cuts, reference, reverse):
    """Check the original ISO1 graph, independently of conversion success."""
    from dataclasses import replace
    finals = {t.key:t for t in tasks if t.key in cuts}
    return [replace(final,key=f'{case}/cut/{final.payload["number"]:05d}/verify/{seed:02d}',
        operation='rebenchmark.shared:verify_seed',dependencies=(key,),
        payload=dict(input=final.payload['input'],cut=key,seed=seed,reference=reference,reversed=reverse))
        for key,final in finals.items() for seed in range(final.payload['seeds'])]


def verify_seed(payload, workspace):
    from mappa import AAMProblem, AAMSearchPlan, AAMResult
    from mappa.mappa_core.domain import AAMSearchMetrics
    from mappa.mappa_core.core.scheduler.search import _session
    from golden_evaluation import evaluate_planned
    from .search import select_contexts
    from mappa.mappa_core.core.scheduler import Blob
    problem, config, _, _ = _session(payload, workspace).context
    reference = workspace.load(Blob(**payload['reference']))
    if not reference.get('input_reference_index_agrees', True):
        return dict(reference_recovery='unknown',reason='input/reference coordinate systems differ',
                    scope='original_iso1_aam',reference=payload['reference'])
    original = AAMProblem(problem.product,problem.reactant,problem.name) if payload['reversed'] else problem
    plan = AAMSearchPlan(original,problem,config,payload['reversed'])
    graph = select_contexts(workspace.dependency(payload['cut'])['graph'],[payload['seed']])
    result = AAMResult(problem,config,graph,AAMSearchMetrics.from_record({},0))
    return dict(evaluate_planned(result,plan,reference['features'],reference['mapping'],seconds=120),
                scope='original_iso1_aam',reference=payload['reference'])


def prepare_campaign(out):
    """Freeze a dataset's explicit cut/seed grid, with parallel input planning."""
    import csv
    import gzip
    import hashlib
    import json
    from collections import Counter
    from concurrent.futures import ProcessPoolExecutor
    from dataclasses import asdict
    from mappa import AAMSearchConfig
    from .decoder_campaign import save
    plan = json.loads((out/'plan.json').read_text())
    manifest_path = Path(plan['input_manifest'])
    if hashlib.sha256(manifest_path.read_bytes()).hexdigest() != plan['input_manifest_sha256']:
        raise ValueError('Input manifest changed')
    manifest = json.loads(manifest_path.read_text())
    plan.setdefault('dataset', 'Golden')
    plan.setdefault('directions', [False, True])
    if plan['dataset'] == 'EnzymeMap':
        dataset = Path(plan['dataset_path'])
        if hashlib.sha256(dataset.read_bytes()).hexdigest() != manifest['dataset_sha256']:
            raise ValueError('EnzymeMap dataset changed')
        with gzip.open(dataset,'rt') as stream:
            items = [dict(row,index=int(row['r_id'])) for row in csv.DictReader(stream)]
        if len(items) != manifest['cases'] or [r['index'] for r in items] != list(range(manifest['cases'])):
            raise ValueError('EnzymeMap row identities differ')
    elif plan['dataset'] == 'Golden':
        items = manifest['inputs']
        if len(items) != 1851 or len(manifest['tasks']) != 3702:
            raise ValueError('Golden input grid is incomplete')
    else:
        raise ValueError('Unknown benchmark dataset')
    config = AAMSearchConfig(**plan['config'])
    if (out/'index.json').exists():
        raise FileExistsError('The frozen plan already exists')
    if 'selected_cases' in plan:
        selected = set(plan['selected_cases'])
        items = [r for r in items if r['index'] in selected]
        if len(items) != len(selected):
            raise ValueError('Selected cases are missing')
    rows, totals = [], Counter()
    import os
    with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
        with ProcessPoolExecutor(max_workers=plan.get('preparation_workers',1),
                initializer=_initialize_planner,initargs=(str(out),scratch)) as pool:
            for done, prepared in enumerate(pool.map(_prepare_record,
                    ((plan,item,str(out)) for item in items), chunksize=1),1):
                rows.extend(prepared)
                for row in prepared:
                    totals.update(row['counts'])
                if done % 100 == 0:
                    save(out/'preparation.json',dict(cases=done,tasks=sum(totals.values())))
    save(out/'index.json',dict(schema='mappa.fresh_benchmark/v1',dataset=plan['dataset'],
        config=asdict(config),directions=rows,counts=dict(totals),tasks=sum(totals.values()),
        cases=len(items),input_manifest_sha256=plan['input_manifest_sha256']))


def prepare_queue(plan, out, scratch):
    """Reuse only verified pilot tasks from the identical frozen task grid."""
    from dataclasses import asdict
    import hashlib
    import json
    import pickle
    import sqlite3
    from mappa.mappa_core.core.scheduler import Blob, TaskQueue, read_blob
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint
    from .decoder_campaign import save
    path = Path(plan['index'])
    if hashlib.sha256(path.read_bytes()).hexdigest() != plan['index_sha256']:
        raise ValueError('Frozen task index changed')
    index = json.loads(path.read_text())
    queue = TaskQueue(Path(scratch)/'queue.sqlite', plan['identity'])
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
    reused, charged = [], {}
    for root in map(Path, plan['pilots']):
        prior_plan = json.loads((root/'plan.json').read_text())
        if prior_plan['commit'] not in plan['pilot_commits']:
            raise ValueError('Pilot commit is not explicitly accepted')
        for key in ('index_sha256', 'memory_bytes'):
            if prior_plan[key] != plan[key]:
                raise ValueError(f'Pilot {key} differs')
        if not json.loads((root/'progress.json').read_text()).get('finished'):
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


def prepare_memory_stage(plan, out, scratch):
    """Continue only authorized memory failures from a retired first stage."""
    import hashlib
    import json
    import sqlite3
    from mappa.mappa_core.core.scheduler import TaskQueue
    from .decoder_campaign import save
    source = Path(plan['previous_stage'])
    progress = json.loads((source/'progress.json').read_text())
    if not progress.get('finished') or any(progress['states'].get(k, 0) for k in ('pending', 'running')):
        raise ValueError('Previous memory stage still owns unfinished work')
    path = source/'queue.sqlite'
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
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
    import json
    import os
    import signal
    from mappa.mappa_core.core.scheduler import connect, run
    from .decoder_campaign import serve_queue
    out = args.output
    if args.mode == 'prepare':
        prepare_campaign(out)
        return
    plan = json.loads((out/'plan.json').read_text())
    if args.mode == 'worker':
        endpoint = json.loads((out/'endpoint.json').read_text())
        queue = connect((endpoint['host'],endpoint['port']),(out/'worker.key').read_bytes())
        claim_until = None
        if 'SLURM_JOB_ID' in os.environ:
            import datetime, re, subprocess, time
            job = subprocess.check_output(['scontrol','show','job',os.environ['SLURM_JOB_ID'],'-o'],text=True)
            end = datetime.datetime.fromisoformat(re.search(r'EndTime=(\S+)',job)[1]).timestamp()
            claim_until = min(end - plan['case_seconds'] - 180, plan['claim_until'])
        states = run(queue, workers=args.workers, identity=plan['identity'], scratch=os.environ['TMPDIR'],
            bundles=out/'bundles', memory_bytes=args.memory_bytes or plan['memory_bytes'], claim_until=claim_until)
        if claim_until is not None:
            from .decoder_campaign import save
            import socket
            folder = out/'retirements'; folder.mkdir(exist_ok=True)
            save(folder/f'{os.environ["SLURM_JOB_ID"]}-{socket.gethostname()}-{os.getpid()}.json',
                dict(states=states, claim_until=claim_until, utc=time.time(), allocation_end=end,
                     memory_bytes=args.memory_bytes or plan['memory_bytes'], workers=args.workers))
        return
    def terminate(signum, frame):
        raise SystemExit(128+signum)
    signal.signal(signal.SIGTERM, terminate)
    with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
        prepare = {'memory-stage': prepare_memory_stage, 'coordinator': prepare_queue}[args.mode]
        queue = prepare(plan,out,scratch)
        try:
            serve_queue(queue,out,dict(campaign='fresh matching and conversion'))
        finally:
            queue.close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('mode',choices=('prepare','coordinator','worker','memory-stage'))
    parser.add_argument('output',type=Path)
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--memory-bytes',type=int,
                        help='Per-worker guard for an explicitly recorded resource stage')
    execute(parser.parse_args())
