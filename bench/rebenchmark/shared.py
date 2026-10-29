"""Golden adapter for the portable execution layer; scoring remains in analysis."""
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
                task_id=int(case), seed=seed, slack=slack, branch_limit=100, random_seed=42)
            yield replace(final, key=f'{case}/cut/{payload["number"]:05d}/seed/{seed:02d}',
                operation='rebenchmark.shared:convert_seed', payload=payload, dependencies=(key,))


def convert_seed(payload, workspace):
    from mappa.mappa_core.core.scheduler.search import _session
    from .decoding import decode_graph
    import time
    start, cpu = time.perf_counter(), time.process_time()
    session = _session(payload, workspace)
    finalized = workspace.dependency(payload['cut'])
    loading = dict(wall_seconds=time.perf_counter()-start, cpu_seconds=time.process_time()-cpu)
    value = decode_graph(session.context[0], finalized['graph'], payload)
    value['timing']['loading'] = loading
    blob = workspace.dependencies[payload['cut']]['blob']
    return dict(value, task_id=payload['task_id'], cut=payload['number'], input=payload['input'],
        finalized=blob, finalized_sha256=blob.sha256, search_rerun=True)


def prepare_campaign(out):
    """Freeze the full task grid in bundles once, before pilot or full execution."""
    from collections import Counter
    from dataclasses import asdict
    import hashlib
    import json
    import os
    from mappa import AAMSearchConfig, AAMSearchPlan, AAMProblem
    from mappa.mappa_core.utils.io.records import problem_from_record
    from mappa.mappa_core.core.scheduler import BundleWriter
    from .decoder_campaign import save
    plan = json.loads((out/'plan.json').read_text())
    manifest_path = Path(plan['input_manifest'])
    if hashlib.sha256(manifest_path.read_bytes()).hexdigest() != plan['input_manifest_sha256']:
        raise ValueError('Input manifest changed')
    manifest = json.loads(manifest_path.read_text())
    if len(manifest['inputs']) != 1851 or len(manifest['tasks']) != 3702:
        raise ValueError('Golden input grid is incomplete')
    config = AAMSearchConfig(**plan['config'])
    if (config.iso_tolerance, config.branch_limit, config.seed_count, config.random_seed, config.sweep_cuts) != (1,100,10,42,True):
        raise ValueError('Unexpected Golden search settings')
    if (out/'index.json').exists():
        raise FileExistsError('The frozen plan already exists')
    with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
        writer = BundleWriter(out/'inputs', scratch)
        rows, totals = [], Counter()
        for item in manifest['inputs']:
            case = item['index']
            folder = manifest_path.parent/'inputs'/str(case)
            raw = {}
            for name in ('input', 'reference'):
                data = (folder/f'{name}.json').read_bytes()
                if hashlib.sha256(data).hexdigest() != item[f'{name}_sha256']:
                    raise ValueError(f'Golden {case} {name} hash differs')
                raw[name] = json.loads(data)
            original = problem_from_record(raw['input'])
            for reverse in (False, True):
                direction = 2*case+int(reverse)
                problem = AAMProblem(original.product,original.reactant,original.name) if reverse else original
                tasks, cuts = plan_search(f'{direction:05d}', problem, config, writer,
                    seed_batch=plan['seed_batch'], case_seconds=plan['case_seconds'], memory_bytes=plan['memory_bytes'])
                converted = list(conversion_tasks(f'{direction:05d}', tasks, cuts, slack=plan['slack']))
                # Keep dependent work close to its producer, so queued later cuts
                # do not consume the direction's deadline before conversion starts.
                tasks = sorted(tasks+converted, key=lambda t: (t.key.split('/')[2],
                    {'search_batch':0,'finalize_cut':1,'convert_seed':2}[t.operation.rsplit(':',1)[1]], t.key))
                counts = Counter(t.operation.rsplit(':',1)[1] for t in tasks)
                totals.update(counts)
                rows.append(dict(direction=direction, case=case, reversed=reverse, cuts=len(cuts),
                    tasks=len(tasks), counts=dict(counts), blob=asdict(writer.put(tasks)),
                    reference=asdict(writer.put(raw['reference'])),
                    input_sha256=item['input_sha256'], reference_sha256=item['reference_sha256']))
            if case % 100 == 0:
                save(out/'preparation.json', dict(cases=case+1, tasks=sum(totals.values())))
        save(out/'index.json', dict(schema='mappa.fresh_golden/v1', config=asdict(config),
            directions=rows, counts=dict(totals), tasks=sum(totals.values()),
            input_manifest_sha256=plan['input_manifest_sha256']))


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
            for name, seconds, deadline in prior.execute('SELECT name,seconds,deadline FROM cases'):
                end = prior.execute('SELECT max(a.finished) FROM attempts a JOIN tasks t ON a.task=t.key WHERE t.case_name=?',(name,)).fetchone()[0]
                charged[name] = charged.get(name,0.) + max(0.,end-(deadline-seconds))
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
            save(folder/f'{os.environ["SLURM_JOB_ID"]}-{socket.gethostname()}.json',
                dict(states=states, claim_until=claim_until, utc=time.time(), allocation_end=end,
                     memory_bytes=args.memory_bytes or plan['memory_bytes'], workers=args.workers))
        return
    def terminate(signum, frame):
        raise SystemExit(128+signum)
    signal.signal(signal.SIGTERM, terminate)
    with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
        queue = prepare_queue(plan,out,scratch)
        try:
            serve_queue(queue,out,dict(campaign='fresh Golden matching and conversion'))
        finally:
            queue.close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('mode',choices=('prepare','coordinator','worker'))
    parser.add_argument('output',type=Path)
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--memory-bytes',type=int,
                        help='Per-worker guard for an explicitly recorded resource stage')
    execute(parser.parse_args())
