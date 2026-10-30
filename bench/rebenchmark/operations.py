"""Scientific task adapters; the public MAPPA APIs own the algorithms."""
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
            tasks.append(Task(evaluation, 'rebenchmark.operations:evaluate_seed',
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
                operation='rebenchmark.operations:convert_seed', payload=payload, dependencies=(key,))


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


def reference_tasks(case, tasks, cuts, reference, reverse):
    """Check the original ISO1 graph, independently of conversion success."""
    from dataclasses import replace
    finals = {t.key:t for t in tasks if t.key in cuts}
    return [replace(final,key=f'{case}/cut/{final.payload["number"]:05d}/verify/{seed:02d}',
        operation='rebenchmark.operations:verify_seed',dependencies=(key,),
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
