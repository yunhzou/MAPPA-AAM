"""Golden adapter for the portable execution layer; scoring remains in analysis."""
from pathlib import Path
import tarfile
import tempfile

from mappa.execution import Task
from mappa.execution.search import plan_search


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
    from mappa.execution.search import _session
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
