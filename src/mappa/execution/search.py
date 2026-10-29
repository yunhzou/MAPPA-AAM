"""Deterministic cut/seed tasks calling the same prepared search as public APIs."""
from dataclasses import asdict

from .queue import Task
from mappa.mappa_core.core.search.session import SearchSession
from mappa.mappa_core.core.search.planning import cut_seed
from mappa.mappa_core.core.search.planning import _generate_seed_orders
from mappa.mappa_core.core.search.planning import cut_sweep_items
from ..domain import AAMProblem, AAMSearchConfig, MolecularEndpoint
from ..frag import build_graph
from mappa.mappa_core.core.search.graph import AAMSearchGraph


def _endpoint(endpoint):
    return dict(elements=endpoint.elements, coordinates=endpoint.coordinates, wbo=endpoint.wbo,
                label=endpoint.label, energy=endpoint.energy, metadata=dict(endpoint.metadata))


def _restore_problem(record):
    return AAMProblem(MolecularEndpoint(**record['reactant']), MolecularEndpoint(**record['product']), record['name'])


def _session(payload, workspace):
    reference = payload['input']
    def prepare():
        record = workspace.load(reference)
        problem = _restore_problem(record)
        return SearchSession(problem, AAMSearchConfig(**record['config']), record['execution'])
    return workspace.prepared(reference, prepare)


def plan_search(case, problem, config, writer, *, seed_batch=10, case_seconds=3600., memory_bytes=2**30):
    """Plan a direction once; combine several plans in the same TaskQueue.

    Inputs occur once in a bundle. Every task names its original seed indices
    and orders; the scientific seed_count remains the full requested count.
    The returned cut keys are in canonical order regardless of task completion.
    """
    if seed_batch < 1:
        raise ValueError('seed_batch must be positive')
    reference = writer.put(dict(reactant=_endpoint(problem.reactant), product=_endpoint(problem.product),
        name=problem.name, config=asdict(config), execution='reused_native'))
    cuts = cut_sweep_items(problem.reactant.wbo, config.cut_floor) if config.sweep_cuts else [()]
    tasks, final_keys = [], []
    for number, cut in enumerate(cuts):
        source = build_graph(problem.reactant.elements, problem.reactant.wbo, bond_cut=config.graph_floor)
        source.remove_edges_from(cut)
        orders = tuple(map(tuple, _generate_seed_orders(source, config.seed_count,
            rng_seed=cut_seed(cut, config.random_seed), seed_selection=config.seed_selection)))
        if not orders:
            raise ValueError('Search plan has no seed orders')
        dependencies = []
        for first in range(0, len(orders), seed_batch):
            indices = tuple(range(first, min(first+seed_batch, len(orders))))
            key = f'{case}/cut/{number:05d}/seeds/{first:05d}'
            tasks.append(Task(key, 'mappa.execution.search:search_batch',
                dict(input=reference, cut=cut, indices=indices, orders=orders[first:first+seed_batch]),
                case=case, case_seconds=case_seconds, memory_bytes=memory_bytes))
            dependencies.append(key)
        key = f'{case}/cut/{number:05d}/finalize'
        tasks.append(Task(key, 'mappa.execution.search:finalize_cut',
            dict(input=reference, cut=cut, number=number, seeds=len(orders), batches=tuple(dependencies)),
            dependencies=tuple(dependencies), case=case, case_seconds=case_seconds, memory_bytes=memory_bytes))
        final_keys.append(key)
    return tasks, tuple(final_keys)


def search_batch(payload, workspace):
    session = _session(payload, workspace)
    graph, counts = session.search_cut(payload['cut'], seed_indices=payload['indices'])
    if tuple(context.seed_order for context in graph.contexts) != payload['orders']:
        raise ValueError('Executed seed orders differ from the planned orders')
    return dict(graph=graph, indices=payload['indices'], counts=counts)


def finalize_cut(payload, workspace):
    from mappa.mappa_core.core.matching.symmetry import ConditionedSymmetryWorkspace
    from mappa.mappa_core.core.search.symmetry import finalize_graph_symmetry
    session = _session(payload, workspace)
    graphs, indices, metrics = [], [], []
    for key in payload['batches']:
        result = workspace.dependency(key)
        graphs.append(result['graph'])
        indices.extend(result['indices'])
        metrics.append(result['counts'])
    if indices != list(range(payload['seeds'])):
        raise ValueError('Missing, repeated or reordered seed batches')
    raw = AAMSearchGraph.combine(graphs)
    del graphs
    if any(tuple(context.cuts) != tuple(payload['cut']) for context in raw.contexts):
        raise ValueError('Seed batch belongs to a different cut')
    _problem, config, target, _orbits = session.context
    graph, statistics = finalize_graph_symmetry(raw, target, iso_tolerance=config.iso_tolerance,
        workspace=ConditionedSymmetryWorkspace(target, config.iso_tolerance))
    return dict(graph=graph, cut=payload['number'], seeds=payload['seeds'],
                search_counts=metrics, finalization_counts=statistics)


def iter_cuts(queue, keys, scratch):
    """Read canonical cuts one at a time; refuse incomplete campaign results."""
    from .bundles import read_blob
    rows = queue.results(keys)
    for key in keys:
        row = rows[key]
        if row['state'] != 'complete':
            raise RuntimeError(f'Incomplete cut {key}: {row["state"]}')
        yield read_blob(row['result']['blob'], scratch)
