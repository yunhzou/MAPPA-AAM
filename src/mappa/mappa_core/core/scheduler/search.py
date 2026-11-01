"""Deterministic cut/seed tasks calling the same prepared search as public APIs."""
from __future__ import annotations
import gc
import multiprocessing as mp
from pathlib import Path
import time

from mappa.mappa_core.core.scheduler.queue import Task
from mappa.mappa_core.core.search.session import SearchSession
from mappa.mappa_core.core.search.planning import cut_seed
from mappa.mappa_core.core.search.planning import _generate_seed_orders
from mappa.mappa_core.core.search.planning import cut_sweep_items
from mappa.mappa_core.core.matching.graphs import build_graph
from mappa.mappa_core.core.search.graph import AAMSearchGraph
from mappa.mappa_core.core.search.branching.symmetry import finalize_graph_symmetry
from .results import SearchInput, SeedBatchResult, CutResult

def _session(payload, workspace):
    reference = payload['input']
    def prepare():
        record = workspace.load(reference)
        return SearchSession(record.problem, record.config, record.execution)
    return workspace.prepared(reference, prepare)


def plan_search(case, problem, config, writer, *, seed_batch=1, case_seconds=3600., memory_bytes=2**30, execution='reused_native', shared_cache=False):
    """Plan a direction once; combine several plans in the same TaskQueue.

    Inputs occur once in the store. Every task names its original seed indices
    and orders; the scientific seed_count remains the full requested count.
    The returned cut keys are in canonical order regardless of task completion.
    """
    if seed_batch < 1:
        raise ValueError('seed_batch must be positive')
    if execution not in ('reference', 'reused_native'):
        raise ValueError('Seed scheduling requires reference or reused_native execution')
    reference = writer.put(SearchInput(problem, config, execution))
    cuts = cut_sweep_items(problem.reactant.wbo, config.cut_floor) if config.sweep_cuts else [()]
    tasks, final_keys = [], []
    original = build_graph(problem.reactant.elements, problem.reactant.wbo, bond_cut=config.graph_floor)
    for number, cut in enumerate(cuts):
        source = original.copy()
        source.remove_edges_from(cut)
        orders = tuple(map(tuple, _generate_seed_orders(source, config.seed_count,
            rng_seed=cut_seed(cut, config.random_seed), seed_selection=config.seed_selection)))
        if not orders:
            raise ValueError('Search plan has no seed orders')
        if len(orders) <= seed_batch:
            key = f'{case}/cut/{number:05d}/finalize'
            tasks.append(Task(key, 'mappa.mappa_core.core.scheduler.search:search_complete_cut',
                dict(input=reference, cut=cut, number=number, seeds=len(orders),
                     indices=tuple(range(len(orders))), orders=orders),
                case=case, case_seconds=case_seconds, memory_bytes=memory_bytes))
            final_keys.append(key)
            continue
        dependencies = []
        for first in range(0, len(orders), seed_batch):
            indices = tuple(range(first, min(first+seed_batch, len(orders))))
            key = f'{case}/cut/{number:05d}/seeds/{first:05d}'
            tasks.append(Task(key, 'mappa.mappa_core.core.scheduler.search:search_batch',
                dict(input=reference, cut=cut, indices=indices, orders=orders[first:first+seed_batch]),
                case=case, case_seconds=case_seconds, memory_bytes=memory_bytes))
            dependencies.append(key)
        key = f'{case}/cut/{number:05d}/finalize'
        tasks.append(Task(key, 'mappa.mappa_core.core.scheduler.search:finalize_cut',
            dict(input=reference, cut=cut, number=number, seeds=len(orders), batches=tuple(dependencies)),
            dependencies=tuple(dependencies), case=case, case_seconds=case_seconds, memory_bytes=memory_bytes))
        final_keys.append(key)
    if shared_cache and execution == 'reused_native':
        from mappa.execution.identity import code_identity
        from .cache import cache_identity
        identity = cache_identity(problem, config, code_identity())
        for task in tasks:
            if task.operation.endswith((':search_batch', ':search_complete_cut')):
                task.payload['cache'] = identity
    return tasks, tuple(final_keys)


def search_batch(payload, workspace):
    session = _session(payload, workspace)
    if workspace.cache.descriptor is not None:
        from mappa.execution.identity import code_identity
        from .cache import cache_identity
        problem, config, *_ = session.context
        identity = workspace.resource('cache_identity', workspace.prepared_key,
            lambda: cache_identity(problem, config, code_identity()))
        descriptor = workspace.cache.descriptor
        if identity != (descriptor['key'], descriptor['direction']):
            raise ValueError('Task cache input, settings or code identity differs')
    workspace.cache.use(session.repair)
    graph, counts = session.search_cut(payload['cut'], seed_indices=payload['indices'])
    if tuple(context.seed_order for context in graph.contexts) != payload['orders']:
        raise ValueError('Executed seed orders differ from the planned orders')
    return SeedBatchResult(graph, payload['indices'], counts)


def finalize_cut(payload, workspace):
    graphs, indices, metrics = [], [], []
    for key in payload['batches']:
        result = workspace.dependency(key)
        graphs.append(result.graph)
        indices.extend(result.indices)
        metrics.append(result.counts)
    if indices != list(range(payload['seeds'])):
        raise ValueError('Missing, repeated or reordered seed batches')
    raw = AAMSearchGraph.combine(graphs)
    del graphs
    if any(tuple(context.cuts) != tuple(payload['cut']) for context in raw.contexts):
        raise ValueError('Seed batch belongs to a different cut')
    return _finalize(raw, metrics, payload, workspace)


def search_complete_cut(payload, workspace):
    """No intermediate transport when one task already owns every seed of a cut."""
    batch = search_batch(payload, workspace)
    return _finalize(batch.graph, [batch.counts], payload, workspace)


def _finalize(raw, metrics, payload, workspace):
    from mappa.mappa_core.core.matching.symmetry.workspace import ConditionedSymmetryWorkspace
    session = _session(payload, workspace)
    _problem, config, target, _orbits = session.context
    symmetry = workspace.resource('finalization', config.iso_tolerance,
        lambda: ConditionedSymmetryWorkspace(target, config.iso_tolerance))
    graph, statistics = finalize_graph_symmetry(raw, target, iso_tolerance=config.iso_tolerance,
        workspace=symmetry)
    return CutResult(graph, payload['number'], payload['seeds'], tuple(metrics), statistics)


def iter_cuts(queue, keys, scratch):
    """Read canonical cuts one at a time; refuse incomplete campaign results."""
    from mappa.mappa_core.core.scheduler.bundles import read_blob
    rows = queue.results(keys)
    for key in keys:
        row = rows[key]
        if row['state'] != 'complete':
            raise RuntimeError(f'Incomplete cut {key}: {row["state"]}')
        yield read_blob(row['result']['blob'], scratch)


def _initialize_search(problem, config, execution='reference'):
    global _SEARCH_SESSION, _SEARCH_CONTEXT, _SEARCH_REPAIR, _FINALIZATION_WORKSPACE, _CHECKPOINT_STORE
    _SEARCH_SESSION = SearchSession(problem, config, execution)
    # Profiling tools inspect these views; computation is owned by the session.
    _SEARCH_CONTEXT, _SEARCH_REPAIR = _SEARCH_SESSION.context, _SEARCH_SESSION.repair
    _FINALIZATION_WORKSPACE = None
    _CHECKPOINT_STORE = None


def _search_cut(cut):
    return _SEARCH_SESSION.search_cut(cut)


def _search_cut_task(payload, *, in_process=False):
    """Checkpoint at the producer; a slow sibling cannot block persistence."""
    index, cut, checkpoint = payload
    graph, counts = _search_cut(cut)
    counts['checkpoint_seconds'] = 0.0
    if checkpoint is not None:
        started = time.perf_counter()
        reference = _checkpoint_store(checkpoint).put(graph)
        counts['checkpoint'] = reference
        counts['checkpoint_seconds'] = time.perf_counter()-started
        # No graph is pickled through IPC or buffered in the parent.
        return index, graph if in_process else reference, counts
    return index, graph, counts


def _initialize_finalization(problem,config,execution='reference'):
    global _FINALIZATION_WORKSPACE
    # Workers own acyclic archive data. Do not scan/copy a fork-inherited heap.
    gc.disable()
    _initialize_search(problem,config)
    if execution in ('reused_native', 'shared_policies'):
        from mappa.mappa_core.core.matching.symmetry.workspace import ConditionedSymmetryWorkspace
        _FINALIZATION_WORKSPACE = ConditionedSymmetryWorkspace(_SEARCH_CONTEXT[2], config.iso_tolerance)


def _restore_finalized_cut(payload):
    from .bundles import Blob, read_blob
    index,data,checkpoint=payload
    if isinstance(checkpoint, Blob):
        return index,read_blob(checkpoint, None),{}
    graph=(data if isinstance(data,AAMSearchGraph) else
           read_blob(data, None))
    _problem,config,target,_orbits=_SEARCH_CONTEXT
    graph,counts=finalize_graph_symmetry(graph,target,iso_tolerance=config.iso_tolerance,
                                       workspace=_FINALIZATION_WORKSPACE)
    if checkpoint is not None:
        counts['checkpoint'] = _checkpoint_store(checkpoint).put(graph)
    return index,graph,counts


_SEARCH_CONTEXT = None


_SEARCH_REPAIR = None


_FINALIZATION_WORKSPACE = None


_SEARCH_SESSION = None

_CHECKPOINT_STORE = None


def _checkpoint_store(directory):
    from .storage import BundleStore
    global _CHECKPOINT_STORE
    if _CHECKPOINT_STORE is None:
        _CHECKPOINT_STORE = BundleStore(directory)
    if _CHECKPOINT_STORE.directory != Path(directory).resolve():
        raise ValueError('Worker checkpoint directory changed')
    return _CHECKPOINT_STORE


def iter_cut_tasks(problem, config, execution, tasks, workers, *, in_process=False):
    """Run cut tasks with the existing worker initialization and checkpoint policy."""
    if workers == 1:
        _initialize_search(problem, config, execution)
        yield from (_search_cut_task(task, in_process=in_process) for task in tasks)
    else:
        with mp.get_context('fork').Pool(workers, initializer=_initialize_search,
                initargs=(problem, config, execution)) as pool:
            yield from pool.imap_unordered(_search_cut_task, tasks,
                                           chunksize=config.task_chunksize)


def iter_finalized_tasks(problem, config, execution, payloads, workers):
    """Restore and finalize saved cuts using the existing worker lifecycle."""
    with mp.get_context('fork').Pool(workers, initializer=_initialize_finalization,
            initargs=(problem, config, execution)) as pool:
        yield from pool.imap_unordered(_restore_finalized_cut, payloads, chunksize=1)

