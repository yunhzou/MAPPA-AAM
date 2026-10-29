"""Local worker lifecycle for search and saved-cut finalization."""

from __future__ import annotations
from mappa.mappa_core.core.search.session import SearchSession
import gc, multiprocessing as mp, time
from pathlib import Path
from mappa.mappa_core.core.search.graph import AAMSearchGraph
from mappa.mappa_core.core.search.symmetry import finalize_graph_symmetry


def _initialize_search(problem, config, execution='reference'):
    global _SEARCH_SESSION, _SEARCH_CONTEXT, _SEARCH_REPAIR, _FINALIZATION_WORKSPACE
    _SEARCH_SESSION = SearchSession(problem, config, execution)
    # Profiling tools inspect these views; computation is owned by the session.
    _SEARCH_CONTEXT, _SEARCH_REPAIR = _SEARCH_SESSION.context, _SEARCH_SESSION.repair
    _FINALIZATION_WORKSPACE = None


def _search_cut(cut):
    return _SEARCH_SESSION.search_cut(cut)


def _search_cut_task(payload, *, in_process=False):
    """Checkpoint at the producer; a slow sibling cannot block persistence."""
    index, cut, checkpoint = payload
    graph, counts = _search_cut(cut)
    counts['checkpoint_seconds'] = 0.0
    if checkpoint is not None:
        from mappa.artifacts import write_raw_cut
        started = time.perf_counter()
        path = Path(checkpoint)
        write_raw_cut(graph,path)
        counts['checkpoint_seconds'] = time.perf_counter()-started
        # No graph is pickled through IPC or buffered in the parent.
        return index, graph if in_process else str(path), counts
    return index, graph, counts


def _initialize_finalization(problem,config,execution='reference'):
    global _FINALIZATION_WORKSPACE
    # Workers own acyclic archive data. Do not scan/copy a fork-inherited heap.
    gc.disable()
    _initialize_search(problem,config)
    if execution in ('reused_native', 'shared_policies'):
        from mappa.mappa_core.core.matching.symmetry import ConditionedSymmetryWorkspace
        _FINALIZATION_WORKSPACE = ConditionedSymmetryWorkspace(_SEARCH_CONTEXT[2], config.iso_tolerance)


def _restore_finalized_cut(payload):
    from mappa.artifacts import read_graph_checkpoint, write_graph_checkpoint, read_raw_cut
    index,data,checkpoint=payload
    if checkpoint is not None and Path(checkpoint).exists():
        return index,read_graph_checkpoint(checkpoint),{}
    graph=(data if isinstance(data,AAMSearchGraph) else
           read_raw_cut(data))
    _problem,config,target,_orbits=_SEARCH_CONTEXT
    graph,counts=finalize_graph_symmetry(graph,target,iso_tolerance=config.iso_tolerance,
                                       workspace=_FINALIZATION_WORKSPACE)
    if checkpoint is not None:write_graph_checkpoint(graph,checkpoint)
    return index,graph,counts


_SEARCH_CONTEXT = None


_SEARCH_REPAIR = None


_FINALIZATION_WORKSPACE = None


_SEARCH_SESSION = None


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
