"""Public AAM search and checkpoint orchestration."""

from __future__ import annotations
from mappa.mappa_core.core.search.planning import checkpoint_manifest, cut_sweep_items
from mappa.execution.local_search import iter_cut_tasks, iter_finalized_tasks
import json, gc, time
from dataclasses import dataclass
from pathlib import Path
from mappa.domain import AAMProblem, AAMResult, AAMSearchConfig, AAMSearchMetrics
from mappa.frag import build_graph
from mappa.mappa_core.core.search.graph import AAMSearchGraph
from mappa.mappa_core.core.search.symmetry import finalize_graph_symmetry


def search_aam(problem: AAMProblem, config: AAMSearchConfig | None = None,
               *, workers: int = 1, intermediate_dir=None, resume=False,
               archive_format='json', execution='reference') -> AAMResult:
    """Return raw matching histories; repair/grouping/ranking are separate calls.

    When supplied, intermediate_dir receives each completed cut graph before
    group finalization and a final reusable graph record. No path or bijection
    enumeration is needed to collect worker results.

    archive_format selects both raw-cut and final persistence: JSON interchange
    records, or trusted compact checkpoints that preserve shared storage.
    Completed cuts are persisted by workers and collected out of order, but
    the returned graph always retains canonical cut/context ordering.

    resume=True reuses cuts only after verifying their input/configuration
    manifest. Invocation timings exclude the earlier checkpointed work;
    resulting graph and cut counts still describe the whole search.

    execution='reused_native' explicitly selects the native sequential scheduler,
    dependency-validated fragment reuse and conditioned symmetry workspace. It
    preserves the seed/cut policy and output graph; it requires the native engine
    and never switches to a different search implementation on failure.

    execution='shared_policies' preserves those same seed orders and frontier
    admissions while sharing equal conditional decisions and fragment DAG
    records across policies within each cut. Independent cuts use the same
    producer-side checkpoint/worker pipeline; completed runtime caches are
    released at each cut boundary. The native fragment-growth kernel is unchanged.
    """
    if not isinstance(problem, AAMProblem):
        raise TypeError('search_aam requires an AAMProblem')
    if execution not in ('reference', 'reused_native', 'shared_policies'):
        raise ValueError('unknown AAM execution backend')
    if execution in ('reused_native', 'shared_policies'):
        from mappa.mappa_core.core.matching.growth.native import available
        if not available():
            raise ValueError('reused_native execution requires the built native engine')
    from mappa.artifacts import raw_cut_paths, read_raw_cut
    config = config or AAMSearchConfig()
    requested_workers=max(1,int(workers))
    if archive_format not in ('json','checkpoint'):
        raise ValueError('archive_format must be json or checkpoint')
    started = time.perf_counter()
    cuts = cut_sweep_items(problem.reactant.wbo, config.cut_floor) if config.sweep_cuts else [()]
    directory = None if intermediate_dir is None else Path(intermediate_dir)
    if resume and directory is None:
        raise ValueError('Resuming requires an intermediate directory')
    if directory is not None:
        directory.mkdir(parents=True, exist_ok=True)
        manifest_path = directory / 'manifest.json'
        identity = checkpoint_manifest(problem,config)
        if execution == 'shared_policies':
            identity['execution'] = execution
        if resume:
            if json.loads(manifest_path.read_text()) != identity:
                raise ValueError('Checkpoint input or configuration differs from this search')
        else:
            if manifest_path.exists() or raw_cut_paths(directory):
                raise ValueError('Existing checkpoints require explicit resume=True')
            manifest_path.write_text(json.dumps(identity)+'\n')
    metrics = {'max_live_branches': 0, 'max_growth_candidates': 0}

    checkpoint_seconds = worker_search_seconds = 0.0
    def collect(payloads):
        nonlocal checkpoint_seconds, worker_search_seconds
        for index, graph, counts in payloads:
            worker_search_seconds += counts['search_seconds']
            checkpoint_seconds += counts['checkpoint_seconds']
            if isinstance(graph, str):
                saved[index] = Path(graph)
            else:
                graphs[index] = graph
            for key in ('max_live_branches', 'max_growth_candidates'):
                metrics[key] = max(metrics[key], counts[key])

    graphs = [None] * len(cuts)
    saved = {int(path.name.split('_')[1].split('.')[0]):path for path in raw_cut_paths(directory)} if resume else {}
    missing = [index for index in range(len(cuts)) if index not in saved]
    raw_suffix = '.raw.pkl.gz' if archive_format=='checkpoint' else '.json'
    tasks = [(index, cuts[index], str(directory/f'cut_{index:05d}{raw_suffix}') if directory else None)
             for index in missing]
    # Finish independent missing cuts before rebuilding the large cached DAG.
    # This also prevents forked cut workers from inheriting that archive heap.
    workers = min(max(1, int(workers)), max(1,len(missing)))
    collect(iter_cut_tasks(problem, config, execution, tasks, workers,
                           in_process=requested_workers == 1))
    # These JSON records are acyclic. Repeated cyclic-GC scans while restoring
    # millions of retained tuples needlessly revisit the entire growing graph.
    # Restore the caller's GC setting before any new matching begins.
    restore_started=time.perf_counter()
    if requested_workers>1:
        payloads=[(index,str(saved[index]) if index in saved else graphs[index],
                   str(directory/f'cut_{index:05d}.finalized.pkl.gz') if directory else None)
                  for index in range(len(cuts))]
        gc_enabled=gc.isenabled()
        try:
            gc.disable()
            for index, graph, counts in iter_finalized_tasks(
                    problem, config, execution, payloads, min(requested_workers, len(cuts))):
                graphs[index] = graph
                for key, value in counts.items():
                    metrics[key] = metrics.get(key, 0) + value
        finally:
            if gc_enabled:gc.enable()
        del payloads
    else:
        gc_enabled = gc.isenabled()
        tuple_pool = {}
        try:
            gc.disable()
            for index in sorted(saved):
                graphs[index] = read_raw_cut(saved[index],tuple_pool=tuple_pool)
        finally:
            if gc_enabled:gc.enable()
        del tuple_pool
    metrics['checkpoint_restore_and_finalize_seconds']=time.perf_counter()-restore_started
    if execution != 'shared_policies':
        for graph in graphs:
            metrics['max_live_branches']=max(metrics['max_live_branches'],
                max((sum(graph.states[t].context==context for t in graph.terminals)
                     for context in range(len(graph.contexts))),default=0))
    merge_started = time.perf_counter()
    gc_enabled = gc.isenabled()
    try:
        gc.disable()
        graph = AAMSearchGraph.combine(graphs)
    finally:
        if gc_enabled:gc.enable()
    metrics['parent_merge_seconds'] = time.perf_counter() - merge_started
    target = build_graph(problem.product.elements, problem.product.wbo,
                         bond_cut=config.graph_floor)
    symmetry_started = time.perf_counter()
    workspace = None
    if execution in ('reused_native', 'shared_policies'):
        from mappa.mappa_core.core.matching.symmetry.workspace import ConditionedSymmetryWorkspace
        workspace = ConditionedSymmetryWorkspace(target, config.iso_tolerance)
    graph, groups = finalize_graph_symmetry(graph, target, iso_tolerance=config.iso_tolerance,
                                          workspace=workspace)
    metrics.update(symmetry_finalization_seconds=time.perf_counter()-symmetry_started,
                   worker_search_seconds=worker_search_seconds, checkpoint_seconds=checkpoint_seconds)
    for key,value in groups.items():metrics[key]=metrics.get(key,0)+value
    metrics.update(cuts=len(cuts), raw_result_count=len(graph.terminals),
        retained_branch_count=len(graph.terminals),
        subtree_branch_cap_count=sum(stop.reason == 'capped' for stop in graph.stops))
    result = AAMResult(problem, config, graph,
                      AAMSearchMetrics.from_record(metrics, time.perf_counter()-started))
    if directory is not None:
        from mappa.artifacts import aam_json, write_aam_checkpoint
        if archive_format=='checkpoint':
            write_aam_checkpoint(result,directory/'aam.pkl.gz')
        else:
            temporary = directory / 'aam.json.tmp'
            temporary.write_text(aam_json(result))
            temporary.replace(directory / 'aam.json')
    return result


@dataclass(frozen=True)
class AAMCheckpointResult:
    """Complete raw cut archives; symmetry is finalized when a cut is consumed."""
    directory: Path
    cut_paths: tuple[Path, ...]
    metrics: AAMSearchMetrics
    capped: bool


def search_aam_checkpoints(problem: AAMProblem, config: AAMSearchConfig | None = None,
                           *, intermediate_dir, resume=False, workers=1,
                           execution='reused_native') -> AAMCheckpointResult:
    """Persist every cut without materializing the combined search history.

    The search and cap semantics match search_aam. Each worker retains only its
    current cut and bounded native reuse cache. Raw archives retain all history;
    consumers finalize and release individual cuts. Memory can still grow with
    one cut's history; this is not a universal memory bound.
    """
    from mappa.artifacts import raw_cut_paths, read_raw_cut
    if not isinstance(problem, AAMProblem):
        raise TypeError('search_aam_checkpoints requires an AAMProblem')
    if execution not in ('reference', 'reused_native', 'shared_policies'):
        raise ValueError('unknown AAM execution backend')
    if execution != 'reference':
        from mappa.mappa_core.core.matching.growth.native import available
        if not available():
            raise ValueError('native execution requires the built native engine')
    config = config or AAMSearchConfig()
    started = time.perf_counter()
    directory = Path(intermediate_dir)
    directory.mkdir(parents=True, exist_ok=True)
    manifest = checkpoint_manifest(problem, config)
    if execution == 'shared_policies':
        manifest['execution'] = execution
    manifest_path = directory / 'manifest.json'
    if resume:
        if json.loads(manifest_path.read_text()) != manifest:
            raise ValueError('Checkpoint input or configuration differs from this search')
    else:
        if manifest_path.exists() or raw_cut_paths(directory):
            raise ValueError('Existing checkpoints require explicit resume=True')
        manifest_path.write_text(json.dumps(manifest) + '\n')
    cuts = cut_sweep_items(problem.reactant.wbo, config.cut_floor) if config.sweep_cuts else [()]
    saved = {int(p.name.split('_')[1].split('.')[0]): p for p in raw_cut_paths(directory)}
    if set(saved) - set(range(len(cuts))):
        raise ValueError('Unexpected cut checkpoint index')
    tasks = [(i, cut, str(directory / f'cut_{i:05d}.raw.pkl.gz'))
             for i, cut in enumerate(cuts) if i not in saved]
    metrics = dict(cuts=len(cuts), worker_search_seconds=0., checkpoint_seconds=0.,
                   max_growth_candidates=0, max_live_branches=0,
                   raw_result_count=0, retained_branch_count=0, subtree_branch_cap_count=0)
    def collect(results):
        for i, path, counts in results:
            saved[i] = Path(path)
            metrics['worker_search_seconds'] += counts['search_seconds']
            metrics['checkpoint_seconds'] += counts['checkpoint_seconds']
            metrics['max_growth_candidates'] = max(metrics['max_growth_candidates'], counts['max_growth_candidates'])
    worker_count = min(max(1, int(workers)), max(1, len(tasks)))
    if tasks:
        collect(iter_cut_tasks(problem, config, execution, tasks, worker_count))
    # Reduce archive statistics one cut at a time, including resumed cuts.
    paths = tuple(saved[i] for i in range(len(cuts)))
    for path in paths:
        graph = read_raw_cut(path)
        count = len(graph.terminals)
        metrics['raw_result_count'] += count
        metrics['retained_branch_count'] += count
        metrics['subtree_branch_cap_count'] += sum(stop.reason == 'capped' for stop in graph.stops)
        by_context = {}
        for terminal in graph.terminals:
            context = graph.states[terminal].context
            by_context[context] = by_context.get(context, 0) + 1
        metrics['max_live_branches'] = max(metrics['max_live_branches'], max(by_context.values(), default=0))
        del graph
    return AAMCheckpointResult(directory, paths,
        AAMSearchMetrics.from_record(metrics, time.perf_counter() - started),
        bool(metrics['subtree_branch_cap_count']))
