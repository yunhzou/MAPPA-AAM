"""Decode saved Golden seed/cut graphs through the public event decoder.

Each seed/cut converts its saved fragment pairs and counts the resulting branches.
Local minimum+slack windows contain the corresponding prefix-union window.
"""
import time

from mappa.execution import Task
from .saved_analysis import prepare_cut
from .search import select_contexts


def seed_tasks(cut_tasks, *, seeds=10, slack=5, branch_limit=100, random_seed=42, case_seconds=3600., memory_bytes=4*2**30):
    """Reuse the verified cut index; one independent seed/cut per queue task."""
    for cut in cut_tasks:
        for seed in range(seeds):
            saved = {name: cut.payload[name] for name in
                     ('task_id', 'cut', 'source', 'raw', 'graph', 'finalized')}
            payload = dict(saved, seed=seed, slack=slack, branch_limit=branch_limit, random_seed=random_seed,
                           decoder='anchored_fragment_conversion')
            yield Task(f'{payload["task_id"]:05d}/cut/{payload["cut"]:05d}/seed/{seed:02d}',
                       'rebenchmark.decoding:decode_seed', payload,
                       case=f'{payload["task_id"]:05d}', case_seconds=case_seconds,
                       memory_bytes=memory_bytes)


def coverage(problem, family):
    """The original Golden heavy-atom/total-atom tier, before event ranking."""
    return (sum(problem.reactant.elements[a] != 'H' for a, _ in family.mapping), len(family.mapping))


def representative_bound(problem, graph, seed, tier, slack):
    """Earlier same-tier witnesses bound this prefix without any optimization.

    Every requested prefix containing this seed also contains these witnesses.
    A larger cut union can only lower the minimum within its chosen tier.
    """
    from mappa.event_patterns import SignedEventIndex
    index = SignedEventIndex(problem)
    vectors, best = [], None
    for terminal in graph.terminals:
        state = graph.states[terminal]
        if state.context <= seed and coverage(problem, state) == tier:
            vectors.append(index.vector(dict(state.mapping)))
            if len(vectors) == 32:
                value = int(index.counts(vectors).sum(axis=1).min())
                best = value if best is None else min(best, value)
                vectors = []
    if vectors:
        value = int(index.counts(vectors).sum(axis=1).min())
        best = value if best is None else min(best, value)
    return None if best is None else best+slack


def decode_seed(payload, workspace):
    from mappa.artifacts import decoded_events_record
    from mappa.final_branches import FinalBranchCatalogue
    from mappa.frag import classify_bonds
    from mappa.postprocessing import EventDecodeConfig, decode_events
    if payload.get('decoder') != 'anchored_fragment_conversion' or 'family_seconds' in payload:
        raise ValueError('Create conversion tasks explicitly; retired solver tasks cannot be resumed here')
    started, cpu = time.perf_counter(), time.process_time()
    plan, _, graph, digest, finalized = prepare_cut(payload, workspace)
    one = select_contexts(graph, [payload['seed']])
    catalogue = FinalBranchCatalogue(plan.problem).add_graph(one)
    all_families = len(catalogue.families)
    tier = max((coverage(plan.problem, f) for f in catalogue.families), default=None)
    selected_ids = {i for i, f in enumerate(catalogue.families) if coverage(plan.problem, f) == tier}
    catalogue.branches = {key: ids & selected_ids for key, ids in catalogue.branches.items()
                          if ids & selected_ids}
    catalogue.coupled_branches = {i: pairs for i, pairs in catalogue.coupled_branches.items()
                                  if i in selected_ids}
    setup = dict(wall_seconds=time.perf_counter()-started, cpu_seconds=time.process_time()-cpu)
    started, cpu = time.perf_counter(), time.process_time()
    result = decode_events(catalogue, EventDecodeConfig(event_slack=payload['slack'],
        branch_limit=payload['branch_limit'], random_seed=payload['random_seed']))
    timing = dict(wall_seconds=time.perf_counter()-started, cpu_seconds=time.process_time()-cpu)
    for candidate in result.candidates:
        broken, formed, *_ = classify_bonds(candidate.mapping, plan.problem.reactant.wbo,
            plan.problem.product.wbo, dwbo_threshold=result.config.threshold,
            elements_R=plan.problem.reactant.elements, elements_P=plan.problem.product.elements,
            metal_dwbo_threshold=result.config.metal_threshold)
        assert candidate.total == len(broken)+len(formed)
        assert candidate.endpoint_events == dict(broken=tuple(sorted((a,b) for a,b,*_ in broken)),
                                                 formed=tuple(sorted((a,b) for a,b,*_ in formed)))
    record = decoded_events_record(result)
    return dict(task_id=payload['task_id'], cut=payload['cut'], seed=payload['seed']+1,
        coverage=tier, conversion_scope='retained_converted_branches',
        all_families=all_families, selected_families=len(selected_ids),
        decoded=record, families=[f.to_record() for f in result.catalogue.families],
        source=payload['source'], raw=payload['raw'], graph=payload['graph'],
        finalized=finalized, finalized_sha256=digest, search_rerun=False,
        timing=dict(setup=setup, decoding=timing))


def merge_windows(rows, *, slack=5):
    """Combine completed local windows at the original best coverage tier.

    Missing rows must be accounted for by the caller's manifest. A positive
    completeness flag here concerns exactly the provided rows, never a dataset.
    """
    rows = list(rows)
    if len({r['decoded'].get('schema', 'mappa.decoded_events/v1') for r in rows}) > 1:
        raise ValueError('Cannot combine solver-decoded windows with converted-branch windows')
    eligible = [row for row in rows if row['coverage'] is not None]
    tier = max((tuple(row['coverage']) for row in eligible), default=None)
    selected = [row for row in eligible if tuple(row['coverage']) == tier]
    minima = [row['decoded']['observed_minimum'] for row in selected if row['decoded']['observed_minimum'] is not None]
    minimum = min(minima, default=None)
    window = minimum+slack if minimum is not None else None
    for row in selected:
        config = row['decoded']['config']
        if config['event_slack'] is not None:
            sufficient = config['event_slack'] >= slack
        else:
            sufficient = config['max_events'] is None or (window is not None and config['max_events'] >= window)
        if not sufficient:
            raise ValueError('Saved windows do not cover this prefix; include its earlier witnesses')
    candidates = {}
    for row in selected:
        for candidate in row['decoded']['candidates']:
            if candidate['total'] <= window:
                candidates.setdefault(candidate['id'], candidate)
    complete = minimum is not None and bool(selected) and all(row['decoded']['complete'] for row in selected)
    return dict(coverage=tier, minimum=minimum, minimum_proven=minimum is not None and complete,
        complete=complete, event_window=window, candidates=sorted(candidates.values(),key=lambda c:(c['total'],c['id'])),
        rows=len(rows), selected_rows=len(selected))
