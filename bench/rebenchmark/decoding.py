"""Decode saved Golden seed/cut graphs through the public event decoder.

Each seed/cut retains its original compressed mappings and local event functions.
Historical converted-branch summaries remain explicitly scoped to their schema.
"""
import time

from mappa.mappa_core.core.scheduler import Task
from .saved_analysis import prepare_cut
from .search import select_contexts


def seed_tasks(cut_tasks, *, seeds=10, slack=5, branch_limit=100, random_seed=42, case_seconds=3600., memory_bytes=4*2**30):
    """Reuse the verified cut index; one independent seed/cut per queue task."""
    for cut in cut_tasks:
        for seed in range(seeds):
            saved = {name: cut.payload[name] for name in
                     ('task_id', 'cut', 'source', 'raw', 'graph', 'finalized')}
            payload = dict(saved, seed=seed, slack=slack, branch_limit=branch_limit, random_seed=random_seed,
                           decoder='compressed_event_interpretation')
            yield Task(f'{payload["task_id"]:05d}/cut/{payload["cut"]:05d}/seed/{seed:02d}',
                       'rebenchmark.decoding:decode_seed', payload,
                       case=f'{payload["task_id"]:05d}', case_seconds=case_seconds,
                       memory_bytes=memory_bytes)


def coverage(problem, family):
    """The original Golden heavy-atom/total-atom tier, before event ranking."""
    return (sum(problem.reactant.elements[a] != 'H' for a, _ in family.mapping), len(family.mapping))


def decode_seed(payload, workspace):
    if payload.get('decoder') != 'compressed_event_interpretation' or 'family_seconds' in payload:
        raise ValueError('Create current event tasks explicitly; retired solver tasks cannot be resumed here')
    started, cpu = time.perf_counter(), time.process_time()
    plan, _, graph, digest, finalized = prepare_cut(payload, workspace)
    loading = dict(wall_seconds=time.perf_counter()-started, cpu_seconds=time.process_time()-cpu)
    value = decode_graph(plan.problem, graph, payload)
    value['timing']['loading'] = loading
    return dict(value, task_id=payload['task_id'], cut=payload['cut'],
        source=payload['source'], raw=payload['raw'], graph=payload['graph'],
        finalized=finalized, finalized_sha256=digest, search_rerun=False)


def dense_event_categories(problem, mapping):
    """Independent original-matrix recount for saved representative validation."""
    import numpy as np
    from mappa.mappa_core.utils.bonds import is_metal_element
    counts = np.zeros(4,dtype=int)
    for reverse in (False,True):
        ep,other,images = ((problem.reactant,problem.product,mapping) if not reverse else
                          (problem.product,problem.reactant,{p:r for r,p in mapping.items()}))
        a,b = np.triu_indices(ep.atom_count,1)
        mapped = np.array([images.get(i,-1) for i in range(ep.atom_count)])
        touched = (mapped[a]>=0)|(mapped[b]>=0)
        both = (mapped[a]>=0)&(mapped[b]>=0)
        opposite = np.zeros(len(a))
        opposite[both] = other.wbo[mapped[a[both]],mapped[b[both]]]
        metal = np.array([is_metal_element(e) for e in ep.elements])
        tau = np.where(metal[a]|metal[b],.3,.5)
        decrease = touched & (ep.wbo[a,b]-opposite>tau)
        persists = opposite>tau
        counts[1 if reverse else 0] += np.count_nonzero(decrease & ~persists)
        counts[2 if reverse else 3] += np.count_nonzero(decrease & persists)
    return dict(zip(('breaking','forming','strengthening','weakening'),map(int,counts)))


def decode_graph(problem, graph, payload):
    """Prepare every retained family's local event functions without rematching."""
    from mappa.mappa_core.utils.io.records import decoded_events_record
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    from mappa.mappa_core.post_processing.decoder.api import decode_events
    started, cpu = time.perf_counter(), time.process_time()
    catalogue = FinalBranchCatalogue(problem).add_graph(select_contexts(graph, [payload['seed']]))
    all_families = len(catalogue.families)
    tier = max((coverage(problem, f) for f in catalogue.families), default=None)
    selected_ids = [i for i,f in enumerate(catalogue.families) if coverage(problem,f) == tier]
    setup = dict(wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu)
    started, cpu = time.perf_counter(), time.process_time()
    result = decode_events(catalogue)
    for fid, partitions in enumerate(result.fragment_partitions):
        for partition in range(len(partitions)):
            result.tables(fid, partition)
    decoding = dict(wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu)
    started, cpu = time.perf_counter(), time.process_time()
    # Independent dense recounts check the sparse factors on each saved branch.
    for fid, partitions in enumerate(result.fragment_partitions):
        vector = result.index.vector(result.realize(fid))
        expected = dense_event_categories(problem,result.realize(fid))
        for partition in range(len(partitions)):
            actual = result.tables(fid,partition).evaluate(vector)['counts']
            if actual != expected:
                raise ValueError(f'event table recount differs for family {fid}')
    verification = dict(wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu)
    started, cpu = time.perf_counter(), time.process_time()
    record = decoded_events_record(result)
    serialization = dict(wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu)
    return dict(seed=payload['seed']+1, coverage=tier,
        decoding_scope='original_compressed_mapping_space', selected_family_ids=selected_ids,
        all_families=all_families,selected_families=len(selected_ids),decoded=record,
        timing=dict(setup=setup,decoding=decoding,verification=verification,record=serialization))


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
