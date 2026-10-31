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


def dense_event_categories(problem, mapping, threshold=.5, metal_threshold=.3):
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
        tau = np.where(metal[a]|metal[b], threshold if metal_threshold is None else metal_threshold, threshold)
        decrease = touched & (ep.wbo[a,b]-opposite>tau)
        persists = opposite>tau
        counts[1 if reverse else 0] += np.count_nonzero(decrease & ~persists)
        counts[2 if reverse else 3] += np.count_nonzero(decrease & persists)
    return dict(zip(('breaking','forming','strengthening','weakening'),map(int,counts)))


def decode_graph(problem, graph, payload):
    """Saved-input adapter; use the same core interpretation as fresh cuts."""
    from mappa.mappa_core.core.scheduler.decoding import decode_graph as interpret
    value = interpret(problem, select_contexts(graph, [payload['seed']]))
    seed = value.pop('seeds')[0]
    verify_event_output(value)
    return dict(value, seed=payload['seed']+1, coverage=seed['coverage'],
        decoding_scope='original_compressed_mapping_space', selected_family_ids=seed['selected_family_ids'],
        all_families=len(seed['family_ids']), selected_families=len(seed['selected_family_ids']))


def verify_event_output(value):
    """Check saved tables and representatives against dense endpoint counts."""
    import numpy as np
    from collections import defaultdict
    from mappa.mappa_core.utils.io.records import problem_from_record
    record = value['decoded']
    problem = problem_from_record(record['problem'])
    representatives = record['representatives']
    families = record['catalogue']['families']
    if len(representatives) != len(families):
        raise ValueError('Incomplete representative catalogue')
    tables = []
    for table in record['tables']:
        groups = defaultdict(list)
        for region in (*table['internal'], *table['boundary']):
            for term in region['terms']:
                groups[term['response']].append(term['atoms'])
        tables.append(({rid: np.array(pairs, dtype=int) for rid, pairs in groups.items()},
                       np.array(table['product_boundary_pairs'], dtype=int).reshape(-1, 2)))
    product = np.pad(problem.product.wbo, ((0, 1), (0, 1)))
    for fid, (representative, family) in enumerate(zip(representatives, families)):
        mapping = {a: representative['mapping'][a] for a in representative['mapped_source']}
        if representative['family_id'] != fid or sorted(mapping.items()) != [tuple(p) for p in family['mapping']]:
            raise ValueError('Representative mapping differs from saved family')
        counts = dense_event_categories(problem, mapping, record['event_policy']['threshold'],
                                        record['event_policy']['metal_threshold'])
        expected = dict(broken=counts['breaking']+counts['weakening'],
                        formed=counts['forming']+counts['strengthening'])
        if representative['counts'] != expected or representative['total'] != sum(expected.values()):
            raise ValueError('Independent representative bond recount differs')
        vector = np.array(representative['mapping'], dtype=int)
        mapped = np.zeros(problem.target_atom_count, dtype=bool)
        mapped[list(mapping.values())] = True
        for tid in record['families'][fid]['tables']:
            groups, boundary = tables[tid]
            actual = np.zeros(5, dtype=int)
            for rid, pairs in groups.items():
                response = record['responses'][rid]
                weights = product[vector[pairs[:, 0]], vector[pairs[:, 1]]]
                allowed = np.array(response['target_weights'])
                indices = np.searchsorted(allowed, weights)
                if np.any(indices >= len(allowed)) or not np.array_equal(allowed[indices], weights):
                    raise ValueError('Saved table omits a representative bond weight')
                codes = np.array(response['categories'])[indices]
                actual += np.bincount(codes, minlength=5)
            actual[2] += np.count_nonzero(mapped[boundary[:, 0]] ^ mapped[boundary[:, 1]])
            if tuple(actual[1:]) != tuple(counts.values()):
                raise ValueError('Saved event table differs from dense representative recount')
    return len(representatives)


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
