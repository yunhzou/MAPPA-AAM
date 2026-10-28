"""Decode saved Golden seed/cut graphs through the public event decoder.

Each seed/cut gets its own minimum+slack window. Their union therefore contains
that window for every seed prefix and any union of cuts, without a minimum stage.
"""
from dataclasses import asdict
import hashlib
import time

from mappa.execution import Task
from .saved_analysis import prepare_cut
from .search import select_contexts


def seed_tasks(cut_tasks, *, seeds=10, slack=5, family_seconds=10., case_seconds=3600., memory_bytes=4*2**30):
    """Reuse the verified cut index; one independent seed/cut per queue task."""
    for cut in cut_tasks:
        for seed in range(seeds):
            payload = dict(cut.payload, seed=seed, slack=slack, family_seconds=family_seconds)
            yield Task(f'{payload["task_id"]:05d}/cut/{payload["cut"]:05d}/seed/{seed:02d}',
                       'rebenchmark.decoding:decode_seed', payload,
                       case=f'{payload["task_id"]:05d}', case_seconds=case_seconds,
                       memory_bytes=memory_bytes)


def coverage(problem, family):
    """The original Golden heavy-atom/total-atom tier, before event ranking."""
    return (sum(problem.reactant.elements[a] != 'H' for a, _ in family.mapping), len(family.mapping))


def decode_seed(payload, workspace):
    from mappa.artifacts import decoded_events_record
    from mappa.final_branches import FinalBranchCatalogue
    from mappa.frag import classify_bonds
    from mappa.postprocessing import EventDecodeConfig, decode_events
    started, cpu = time.perf_counter(), time.process_time()
    plan, _, graph, digest, finalized = prepare_cut(payload, workspace)
    one = select_contexts(graph, [payload['seed']])
    catalogue = FinalBranchCatalogue(plan.problem).add_graph(one)
    all_families = len(catalogue.families)
    tier = max((coverage(plan.problem, f) for f in catalogue.families), default=None)
    catalogue.families = [f for f in catalogue.families if coverage(plan.problem, f) == tier]
    setup = dict(wall_seconds=time.perf_counter()-started, cpu_seconds=time.process_time()-cpu)
    def query_log(text, status, seconds):
        return asdict(workspace.writer.put(dict(smt2=text, status=status, cpu_seconds=seconds,
            sha256=hashlib.sha256(text.encode()).hexdigest())))
    started, cpu = time.perf_counter(), time.process_time()
    result = decode_events(catalogue, EventDecodeConfig(event_slack=payload['slack'],
        seconds_per_family=payload['family_seconds']), query_log=query_log)
    timing = dict(wall_seconds=time.perf_counter()-started, cpu_seconds=time.process_time()-cpu)
    for candidate in result.candidates:
        broken, formed, *_ = classify_bonds(candidate.mapping, plan.problem.reactant.wbo,
            plan.problem.product.wbo, dwbo_threshold=result.config.threshold,
            elements_R=plan.problem.reactant.elements, elements_P=plan.problem.product.elements,
            metal_dwbo_threshold=result.config.metal_threshold)
        assert candidate.total == len(broken)+len(formed)
        assert candidate.endpoint_events == dict(broken=tuple(sorted(broken)), formed=tuple(sorted(formed)))
    record = decoded_events_record(result)
    return dict(task_id=payload['task_id'], cut=payload['cut'], seed=payload['seed']+1,
        coverage=tier, all_families=all_families, selected_families=len(catalogue.families),
        decoded=record, families=[f.to_record() for f in catalogue.families],
        source=payload['source'], raw=payload['raw'], graph=payload['graph'],
        finalized=finalized, finalized_sha256=digest, search_rerun=False,
        timing=dict(setup=setup, decoding=timing))


def merge_windows(rows, *, slack=5):
    """Combine completed local windows at the original best coverage tier.

    Missing rows must be accounted for by the caller's manifest. A positive
    completeness flag here concerns exactly the provided rows, never a dataset.
    """
    rows = list(rows)
    eligible = [row for row in rows if row['coverage'] is not None]
    tier = max((tuple(row['coverage']) for row in eligible), default=None)
    selected = [row for row in eligible if tuple(row['coverage']) == tier]
    if any(row['decoded']['config']['event_slack'] < slack for row in selected):
        raise ValueError('Saved local windows are too narrow')
    minima = [row['decoded']['observed_minimum'] for row in selected if row['decoded']['observed_minimum'] is not None]
    minimum = min(minima, default=None)
    window = minimum+slack if minimum is not None else None
    candidates = {}
    for row in selected:
        for candidate in row['decoded']['candidates']:
            if candidate['total'] <= window:
                candidates.setdefault(candidate['id'], candidate)
    complete = bool(selected) and all(row['decoded']['complete'] for row in selected)
    return dict(coverage=tier, minimum=minimum, minimum_proven=minimum is not None and (minimum == 0 or complete),
        complete=complete, event_window=window, candidates=sorted(candidates.values(),key=lambda c:(c['total'],c['id'])),
        rows=len(rows), selected_rows=len(selected))
