"""Time-bounded improvement of a Fragment's starting correspondence.

The complete saved relation stays intact. Search improves an actual witness;
attaining a lower bound or certifying complete symmetry coverage proves a minimum.
"""
from dataclasses import fields, replace
from math import isfinite
from time import perf_counter

import numpy as np

from .minimum import (FragmentMinimum, FragmentSymmetryEvents, SearchDeadlineReached,
                      SymmetryWork, check_deadline)
from .objective import FragmentObjective


class _StartingMappings:
    """Native bounded discoveries; only improvements cross the Python boundary."""
    def __init__(self, fragment, objective, deadline, ceiling, work):
        from mappa._group_ops import AssignmentSearch
        self.relation = fragment.relation
        self.deadline, self.ceiling, self.work = deadline, ceiling, work
        sources = {atom: i for i, atom in enumerate(self.relation.source_atoms)}
        targets = objective.target_positions
        domains = {r: tuple(sorted(targets[p] for p in values))
                   for r, values in self.relation.domains}
        domains.update((r, (targets[p],)) for r, p in self.relation.anchors)
        edges = [(sources[edge.atoms[0]], sources[edge.atoms[1]],
                  [(targets[p], targets[q]) for p, q in self.relation.discovery_pairs(edge)])
                 for edge in self.relation.accepted_edges]
        source_rows, target_rows = np.sort(objective.source, axis=1), np.sort(objective.target, axis=1)
        # Keep the existing floating-point proposal order, including its ties.
        # Only one row-sized temporary is needed, rather than an n-cubed array.
        priorities = [np.sum(np.abs(row - target_rows), axis=1).tolist() for row in source_rows]
        degrees = np.count_nonzero(fragment.source_symmetry.weights >
                                  fragment.source_symmetry.bond_threshold, axis=1)
        self.search = AssignmentSearch(objective.source.tolist(), objective.target.tolist(),
            objective.thresholds.tolist(), priorities,
            [domains[r] for r in self.relation.source_atoms], degrees.tolist(), edges)

    def __iter__(self):
        while True:
            images = self.search.next(self.ceiling, max(0., self.deadline - perf_counter()))
            self.work.discovery_nodes = self.search.nodes
            self.work.discovery_improvements = self.search.improvements
            if self.search.expired:
                raise SearchDeadlineReached
            if images is None:
                return
            witness = tuple((r, self.relation.target_atoms[p])
                            for r, p in zip(self.relation.source_atoms, images))
            self.relation._validate_discovery(witness)
            yield witness


def _accumulate(total, part):
    for item in fields(total):
        name = item.name
        before, value = getattr(total, name), getattr(part, name)
        setattr(total, name, (before or value) if isinstance(value, bool) else
                max(before, value) if name.startswith('peak_') else before + value)


def optimize_fragment(fragment, *, time_limit=5., event_problem=None,
                      threshold=.5, metal_threshold=.3, cache_entries=2048):
    """Improve one valid event witness under a single wall-clock search budget.

    Input is the saved Fragment; atom sets, constraints and candidate count are
    unchanged. A full-relation proof requires attaining the weight lower bound
    or solving symmetry that covers every discovery. Otherwise the witness is
    best found. The same deadline covers preparation, new discoveries
    and all symmetry searches. Checks are cooperative; an opaque native call
    can overrun until it returns. Execution-layer watchdogs enforce hard limits.
    """
    if not isfinite(time_limit) or time_limit < 0:
        raise ValueError('time_limit must be finite and nonnegative')
    if cache_entries < 0:
        raise ValueError('cache_entries must be nonnegative')
    started = perf_counter()
    deadline = started + time_limit
    objective = FragmentObjective(fragment, event_problem, threshold, metal_threshold)
    identity = tuple(range(objective.size))
    work = SymmetryWork()
    best = FragmentMinimum(fragment.witness, objective.count(objective.images(fragment.witness)),
        False, 'saved witness', 0., work, identity, identity,
        scope='fragment_relation', discovery_witness=fragment.witness)
    lower_bound = 0
    covered_symmetry = False

    def optimize_start(witness):
        nonlocal best, lower_bound, covered_symmetry
        # Publish the valid discovery before any potentially interrupted work.
        count = objective.count(objective.images(witness))
        if count < best.event_count:
            best = replace(best, witness=witness, event_count=count,
                discovery_witness=witness, source_operation=identity, target_operation=identity)
        if best.event_count == lower_bound:
            return
        check_deadline(deadline)
        current = replace(fragment, relation=replace(fragment.relation, witness=witness))
        prepared = FragmentSymmetryEvents(current, objective=objective, deadline=deadline)
        answer = prepared.solve(max_nodes=None, cache_entries=cache_entries)
        _accumulate(work, answer.work)
        if answer.event_count < best.event_count:
            best = replace(answer, work=work, scope='fragment_relation')
        if answer.complete and fragment.relation.discovery_in_witness_orbit:
            lower_bound = answer.event_count
            covered_symmetry = True

    try:
        check_deadline(deadline)
        lower_bound = objective.lower_bound()
        if best.event_count > lower_bound:
            lower_bound = max(lower_bound, objective.incidence_lower_bound(
                max(0., deadline - perf_counter())))
        check_deadline(deadline)
        optimize_start(fragment.witness)
        if best.event_count > lower_bound:
            starts = _StartingMappings(fragment, objective, deadline, best.event_count, work)
            for witness in starts:
                optimize_start(witness)
                starts.ceiling = best.event_count
                if best.event_count == lower_bound:
                    break
    except SearchDeadlineReached:
        work.time_limit_reached = True
    complete = best.event_count == lower_bound
    elapsed = perf_counter() - started
    return replace(best, complete=complete, lower_bound=lower_bound,
        scope='fragment_relation', elapsed_seconds=elapsed, search_seconds=elapsed,
        reason=('complete discovery covered by solved symmetry' if covered_symmetry else
                'weight lower bound attained' if complete else
                'time limit reached; best found' if work.time_limit_reached else
                'no better discovery found; full-relation minimum unproved'))
