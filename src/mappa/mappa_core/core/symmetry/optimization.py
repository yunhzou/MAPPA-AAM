"""Time-bounded improvement of a Fragment's starting correspondence.

The complete saved relation stays intact. Search improves an actual witness;
attaining a lower bound or certifying complete symmetry coverage proves a minimum.
"""
from dataclasses import fields, replace
from math import isfinite
from time import perf_counter

import numpy as np

from .constraints import FiniteConstraints, compile_discovery
from .minimum import (FragmentMinimum, FragmentSymmetryEvents, SearchDeadlineReached,
                      SymmetryWork, check_deadline)
from .objective import FragmentObjective


class _StartingMappings:
    """Search for better discoveries with the shared growth-constraint compiler.

    Partial event counts bound the unshuffled discovery objective only. They do
    not certify unexplored symmetry classes: this is deliberate bounded witness
    improvement, not exhaustive event decoding or a replacement growth matcher.
    """
    def __init__(self, fragment, objective, deadline, ceiling, work):
        self.relation = fragment.relation
        self.objective = objective
        self.deadline = deadline
        self.ceiling = ceiling
        self.work = work
        self.system = FiniteConstraints()
        self.values = compile_discovery(self.system, self.relation, self.check)
        self.source_rows = np.sort(objective.source, axis=1)
        self.target_rows = np.sort(objective.target, axis=1)
        self.degrees = np.count_nonzero(
            fragment.source_symmetry.weights > fragment.source_symmetry.bond_threshold, axis=1)

    def check(self):
        check_deadline(self.deadline)

    def _domains(self, system):
        positions = self.objective.target_positions
        return [tuple(sorted(positions[p] for p in system.domain(self.values[r][0])))
                for r in self.relation.source_atoms]

    def _bound(self, domains):
        cost = 0
        data = self.objective
        for i, j in data.pairs:
            # Unfixed pairs contribute zero. Fixed/unfixed terms use their best
            # remaining endpoint, without assuming those separate choices agree.
            if len(domains[i]) == 1 or len(domains[j]) == 1:
                values = data.target[np.ix_(domains[i], domains[j])]
                cost += int(np.all(np.abs(data.source[i, j] - values) > data.thresholds[i, j]))
                if cost >= self.ceiling:
                    break
        return cost

    def _search(self, system):
        self.check()
        self.work.discovery_nodes += 1
        system.propagate(self.check)
        if system.failed:
            return
        domains = self._domains(system)
        if self._bound(domains) >= self.ceiling:
            return
        moving = [i for i, values in enumerate(domains) if len(values) > 1]
        if not moving:
            images = tuple(values[0] for values in domains)
            witness = tuple((r, self.relation.target_atoms[p])
                            for r, p in zip(self.relation.source_atoms, images))
            self.relation._validate_discovery(witness)
            self.work.discovery_improvements += 1
            yield witness
            return
        atom = min(moving, key=lambda i: (len(domains[i]), -int(self.degrees[i]), i))
        # Weight profiles order proposals only; all pruning uses valid bounds.
        choices = sorted(domains[atom], key=lambda p: (
            float(np.sum(np.abs(self.source_rows[atom] - self.target_rows[p]))), p))
        value = self.values[self.relation.source_atoms[atom]][0]
        for image in choices:
            self.check()
            child = FiniteConstraints()
            child.domains = [set(domain) for domain in system.domains]
            child.constraints = list(system.constraints)
            child.restrict(value, {self.relation.target_atoms[image]})
            yield from self._search(child)

    def __iter__(self):
        return self._search(self.system)


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
