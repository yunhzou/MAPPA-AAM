"""Target event queries over the shared binary-symmetry response.

Selections retain the full relationship. Realization finds one requested member
with the same stabilizer and independence operations used for minimization.
"""
from dataclasses import dataclass
from math import isfinite
from time import perf_counter

import numpy as np

from .groups import event_components, support
from .minimum import _Search, SymmetryWork, SearchDeadlineReached


@dataclass(frozen=True)
class EventQuery:
    status: str
    witness: tuple | None
    choices: object
    event_count: int | None
    seconds: float
    work: SymmetryWork


@dataclass(frozen=True)
class ImageCondition:
    atom: int
    image: int
    scope: tuple

    def __iter__(self):
        return iter(self.scope)


@dataclass(frozen=True)
class EventSelection:
    response: object
    lower: int
    upper: int

    def __post_init__(self):
        if any(isinstance(v, bool) or not isinstance(v, int) for v in (self.lower, self.upper)):
            raise ValueError('event bounds must be integers')
        if not 0 <= self.lower <= self.upper:
            raise ValueError('event bounds must be ordered and nonnegative')

    def find(self, *, assignments=None, seconds=5., max_nodes=10000, cache_entries=2048):
        if not isfinite(seconds) or seconds <= 0 or max_nodes < 0 or cache_entries < 0:
            raise ValueError('a positive finite deadline and nonnegative budgets are required')
        data = self.response
        started = perf_counter()
        conditions = tuple(ImageCondition(atom, image, data.condition_scope(atom))
                           for atom, image in sorted((assignments or {}).items()))
        search = _EventSearch(data, max_nodes, cache_entries)
        search.deadline = started + seconds
        try:
            prefix = np.arange(data.size)
            count = data.fixed_events + int(data.score(prefix, data.pairs))
            mapping = data.mapping(prefix) if conditions else {}
            if self.lower <= count <= self.upper and all(mapping[c.atom] == c.image for c in conditions):
                pass
            elif self.upper < data.fixed_events or self.lower > data.fixed_events + len(data.pairs):
                prefix = None
            else:
                prefix = search.find(data.group, prefix, data.pairs,
                    self.lower - data.fixed_events, self.upper - data.fixed_events, conditions)
            status = 'absent' if prefix is None else 'found'
        except SearchDeadlineReached:
            prefix, status = None, 'unknown'
        witness, choices, count = None, None, None
        if prefix is not None:
            count = data.fixed_events + int(data.score(prefix, data.pairs))
            witness, choices = data.query_result(prefix, count, perf_counter() - started, search.work)
            assert all(dict(witness)[condition.atom] == condition.image for condition in conditions)
        return EventQuery(status, witness, choices, count, perf_counter() - started, search.work)


class _EventSearch(_Search):
    def find(self, group, prefix, pairs, lower, upper, conditions=()):
        lower, upper = max(0, lower), min(len(pairs), upper)
        if lower > upper:
            return None
        if perf_counter() >= self.deadline:
            self.work.time_limit_reached = True
            raise SearchDeadlineReached
        initial = int(self._score(prefix, pairs))
        images = self.data.mapping(prefix) if conditions else {}
        if lower <= initial <= upper and all(images[c.atom] == c.image for c in conditions):
            return prefix.copy()
        key, (moved, factors, point) = self._metadata(group)
        if not moved:
            return None
        moved_set = set(moved)
        active_conditions = tuple(c for c in conditions if moved_set.intersection(c))
        if any(images[c.atom] != c.image for c in conditions if c not in active_conditions):
            return None
        relevant = sorted(moved_set.union(atom for term in (*pairs, *conditions) for atom in term))
        cache_key = key, pairs, conditions, tuple(int(prefix[i]) for i in relevant), lower, upper
        if cache_key in self.answers:
            self.work.cache_hits += 1
            self.answers.move_to_end(cache_key)
            changes = self.answers[cache_key]
            if changes is None:
                return None
            answer = prefix.copy()
            for atom, image in changes:
                answer[atom] = image
            return answer
        if self.work.nodes >= self.max_nodes:
            self.work.cap_reached = True
            raise SearchDeadlineReached
        self.work.nodes += 1
        active = tuple(pair for pair in pairs if moved_set.intersection(pair))
        fixed = initial - int(self._score(prefix, active))
        low, high = max(0, lower - fixed), min(len(active), upper - fixed)
        answer = None
        if low <= high:
            components = event_components(factors, (*active, *active_conditions))
            self.work.event_coupling_merges += len(factors) - len(components)
            if len(components) > 1:
                self.work.independent_splits += 1
                supports = [set(support(component)) for component in components]
                terms = [tuple(pair for pair in active if atoms.intersection(pair)) for atoms in supports]
                constraints = [tuple(c for c in active_conditions if atoms.intersection(c)) for atoms in supports]
                answer = self._join(components, terms, constraints, prefix, low, high)
            else:
                affected = moved_set.intersection(atom for term in (*active, *active_conditions) for atom in term)
                point = min(affected) if affected else point
                step_key = key, point
                if step_key not in self.steps:
                    representatives = sorted(group.transversal(point),
                        key=lambda item: (item[0] != point, item[0]))
                    self.steps[step_key] = group.stabilizer((point,)), representatives
                    self.work.group_steps_prepared += 1
                subgroup, representatives = self.steps[step_key]
                # Count invariance alone does not preserve a requested endpoint.
                if not active_conditions:
                    representatives = self._representatives(prefix, representatives, moved)
                for _, operation in representatives:
                    image = prefix[np.asarray(operation[:self.data.size])]
                    answer = self.find(subgroup, image, active, low, high, active_conditions)
                    if answer is not None:
                        break
        if self.cache_entries:
            self.answers[cache_key] = (None if answer is None else
                tuple((atom, int(answer[atom])) for atom in moved))
            if len(self.answers) > self.cache_entries:
                self.answers.popitem(last=False)
            self.work.peak_cache_entries = max(self.work.peak_cache_entries, len(self.answers))
        return answer

    def _join(self, groups, terms, conditions, prefix, lower, upper):
        if not groups:
            return prefix if lower <= 0 <= upper else None
        if len(groups) == 1:
            return self.find(groups[0], prefix, terms[0], lower, upper, conditions[0])
        remaining = sum(map(len, terms[1:]))
        for count in range(max(0, lower - remaining), min(len(terms[0]), upper) + 1):
            image = self.find(groups[0], prefix, terms[0], count, count, conditions[0])
            if image is not None:
                answer = self._join(groups[1:], terms[1:], conditions[1:], image, lower - count, upper - count)
                if answer is not None:
                    return answer
        return None
