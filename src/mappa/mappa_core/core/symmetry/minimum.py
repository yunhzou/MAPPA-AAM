"""Capped minimum-event witnesses over a fragment's two-sided binary symmetry.

The saved discovery witness fixes this orbit. Other discovery mappings remain
part of the Fragment relation and are not silently replaced by its orbit.
"""
from collections import OrderedDict
from dataclasses import dataclass
from functools import cached_property
from time import perf_counter

import numpy as np

from .groups import (binary_group, constant_response, event_components, group_key,
                     independent_factors, pair_orbits, permutation_group,
                     product_group, support)
from .objective import FragmentObjective


class SearchDeadlineReached(Exception):
    """An explicit search limit; the caller retains its valid incumbent."""


def check_deadline(deadline):
    if deadline is not None and perf_counter() >= deadline:
        raise SearchDeadlineReached


@dataclass
class SymmetryWork:
    nodes: int = 0
    score_evaluations: int = 0
    bond_comparisons: int = 0
    cache_hits: int = 0
    peak_cache_entries: int = 0
    independent_splits: int = 0
    event_coupling_merges: int = 0
    equivalent_positions_skipped: int = 0
    cap_reached: bool = False
    restriction_cache_hits: int = 0
    restriction_parent_reuses: int = 0
    peak_restriction_cache_entries: int = 0
    group_steps_prepared: int = 0
    discovery_nodes: int = 0
    discovery_improvements: int = 0
    time_limit_reached: bool = False


@dataclass(frozen=True)
class FragmentMinimum:
    witness: tuple[tuple[int, int], ...]
    event_count: int
    complete: bool
    reason: str
    search_seconds: float
    work: SymmetryWork
    source_operation: tuple[int, ...]
    target_operation: tuple[int, ...]
    scope: str = 'saved_witness_binary_symmetry'
    discovery_witness: tuple[tuple[int, int], ...] = ()
    lower_bound: int | None = None
    elapsed_seconds: float = 0.


@dataclass(frozen=True)
class _Answer:
    cost: int
    images: tuple[tuple[int, int], ...]
    complete: bool


class FragmentSymmetryEvents:
    """Prepare the actual Fragment once, retaining its anchors and atom frames.

    Minimum and target-count queries use the saved witness's binary-symmetry
    orbit. They do not exhaust the Fragment's other discovery alternatives.
    """

    def __init__(self, fragment, *, event_problem=None, threshold=.5, metal_threshold=.3,
                 deadline=None, objective=None):
        started = perf_counter()
        self.deadline = deadline
        check_deadline(deadline)

        self.fragment = fragment
        source, target = fragment.source_symmetry, fragment.target_symmetry
        positions = {atom: i for i, atom in enumerate(target.atoms)}
        images = dict(fragment.witness)
        self.witness = tuple(positions[images[a]] for a in source.atoms)
        self.size = len(self.witness)
        inverse = np.argsort(self.witness)
        self.source_binary = source.weights > source.bond_threshold
        self.target_binary = (target.weights > target.bond_threshold)[np.ix_(self.witness, self.witness)]
        objective = objective or FragmentObjective(fragment, event_problem, threshold, metal_threshold)
        self.source = objective.source
        self.target = objective.target[np.ix_(self.witness, self.witness)]
        self.thresholds = objective.thresholds
        self.binary = self.target_binary
        self.edges = tuple(zip(*np.where(np.triu(self.binary, 1))))
        anchored = {i for i, atom in enumerate(source.atoms) if atom in source.anchors}
        self.colors = [('atom', element, i if i in anchored else None)
                       for i, element in enumerate(source.elements)]
        self.target_colors = [('atom', target.elements[p], i if i in anchored else None)
                              for i, p in enumerate(self.witness)]
        self.group = permutation_group(
            [tuple(int(inverse[g[p]]) for p in self.witness) for g in target.generators], self.size)
        self.source_group = permutation_group(source.generators, self.size)
        check_deadline(deadline)
        self.pair_orbit_states = 0
        self.two_sided = (not np.array_equal(self.source_binary, self.binary)
                          or self.colors != self.target_colors)
        if self.two_sided:
            self._prepare_two_sided()
        else:
            self.fixed_events, self.pairs = self._event_terms()
            self.invariant_group = self._invariant_group(self.target, self.source,
                self.target_binary, self.target_colors) if self.pairs else self.group
        self.preparation_seconds = perf_counter() - started
        check_deadline(deadline)

    def __getstate__(self):
        state = dict(self.__dict__, deadline=None)
        for name in ('group', 'source_group', 'invariant_group'):
            group = state[name]
            state[name] = group.degree, tuple(tuple(g.array_form) for g in group.generators)
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        for name in ('group', 'source_group', 'invariant_group'):
            degree, generators = state[name]
            setattr(self, name, permutation_group(generators, degree))
        for array in (self.source, self.target, self.thresholds):
            array.setflags(write=False)

    def select(self, lower, upper=None):
        from .queries import EventSelection
        return EventSelection(self, lower, lower if upper is None else upper)

    def to_record(self):
        """Store decisions and costs; endpoint matrices stay with the reaction."""
        return dict(schema='mappa.fragment_events/1', fixed_events=int(self.fixed_events),
            pairs=[tuple(map(int, pair)) for pair in self.pairs], two_sided=self.two_sided,
            generators=[tuple(g.array_form) for g in self.group.generators],
            invariant_generators=[tuple(g.array_form) for g in self.invariant_group.generators],
            pair_orbit_states=self.pair_orbit_states)

    @classmethod
    def from_record(cls, fragment, record, *, event_problem=None, threshold=.5, metal_threshold=.3):
        if record['schema'] != 'mappa.fragment_events/1':
            raise ValueError('unsupported fragment event response')
        result = cls.__new__(cls)
        result.fragment, result.deadline = fragment, None
        positions = {atom:i for i,atom in enumerate(fragment.target_atoms)}
        mapping = dict(fragment.witness)
        result.witness = tuple(positions[mapping[a]] for a in fragment.source_atoms)
        objective = FragmentObjective(fragment, event_problem, threshold, metal_threshold)
        result.source, result.thresholds = objective.source, objective.thresholds
        result.target = objective.target[np.ix_(result.witness, result.witness)]
        result.source_binary = fragment.source_symmetry.weights > fragment.source_symmetry.bond_threshold
        result.target_binary = (fragment.target_symmetry.weights > fragment.target_symmetry.bond_threshold)[
            np.ix_(result.witness, result.witness)]
        n = len(result.witness)
        result.two_sided = record['two_sided']
        result.size = n * (2 if result.two_sided else 1)
        result.binary = result.target_binary
        if result.two_sided:
            result.binary = np.zeros((2*n,2*n), dtype=bool)
            result.binary[:n,:n], result.binary[n:,n:] = result.source_binary, result.target_binary
        result.edges = tuple(zip(*np.where(np.triu(result.target_binary,1))))
        result.source_group = permutation_group(fragment.source_symmetry.generators, n)
        result.group = permutation_group(record['generators'], result.size)
        result.invariant_group = permutation_group(record['invariant_generators'], result.size)
        result.fixed_events = int(record['fixed_events'])
        result.pairs = tuple(tuple(pair) for pair in record['pairs'])
        result.pair_orbit_states = record['pair_orbit_states']
        result.preparation_seconds = 0.
        return result

    def score(self, prefix, pairs):
        if self.two_sided:
            n = len(self.witness)
            return sum(abs(self.source[prefix[i], prefix[j]]
                           - self.target[prefix[k]-n, prefix[l]-n]) > self.thresholds[i, j]
                       for i, j, k, l in pairs)
        return sum(abs(self.source[i, j] - self.target[prefix[i], prefix[j]])
                   > self.thresholds[i, j] for i, j in pairs)

    def mapping(self, prefix):
        images, _ = self.correspondence(prefix)
        return {r: self.fragment.target_atoms[p] for r, p in zip(self.fragment.source_atoms, images)}

    def condition_scope(self, atom):
        position = self.fragment.source_atoms.index(atom)
        if not self.two_sided:
            return (position,)
        orbit = self.source_group.orbit(position)
        n = len(self.witness)
        return tuple(sorted(orbit)) + tuple(n + i for i in sorted(orbit))

    def query_result(self, prefix, count, seconds, work):
        result = self.result(prefix, count, True, seconds, work)
        choices = dict(discovery=dict(result.discovery_witness),
            source_permutation=result.source_operation, target_permutation=result.target_operation)
        return result.witness, choices

    @cached_property
    def group_order(self):
        return int(self.group.order())

    @cached_property
    def invariant_order(self):
        return int(self.invariant_group.order())

    def _event_terms(self):
        remaining = set(zip(*np.triu_indices(self.size, 1)))
        constant = 0
        variable = []
        # A whole structural edge/nonedge class can sometimes be settled immediately.
        for present in (False, True):
            pairs = sorted(pair for pair in remaining if bool(self.binary[pair]) == present)
            if not pairs:
                continue
            target_values = np.unique([self.target[pair] for pair in pairs])
            for pair in pairs:
                response = np.abs(self.source[pair] - target_values) > self.thresholds[pair]
                if np.all(response == response[0]):
                    constant += int(response[0])
                    remaining.remove(pair)
        # Orbit aggregation needs all pairs in the orbit, including any individually
        # constant ones, so process remaining orbits against their complete membership.
        for ordered in pair_orbits(self.group, self.size, remaining):
            check_deadline(self.deadline)
            self.pair_orbit_states += len(ordered)
            targets = np.unique([self.target[pair] for pair in ordered])
            signatures = {tuple(np.abs(self.source[pair] - targets) > self.thresholds[pair]) for pair in ordered}
            if len(signatures) == 1:
                # Already settled terms are counted only once.
                constant += sum(int(abs(self.source[pair] - self.target[pair]) > self.thresholds[pair])
                                for pair in ordered if pair in remaining)
            else:
                for pair in ordered:
                    if pair not in remaining:
                        continue
                    response = np.abs(self.source[pair] - targets) > self.thresholds[pair]
                    if np.all(response == response[0]):
                        constant += int(response[0])
                    else:
                        variable.append(pair)
        return constant, tuple(sorted(variable))

    def _invariant_group(self, weights, other, binary, colors):
        check_deadline(self.deadline)
        n = len(weights)
        values = np.unique(other[np.triu_indices(n, 1)])
        profiles = {(i, j): (bool(binary[i, j]),
                    tuple(np.abs(values - weights[i, j]) > self.thresholds[i, j]))
                    for i in range(n) for j in range(i + 1, n)}
        nonedge_profiles = {profile for pair, profile in profiles.items() if not binary[pair]}
        graph_edges = []
        colors = list(colors)
        for (i, j), profile in profiles.items():
            if not profile[0] and len(nonedge_profiles) <= 1:
                continue
            edge_vertex = len(colors)
            colors.append(('pair', profile))
            graph_edges.extend(((i, edge_vertex), (j, edge_vertex)))
        encoded = binary_group(len(colors), graph_edges, colors)
        check_deadline(self.deadline)
        return permutation_group([g.array_form[:n] for g in encoded.generators], n)

    def _prepare_two_sided(self):
        n = self.size
        source_binary = self.source_binary
        source_group = self.source_group
        target_group = self.group
        values, orbits = [], []
        for weights, group in ((self.source, source_group), (self.target, target_group)):
            lookup, reachable = {}, []
            for orbit in pair_orbits(group, n):
                check_deadline(self.deadline)
                self.pair_orbit_states += len(orbit)
                lookup.update((pair, len(reachable)) for pair in orbit)
                reachable.append(np.unique([weights[pair] for pair in orbit]))
            values.append(reachable)
            orbits.append(lookup)
        self.fixed_events, terms, responses = 0, [], {}
        for i, j in zip(*np.triu_indices(n, 1)):
            key = orbits[0][i, j], orbits[1][i, j], self.thresholds[i, j]
            if key not in responses:
                responses[key] = constant_response(values[0][key[0]], values[1][key[1]], key[2])
            response = responses[key]
            if response is not None:
                self.fixed_events += int(response)
            else:
                terms.append((i, j, n+i, n+j))
        self.pairs = tuple(terms)
        self.invariant_group = product_group(
            self._invariant_group(self.source, self.target, source_binary, self.colors),
            self._invariant_group(self.target, self.source, self.target_binary, self.target_colors)) if terms else None
        self.group = product_group(source_group, target_group)
        if not terms:
            self.invariant_group = self.group
        self.binary = np.zeros((2*n, 2*n), dtype=bool)
        self.binary[:n, :n], self.binary[n:, n:] = source_binary, self.target_binary
        self.size = 2*n

    @property
    def factor_graph(self):
        return self.binary

    def solve(self, *, max_nodes=10000, cache_entries=2048):
        if (max_nodes is not None and max_nodes < 0) or cache_entries < 0:
            raise ValueError('budgets must be nonnegative')
        return _Search(self, max_nodes, cache_entries).run()

    def correspondence(self, prefix):
        n = len(self.witness)
        if not self.two_sided:
            return tuple(self.witness[int(p)] for p in prefix), tuple(range(n))
        source_operation = tuple(map(int, prefix[:n]))
        mapping = [None] * n
        for i, atom in enumerate(source_operation):
            mapping[atom] = self.witness[int(prefix[n+i])-n]
        return tuple(mapping), source_operation

    def result(self, prefix, count, complete, seconds, work):
        images, source_operation = self.correspondence(prefix)
        target_operation = [0] * len(images)
        for i, p in enumerate(self.witness):
            target_operation[p] = images[source_operation[i]]
        mapping = self.fragment.relation.realize_permutations(
            source_permutation=source_operation, target_permutation=target_operation)
        aligned = np.argsort(self.witness)[list(images)]
        observed = np.count_nonzero(np.triu(
            np.abs(self.source - self.target[np.ix_(aligned, aligned)]) > self.thresholds, 1))
        assert observed == count, 'decomposition lost an event dependency'
        return FragmentMinimum(tuple(sorted(mapping.items())), count, complete,
            ('minimum established in witness symmetry' if complete else
             'time limit reached' if work.time_limit_reached else 'node cap reached'),
            seconds, work, source_operation, tuple(target_operation),
            discovery_witness=self.fragment.witness)


class _Search:
    def __init__(self, prepared, max_nodes, cache_entries):
        self.data = prepared
        self.max_nodes = max_nodes
        self.cache_entries = cache_entries
        self.work = SymmetryWork()
        self.answers = OrderedDict()
        self.groups = {}
        self.steps = {}
        self.restrictions = OrderedDict()

    def stopped(self):
        if self.data.deadline is not None and perf_counter() >= self.data.deadline:
            self.work.time_limit_reached = True
        if self.max_nodes is not None and self.work.nodes >= self.max_nodes:
            self.work.cap_reached = True
        return self.work.time_limit_reached or self.work.cap_reached

    def _score(self, prefix, pairs):
        self.work.score_evaluations += 1
        self.work.bond_comparisons += len(pairs)
        return self.data.score(prefix, pairs)

    def _metadata(self, group):
        key = group_key(group)
        if key not in self.groups:
            moved = support(group)
            factors = independent_factors(group, self.data.factor_graph)
            orbits = [orbit for orbit in group.orbits() if len(orbit) > 1]
            point = min((i for orbit in orbits for i in orbit),
                key=lambda i: (next(len(orbit) for orbit in orbits if i in orbit),
                               -np.count_nonzero(self.data.binary[i]) if self.data.binary is not None else 0,
                               i)) if moved else None
            self.groups[key] = moved, factors, point
        return key, self.groups[key]

    def _remember(self, key, answer):
        if answer.complete and self.cache_entries:
            self.answers[key] = answer
            self.answers.move_to_end(key)
            if len(self.answers) > self.cache_entries:
                self.answers.popitem(last=False)
            self.work.peak_cache_entries = max(self.work.peak_cache_entries, len(self.answers))
        return answer

    def _representatives(self, prefix, representatives, moved):
        if all(g.is_Identity for g in self.data.invariant_group.generators):
            return representatives
        outside = tuple(sorted(int(prefix[i]) for i in range(self.data.size) if i not in moved))
        fixed = sum(1 << atom for atom in outside)
        if fixed in self.restrictions:
            self.work.restriction_cache_hits += 1
            self.restrictions.move_to_end(fixed)
            invariant = self.restrictions[fixed]
        else:
            inherited = 0
            invariant = self.data.invariant_group
            for previous, restricted in self.restrictions.items():
                if (previous & fixed) == previous and previous.bit_count() > inherited.bit_count():
                    inherited, invariant = previous, restricted
            if inherited:
                self.work.restriction_parent_reuses += 1
            moving = set(support(invariant))
            additional = [atom for atom in outside if atom in moving and not (inherited & (1 << atom))]
            if additional:
                invariant = invariant.pointwise_stabilizer(additional, incremental=False)
            if self.cache_entries:
                self.restrictions[fixed] = invariant
                if len(self.restrictions) > self.cache_entries:
                    self.restrictions.popitem(last=False)
                self.work.peak_restriction_cache_entries = max(
                    self.work.peak_restriction_cache_entries, len(self.restrictions))
        # Every search group is the full binary symmetry fixing its outside
        # atoms. The event-invariant subgroup, restricted to those target images,
        # therefore already satisfies the current search constraints.
        visited = set()
        retained = []
        for image, operation in representatives:
            target = int(prefix[image])
            if target in visited:
                self.work.equivalent_positions_skipped += 1
                continue
            retained.append((image, operation))
            visited.update(invariant.orbit(target))
        return retained

    def _solve(self, group, prefix, pairs):
        if not pairs:
            return _Answer(0, (), True)
        initial = int(self._score(prefix, pairs))
        if initial == 0:
            return _Answer(0, (), True)
        if self.stopped():
            return _Answer(initial,
                tuple((i, int(prefix[i])) for i in support(group)), False)
        key, (moved, factors, point) = self._metadata(group)
        relevant = sorted(set(moved).union(atom for pair in pairs for atom in pair))
        answer_key = key, pairs, tuple(int(prefix[i]) for i in relevant)
        if answer_key in self.answers:
            self.work.cache_hits += 1
            self.answers.move_to_end(answer_key)
            return self.answers[answer_key]
        best = _Answer(initial, tuple((i, int(prefix[i])) for i in moved), initial == 0)
        if initial == 0 or not moved:
            return self._remember(answer_key, _Answer(initial, best.images, True))
        if self.stopped():
            return best
        self.work.nodes += 1
        moved_set = set(moved)
        active = tuple(pair for pair in pairs if moved_set.intersection(pair))
        fixed = initial - int(self._score(prefix, active))
        if not active:
            return self._remember(answer_key, _Answer(fixed, best.images, True))
        components = event_components(factors, active)
        self.work.event_coupling_merges += len(factors) - len(components)
        if len(components) > 1:
            self.work.independent_splits += 1
            images = dict(best.images)
            cost, complete = fixed, True
            for component in components:
                atoms = set(support(component))
                local = tuple(pair for pair in active if atoms.intersection(pair))
                answer = self._solve(component, prefix, local)
                images.update(answer.images)
                cost += answer.cost
                complete = complete and answer.complete
            return self._remember(answer_key, _Answer(cost, tuple(sorted(images.items())), complete))
        affected = moved_set.intersection(atom for pair in active for atom in pair)
        if affected and point not in affected:
            point = min(affected)
        step_key = key, point
        if step_key not in self.steps:
            representatives = group.orbit_transversal(point, pairs=True)
            representatives.sort(key=lambda item: (item[0] != point, item[0]))
            self.steps[step_key] = group.stabilizer(point), representatives
            self.work.group_steps_prepared += 1
        subgroup, representatives = self.steps[step_key]
        complete = True
        for _, operation in self._representatives(prefix, representatives, moved):
            next_prefix = prefix[np.asarray(operation.array_form[:self.data.size])]
            answer = self._solve(subgroup, next_prefix, active)
            complete = complete and answer.complete
            cost = fixed + answer.cost
            if cost < best.cost:
                image = next_prefix.copy()
                for i, p in answer.images:
                    image[i] = p
                best = _Answer(cost, tuple((i, int(image[i])) for i in moved), False)
            if best.cost == fixed:
                return self._remember(answer_key, _Answer(best.cost, best.images, True))
            if self.stopped() and not answer.complete:
                complete = False
                break
        return self._remember(answer_key, _Answer(best.cost, best.images, complete))

    def run(self):
        started = perf_counter()
        prefix = np.arange(self.data.size)
        answer = self._solve(self.data.group, prefix, self.data.pairs)
        for i, p in answer.images:
            prefix[i] = p
        count = self.data.fixed_events + answer.cost
        return self.data.result(prefix, count, answer.complete, perf_counter() - started, self.work)


def fragment_symmetry_minimum(fragment, *, event_problem=None, threshold=.5,
                              metal_threshold=.3, max_nodes=10000, cache_entries=2048):
    """Return one minimum-event witness in the actual Fragment's saved orbit.

    ``complete`` certifies this orbit only. The original Fragment retains its
    other discovery alternatives. A cap returns an actual witness and explicitly
    incomplete minimization; it never selects a different algorithm.
    """
    return FragmentSymmetryEvents(fragment, event_problem=event_problem,
        threshold=threshold, metal_threshold=metal_threshold).solve(
            max_nodes=max_nodes, cache_entries=cache_entries)
