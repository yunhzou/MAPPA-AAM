"""Shared fragment event counts and endpoint-conditioned sweep responses.

Counts are contracted from the original joint choices. The window consumes
these responses; it never constructs feasibility queries or probes event totals.
"""
from collections import defaultdict
from functools import cached_property

import numpy as np

from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints, compile_relaxed_relation
from mappa.mappa_core.core.symmetry.event_counts import EventCountProgram
from mappa.mappa_core.utils.bonds import bond_event_threshold


def count_bounds(index, response):
    """Independent-term lower bound and the saved valid witness's count."""
    mapping = dict(response.mapping)
    lower = sum(bool(np.all(term.response.categories != 0)) for term in response.terms)
    upper = sum(abs(float(index.p[mapping[a], mapping[b]]) - term.response.source_weight)
                > term.response.threshold for term in response.terms for a, b in (term.atoms,))
    return int(lower), int(upper)


class FragmentCountIndex:
    """Reaction-owned certificates and weight bags, shared across fragment uses."""
    def __init__(self, index):
        self.index = index
        self.groups = {}
        self.bags = {}

    def preserves_weights(self, action):
        if action not in self.groups:
            kind, generators = action
            elements = np.asarray(self.index.problem.product.elements)
            self.groups[action] = kind == 'group' and all(
                np.array_equal(elements, elements[np.asarray(g)]) and
                np.array_equal(self.index.p, self.index.p[np.ix_(g, g)]) for g in generators)
        return self.groups[action]

    def fixed_pair_bag(self, response):
        """Set permutations and subsequent weighted placements preserve pair bags."""
        selected = frozenset(p for _, p in response.mapping)
        groups_started = False
        for action in response.actions:
            kind, value = action
            if kind == 'group':
                if not self.preserves_weights(action):
                    return False
                groups_started = True
            elif kind == 'relation':
                support = frozenset(value.target_atoms)
                if groups_started or (selected & support and not support <= selected):
                    return False
            else:
                raise ValueError(f'Unknown fragment action: {kind}')
        return True

    def bag(self, side, atoms, anchor=None):
        key = side, atoms, anchor
        if key not in self.bags:
            endpoint = getattr(self.index.problem, side)
            grouped = defaultdict(list)
            for i, a in enumerate(atoms):
                for b in atoms[i + 1:]:
                    elements = tuple(sorted((endpoint.elements[a], endpoint.elements[b])))
                    grouped[elements, anchor in (a, b)].append(float(endpoint.wbo[a, b]))
            self.bags[key] = {elements: np.sort(values) for elements, values in grouped.items()}
        return self.bags[key]

    def distance(self, source, target):
        representatives = {element: a for a, element in enumerate(self.index.elements)}
        lower = 0
        for key, xs in source.items():
            elements, _ = key
            ys = target[key]
            a, b = (representatives[e] for e in elements)
            threshold = bond_event_threshold(self.index.elements, a, b,
                default_threshold=self.index.threshold, metal_threshold=self.index.metal_threshold)
            i = j = matches = 0
            while i < len(xs) and j < len(ys):
                delta = float(xs[i]) - float(ys[j])
                if delta > threshold:
                    j += 1
                elif delta < -threshold:
                    i += 1
                else:
                    i += 1
                    j += 1
                    matches += 1
            lower += len(xs) - matches
        return lower

    def aggregate_bound(self, response):
        if not self.fixed_pair_bag(response):
            return 0
        source_atoms = tuple(a for a, _ in response.mapping)
        target_atoms = tuple(sorted(p for _, p in response.mapping))
        return self.distance(self.bag('reactant', source_atoms), self.bag('product', target_atoms))

    def endpoint_bounds(self, response, atom, lower):
        """Keep incident and nonincident pair costs separate for a sweep endpoint."""
        if not self.fixed_pair_bag(response):
            return {p: lower for p in response.atom_domains[atom]}
        source_atoms = tuple(a for a, _ in response.mapping)
        target_atoms = tuple(sorted(p for _, p in response.mapping))
        source = self.bag('reactant', source_atoms, atom)
        target_elements = self.index.problem.product.elements
        bounds = {p: max(lower, self.distance(source, self.bag('product', target_atoms, p)))
                  for p in target_atoms if target_elements[p] == self.index.elements[atom]}
        for action in response.actions:
            if action[0] == 'group':
                domains = self.index.action_domains(action)
                advanced = {}
                for p, count in bounds.items():
                    for q in domains[p]:
                        advanced[q] = min(advanced.get(q, count), count)
                bounds = advanced
        return {p: bounds[p] for p in response.atom_domains[atom] if p in bounds}


class FragmentEventCounts:
    """One local correlated response, reused for counts and sweep endpoints."""
    def __init__(self, owner, response):
        self.owner, self.response = owner, response
        self.bounds = count_bounds(owner.index, response)
        self.lower = self.bounds[0]
        if self.lower != self.bounds[1]:
            self.lower = max(self.lower, owner.aggregate_bound(response))
        self._endpoints = {}
        self._endpoint_bounds = {}

    @cached_property
    def minimum(self):
        if self.lower == self.bounds[1]:
            return self.lower
        totals = self.project(())
        if not totals:
            raise ValueError('Saved fragment relationship has no attainable event count')
        return totals[()]

    @cached_property
    def program(self):
        response, index = self.response, self.owner.index
        system = FiniteConstraints()
        current = {a: system.constant(p) for a, p in response.mapping}
        for position, (kind, data) in enumerate(response.actions):
            if kind == 'relation':
                compiled = compile_relaxed_relation(system, data.relation, position, images_only=True)
                transport = {p: compiled.images[r] for r, p in data.witness}
                current = {a: system.lookup(value, transport) for a, value in current.items()}
            elif kind == 'group':
                current = dict(zip(current, system.act(list(current.values()), data, position)))
            else:
                raise ValueError(f'Unknown fragment action: {kind}')
        system.propagate()
        costs = []
        for term in response.terms:
            a, b = term.atoms
            left, right = current[a][0], current[b][0]
            rows = {(x, y): int(abs(float(index.p[x, y]) - term.response.source_weight)
                                > term.response.threshold)
                    for x in system.domain(left) for y in system.domain(right)}
            costs.append((left, right, rows))
        return EventCountProgram(system, costs), {a: value[0] for a, value in current.items()}

    def project(self, atoms):
        program, images = self.program
        return program.project(tuple(images[a] for a in atoms))

    def endpoint(self, atom):
        if atom not in self._endpoints:
            self._endpoints[atom] = {images[0]: count for images, count in self.project((atom,)).items()}
        return self._endpoints[atom]

    def endpoint_bounds(self, atom):
        if atom not in self._endpoint_bounds:
            self._endpoint_bounds[atom] = self.owner.endpoint_bounds(self.response, atom, self.lower)
        return self._endpoint_bounds[atom]

    def witness_endpoints(self, atom):
        """Attained counts from the saved choice and its exact placement orbits."""
        reached = {dict(self.response.mapping)[atom]}
        for action in self.response.actions:
            if action[0] != 'group' or not self.owner.preserves_weights(action):
                continue
            domains = self.owner.index.action_domains(action)
            reached = {q for p in reached for q in domains[p]}
        return {image: self.bounds[1] for image in reached}
