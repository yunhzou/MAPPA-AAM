"""Shared binary symmetry responses and their endpoint-conditioned queries."""
from collections import defaultdict
from dataclasses import replace
from functools import cached_property

import numpy as np

from mappa.mappa_core.core.symmetry.program import OrderedSymmetryEvents
from mappa.mappa_core.core.symmetry.minimum import SymmetryWork
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
    def __init__(self, index, decoder):
        self.index = index
        self.decoder = decoder
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
    """An event projection shared by all uses of one contextual fragment.

    A simple completed fragment uses its prepared response directly. Historical
    overlapping actions retain their order in the same core symmetry search.
    Neither representation builds a full attainable-count table.
    """
    def __init__(self, owner, response):
        self.owner, self.response = owner, response
        self.lower = max(count_bounds(owner.index, response)[0], owner.aggregate_bound(response))
        self._endpoints = {}

    @cached_property
    def direct(self):
        actions = self.response.actions
        if len(actions) != 1 or actions[0][0] != 'relation':
            return False
        fragment = actions[0][1]
        return (self.response.mapping == fragment.witness
                and self.response.piece.source_atoms == fragment.source_atoms)

    @cached_property
    def program(self):
        owner, response = self.owner, self.response
        if self.direct:
            return owner.decoder.prepare(response.actions[0][1]).response
        actions = tuple((action, owner.decoder.prepare(action[1]).response if action[0] == 'relation' else None)
                        for action in response.actions)
        # A whole fragment followed only by weight-preserving placements has
        # exactly its prepared internal costs. Placements still move endpoints.
        if actions and actions[0][0][0] == 'relation':
            fragment = actions[0][0][1]
            if (response.mapping == fragment.witness
                    and response.piece.source_atoms == fragment.source_atoms
                    and all(owner.preserves_weights(action) for action in response.actions[1:])):
                return OrderedSymmetryEvents(response.mapping, actions, owner.index.p, (),
                    projections=((actions[0][1], (0,)),))
        terms = tuple((term.atoms, term.response.source_weight, term.response.threshold,
                       term.response.constant) for term in response.terms)
        return OrderedSymmetryEvents(response.mapping, actions, owner.index.p, terms)

    @cached_property
    def bounds(self):
        data = self.program
        prefix = np.arange(data.size)
        return self.lower, data.fixed_events + int(data.score(prefix, data.pairs))

    @cached_property
    def minimum_result(self):
        if self.direct:
            prepared = self.owner.decoder.prepare(self.response.actions[0][1])
            return prepared.optimization
        lower, attained = self.bounds
        if lower == attained:
            return self.program.result(np.arange(self.program.size), attained, True, 0., SymmetryWork())
        return self.program.solve(seconds=self.owner.decoder.fragment_seconds)

    @property
    def minimum(self):
        return self.minimum_result.event_count

    @property
    def minimum_proven(self):
        return self.minimum_result.complete

    @property
    def minimum_choices(self):
        result = self.minimum_result
        if self.direct:
            return (dict(discovery=dict(result.discovery_witness),
                source_permutation=result.source_operation, target_permutation=result.target_operation),)
        return result.choices

    @property
    def constant(self):
        return self.program.fixed_events if not self.program.pairs else None

    def find(self, lower, upper=None, *, assignments=None, seconds=None):
        answer = self.program.select(lower, upper).find(assignments=assignments,
            seconds=self.owner.decoder.fragment_seconds if seconds is None else seconds)
        if self.direct and answer.status == 'found':
            answer = replace(answer, choices=(answer.choices,))
        return answer

    def endpoint_bounds(self, atom):
        return self.owner.endpoint_bounds(self.response, atom, self.lower)

    def witness_endpoints(self, atom):
        answer = self.minimum_result
        return {dict(answer.witness)[atom]: answer.event_count}
