"""Ordered binary symmetry actions with shared, sparse original-weight costs.

Each action retains its own choice. Their disjoint permutation frames are a
direct product; applying them still follows the saved order. Taking the closure
of all target generators together would instead invent additional mappings.
"""
from dataclasses import dataclass
from functools import cached_property
from time import perf_counter

import numpy as np

from .groups import permutation_group
from .minimum import _Search, SymmetryWork
from .queries import EventSelection


@dataclass(frozen=True)
class EventTerm:
    atoms: tuple
    weight: float
    threshold: float
    scope: tuple
    excluded_pairs: frozenset | None = None

    def __iter__(self):
        return iter(self.scope)


@dataclass(frozen=True)
class BoundEventTerm:
    """One prepared cost, referencing the same choices in a larger program."""
    response: object
    term: object
    positions: tuple
    scope: tuple

    def __iter__(self):
        return iter(self.scope)


@dataclass(frozen=True)
class ProgramMinimum:
    witness: tuple
    choices: tuple
    event_count: int
    lower_bound: int
    complete: bool
    search_seconds: float
    work: SymmetryWork


class OrderedSymmetryEvents:
    """An event projection of shared fragment responses and placement actions.

    ``actions`` contains (saved action, prepared fragment response or None).
    Only variable pair terms enter search. No discovery variables or complete
    atom-assignment tables are reconstructed here.
    """
    def __init__(self, mapping, actions, target_weights, terms,
                 *, fixed_events=0, sweep=None, projections=()):
        self.base = tuple(mapping)
        self.actions = tuple(actions)
        self.target = target_weights
        self.deadline = None
        self.binary = None
        frames, size = [], 0
        for (kind, data), prepared in actions:
            if kind == 'relation':
                frame = (), prepared.size
            elif kind == 'group':
                atoms = tuple(sorted({a for g in data for a, b in enumerate(g) if a != b}))
                frame = atoms, len(atoms)
            else:
                raise ValueError(f'Unsupported symmetry action: {kind}')
            frames.append((size, *frame))
            size += frame[1]
        self.frames, self.size = tuple(frames), size
        scopes = self._endpoint_scopes()
        self.scopes = scopes
        self.fixed_events = fixed_events
        pairs = []
        for atoms, weight, threshold, constant in terms:
            if constant is not None:
                self.fixed_events += int(bool(constant))
            else:
                scope = tuple(sorted(set(scopes[atoms[0]]) | set(scopes[atoms[1]])))
                pairs.append(EventTerm(atoms, weight, threshold, scope))
        self.projections = []
        for response, binding in projections:
            local_frames = (response.frames if isinstance(response, OrderedSymmetryEvents)
                            else ((0, (), response.size),))
            positions = tuple(self.frames[outer][0] + i for outer, (_, _, degree)
                              in zip(binding, local_frames) for i in range(degree))
            self.fixed_events += response.fixed_events
            self.projections.append((response, positions))
            pairs.extend(BoundEventTerm(response, term, positions,
                                       tuple(positions[i] for i in term)) for term in response.pairs)
        if sweep is not None:
            atoms, weight, threshold, reserve, boundary_pairs = sweep
            scope = tuple(sorted(set(scopes[atoms[0]]) | set(scopes[atoms[1]])))
            pairs.append(EventTerm(atoms, weight, threshold, scope))
            if reserve:
                pairs.append(EventTerm(atoms, 0., 0., scope, boundary_pairs))
        self.pairs = tuple(pairs)

    def _endpoint_scopes(self):
        """Trace just the positions that can determine each endpoint image."""
        from .actions import _action_orbits
        rows = []
        for ((kind, data), prepared), (start, atoms, _) in zip(self.actions, self.frames):
            if kind == 'group':
                orbits = _action_orbits((kind, data), len(self.target))
                rows.append({atom: ((start + i,), orbits[atom]) for i, atom in enumerate(atoms)})
            else:
                n = len(prepared.witness)
                row = {}
                for atom, image in data.witness:
                    i = prepared.fragment.source_atoms.index(atom)
                    origins = prepared.source_group.orbit(i) if prepared.two_sided else (i,)
                    offset = n if prepared.two_sided else 0
                    reached = {p - offset for j in origins for p in prepared.group.orbit(offset + j)}
                    targets = {prepared.fragment.target_atoms[prepared.witness[p]] for p in reached}
                    scope = tuple(start + p for p in prepared.condition_scope(atom))
                    row[image] = scope, targets
                rows.append(row)
        scopes = {}
        for atom, image in self.base:
            domain, scope = {image}, set()
            for row in rows:
                reached = set()
                for p in domain:
                    positions, images = row.get(p, ((), (p,)))
                    scope.update(positions)
                    reached.update(images)
                domain = reached
            scopes[atom] = tuple(sorted(scope))
        return scopes

    @cached_property
    def group(self):
        generators = []
        for ((kind, data), prepared), (start, atoms, degree) in zip(self.actions, self.frames):
            if kind == 'relation':
                local = (g.array_form for g in prepared.group.generators)
            else:
                positions = {atom: i for i, atom in enumerate(atoms)}
                local = ([positions[g[atom]] for atom in atoms] for g in data)
            for images in local:
                operation = list(range(self.size))
                operation[start:start + degree] = [start + int(p) for p in images]
                generators.append(operation)
        return permutation_group(generators, self.size)

    @cached_property
    def invariant_group(self):
        # Reuse a local certificate only after fixing every possible image
        # read by other costs. This stays valid after any earlier choice.
        generators = []
        orbits = {i: orbit for orbit in self.group.orbits() for i in orbit}
        for response, positions in self.projections:
            other_scope = {i for term in self.pairs if not (
                isinstance(term, BoundEventTerm) and term.response is response
                and term.positions == positions) for i in term}
            protected = set().union(*(orbits[i] for i in other_scope))
            fixed = [i for i, position in enumerate(positions) if position in protected]
            invariant = response.invariant_group
            if fixed:
                invariant = invariant.pointwise_stabilizer(fixed, incremental=False)
            for g in invariant.generators:
                if g.is_Identity:
                    continue
                operation = list(range(self.size))
                for i, position in enumerate(positions):
                    operation[position] = positions[g(i)]
                generators.append(operation)
        return permutation_group(generators, self.size)

    def __getstate__(self):
        return {k:v for k,v in self.__dict__.items() if k not in ('group','invariant_group','deadline')}

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.deadline = None

    def _apply(self, prefix, *, proofs=False):
        images = np.array([p for _, p in self.base])
        choices = []
        for ((kind, data), prepared), (start, atoms, degree) in zip(self.actions, self.frames):
            operation = np.asarray(prefix[start:start + degree]) - start
            transport = np.arange(len(self.target))
            if kind == 'relation':
                selected = prepared.mapping(operation)
                for atom, image in data.witness:
                    transport[image] = selected[atom]
                if proofs:
                    count = prepared.fixed_events + int(prepared.score(operation, prepared.pairs))
                    _, choice = prepared.query_result(operation, count, 0., SymmetryWork())
            else:
                transport[list(atoms)] = [atoms[p] for p in operation]
                if proofs:
                    choice = dict(permutation=tuple(map(int, transport)))
            images = transport[images]
            if proofs:
                choices.append(choice)
        return dict(zip((a for a, _ in self.base), map(int, images))), tuple(choices)

    def mapping(self, prefix):
        return self._apply(prefix)[0]

    def score(self, prefix, terms):
        if not terms:
            return 0
        raw, projected = [], {}
        for term in terms:
            if isinstance(term, BoundEventTerm):
                projected.setdefault((term.response, term.positions), []).append(term.term)
            else:
                raw.append(term)
        count = self.count_terms(self.mapping(prefix), raw) if raw else 0
        for (response, positions), costs in projected.items():
            inverse = {position: i for i, position in enumerate(positions)}
            local = np.asarray([inverse[int(prefix[i])] for i in positions])
            count += response.score(local, costs)
        return int(count)

    def count_terms(self, images, terms):
        count = 0
        for term in terms:
            pair = tuple(sorted(images[a] for a in term.atoms))
            count += (pair not in term.excluded_pairs if term.excluded_pairs is not None else
                      abs(term.weight - self.target[pair]) > term.threshold)
        return int(count)

    def condition_scope(self, atom):
        return self.scopes[atom]

    def query_result(self, prefix, count, seconds, work):
        mapping, choices = self._apply(prefix, proofs=True)
        return tuple(sorted(mapping.items())), choices

    def result(self, prefix, count, complete, seconds, work):
        witness, choices = self.query_result(prefix, count, seconds, work)
        return ProgramMinimum(witness, choices, count, count if complete else self.fixed_events,
                              complete, seconds, work)

    def solve(self, *, seconds=5., max_nodes=10000, cache_entries=2048):
        started = perf_counter()
        self.deadline = started + seconds
        try:
            if not self.pairs:
                return self.result(np.arange(self.size), self.fixed_events, True,
                                   perf_counter() - started, SymmetryWork())
            return _Search(self, max_nodes, cache_entries).run()
        finally:
            self.deadline = None

    def select(self, lower, upper=None):
        return EventSelection(self, lower, lower if upper is None else upper)
