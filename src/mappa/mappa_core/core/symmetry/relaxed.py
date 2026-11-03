"""Broad structural shuffles constrained against original internal bond weights.

The permitted shuffles need not form a group. Structural generators describe
the broad space; shared bond responses constrain their final joint action.
"""
from collections import defaultdict
from dataclasses import dataclass, field
from functools import cached_property
import math
from numbers import Integral

import numpy as np

from .actions import apply_action
from .responses import BondResponse, BondResponseIndex


@dataclass(frozen=True)
class BondConstraint:
    """Pairs sharing a source weight and one correlated pair-orbit response."""
    pairs: tuple[tuple[int, int], ...]
    response: BondResponse


@dataclass(frozen=True, eq=False)
class RelaxedAutomorphism:
    """A global or anchored fragment self-relationship, in a local atom frame.

    ``atoms`` maps local positions to original endpoint indices. Generators use
    local positions. Every constraint evaluates the same final permutation, so
    separately permitted shuffles cannot silently become unrestricted compositions.
    """
    atoms: tuple[int, ...]
    elements: tuple[str, ...]
    anchors: tuple[int, ...]
    bond_threshold: float
    tolerance: float
    generators: tuple[tuple[int, ...], ...]
    constraints: tuple[BondConstraint, ...]
    weights: np.ndarray = field(repr=False)
    edges: tuple[tuple[int, int], ...] = field(repr=False)
    projected_pair_states: int = field(repr=False)

    @cached_property
    def _evaluation(self):
        return tuple((np.asarray(term.pairs, dtype=int), term.response)
                     for term in self.constraints)

    def _accepts(self, permutation):
        for pairs, response in self._evaluation:
            images = permutation[pairs]
            weights = self.weights[images[:, 0], images[:, 1]]
            if np.any(response.lookup(weights) != 0):
                return False
        return True

    def _permutation(self, word):
        return apply_action(np.arange(len(self.atoms)),
                            ('group', self.generators), word, len(self.atoms))

    def allows(self, word=()):
        """Check a signed generator word against the ORIGINAL self-reference."""
        return self._accepts(self._permutation(word))

    def realize(self, word=()):
        """Return one explicitly requested, permitted original-index mapping."""
        permutation = self._permutation(word)
        if not self._accepts(permutation):
            raise ValueError('rearrangement violates relaxed self tolerance')
        return {atom: self.atoms[int(image)]
                for atom, image in zip(self.atoms, permutation)}

    def contains(self, permutation):
        """Check a supplied local-index permutation without group enumeration."""
        permutation = np.asarray(permutation)
        n = len(self.atoms)
        if permutation.shape != (n,) or (n and permutation.dtype.kind not in 'iu'):
            return False
        permutation = permutation.astype(int, copy=False)
        if not np.array_equal(np.sort(permutation), np.arange(n)):
            return False
        elements = np.asarray(self.elements)
        if not np.array_equal(elements[permutation], elements):
            return False
        positions = {atom: i for i, atom in enumerate(self.atoms)}
        if any(permutation[positions[atom]] != positions[atom] for atom in self.anchors):
            return False
        if self.edges:
            pairs = permutation[np.asarray(self.edges, dtype=int)]
            if np.any(self.weights[pairs[:, 0], pairs[:, 1]] <= self.bond_threshold):
                return False
        # A bijection preserving every edge preserves the entire binary graph.
        return self._accepts(permutation)

    @property
    def statistics(self):
        return dict(atoms=len(self.atoms), generators=len(self.generators),
                    response_terms=len(self.constraints),
                    constrained_pairs=sum(len(term.pairs) for term in self.constraints),
                    projected_pair_states=self.projected_pair_states)


def _constraints(weights, tolerance, generators):
    """Prepare one response per weight and pair orbit, never full shuffles.

Nonnegative weights can differ by more than tolerance only when at least one
exceeds it. Closing those pair orbits also includes their possible lower-weight
partners, including subthreshold contacts. Each orbit is visited only once.
"""
    active = set(zip(*np.nonzero(np.triu(weights > tolerance, 1))))
    if not generators or not active:
        return (), 0
    index = BondResponseIndex(weights)
    actions = (('group', generators),)
    constraints, projected = [], 0
    for pair in sorted(active):
        if pair not in active:
            continue
        reached = index.endpoint_pairs(tuple(map(int, pair)), actions)
        projected += len(reached)
        orbit = sorted({(min(a, b), max(a, b)) for a, b in reached})
        active.difference_update(orbit)
        pairs = np.asarray(orbit, dtype=int)
        values = weights[pairs[:, 0], pairs[:, 1]]
        distinct = np.unique(values)
        if distinct[-1] - distinct[0] <= tolerance:
            continue
        by_weight = defaultdict(list)
        for atoms, weight in zip(orbit, values):
            by_weight[float(weight)].append(atoms)
        for weight, atoms in by_weight.items():
            response = index.response_from_weights(weight, tolerance, distinct)
            if response.constant != 0:
                constraints.append(BondConstraint(tuple(atoms), response))
    return tuple(constraints), projected


def relaxed_automorphism(elements, weights, *, tolerance, bond_threshold=0.2,
                         atoms=None, anchors=()):
    """Construct a relaxed global or anchor-aware fragment self-automorph.

Omit ``atoms`` for the whole endpoint; select original atom indices for a
fragment. Only internal weights participate. The broad binary graph has an edge
iff its weight is strictly above ``bond_threshold``. The final rearrangement
must preserve every internal pair's weight within inclusive ``tolerance``.
Neither tolerance is inferred from source-to-target discovery settings.
"""
    elements = tuple(elements)
    weights = np.asarray(weights, dtype=float)
    n = len(elements)
    if (weights.shape != (n, n) or not np.all(np.isfinite(weights))
            or np.any(weights < 0) or not np.array_equal(weights, weights.T)):
        raise ValueError('weights must be a finite nonnegative symmetric endpoint matrix')
    if any(not math.isfinite(value) or value < 0 for value in (tolerance, bond_threshold)):
        raise ValueError('tolerance and bond threshold must be finite and nonnegative')
    atoms = tuple(range(n)) if atoms is None else tuple(atoms)
    anchors = tuple(anchors)
    if (any(not isinstance(a, Integral) or not 0 <= a < n for a in atoms)
            or len(set(atoms)) != len(atoms)):
        raise ValueError('fragment atoms must be distinct endpoint indices')
    if (any(not isinstance(a, Integral) or a not in atoms for a in anchors)
            or len(set(anchors)) != len(anchors)):
        raise ValueError('anchors must be distinct selected atom indices')
    atoms, anchors = tuple(sorted(map(int, atoms))), tuple(sorted(map(int, anchors)))
    local = weights[np.ix_(atoms, atoms)].copy()
    local.setflags(write=False)
    labels = tuple(elements[a] for a in atoms)
    edges = tuple((int(a), int(b)) for a, b in
                  zip(*np.nonzero(np.triu(local > bond_threshold, 1))))
    colors = defaultdict(list)
    for position, (atom, element) in enumerate(zip(atoms, labels)):
        colors[element, atom if atom in anchors else None].append(position)
    if atoms:
        from mappa._engine import AutGraph
        generators = tuple(tuple(g) for g in
                           AutGraph(len(atoms), edges).generators(list(colors.values())))
    else:
        generators = ()
    constraints, projected = _constraints(local, float(tolerance), generators)
    return RelaxedAutomorphism(atoms, labels, anchors, float(bond_threshold),
        float(tolerance), generators, constraints, local, edges, projected)
