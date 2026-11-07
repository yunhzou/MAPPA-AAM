"""Element- and anchor-preserving symmetry of a binary molecular graph."""
from collections import defaultdict
from dataclasses import dataclass, field, fields
from functools import cached_property, total_ordering
import math
from numbers import Integral

import numpy as np

from .actions import apply_action


@total_ordering
class _MatrixKey:
    """Exact byte comparison without retaining a second immutable matrix."""
    __slots__ = ('data', '_hash')

    def __init__(self, weights):
        self.data = weights.view(np.uint8).reshape(-1)
        self._hash = hash(self.data.tobytes())

    def __hash__(self):
        return self._hash

    def __reduce__(self):
        return type(self), (self.data,)

    def __eq__(self, other):
        if not isinstance(other, _MatrixKey):
            return NotImplemented
        return self is other or (self._hash == other._hash
                                 and np.array_equal(self.data, other.data))

    def __lt__(self, other):
        if not isinstance(other, _MatrixKey):
            return NotImplemented
        return self.data.tobytes() < other.data.tobytes()


def structural_generators(size, edges, colors):
    """Use the same nauty owner as the native matcher; retain generators only."""
    if not size:
        return ()
    from mappa._engine import AutGraph
    cells = defaultdict(list)
    for atom, color in enumerate(colors):
        cells[color].append(atom)
    return tuple(tuple(g) for g in AutGraph(size, list(edges)).generators(list(cells.values())))


@dataclass(frozen=True, eq=False)
class BinarySymmetry:
    """Compact correlated shuffles in a local atom frame.

    ``atoms`` relates local positions to endpoint indices. Original weights are
    shared with discovery and event counting; they do not filter self shuffles.
    """
    atoms: tuple[int, ...]
    elements: tuple[str, ...]
    anchors: tuple[int, ...]
    bond_threshold: float
    _generators: tuple[tuple[int, ...], ...] | None
    weights: np.ndarray = field(repr=False)
    edges: tuple[tuple[int, int], ...] = field(repr=False)

    @cached_property
    def key(self):
        return (self.atoms, self.elements, self.anchors, self.bond_threshold,
                _MatrixKey(self.weights))

    @property
    def generators(self):
        if self._generators is None:
            colors = [(element, atom if atom in self.anchors else None)
                      for atom, element in zip(self.atoms, self.elements)]
            object.__setattr__(self, '_generators',
                structural_generators(len(self.atoms), self.edges, colors))
        return self._generators

    def share_weights(self, pool):
        """Intern identical storage without changing this record or its keys."""
        key = self.key[-1]
        weights = pool.setdefault((self.weights.shape, key), self.weights)
        if weights is not self.weights:
            object.__setattr__(self, 'weights', weights)
            # Existing fragment keys borrow this key. Repoint its equal bytes
            # too, so those keys do not keep the discarded allocation alive.
            key.data = weights.view(np.uint8).reshape(-1)

    def __reduce__(self):
        return type(self)._from_binary, tuple(getattr(self, item.name) for item in fields(self))

    @classmethod
    def _from_binary(cls, *values):
        result = cls(*values)
        result.weights.setflags(write=False)
        return result

    def to_record(self, *, pool=None, array_lists=True):
        if pool is not None and self.key in pool:
            return pool[self.key]
        record = dict(atoms=self.atoms, elements=self.elements, anchors=self.anchors,
                      bond_threshold=self.bond_threshold, generators=self.generators,
                      weights=self.weights.tolist() if array_lists else self.weights)
        if pool is not None:
            pool[self.key] = record
        return record

    @classmethod
    def from_record(cls, record, *, pool=None):
        """Restore saved preparation without repeating graph matching."""
        if 'tolerance' in record:
            raise ValueError('historical self-tolerance records require their producing revision')
        weights = np.asarray(record['weights'], dtype=float)
        weights.setflags(write=False)
        generators = tuple(tuple(g) for g in record['generators'])
        threshold = float(record['bond_threshold'])
        key = (tuple(record['atoms']), tuple(record['elements']), tuple(record['anchors']),
               threshold, _MatrixKey(weights))
        if pool is not None and key in pool:
            return pool[key]
        edges = tuple((int(a), int(b)) for a, b in zip(*np.nonzero(np.triu(weights > threshold, 1))))
        result = cls(*key[:4], generators, weights, edges)
        if pool is not None:
            pool[key] = result
        return result

    def _permutation(self, word):
        return apply_action(np.arange(len(self.atoms)),
                            ('group', self.generators), word, len(self.atoms))

    def realize(self, word=()):
        """Apply one requested generator word, without expanding its group."""
        return {atom: self.atoms[int(image)]
                for atom, image in zip(self.atoms, self._permutation(word))}

    def contains(self, permutation):
        """Check one permutation against connectivity, elements and anchors."""
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
        # A bijection preserving every edge also preserves all nonedges.
        return True


def binary_symmetry(elements, weights, *, bond_threshold=0.2, atoms=None, anchors=()):
    """Prepare whole-endpoint or fragment-local symmetry with fixed anchors."""
    elements = tuple(elements)
    weights = np.asarray(weights, dtype=float)
    n = len(elements)
    if (weights.shape != (n, n) or not np.all(np.isfinite(weights))
            or np.any(weights < 0) or not np.array_equal(weights, weights.T)):
        raise ValueError('weights must be a finite nonnegative symmetric endpoint matrix')
    if not math.isfinite(bond_threshold) or bond_threshold < 0:
        raise ValueError('bond threshold must be finite and nonnegative')
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
    edges = tuple((int(a), int(b)) for a, b in zip(*np.nonzero(np.triu(local > bond_threshold, 1))))
    return BinarySymmetry(atoms, labels, anchors, float(bond_threshold), None, local, edges)
