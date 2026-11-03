"""Complete discovery constraints with one permitted self action on each side."""
from collections.abc import Mapping
from dataclasses import dataclass, fields
from functools import cached_property
import math
from numbers import Integral

from .relaxed import RelaxedAutomorphism


def _pairs(values):
    return tuple(values.items()) if isinstance(values, Mapping) else tuple(values)


@dataclass(frozen=True, slots=True)
class DiscoveryEdge:
    """One accepted source edge and its original discovery predicate."""
    atoms: tuple[int, int]
    weight: float
    tolerance: float
    inclusive: bool = True

    def __post_init__(self):
        atoms = tuple(self.atoms)
        if (len(atoms) != 2 or atoms[0] == atoms[1]
                or any(not isinstance(atom, Integral) or atom < 0 for atom in atoms)):
            raise ValueError('discovery edge requires two distinct source atoms')
        if any(not math.isfinite(value) or value < 0 for value in (self.weight, self.tolerance)):
            raise ValueError('discovery weight and tolerance must be finite and nonnegative')
        object.__setattr__(self, 'atoms', tuple(sorted(map(int, atoms))))
        object.__setattr__(self, 'weight', float(self.weight))
        object.__setattr__(self, 'tolerance', float(self.tolerance))

    def allows(self, target_weight):
        difference = abs(self.weight - target_weight)
        return difference <= self.tolerance if self.inclusive else difference < self.tolerance


@dataclass(frozen=True, eq=False)
class RelaxedIsomorphism:
    """All valid discovery bijections, with independently constrained self actions.

    The complete relation is beta o discovery o alpha^-1. Discovery is described
    by domains, injectivity, anchors and accepted edges, never by a list of maps.
    Each self action is checked against its own original endpoint; the resulting
    relation is not closed repeatedly. Atom indices in this record are original
    endpoint indices; self generators use their existing local frames.
    """
    source: RelaxedAutomorphism
    target: RelaxedAutomorphism
    witness: tuple[tuple[int, int], ...]
    anchors: tuple[tuple[int, int], ...]
    domains: tuple[tuple[int, tuple[int, ...]], ...]
    accepted_edges: tuple[DiscoveryEdge, ...]
    target_edges: tuple[tuple[int, int], ...]
    graph_floor: float

    def __reduce__(self):
        return type(self), tuple(getattr(self, item.name) for item in fields(self))

    def __post_init__(self):
        witness, anchors = _pairs(self.witness), _pairs(self.anchors)
        domains = tuple((atom, tuple(values)) for atom, values in _pairs(self.domains))
        target_edges = tuple(tuple(edge) for edge in self.target_edges)
        if len(self.source_atoms) != len(self.target_atoms):
            raise ValueError('a fragment relation requires equal source and target atom counts')
        if not math.isfinite(self.graph_floor) or self.graph_floor < 0:
            raise ValueError('discovery graph floor must be finite and nonnegative')
        indices = (atom for pair in anchors for atom in pair)
        domain_indices = (atom for r, values in domains for atom in (r, *values))
        edge_indices = (atom for edge in target_edges for atom in edge)
        if any(not isinstance(atom, Integral) for values in (indices, domain_indices, edge_indices)
               for atom in values):
            raise ValueError('relation constraints require integer atom indices')
        if (len(anchors) != len(dict(anchors))
                or {r for r, _ in anchors} != set(self.source.anchors)
                or len(anchors) != len({p for _, p in anchors})
                or {p for _, p in anchors} != set(self.target.anchors)):
            raise ValueError('relation anchors must equal both self anchor projections')
        source_atoms, target_atoms = set(self.source_atoms), set(self.target_atoms)
        if (len(domains) != len(source_atoms) or {r for r, _ in domains} != source_atoms
                or any(not values or len(values) != len(set(values))
                       or not set(values) <= target_atoms for _, values in domains)):
            raise ValueError('discovery domains must cover each source atom within the fixed target set')
        edges = tuple(self.accepted_edges)
        if (len({edge.atoms for edge in edges}) != len(edges)
                or any(not set(edge.atoms) <= source_atoms for edge in edges)):
            raise ValueError('accepted edges must be distinct pairs of source atoms')
        if any(len(edge) != 2 or edge[0] == edge[1]
               or not set(edge) <= target_atoms for edge in target_edges):
            raise ValueError('active target edges must stay within the fixed target set')
        object.__setattr__(self, 'witness', tuple(sorted(witness)))
        object.__setattr__(self, 'anchors', tuple(sorted(anchors)))
        object.__setattr__(self, 'domains', tuple(sorted((r, tuple(sorted(values))) for r, values in domains)))
        object.__setattr__(self, 'accepted_edges', tuple(sorted(edges, key=lambda edge: edge.atoms)))
        object.__setattr__(self, 'target_edges', tuple(sorted({tuple(sorted(edge)) for edge in target_edges})))
        object.__setattr__(self, 'graph_floor', float(self.graph_floor))
        self._validate_discovery(self.witness)

    @property
    def source_atoms(self):
        return self.source.atoms

    @property
    def target_atoms(self):
        return self.target.atoms

    @cached_property
    def _target_positions(self):
        return {atom: index for index, atom in enumerate(self.target_atoms)}

    @cached_property
    def _active_edges(self):
        return frozenset(self.target_edges)

    def _allows_edge(self, edge, p, q):
        positions = self._target_positions
        weight = float(self.target.weights[positions[p], positions[q]])
        return (tuple(sorted((p, q))) in self._active_edges
                and weight >= self.graph_floor and edge.allows(weight))

    def discovery_pairs(self, edge):
        """Allowed endpoint pairs for one edge; global injectivity stays joint."""
        domains = dict(self.domains)
        left, right = (frozenset(domains[r]) for r in edge.atoms)
        return tuple((p, q) for a, b in self.target_edges for p, q in ((a, b), (b, a))
                     if p in left and q in right and self._allows_edge(edge, p, q))

    def _validate_discovery(self, discovery):
        pairs = _pairs(discovery)
        if any(not isinstance(atom, Integral) for pair in pairs for atom in pair):
            raise ValueError('discovery assignments require integer atom indices')
        mapping = dict(pairs)
        if (len(pairs) != len(self.source_atoms) or set(mapping) != set(self.source_atoms)
                or set(mapping.values()) != set(self.target_atoms)):
            raise ValueError('discovery must biject exactly the fixed source and target atom sets')
        if any(mapping[r] not in values for r, values in self.domains):
            raise ValueError('discovery violates an atom domain')
        if any(mapping[r] != p for r, p in self.anchors):
            raise ValueError('discovery violates an anchor')
        for edge in self.accepted_edges:
            p, q = (mapping[r] for r in edge.atoms)
            if not self._allows_edge(edge, p, q):
                raise ValueError('discovery violates an accepted edge')
        return mapping

    def allows_discovery(self, discovery):
        """Check one supplied discovery bijection without searching for another."""
        try:
            self._validate_discovery(discovery)
        except (TypeError, ValueError):
            return False
        return True

    def realize(self, discovery=None, *, source_word=(), target_word=()):
        """Return beta o discovery o alpha^-1 after all joint checks pass."""
        return self.realize_permutations(discovery,
            source_permutation=self.source._permutation(source_word),
            target_permutation=self.target._permutation(target_word))

    def realize_permutations(self, discovery=None, *, source_permutation=None,
                             target_permutation=None):
        """Validate explicit local-index self choices, including solver models."""
        discovery = self._validate_discovery(self.witness if discovery is None else discovery)
        alpha = tuple(range(len(self.source_atoms))) if source_permutation is None else source_permutation
        beta = tuple(range(len(self.target_atoms))) if target_permutation is None else target_permutation
        if not self.source.contains(alpha):
            raise ValueError('source rearrangement violates relaxed self tolerance or structural constraints')
        if not self.target.contains(beta):
            raise ValueError('target rearrangement violates relaxed self tolerance or structural constraints')
        inverse = {self.source_atoms[int(image)]: atom for atom, image in zip(self.source_atoms, alpha)}
        beta = {atom: self.target_atoms[int(image)] for atom, image in zip(self.target_atoms, beta)}
        return {atom: beta[discovery[inverse[atom]]] for atom in self.source_atoms}
