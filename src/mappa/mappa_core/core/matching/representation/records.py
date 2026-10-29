"""Fragment matching inputs and compressed growth results."""
from __future__ import annotations

from dataclasses import dataclass, field
from collections.abc import Mapping
from mappa.mappa_core.core.matching.policy import BondTolerance


class _IsoResult(dict):
    __slots__ = ('deferred_edges', 'fragment', 'symmetry')

    def __init__(self, mapping=None, deferred_edges=(), fragment=(),
                 symmetry=None):
        super().__init__(mapping or {})
        self.deferred_edges = frozenset(tuple(sorted(e)) for e in deferred_edges)
        self.fragment = frozenset(fragment)
        self.symmetry = symmetry or {}


class IslandBranchLimitExceeded(RuntimeError):
    """A completed fragment exceeded its distinct-placement branch limit."""

    def __init__(self, count, limit, *, seed=None):
        self.count = int(count)
        self.limit = int(limit)
        self.seed = seed
        message = f"fragment subtree branch cap hit: {count}>{limit}"
        if seed is not None:
            message += f" seed={seed}"
        super().__init__(message)


@dataclass(frozen=True)
class FragmentPlacement(Mapping):
    """Matched R–P fragment: witness plus correlated compressed alternatives."""
    mapping: dict
    fragment: frozenset
    deferred_edges: frozenset
    symmetry: dict
    preserved_bonds: tuple

    def __iter__(self):
        return iter(self.mapping)

    def __len__(self):
        return len(self.mapping)

    def __getitem__(self, atom):
        return self.mapping[atom]

    def items(self):
        return self.mapping.items()

    @classmethod
    def from_match(cls, match, source):
        """Attach the induced bonds once, without scanning unrelated atoms."""
        atoms, cuts = match.fragment, match.deferred_edges
        bonds = tuple(sorted((a, b) for a in atoms for b in source[a]
                             if a <= b and b in atoms and (a, b) not in cuts))
        return cls(dict(match), atoms, cuts, match.symmetry, bonds)


@dataclass(frozen=True)
class FragmentMatchContext:
    locked_mapping: dict = field(default_factory=dict)
    islands: dict | None = None
    deferred_edges: tuple = ()
    source_orbits: object = None
    target_orbits: object = None
    growth_replay: object = None

    @classmethod
    def from_saved_fragment(cls, fragment):
        """Use recorded fixed assignments as anchors in the fragment's frame.

        Within the supplied R/P fragment graphs, the matcher rediscovers free
        assignments and symmetry under its new bond tolerance.
        """
        anchors = dict(fragment.fixed_assignments)
        return cls(locked_mapping=anchors, islands={r: r for r in anchors},
                   deferred_edges=fragment.deferred_edges)


@dataclass(frozen=True)
class FragmentMatchConfig:
    graph_floor: float = 0.2
    iso_tolerance: float = 0.5
    minimum_size: int = 1
    branch_limit: int = 100  # Distinct placements after growth and deduplication.
    node_policy: object = None
    allow_mapped_seed: bool = False
    orbit_dedup: bool = True
    random_seed: int = 42
    bond_tolerance: BondTolerance | None = None


@dataclass(frozen=True)
class FragmentMatchResult:
    matches: tuple
    capped: bool
    branch_count: int
    branch_limit: int
    elapsed_seconds: float
