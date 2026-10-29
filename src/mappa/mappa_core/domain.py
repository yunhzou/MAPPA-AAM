"""Typed symmetry objects shared by mapping, conversion and TS processing."""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class AtomPermutation:
    """One exact permutation written as ``atom -> image``."""

    images: tuple[int, ...]

    def __post_init__(self):
        images = tuple(int(value) for value in self.images)
        if set(images) != set(range(len(images))):
            raise ValueError("atom permutation is not bijective")
        object.__setattr__(self, "images", images)

    @classmethod
    def identity(cls, degree: int):
        return cls(tuple(range(int(degree))))

    @classmethod
    def from_mapping(cls, mapping: Mapping[int, int], degree: int | None = None):
        raw = {int(atom): int(image) for atom, image in mapping.items()}
        degree = len(raw) if degree is None else int(degree)
        if set(raw) != set(range(degree)):
            raise ValueError("permutation mapping has an incomplete domain")
        return cls(tuple(raw[atom] for atom in range(degree)))

    @property
    def degree(self):
        return len(self.images)

    @property
    def is_identity(self):
        return all(atom == image for atom, image in enumerate(self.images))

    def inverse(self):
        inverse = [0] * self.degree
        for atom, image in enumerate(self.images):
            inverse[image] = atom
        return AtomPermutation(tuple(inverse))

    def then(self, other: "AtomPermutation"):
        """Return ``other(self(atom))``."""
        if self.degree != other.degree:
            raise ValueError("cannot compose permutations of different degree")
        return AtomPermutation(tuple(other.images[image] for image in self.images))

    def cycles(self):
        cycles, unseen = [], set(range(self.degree))
        while unseen:
            start = min(unseen)
            cycle, atom = [], start
            while atom in unseen:
                unseen.remove(atom)
                cycle.append(atom)
                atom = self.images[atom]
            if len(cycle) > 1:
                cycles.append(tuple(cycle))
        return tuple(cycles)


@dataclass(frozen=True)
class AtomBijection:
    """One complete R-to-P assignment."""

    images: tuple[int, ...]

    def __post_init__(self):
        images = tuple(int(value) for value in self.images)
        if set(images) != set(range(len(images))):
            raise ValueError("atom mapping is not a complete bijection")
        object.__setattr__(self, "images", images)

    @classmethod
    def from_mapping(cls, mapping: Mapping[int, int], degree: int | None = None):
        raw = {int(atom): int(image) for atom, image in mapping.items()}
        degree = len(raw) if degree is None else int(degree)
        if set(raw) != set(range(degree)):
            raise ValueError("atom mapping has an incomplete R domain")
        return cls(tuple(raw[atom] for atom in range(degree)))

    @property
    def degree(self):
        return len(self.images)

    def as_dict(self):
        return {atom: image for atom, image in enumerate(self.images)}

    def inverse(self):
        return AtomBijection(AtomPermutation(self.images).inverse().images)

    def act(self, source: AtomPermutation | None = None,
            target: AtomPermutation | None = None):
        """Return ``target o self o source^-1``.

        Source and target actions are explicit group elements.  No branch
        witness or correspondence search participates in this operation.
        """
        source = source or AtomPermutation.identity(self.degree)
        target = target or AtomPermutation.identity(self.degree)
        if source.degree != self.degree or target.degree != self.degree:
            raise ValueError("group action degree differs from atom mapping")
        source_inverse = source.inverse()
        return AtomBijection(tuple(
            target.images[self.images[source_inverse.images[atom]]]
            for atom in range(self.degree)))

    def product_in_reactant_order(self, product_coords):
        product = np.asarray(product_coords, dtype=float)
        if product.shape != (self.degree, 3):
            raise ValueError("product coordinates disagree with mapping degree")
        return product[np.asarray(self.images, dtype=int)]


@dataclass(frozen=True)
class PermutationGroup:
    """Generator representation of a finite atom-permutation group."""

    degree: int
    generators: tuple[AtomPermutation, ...] = ()

    def __post_init__(self):
        degree = int(self.degree)
        unique = {}
        for generator in self.generators:
            if generator.degree != degree:
                raise ValueError("group generator has the wrong degree")
            if not generator.is_identity:
                unique.setdefault(generator.images, generator)
        object.__setattr__(self, "degree", degree)
        object.__setattr__(self, "generators", tuple(unique.values()))

    @classmethod
    def trivial(cls, degree: int):
        return cls(int(degree), ())

    @classmethod
    def from_generator_mappings(cls, degree: int, generators):
        permutations = []
        for generator in generators:
            if isinstance(generator, Mapping):
                permutations.append(
                    AtomPermutation.from_mapping(generator, degree))
            else:
                permutation = AtomPermutation(tuple(map(int, generator)))
                if permutation.degree != int(degree):
                    raise ValueError("group generator has the wrong degree")
                permutations.append(permutation)
        return cls(int(degree), tuple(permutations))

    def orbits(self):
        parent = list(range(self.degree))

        def find(atom):
            while parent[atom] != atom:
                parent[atom] = parent[parent[atom]]
                atom = parent[atom]
            return atom

        def union(left, right):
            left, right = find(left), find(right)
            if left != right:
                parent[right] = left

        for generator in self.generators:
            for atom, image in enumerate(generator.images):
                union(atom, image)
        groups = {}
        for atom in range(self.degree):
            groups.setdefault(find(atom), []).append(atom)
        return tuple(tuple(members) for members in groups.values())


@dataclass(frozen=True)
class SymmetryDomain:
    """One compressed assignment domain carried by an AAM fragment."""

    r_atoms: tuple[int, ...]
    p_atoms: tuple[int, ...]
    source: str
    extendable: bool = False

    def __post_init__(self):
        object.__setattr__(self, "r_atoms", tuple(sorted(map(int, self.r_atoms))))
        object.__setattr__(self, "p_atoms", tuple(sorted(map(int, self.p_atoms))))


@dataclass(frozen=True)
class FragmentMatch:
    fragment_index: int
    island_index: int
    r_atoms: tuple[int, ...]
    deferred_edges: tuple[tuple[int, int], ...] = ()
    symmetry_domains: tuple[SymmetryDomain, ...] = ()
    representative_assignments: tuple[tuple[int, int], ...] = ()
    exact_fixed: tuple[int, ...] = ()
    multiplicity: int = 1
    automorph_domains: tuple[SymmetryDomain, ...] = ()
    target_generators: tuple[AtomPermutation, ...] | None = None
    compression_version: int = 1

    def __post_init__(self):
        object.__setattr__(self, "r_atoms", tuple(sorted(map(int, self.r_atoms))))
        object.__setattr__(self, "deferred_edges", tuple(sorted(
            tuple(sorted(map(int, edge))) for edge in self.deferred_edges)))
        object.__setattr__(self, "symmetry_domains", tuple(self.symmetry_domains))
        object.__setattr__(self, "representative_assignments", tuple(sorted(
            (int(source), int(target))
            for source, target in self.representative_assignments)))
        object.__setattr__(self, "exact_fixed", tuple(sorted(
            map(int, self.exact_fixed))))
        object.__setattr__(self, "multiplicity", int(self.multiplicity))
        object.__setattr__(self, "automorph_domains", tuple(
            self.automorph_domains))
        if self.target_generators is not None:
            object.__setattr__(self, "target_generators", tuple(
                self.target_generators))

    @property
    def has_exact_target_group(self):
        """Whether AAM supplied the exact group (including a trivial one)."""
        return self.target_generators is not None

    @property
    def fixed_assignments(self):
        """Assignments fixed by both saved pools and correlated target actions.

        ``exact_fixed`` locks an assignment inside the local pools. A recorded
        target action can still transport it together with the whole fragment.
        This is local to this fragment's recorded frame, not a whole-path lock.
        """
        if self.target_generators is None:
            raise ValueError('Finalize saved fragment symmetry before reading anchors')
        pools = tuple(d for d in self.symmetry_domains
                      if d.source != 'exact_automorph_group')
        pooled = {r for d in pools for r in d.r_atoms}
        fixed_targets = {p for r, p in self.representative_assignments
                         if r in self.exact_fixed or r not in pooled}
        moving = set()
        for domain in pools:
            free = set(domain.p_atoms) - fixed_targets
            if len(free) > 1:
                moving.update(free)
        for generator in self.target_generators:
            moving.update(p for p, q in enumerate(generator.images) if p != q)
        return tuple((r, p) for r, p in self.representative_assignments
                     if p not in moving)


@dataclass(frozen=True)
class AAMHierarchy:
    fragments: tuple[FragmentMatch, ...]

    @property
    def segments(self):
        """Shared base fragments and their correlated coordinate transforms."""
        return ((self, ()),)

    def relabel_target(self, action):
        """Keep a correlated transform as a shared view until inspected."""
        action = tuple(sorted(dict(action).items()))
        return AAMHierarchyView(self, action) if action else self

    def _materialize_target(self, action):
        """Transport a hierarchy by one correlated target permutation.

        Unmentioned atoms (including augmentation copies) are fixed. Exact
        generators are conjugated, not copied into the old target frame.
        """
        from dataclasses import replace
        from mappa._group_ops import conjugate_generators
        action = dict(action)
        extent = max(max(action, default=-1), max(action.values(), default=-1)) + 1
        generators = tuple(g for f in self.fragments for g in (f.target_generators or ()))
        if all(a == b for a, b in action.items()) and all(g.degree >= extent for g in generators):
            return self
        unique = {g.images: g for g in generators}
        transformed = conjugate_generators(tuple(unique), tuple(action.get(a, a) for a in range(extent)))
        transported = {old: unique[old] if old == moved else AtomPermutation(moved)
                       for old, moved in zip(unique, transformed, strict=True)}
        def image(atom):
            return int(action.get(atom, atom))
        def domain(item):
            return replace(item, p_atoms=tuple(image(a) for a in item.p_atoms))
        def generator(item):
            return transported[item.images]
        return AAMHierarchy(tuple(replace(fragment,
            representative_assignments=tuple((a, image(b)) for a, b in
                                               fragment.representative_assignments),
            symmetry_domains=tuple(domain(d) for d in fragment.symmetry_domains),
            automorph_domains=tuple(domain(d) for d in fragment.automorph_domains),
            target_generators=(None if fragment.target_generators is None else
                               tuple(generator(g) for g in fragment.target_generators)))
            for fragment in self.fragments))

    @classmethod
    def from_record(cls, branch_symmetry, *, generator_pool=None):
        # The owner controls the pool lifetime; never retain graphs globally.
        generator_pool = {} if generator_pool is None else generator_pool
        def permutation(raw):
            images = tuple(map(int, raw))
            if images not in generator_pool:
                generator_pool[images] = AtomPermutation(images)
            return generator_pool[images]
        fragments = []
        for position, raw in enumerate(
                dict(branch_symmetry or {}).get("fragments") or ()):
            symmetry = raw.get("symmetry") or {}
            raw_generators = symmetry.get("automorph_generators")
            domains = []
            for block in symmetry.get("blocks") or ():
                domains.append(SymmetryDomain(
                    tuple(block.get("r_atoms") or ()),
                    tuple(block.get("p_atoms") or ()),
                    str(block.get("source") or "sym_block"),
                    bool(block.get("extendable", False))))
            automorph_domains = tuple(SymmetryDomain(
                tuple(block.get("r_atoms") or ()),
                tuple(block.get("p_atoms") or ()),
                str(block.get("source") or "exact_automorph_group"),
                bool(block.get("extendable", False)))
                for block in symmetry.get("automorph_blocks") or ())
            fragments.append(FragmentMatch(
                fragment_index=int(raw.get("fragment_index", position)),
                island_index=int(raw.get("island_idx", position)),
                r_atoms=tuple(sorted(map(int, raw.get("fragment") or ()))),
                deferred_edges=tuple(
                    tuple(sorted(map(int, edge)))
                    for edge in raw.get("deferred_edges") or ()),
                symmetry_domains=tuple(domains),
                representative_assignments=tuple(
                    (int(source), int(target))
                    for source, target in dict(
                        symmetry.get("witness") or {}).items()),
                exact_fixed=tuple(map(
                    int, symmetry.get("exact_fixed") or ())),
                multiplicity=int(symmetry.get("multiplicity", 1)),
                compression_version=int(symmetry.get("compression_version", 1)),
                automorph_domains=automorph_domains,
                target_generators=(
                    None if raw_generators is None else tuple(
                        permutation(generator)
                        for generator in raw_generators))))
        return cls(tuple(fragments))

    @property
    def has_complete_exact_target_groups(self):
        return bool(self.fragments) and all(
            fragment.has_exact_target_group for fragment in self.fragments)

    def to_record(self):
        def domain_record(domain):
            return {
                "r_atoms": list(domain.r_atoms),
                "p_atoms": list(domain.p_atoms),
                "source": domain.source,
                "extendable": bool(domain.extendable),
            }

        fragments = []
        for fragment in self.fragments:
            symmetry = {
                "witness": dict(fragment.representative_assignments),
                "blocks": [
                    domain_record(domain)
                    for domain in fragment.symmetry_domains
                ],
                "exact_fixed": list(fragment.exact_fixed),
                "multiplicity": int(fragment.multiplicity),
                "compression_version": fragment.compression_version,
                "automorph_blocks": [
                    domain_record(domain)
                    for domain in fragment.automorph_domains
                ],
            }
            if fragment.target_generators is not None:
                symmetry["automorph_generators"] = [
                    list(generator.images)
                    for generator in fragment.target_generators
                ]
            fragments.append({
                "fragment_index": int(fragment.fragment_index),
                "island_idx": int(fragment.island_index),
                "fragment": list(fragment.r_atoms),
                "deferred_edges": [
                    list(edge) for edge in fragment.deferred_edges],
                "symmetry": symmetry,
            })
        return {
            "rule": "typed_aam_hierarchy",
            "fragments": fragments,
            "blocks": [],
        }


@dataclass(frozen=True)
class AAMHierarchyView:
    """A base hierarchy plus one coordinate transform, without copied groups."""

    base: AAMHierarchy
    target_action: tuple[tuple[int, int], ...]

    @property
    def segments(self):
        return ((self.base, self.target_action),)

    @cached_property
    def materialized(self):
        return self.base._materialize_target(dict(self.target_action))

    @property
    def fragments(self):
        return self.materialized.fragments

    @property
    def has_complete_exact_target_groups(self):
        return self.base.has_complete_exact_target_groups

    def relabel_target(self, action):
        action, prior = dict(action), dict(self.target_action)
        if not action:
            return self
        combined = tuple(sorted((a, action.get(prior.get(a, a), prior.get(a, a)))
                                for a in set(prior) | set(action)))
        return AAMHierarchyView(self.base, combined)

    def to_record(self):
        """Explicitly requested, fully materialized legacy hierarchy record."""
        return self.materialized.to_record()


@dataclass(frozen=True)
class AAMHierarchyChain:
    """Concatenated fragment chains, each retaining its own coordinate frame.

    In particular, an augmented residual action must not transform the locked
    initial fragment. Composition preserves those independent frames without
    copying the fragment objects or conjugating their generators eagerly.
    """

    parts: tuple[AAMHierarchy | AAMHierarchyView, ...]

    @property
    def segments(self):
        return tuple(segment for part in self.parts for segment in part.segments)

    @cached_property
    def fragments(self):
        """Materialize only when a consumer explicitly requests full fragments."""
        return tuple(fragment for part in self.parts for fragment in part.fragments)

    @property
    def has_complete_exact_target_groups(self):
        return bool(self.parts) and all(part.has_complete_exact_target_groups for part in self.parts)

    def relabel_target(self, action):
        return AAMHierarchyChain(tuple(part.relabel_target(action) for part in self.parts))

    def to_record(self):
        return AAMHierarchy(self.fragments).to_record()


from collections import Counter
from dataclasses import dataclass, field, replace
from types import MappingProxyType
from typing import Any, Mapping

import numpy as np
import math

from mappa.mappa_core.core.search.graph import AAMSearchGraph, SearchBranch


def _readonly_array(value, *, shape=None, ndim=None):
    array = np.array(value, dtype=float, copy=True)
    if ndim is not None and array.ndim != int(ndim):
        raise ValueError(f"expected a {ndim}-dimensional array")
    if shape is not None and array.shape != tuple(shape):
        raise ValueError(f"expected array shape {tuple(shape)}, got {array.shape}")
    array.setflags(write=False)
    return array


def _frozen_mapping(value: Mapping[str, Any] | None):
    return MappingProxyType(dict(value or {}))


@dataclass(frozen=True)
class MolecularEndpoint:
    """One fully materialized molecular endpoint.

    Atom indices are local to the endpoint.  Coordinates and the complete WBO
    matrix are owned immutable copies, so later stages cannot silently mutate
    an AAM problem after search.
    """

    elements: tuple[str, ...]
    coordinates: np.ndarray
    wbo: np.ndarray
    label: str = ""
    energy: float | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        elements = tuple(str(element) for element in self.elements)
        if not elements:
            raise ValueError("a molecular endpoint must contain atoms")
        coordinates = _readonly_array(
            self.coordinates, shape=(len(elements), 3))
        wbo = _readonly_array(
            self.wbo, shape=(len(elements), len(elements)))
        if not np.allclose(wbo, wbo.T, atol=1e-12, rtol=0.0):
            raise ValueError("WBO matrix must be symmetric")
        if not np.all(np.isfinite(coordinates)) or not np.all(np.isfinite(wbo)):
            raise ValueError("endpoint arrays must contain finite values")
        object.__setattr__(self, "elements", elements)
        object.__setattr__(self, "coordinates", coordinates)
        object.__setattr__(self, "wbo", wbo)
        object.__setattr__(self, "label", str(self.label))
        object.__setattr__(
            self, "energy", None if self.energy is None else float(self.energy))
        object.__setattr__(self, "metadata", _frozen_mapping(self.metadata))

    @property
    def atom_count(self):
        return len(self.elements)


@dataclass(frozen=True)
class AAMProblem:
    """Reactant/product endpoint pair supplied to AAM search."""

    reactant: MolecularEndpoint
    product: MolecularEndpoint
    name: str = ""

    def __post_init__(self):
        object.__setattr__(self, "name", str(self.name))

    @property
    def balanced(self):
        return Counter(self.reactant.elements) == Counter(self.product.elements)

    @property
    def source_atom_count(self):
        return self.reactant.atom_count

    @property
    def target_atom_count(self):
        return self.product.atom_count

    @property
    def atom_count(self):
        """Source endpoint size; use target_atom_count for the other side."""
        return self.reactant.atom_count


@dataclass(frozen=True)
class AAMSearchConfig:
    """All hypotheses that affect AAM search and mechanism classification."""

    cut_floor: float = 0.2
    graph_floor: float = 0.2
    iso_tolerance: float = 1.0
    event_threshold: float = 0.5
    metal_event_threshold: float | None = 0.3
    seed_count: int = 1
    branch_limit: int = 100
    task_chunksize: int = 1
    anchors: tuple[tuple[int, int], ...] = ()
    seed_selection: str = 'random'
    random_seed: int = 42
    sweep_cuts: bool = True

    def __post_init__(self):
        for name in ('cut_floor', 'graph_floor', 'iso_tolerance', 'event_threshold'):
            if not math.isfinite(getattr(self, name)):
                raise ValueError(f'{name} must be finite')
        if self.metal_event_threshold is not None and (
                not math.isfinite(self.metal_event_threshold) or self.metal_event_threshold <= 0):
            raise ValueError('metal event threshold must be finite and positive')
        if not isinstance(self.random_seed, int):
            raise ValueError('random_seed must be an integer')
        if not isinstance(self.sweep_cuts, bool):
            raise ValueError('sweep_cuts must be boolean')
        if self.seed_selection not in ('random', 'distance'):
            raise ValueError('unknown seed selection policy')
        if self.cut_floor <= 0 or self.graph_floor <= 0:
            raise ValueError("graph and cut floors must be positive")
        if self.iso_tolerance <= 0:
            raise ValueError("isomorphism tolerance must be positive")
        if self.event_threshold <= 0:
            raise ValueError("event threshold must be positive")
        if self.seed_count < 1 or self.branch_limit < 1:
            raise ValueError("seed count and branch limit must be positive")
        if self.task_chunksize < 1:
            raise ValueError("task chunksize must be positive")
        anchors = tuple(sorted(
            (int(source), int(target)) for source, target in self.anchors))
        if len({source for source, _target in anchors}) != len(anchors):
            raise ValueError("an AAM source atom cannot have multiple anchors")
        if len({target for _source, target in anchors}) != len(anchors):
            raise ValueError("an AAM target atom cannot have multiple anchors")
        object.__setattr__(self, "anchors", anchors)


@dataclass(frozen=True)
class AAMSearchMetrics:
    elapsed_seconds: float
    cut_count: int
    raw_result_count: int
    retained_branch_count: int
    max_live_branches: int
    max_growth_candidates: int
    subtree_branch_cap_count: int
    worker_returned_branch_count: int
    parent_merge_seconds: float
    completed_group_requests: int = 0
    completed_group_calculations: int = 0
    completed_group_cache_hits: int = 0
    worker_search_seconds: float | None = None
    checkpoint_seconds: float | None = None
    symmetry_finalization_seconds: float | None = None
    checkpoint_restore_and_finalize_seconds: float | None = None

    @classmethod
    def from_record(cls, record, elapsed_seconds):
        record = dict(record or {})
        return cls(
            elapsed_seconds=float(elapsed_seconds),
            cut_count=int(record.get("cuts", 0)),
            raw_result_count=int(record.get("raw_result_count", 0)),
            retained_branch_count=int(record.get("retained_branch_count", 0)),
            max_live_branches=int(record.get("max_live_branches", 0)),
            max_growth_candidates=int(record.get("max_growth_candidates", 0)),
            subtree_branch_cap_count=int(
                record.get("subtree_branch_cap_count", 0)),
            worker_returned_branch_count=int(
                record.get("worker_returned_branch_count", 0)),
            parent_merge_seconds=float(record.get("parent_merge_seconds", 0.0)),
            completed_group_requests=int(record.get(
                "completed_candidate_group_requests", 0)),
            completed_group_calculations=int(record.get(
                "completed_candidate_group_calculations", 0)),
            completed_group_cache_hits=int(record.get(
                "completed_candidate_group_cache_hits", 0)),
            worker_search_seconds=record.get('worker_search_seconds'),
            checkpoint_seconds=record.get('checkpoint_seconds'),
            symmetry_finalization_seconds=record.get('symmetry_finalization_seconds'),
            checkpoint_restore_and_finalize_seconds=record.get('checkpoint_restore_and_finalize_seconds'),
        )


@dataclass(frozen=True)
class AAMResult:
    """Grouped final fragment sets and their unmodified search alternatives."""

    problem: AAMProblem
    config: AAMSearchConfig
    graph: AAMSearchGraph
    metrics: AAMSearchMetrics

    branches: tuple[SearchBranch, ...] = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        branches = self.graph.branches()
        object.__setattr__(self, 'branches', branches)
        object.__setattr__(self, 'metrics', replace(self.metrics,
            retained_branch_count=len(branches)))

    def final_catalogue(self):
        """Unordered branches and unique flat families for exact postprocessing."""
        from mappa.mappa_core.post_processing.final_branches import deduplicate_final_branches
        return deduplicate_final_branches(self)


@dataclass(frozen=True)
class ResolvedMechanism:
    """Explicit R/P mapping and events consumed by downstream TS processing.

    This is the exact interface required by downstream TS processing. The
    caller supplies the mapping and events explicitly; this object does not
    select an orientation or repeat the R/P search.
    """

    mapping: AtomBijection
    broken_bonds: tuple[tuple[int, int], ...]
    formed_bonds: tuple[tuple[int, int], ...]
    core_atoms: tuple[int, ...]

    def __post_init__(self):
        object.__setattr__(self, "broken_bonds", tuple(sorted(
            tuple(sorted(map(int, bond))) for bond in self.broken_bonds)))
        object.__setattr__(self, "formed_bonds", tuple(sorted(
            tuple(sorted(map(int, bond))) for bond in self.formed_bonds)))
        expected_core = tuple(sorted({
            atom for bond in (*self.broken_bonds, *self.formed_bonds)
            for atom in bond
        }))
        core = tuple(sorted(map(int, self.core_atoms)))
        if core != expected_core:
            raise ValueError("resolved mechanism core differs from its events")
        object.__setattr__(self, "core_atoms", core)


@dataclass(frozen=True)
class ReactionContext:
    """Materialized endpoint pair and selected mechanisms consumed by TS."""

    problem: AAMProblem
    config: AAMSearchConfig
    mechanisms: tuple[ResolvedMechanism, ...]

    def __post_init__(self):
        mechanisms = tuple(self.mechanisms)
        if any(item.mapping.degree != self.problem.atom_count
               for item in mechanisms):
            raise ValueError("resolved mechanism has the wrong mapping degree")
        object.__setattr__(self, "mechanisms", mechanisms)


@dataclass(frozen=True)
class VibrationalModes:
    """Frequencies and Cartesian normal modes in one target atom order."""

    frequencies: np.ndarray
    displacements: np.ndarray

    def __post_init__(self):
        frequencies = _readonly_array(self.frequencies, ndim=1)
        displacements = _readonly_array(self.displacements, ndim=3)
        if displacements.shape[0] != frequencies.shape[0]:
            raise ValueError("mode and frequency counts differ")
        if displacements.shape[2] != 3:
            raise ValueError("normal-mode displacements must be Cartesian")
        object.__setattr__(self, "frequencies", frequencies)
        object.__setattr__(self, "displacements", displacements)




@dataclass(frozen=True)
class AtomAssignment:
    """An injective partial atom assignment with an explicit source domain."""

    pairs: tuple[tuple[int, int], ...]

    def __post_init__(self):
        pairs = tuple(sorted((int(a), int(b)) for a, b in self.pairs))
        if len({a for a, _ in pairs}) != len(pairs):
            raise ValueError("partial assignment repeats a source atom")
        if len({b for _, b in pairs}) != len(pairs):
            raise ValueError("partial assignment repeats a target atom")
        object.__setattr__(self, "pairs", pairs)

    @classmethod
    def from_mapping(cls, mapping):
        return cls(tuple(dict(mapping).items()))

    def as_dict(self):
        return dict(self.pairs)

    @property
    def source_atoms(self):
        return tuple(source for source, _target in self.pairs)


@dataclass(frozen=True)
class CoreAAMBranch:
    """One partial-AAM growth branch and its exact correlated core orbit."""

    representative: AtomAssignment
    hierarchy: AAMHierarchy
    exact_assignments: tuple[AtomAssignment, ...]
    seed_index: int
    search_path: object = None

    def __post_init__(self):
        exact = tuple(self.exact_assignments)
        if not exact:
            raise ValueError("a core branch must realize an exact assignment")
        domain = self.representative.source_atoms
        if any(item.source_atoms != domain for item in exact):
            raise ValueError("core branch assignments disagree on their domain")
        object.__setattr__(self, "exact_assignments", exact)
        object.__setattr__(self, "seed_index", int(self.seed_index))


@dataclass(frozen=True)
class CoreAAMResult:
    """Exact correlated core assignments from endpoint to target."""

    source: MolecularEndpoint
    target: MolecularEndpoint
    core_atoms: tuple[int, ...]
    branches: tuple[CoreAAMBranch, ...]
    assignments: tuple[AtomAssignment, ...]
    elapsed_seconds: float
    capped_seed_count: int = 0

    def __post_init__(self):
        core = tuple(sorted(map(int, self.core_atoms)))
        if any(atom < 0 or atom >= self.source.atom_count for atom in core):
            raise ValueError("core atom lies outside the source molecule")
        assignments = tuple(self.assignments)
        if any(item.source_atoms != core for item in assignments):
            raise ValueError("core assignment has the wrong source domain")
        if any(target < 0 or target >= self.target.atom_count
               for item in assignments for _source, target in item.pairs):
            raise ValueError("core assignment lies outside the target molecule")
        object.__setattr__(self, "core_atoms", core)
        object.__setattr__(self, "branches", tuple(self.branches))
        object.__setattr__(self, "assignments", assignments)
        object.__setattr__(self, "elapsed_seconds", float(self.elapsed_seconds))
        object.__setattr__(self, "capped_seed_count", int(self.capped_seed_count))
