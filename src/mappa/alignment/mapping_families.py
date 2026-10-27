"""Exact structural and event-conditioned mapping families, without 3D selection."""
from __future__ import annotations

from collections import Counter, defaultdict
from typing import Mapping, Sequence

import numpy as np

from ..frag import bond_event_threshold, build_graph
from ..matcher.orbits import _wbo_tolerance_bucket_lookup


class MappingFamilyError(ValueError):
    """Invalid mapping input or selected-fragment metadata."""


class MappingFamilyConflict(MappingFamilyError):
    """The saved mapping-family relation has no endpoint isomorphism."""

    def __init__(self, message, *, diagnostics=None):
        super().__init__(message)
        self.diagnostics = diagnostics


class _RelationalGraph:
    """Two equally shaped colored graphs with shared incidence topology."""

    def __init__(self, atom_colors_A, atom_colors_B):
        if len(atom_colors_A) != len(atom_colors_B):
            raise MappingFamilyError("relational atom color lengths differ")
        self.colors_A = list(atom_colors_A)
        self.colors_B = list(atom_colors_B)
        self.adjacency = [set() for _ in atom_colors_A]
        self.atom_count = len(atom_colors_A)
        self.relation_records_A = []
        self.relation_records_B = []

    def _vertex(self, color_A, color_B=None):
        index = len(self.adjacency)
        self.adjacency.append(set())
        self.colors_A.append(color_A)
        self.colors_B.append(color_A if color_B is None else color_B)
        return index

    def _edge(self, left, right):
        self.adjacency[left].add(right)
        self.adjacency[right].add(left)

    def add_pair(self, left, right, color_A, color_B=None):
        color_B = color_A if color_B is None else color_B
        relation = self._vertex(color_A, color_B)
        self._edge(relation, int(left))
        self._edge(relation, int(right))
        atoms = tuple(sorted((int(left), int(right))))
        self.relation_records_A.append(("pair", repr(color_A), atoms))
        self.relation_records_B.append(("pair", repr(color_B), atoms))


    def clone(self):
        result = _RelationalGraph([], [])
        result.colors_A = list(self.colors_A)
        result.colors_B = list(self.colors_B)
        result.adjacency = [set(neighbors) for neighbors in self.adjacency]
        result.atom_count = self.atom_count
        result.relation_records_A = list(self.relation_records_A)
        result.relation_records_B = list(self.relation_records_B)
        return result

    def graph(self, side, *, individualized=()):
        try:
            import pynauty
        except ImportError as exc:
            raise RuntimeError(
                "pynauty is required for exact mapping-family relations"
            ) from exc
        colors = list(self.colors_A if side == "A" else self.colors_B)
        for rank, vertex in enumerate(individualized):
            vertex = int(vertex)
            colors[vertex] = (colors[vertex], "individualized", rank)
        classes = defaultdict(set)
        for vertex, color in enumerate(colors):
            classes[color].add(vertex)
        coloring = [
            set(vertices) for _, vertices in
            sorted(classes.items(), key=lambda item: repr(item[0]))
        ]
        return pynauty.Graph(
            len(self.adjacency), directed=False,
            adjacency_dict={
                vertex: sorted(neighbors)
                for vertex, neighbors in enumerate(self.adjacency)
            },
            vertex_coloring=coloring,
        )


def _int_mapping(mapping: Mapping[int, int]) -> dict[int, int]:
    return {int(r): int(p) for r, p in dict(mapping).items()}


def validate_mapping(mapping, elements_R, elements_P):
    mapping = _int_mapping(mapping)
    atom_count = len(elements_R)
    if len(elements_P) != atom_count:
        raise MappingFamilyError("endpoint atom counts differ")
    if set(mapping) != set(range(atom_count)):
        raise MappingFamilyError("mapping is not complete")
    if set(mapping.values()) != set(range(atom_count)):
        raise MappingFamilyError("mapping is not bijective")
    for r, p in mapping.items():
        if str(elements_R[r]) != str(elements_P[p]):
            raise MappingFamilyError(f"element mismatch at R{r}->P{p}")
    return mapping


def mapping_event_signature(mapping, wbo_R, wbo_P, elements_R, *,
                            dwbo_threshold=0.5,
                            metal_dwbo_threshold=0.3):
    """Exact broken/formed R-index pairs for one complete mapping."""
    mapping = _int_mapping(mapping)
    broken, formed = [], []
    for left in range(len(elements_R)):
        for right in range(left + 1, len(elements_R)):
            threshold = bond_event_threshold(
                elements_R, left, right,
                default_threshold=float(dwbo_threshold),
                metal_threshold=metal_dwbo_threshold)
            delta = float(
                wbo_R[left, right]
                - wbo_P[mapping[left], mapping[right]])
            if delta >= threshold:
                broken.append((left, right))
            elif -delta >= threshold:
                formed.append((left, right))
    return tuple(broken), tuple(formed)


def _event_class(r_wbo, p_wbo, threshold):
    delta = float(r_wbo) - float(p_wbo)
    if delta >= threshold:
        return 1
    if -delta >= threshold:
        return -1
    return 0


def _selected_fragments(branch_symmetry, source, fixed_singletons=()):
    raw = list(dict(branch_symmetry or {}).get("fragments") or ())
    if not raw:
        raw = [{"fragment_index": 0, "fragment": sorted(source)}]
    fragments, owner = [], {}
    for position, record in enumerate(raw):
        atoms = tuple(sorted(int(atom) for atom in record.get("fragment", ())))
        if not atoms:
            raise MappingFamilyError("selected AAM fragment is empty")
        for atom in atoms:
            if atom not in source:
                raise MappingFamilyError(
                    f"selected fragment contains unknown R atom {atom}")
            if atom in owner:
                raise MappingFamilyError(
                    f"selected fragments overlap at R atom {atom}")
            owner[atom] = position
        fragments.append(atoms)
    missing = sorted(set(source) - set(owner))
    fixed_singletons = set(map(int, fixed_singletons or ()))
    uncovered = sorted(set(missing) - fixed_singletons)
    if uncovered:
        raise MappingFamilyError(
            f"selected AAM fragments do not cover R atoms {uncovered}")
    # Hard anchors are preloaded into the branch before island growth.  If an
    # anchor never participates in a grown fragment, it is still a complete
    # individually fixed analytical fragment rather than missing hierarchy.
    for atom in missing:
        owner[atom] = len(fragments)
        fragments.append((atom,))
    return tuple(fragments), owner


def analytical_family_static_context(
        elements_R, wbo_R, elements_P, wbo_P, *, graph_floor=0.2,
        dwbo_threshold=0.5, metal_dwbo_threshold=0.3):
    """Precompute endpoint-only data shared by every analytical branch."""
    wbo_R = np.asarray(wbo_R)
    wbo_P = np.asarray(wbo_P)
    atom_count = len(elements_R)
    pair_groups_R = defaultdict(list)
    pair_groups_P = defaultdict(list)
    records_R, records_P = {}, {}
    for side, elements, wbo, groups, records in (
            ('R', elements_R, wbo_R, pair_groups_R, records_R),
            ('P', elements_P, wbo_P, pair_groups_P, records_P)):
        for left in range(atom_count):
            for right in range(left + 1, atom_count):
                element_pair = tuple(sorted((
                    str(elements[left]), str(elements[right]))))
                threshold = bond_event_threshold(
                    elements, left, right,
                    default_threshold=float(dwbo_threshold),
                    metal_threshold=metal_dwbo_threshold)
                group_key = (element_pair, float(threshold))
                value = float(wbo[left, right])
                groups[group_key].append(value)
                records[(left, right)] = (group_key, value)
    values_R = {
        key: tuple(sorted(set(values)))
        for key, values in pair_groups_R.items()
    }
    values_P = {
        key: tuple(sorted(set(values)))
        for key, values in pair_groups_P.items()
    }
    behavior_R = {}
    for pair, (group_key, value) in records_R.items():
        threshold = group_key[1]
        behavior_R[pair] = tuple(
            _event_class(value, other, threshold)
            for other in values_P[group_key])
    behavior_P = {}
    for pair, (group_key, value) in records_P.items():
        threshold = group_key[1]
        behavior_P[pair] = tuple(
            _event_class(other, value, threshold)
            for other in values_R[group_key])
    return {
        'graph_R': build_graph(elements_R, wbo_R, bond_cut=graph_floor),
        'graph_P': build_graph(elements_P, wbo_P, bond_cut=graph_floor),
        'pair_records_R': records_R,
        'pair_records_P': records_P,
        'behavior_R': behavior_R,
        'behavior_P': behavior_P,
    }


def _masked_relation_data(source, branch_symmetry, elements_R, wbo_R,
                          elements_P, wbo_P, graph_floor,
                          symmetry_wbo_tol, dwbo_threshold,
                          metal_dwbo_threshold, anchor_map,
                          static_context=None, include_event_relations=True):
    atom_count = len(source)
    inverse = {p: r for r, p in source.items()}
    anchors = {int(r): int(p) for r, p in dict(anchor_map or {}).items()}
    for r, p in anchors.items():
        if r not in source or p not in inverse:
            raise MappingFamilyError(f"invalid anchor R{r}->P{p}")
        if str(elements_R[r]) != str(elements_P[p]):
            raise MappingFamilyError(f"anchor element mismatch R{r}->P{p}")
    fragments, owner_R = _selected_fragments(
        branch_symmetry, source, fixed_singletons=anchors)
    owner_P = {source[r]: owner_R[r] for r in source}

    source_anchor_tags = defaultdict(list)
    target_anchor_tags = defaultdict(list)
    for r, p in sorted(anchors.items()):
        source_anchor_tags[source[r]].append(r)
        target_anchor_tags[p].append(r)
    atom_colors_A, atom_colors_B = [], []
    for p in range(atom_count):
        common = ("atom", str(elements_P[p]), owner_P[p])
        atom_colors_A.append(common + (tuple(source_anchor_tags[p]),))
        atom_colors_B.append(common + (tuple(target_anchor_tags[p]),))
    relation = _RelationalGraph(atom_colors_A, atom_colors_B)

    static_context = static_context or analytical_family_static_context(
        elements_R, wbo_R, elements_P, wbo_P,
        graph_floor=graph_floor,
        dwbo_threshold=dwbo_threshold,
        metal_dwbo_threshold=metal_dwbo_threshold)
    g_R = static_context['graph_R']
    g_P = static_context['graph_P']
    persistent_P = {p: set() for p in range(atom_count)}
    for fragment_id, fragment_R in enumerate(fragments):
        fragment_P = tuple(source[r] for r in fragment_R)
        sub_R = g_R.subgraph(fragment_R).copy()
        sub_P = g_P.subgraph(fragment_P).copy()
        buckets_P, zero_P = _wbo_tolerance_bucket_lookup(
            sub_P, float(symmetry_wbo_tol))
        # A completed AAM branch is a target-group coset acting on its
        # validated R->P mapping.  Reimposing masked R automorphisms here
        # would intersect Aut(R) with Aut(P) and delete valid target shuffles.
        # Put the selected masked P relation identically on both sides; event
        # and anchor relations below take the required stabilizers.
        for (left_P, right_P), bucket in sorted(buckets_P.items()):
            if bucket != zero_P:
                relation.add_pair(
                    left_P, right_P,
                    ("masked_target_wbo", fragment_id, bucket))
        for left_R in fragment_R:
            for right_R in fragment_R:
                if left_R >= right_R or not sub_R.has_edge(left_R, right_R):
                    continue
                left_P, right_P = source[left_R], source[right_R]
                if sub_P.has_edge(left_P, right_P):
                    persistent_P[left_P].add(right_P)
                    persistent_P[right_P].add(left_P)

    if not include_event_relations:
        return relation, persistent_P, inverse, fragments

    # These complete-pair colors are a conservative exact-event stabilizer:
    # they preserve the event classification against every endpoint scalar
    # value, not merely the event realized by this representative.  The
    # quotient between the structural group above and this stabilizer is
    # handled analytically by AnalyticalMappingFamily; otherwise valid
    # same-event cosets can be lost at a threshold boundary.
    #
    # Within this stabilizer, the complete-pair colors are the coarsest
    # endpoint scalar relations needed to make every broken/formed event
    # classification invariant.
    # They split an iso-tolerance class only when crossing an event boundary.
    pair_records = []
    for left_P in range(atom_count):
        for right_P in range(left_P + 1, atom_count):
            left_R, right_R = inverse[left_P], inverse[right_P]
            pair_R = tuple(sorted((left_R, right_R)))
            pair_P = (left_P, right_P)
            group_key = static_context['pair_records_R'][pair_R][0]
            if static_context['pair_records_P'][pair_P][0] != group_key:
                raise MappingFamilyError(
                    "mapping changed an endpoint element-pair class")
            pair_records.append((left_P, right_P, group_key,
                                 static_context['behavior_R'][pair_R],
                                 static_context['behavior_P'][pair_P]))
    colored_pair_records = []
    color_counts = defaultdict(Counter)
    for left_P, right_P, group_key, r_behavior, p_behavior in pair_records:
        color = (
            "event_invariant_pair", group_key, r_behavior, p_behavior)
        colored_pair_records.append((left_P, right_P, group_key, color))
        color_counts[group_key][color] += 1

    # This is a complete edge-colored relation within each element/threshold
    # pair class.  Its most frequent color can be represented by absence:
    # atom colors already preserve the pair class, and any permutation that
    # preserves every exceptional colored pair necessarily maps the remaining
    # (baseline) pairs among themselves.  This is exactly equivalent to
    # materializing O(N^2) relation vertices, while molecular graphs normally
    # retain only O(E) exceptional event relations.
    baseline_color = {
        group_key: min(
            counts,
            key=lambda color: (-counts[color], repr(color)),
        )
        for group_key, counts in color_counts.items()
    }
    for left_P, right_P, group_key, color in colored_pair_records:
        if color == baseline_color[group_key]:
            continue
        relation.add_pair(
            left_P, right_P,
            color)

    return relation, persistent_P, inverse, fragments


def _canonical_isomorphism(graph_A, graph_B):
    import pynauty

    if pynauty.certificate(graph_A) != pynauty.certificate(graph_B):
        return None
    label_A = pynauty.canon_label(graph_A)
    label_B = pynauty.canon_label(graph_B)
    return {int(label_A[i]): int(label_B[i])
            for i in range(len(label_A))}


class AnalyticalMappingFamily:
    """Exact isomorphism coset represented by one completed AAM branch.

    The branch relation defines ``Iso(A, B)``.  Restricting those
    isomorphisms to atom vertices and composing with the branch's R->P source
    mapping gives every concrete bijection represented by the branch.  Group
    equality is proven by generator containment; no elements are enumerated.
    """

    def __init__(self, source_mapping, relation, structural_relation=None):
        import pynauty

        self.source_mapping = dict(source_mapping)
        self.degree = len(self.source_mapping)
        self.relation = relation
        self.graph_A = relation.graph("A")
        self.graph_B = relation.graph("B")
        self.certificate_A = pynauty.certificate(self.graph_A)
        self.certificate_B = pynauty.certificate(self.graph_B)
        isomorphism = _canonical_isomorphism(self.graph_A, self.graph_B)
        if isomorphism is None:
            raise MappingFamilyConflict(
                "analytical AAM branch relation has no endpoint isomorphism")
        self.representative_mapping = {
            r: int(isomorphism[self.source_mapping[r]])
            for r in sorted(self.source_mapping)
        }
        (raw_generators, mantissa, exponent,
         raw_orbits, _orbit_count) = pynauty.autgrp(self.graph_B)
        identity = tuple(range(self.degree))
        self.target_generators = tuple(dict.fromkeys(
            tuple(int(generator[atom]) for atom in range(self.degree))
            for generator in raw_generators
            if tuple(int(generator[atom])
                     for atom in range(self.degree)) != identity
        ))
        orbit_members = {}
        for atom in range(self.degree):
            orbit_members.setdefault(int(raw_orbits[atom]), []).append(atom)
        self.target_orbits = tuple(sorted(
            (tuple(members) for members in orbit_members.values()),
            key=lambda members: members[0]))
        self.group_order = (round(float(mantissa), 12), int(exponent))
        self._membership_cache = {}
        self._relation_record_counts_B = Counter(
            self.relation.relation_records_B)
        self.structural_relation = structural_relation
        self.structural_target_generators = ()
        self.structural_group_order = self.group_order
        if structural_relation is not None:
            structural_A = structural_relation.graph("A")
            structural_B = structural_relation.graph("B")
            if _canonical_isomorphism(structural_A, structural_B) is None:
                raise MappingFamilyConflict(
                    "pre-event AAM branch relation has no endpoint "
                    "isomorphism")
            (structural_generators, structural_mantissa,
             structural_exponent, _structural_orbits,
             _structural_orbit_count) = pynauty.autgrp(structural_B)
            self.structural_target_generators = tuple(dict.fromkeys(
                tuple(int(generator[atom]) for atom in range(self.degree))
                for generator in structural_generators
                if tuple(int(generator[atom])
                         for atom in range(self.degree)) != identity
            ))
            self.structural_group_order = (
                round(float(structural_mantissa), 12),
                int(structural_exponent))

    @property
    def invariant(self):
        return self.degree, self.group_order, self.target_orbits

    @property
    def equivalence_bucket(self):
        """Fast canonical bucket; equality is still proven explicitly."""
        return (self.invariant, self.certificate_A, self.certificate_B)

    def with_coset_representative(self, source_mapping):
        """Reuse this exact target subgroup for another right-coset member.

        If ``K`` is the compiled conservative target subgroup and ``g`` is a
        structural quotient representative, the corresponding mapping family
        is ``K g m``.  Its target automorphism group is still exactly ``K``;
        only the R->P coset representative changes.  Build an A-side copy of
        the compiled B relation so the identity A->B isomorphism represents
        the new source, without rerunning graph construction or pynauty.
        """
        source = _int_mapping(source_mapping)
        if (set(source) != set(range(self.degree))
                or set(source.values()) != set(range(self.degree))):
            raise MappingFamilyError(
                "translated analytical coset mapping is not bijective")
        relation = self.relation.clone()
        relation.colors_A = list(relation.colors_B)
        relation.relation_records_A = list(relation.relation_records_B)

        translated = object.__new__(AnalyticalMappingFamily)
        translated.source_mapping = source
        translated.degree = self.degree
        translated.relation = relation
        translated.graph_A = relation.graph("A")
        translated.graph_B = self.graph_B
        translated.certificate_A = self.certificate_B
        translated.certificate_B = self.certificate_B
        translated.representative_mapping = dict(source)
        translated.target_generators = self.target_generators
        translated.target_orbits = self.target_orbits
        translated.group_order = self.group_order
        translated._membership_cache = {}
        translated._relation_record_counts_B = (
            self._relation_record_counts_B)
        translated.structural_relation = None
        translated.structural_target_generators = ()
        translated.structural_group_order = self.group_order
        return translated

    def contains(self, mapping):
        mapping = validate_mapping(
            mapping,
            ["X"] * self.degree,
            ["X"] * self.degree)
        key = tuple(mapping[r] for r in range(self.degree))
        cached = self._membership_cache.get(key)
        if cached is not None:
            return cached
        source_order = tuple(self.source_mapping[r]
                             for r in range(self.degree))
        target_order = tuple(mapping[r] for r in range(self.degree))
        sigma = {}
        for source_atom, target_atom in zip(source_order, target_order):
            sigma[int(source_atom)] = int(target_atom)
        atom_colors_match = all(
            self.relation.colors_A[atom]
            == self.relation.colors_B[sigma[atom]]
            for atom in range(self.degree)
        )
        transported = Counter()
        if atom_colors_match:
            for kind, color, atoms in self.relation.relation_records_A:
                images = tuple(sigma[atom] for atom in atoms)
                if kind == "pair":
                    images = tuple(sorted(images))
                record = (kind, color, images)
                if record not in self._relation_record_counts_B:
                    atom_colors_match = False
                    break
                transported[record] += 1
        result = (atom_colors_match
                  and transported == self._relation_record_counts_B)
        self._membership_cache[key] = bool(result)
        return bool(result)

    @staticmethod
    def _left_act(generator, mapping):
        return {
            r: int(generator[mapping[r]])
            for r in mapping
        }

    def exact_event_coset_representatives(
            self, wbo_R, wbo_P, elements_R, *,
            dwbo_threshold=0.5, metal_dwbo_threshold=0.3):
        """Return one mapping per same-event structural quotient coset.

        The structural AAM group can be larger than the conservative colored
        relation used to stabilize event thresholds.  We traverse the finite
        left-coset quotient by Schreier-style generator steps and retain only
        cosets whose representative realizes the source mapping's exact
        broken/formed event.  This enumerates event families, never group
        elements or atom bijections.
        """
        from sympy.combinatorics import Permutation, PermutationGroup

        identity = tuple(range(self.degree))
        structural_generators = (
            self.structural_target_generators or (identity,))
        conservative_generators = self.target_generators or (identity,)
        structural_group = PermutationGroup([
            Permutation(generator, size=self.degree)
            for generator in structural_generators
        ])
        conservative_group = PermutationGroup([
            Permutation(generator, size=self.degree)
            for generator in conservative_generators
        ])
        # SymPy's deterministic Schreier--Sims transversal constructs the
        # exact right quotient K\\G directly.  The previous breadth-first
        # traversal compared every new representative with every old one by
        # transporting the complete relation, making a q-element quotient
        # O(q^2 * |relation|).
        representatives = tuple(
            tuple(int(representative(atom))
                  for atom in range(self.degree))
            for representative in structural_group.coset_transversal(
                conservative_group)
        )

        # Filter quotient representatives against the concrete event in
        # vectorized batches.  This is exactly mapping_event_signature for
        # every representative, but avoids q repeated Python O(N^2) loops.
        wbo_R = np.asarray(wbo_R)
        wbo_P = np.asarray(wbo_P)
        left, right, thresholds = [], [], []
        for atom_left in range(self.degree):
            for atom_right in range(atom_left + 1, self.degree):
                left.append(atom_left)
                right.append(atom_right)
                thresholds.append(bond_event_threshold(
                    elements_R, atom_left, atom_right,
                    default_threshold=float(dwbo_threshold),
                    metal_threshold=metal_dwbo_threshold))
        left = np.asarray(left, dtype=int)
        right = np.asarray(right, dtype=int)
        thresholds = np.asarray(thresholds, dtype=float)
        r_values = wbo_R[left, right]
        source_classes = np.zeros(len(left), dtype=np.int8)
        source_delta = r_values - wbo_P[
            np.asarray([self.source_mapping[r] for r in left], dtype=int),
            np.asarray([self.source_mapping[r] for r in right], dtype=int)]
        source_classes[source_delta >= thresholds] = 1
        source_classes[-source_delta >= thresholds] = -1

        source_images = np.asarray([
            self.source_mapping[r] for r in range(self.degree)], dtype=int)
        valid_actions = []
        batch_size = 256
        for start in range(0, len(representatives), batch_size):
            actions = np.asarray(
                representatives[start:start + batch_size], dtype=int)
            candidate_mappings = actions[:, source_images]
            candidate_delta = (
                r_values[None, :]
                - wbo_P[candidate_mappings[:, left],
                        candidate_mappings[:, right]])
            candidate_classes = np.zeros(
                candidate_delta.shape, dtype=np.int8)
            candidate_classes[candidate_delta >= thresholds[None, :]] = 1
            candidate_classes[
                -candidate_delta >= thresholds[None, :]] = -1
            valid = np.all(
                candidate_classes == source_classes[None, :], axis=1)
            valid_actions.extend(
                tuple(map(int, action))
                for action in actions[valid])

        mappings = []
        seen = set()
        for action in valid_actions:
            mapping = self._left_act(action, self.source_mapping)
            key = tuple(mapping[r] for r in range(self.degree))
            if key not in seen:
                seen.add(key)
                mappings.append(mapping)
        return tuple(mappings)

    def equivalent(self, other):
        if not isinstance(other, AnalyticalMappingFamily):
            return False
        return self.is_subset_of(other) and other.is_subset_of(self)

    def is_subset_of(self, other):
        """Prove coset inclusion by representative/generator membership."""
        if not isinstance(other, AnalyticalMappingFamily):
            return False
        if self.degree != other.degree:
            return False
        self_log_order = (np.log10(self.group_order[0])
                          + self.group_order[1])
        other_log_order = (np.log10(other.group_order[0])
                           + other.group_order[1])
        if self_log_order > other_log_order + 1e-10:
            return False
        other_orbit = {
            atom: orbit_index
            for orbit_index, orbit in enumerate(other.target_orbits)
            for atom in orbit
        }
        # A subgroup orbit cannot cross two orbits of the containing group.
        # This is a necessary exact group-inclusion condition and rejects most
        # unrelated family pairs before any dense relation transport.
        if any(len({other_orbit[atom] for atom in orbit}) > 1
               for orbit in self.target_orbits):
            return False
        if any(
            other_orbit[self.representative_mapping[r]]
            != other_orbit[other.representative_mapping[r]]
            for r in range(self.degree)
        ):
            return False
        if not other.contains(self.representative_mapping):
            return False
        return not any(not other.contains(self._left_act(
            generator, self.representative_mapping))
            for generator in self.target_generators)

    def record(self):
        return {
            "representative_mapping": dict(self.representative_mapping),
            "target_generators": [list(generator)
                                  for generator in self.target_generators],
            "target_orbits": [list(orbit) for orbit in self.target_orbits],
            "group_order": {
                "mantissa": self.group_order[0],
                "decimal_exponent": self.group_order[1],
            },
            "pre_event_group_order": {
                "mantissa": self.structural_group_order[0],
                "decimal_exponent": self.structural_group_order[1],
            },
        }


def compile_analytical_mapping_family(
        source_mapping, branch_symmetry,
        elements_R: Sequence[str], wbo_R,
        elements_P: Sequence[str], wbo_P, *,
        graph_floor=0.2, symmetry_wbo_tol=0.2,
        dwbo_threshold=0.5, metal_dwbo_threshold=0.3,
        anchor_map=None, static_context=None, include_event_relations=True):
    """Compile one completed branch into its exact structural mapping coset."""
    source = validate_mapping(source_mapping, elements_R, elements_P)
    relation, _persistent, _inverse, _fragments = _masked_relation_data(
        source, branch_symmetry, elements_R, np.asarray(wbo_R),
        elements_P, np.asarray(wbo_P), graph_floor,
        symmetry_wbo_tol, dwbo_threshold,
        metal_dwbo_threshold, anchor_map or {},
        static_context=static_context, include_event_relations=include_event_relations)
    structural_relation, _persistent, _inverse, _fragments = (
        _masked_relation_data(
            source, branch_symmetry, elements_R, np.asarray(wbo_R),
            elements_P, np.asarray(wbo_P), graph_floor,
            symmetry_wbo_tol, dwbo_threshold,
            metal_dwbo_threshold, anchor_map or {},
            static_context=static_context,
            include_event_relations=False))
    return AnalyticalMappingFamily(
        source, relation, structural_relation=structural_relation)
