import numpy as np
from itertools import permutations
from mappa.coordinate_orientation import _simplex_measure
from mappa.alignment.mapping_families import (
    _masked_relation_data, analytical_family_static_context,
    compile_analytical_mapping_family, mapping_event_signature,
)


def test_four_ligand_simplex_degeneracy_is_permutation_invariant():
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.1, 0.0, 0.0],
        [0.1, 1.7, 0.0],
        [0.2, 0.3, 1.3],
        [3.0, 3.0, 3.0],  # unused center for the affine four-point frame
    ])
    measures = [
        _simplex_measure(
            coords, 4, order, 0.1)
        for order in permutations((0, 1, 2, 3))
    ]

    assert all(measure.sign in {-1, 1} for measure in measures)
    assert len({round(abs(measure.normalized), 14)
                for measure in measures}) == 1
    canonical_sign = measures[0].sign
    for order, measure in zip(permutations((0, 1, 2, 3)), measures):
        inversions = sum(
            order[left] > order[right]
            for left in range(4) for right in range(left + 1, 4))
        assert measure.sign == canonical_sign * (-1 if inversions % 2 else 1)


def test_event_relation_uses_lossless_sparse_baseline_encoding():
    import pynauty

    elements = ["C"] * 4
    wbo_R = np.zeros((4, 4))
    wbo_P = np.zeros((4, 4))
    wbo_P[0, 1] = wbo_P[1, 0] = 0.8
    identity = {index: index for index in range(4)}

    relation, _persistent, _inverse, _fragments = _masked_relation_data(
        identity, {}, elements, wbo_R, elements, wbo_P,
        1.0, 0.2, 0.5, 0.3, None)
    graph = relation.graph("B")
    _generators, _m1, _m2, orbits, _count = pynauty.autgrp(graph)

    # Six complete atom pairs have two event colors.  Absence represents the
    # five-pair baseline exactly, leaving one exceptional relation vertex.
    assert graph.number_of_vertices == 5
    assert orbits[0] == orbits[1]
    assert orbits[2] == orbits[3]
    assert orbits[0] != orbits[2]


def test_actual_event_quotient_recovers_same_event_threshold_cosets():
    elements = ["C"] * 3
    wbo_R = np.zeros((3, 3))
    wbo_P = np.zeros((3, 3))
    for matrix in (wbo_R, wbo_P):
        for (left, right), value in {
                (0, 1): 0.6, (0, 2): 1.0, (1, 2): 0.2}.items():
            matrix[left, right] = matrix[right, left] = value
    identity = {atom: atom for atom in range(3)}
    hierarchy = {"fragments": [{"fragment": [0, 1, 2]}]}

    family = compile_analytical_mapping_family(
        identity, hierarchy,
        elements, wbo_R, elements, wbo_P,
        graph_floor=0.1, symmetry_wbo_tol=0.5)
    representatives = family.exact_event_coset_representatives(
        wbo_R, wbo_P, elements)

    # The conservative all-threshold relation separates the two assignments,
    # while the actual selected event (no broken/formed bonds) does not.  The
    # quotient must retain both analytical cosets without enumerating either
    # group's individual bijections.
    assert family.group_order == (1.0, 0)
    assert family.structural_group_order == (2.0, 0)
    assert len(representatives) == 2
    assert all(mapping_event_signature(
        mapping, wbo_R, wbo_P, elements) == ((), ())
        for mapping in representatives)


def test_cached_event_relation_is_identical_to_uncached_construction():
    rng = np.random.default_rng(73)
    elements = ["C", "C", "O", "H", "H", "H"]
    wbo_R = rng.uniform(0.0, 1.2, size=(6, 6))
    wbo_R = (wbo_R + wbo_R.T) / 2.0
    np.fill_diagonal(wbo_R, 0.0)
    wbo_P = wbo_R.copy()
    mapping = {index: index for index in range(6)}
    hierarchy = {"fragments": [
        {"fragment": [0, 1, 2]},
        {"fragment": [3, 4, 5]},
    ]}
    uncached = _masked_relation_data(
        mapping, hierarchy, elements, wbo_R, elements, wbo_P,
        0.2, 0.2, 0.5, 0.3, {})[0]
    context = analytical_family_static_context(
        elements, wbo_R, elements, wbo_P,
        graph_floor=0.2, dwbo_threshold=0.5,
        metal_dwbo_threshold=0.3)
    cached = _masked_relation_data(
        mapping, hierarchy, elements, wbo_R, elements, wbo_P,
        0.2, 0.2, 0.5, 0.3, {}, static_context=context)[0]

    assert cached.colors_A == uncached.colors_A
    assert cached.colors_B == uncached.colors_B
    assert cached.relation_records_A == uncached.relation_records_A
    assert cached.relation_records_B == uncached.relation_records_B
