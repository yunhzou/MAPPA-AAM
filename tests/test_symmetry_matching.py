import numpy as np
import pytest

from mappa import classify_bonds
from mappa import build_graph

from mappa.mappa_core.core.search.planning import _generate_seed_orders
from mappa.mappa_core.core.search.queries.partial import _branch_symmetry_record
from mappa.mappa_workflow.mapping_analytics.family import project_path_assignments
import mappa.mappa_core.core.search.branching.construction as branch_mod
from mappa.mappa_core.core.matching.symmetry.canonical import _PartialMappingCanonicalizer
from mappa.mappa_core.core.matching.representation.state import _cand_possible_p_atoms
from mappa.mappa_core.core.matching.symmetry.orbits import _nauty_orbits
from mappa.mappa_core.core.matching.symmetry.orbits import _wbo_tolerance_bucket_lookup
from mappa.mappa_core.core.matching.representation.records import _IsoResult
from mappa.mappa_core.core.matching.representation.records import IslandBranchLimitExceeded
from mappa.mappa_core.core.matching.growth.engine import grow_island
from mappa.mappa_core.core.matching.growth.existence import FragmentQuery
from mappa.mappa_core.core.matching.policy import DEFAULT_NODE_POLICY
from test_growth_compression_exactness import complete, explicit, represented


def _tetramethyl_metal_graph():
    elements = ["Pd"]
    edges = []
    for i in range(4):
        c = 1 + 4 * i
        elements.append("C")
        edges.append((0, c, 0.8))
        for j in range(3):
            h = c + 1 + j
            elements.append("H")
            edges.append((c, h, 1.0))

    wbo = np.zeros((len(elements), len(elements)))
    for i, j, w in edges:
        wbo[i, j] = wbo[j, i] = w
    return elements, wbo


def test_fragment_branch_cap_is_loud_and_allows_exact_limit(monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', '0')
    g_r = build_graph(["C"], np.zeros((1, 1)), bond_cut=0.2)
    wbo_two = np.zeros((4, 4))
    wbo_two[0, 2] = wbo_two[2, 0] = 1.0
    wbo_two[1, 3] = wbo_two[3, 1] = 1.0
    g_p_two = build_graph(["C", "C", "O", "N"], wbo_two, bond_cut=0.2)
    wbo_three = np.zeros((6, 6))
    for carbon, ligand in ((0, 3), (1, 4), (2, 5)):
        wbo_three[carbon, ligand] = wbo_three[ligand, carbon] = 1.0
    g_p_three = build_graph(
        ["C", "C", "C", "O", "N", "F"], wbo_three, bond_cut=0.2)

    assert len(grow_island(g_r, g_p_two, 0, {}, max_branches=2)) == 2
    with pytest.raises(IslandBranchLimitExceeded) as exc_info:
        grow_island(g_r, g_p_three, 0, {}, max_branches=2)
    assert exc_info.value.count == 3
    assert exc_info.value.limit == 2


def test_completed_fragment_enforces_cap_on_distinct_choices():
    g_r = build_graph(
        ["Si", "N"], np.array([[0.0, 1.0], [1.0, 0.0]]), bond_cut=0.2)
    wbo = np.zeros((9, 9))
    for nitrogen, marker in ((1, 5), (2, 6), (3, 7), (4, 8)):
        wbo[0, nitrogen] = wbo[nitrogen, 0] = 1.0
        wbo[nitrogen, marker] = wbo[marker, nitrogen] = 1.0
    g_p = build_graph(
        ["Si", "N", "N", "N", "N", "O", "F", "Cl", "Br"],
        wbo,
        bond_cut=0.2,
    )

    with pytest.raises(IslandBranchLimitExceeded) as exc_info:
        grow_island(
            g_r,
            g_p,
            0,
            {},
            max_branches=2,
            p_orbits=_nauty_orbits(g_p),
        )

    assert exc_info.value.count == 4
    assert exc_info.value.limit == 2


def test_partial_mapping_certificate_requires_exact_endpoint_transporter():
    source = build_graph(["C"], np.zeros((1, 1)), bond_cut=0.2)
    target_wbo = np.zeros((3, 3))
    target_wbo[0, 1] = target_wbo[1, 0] = 1.0
    target_wbo[1, 2] = target_wbo[2, 1] = 1.0
    target = build_graph(["C", "C", "C"], target_wbo, bond_cut=0.2)
    canonicalizer = _PartialMappingCanonicalizer(
        source, target, wbo_tol=0.5)

    left_end = canonicalizer.certificate({0: 0})
    right_end = canonicalizer.certificate({0: 2})
    center = canonicalizer.certificate({0: 1})

    assert left_end == right_end
    assert left_end != center
    assert canonicalizer.equivalent({0: 0}, {0: 2})
    assert not canonicalizer.equivalent({0: 0}, {0: 1})


def test_heavy_only_cut_items_keep_hydrogens_but_do_not_cut_xh_edges():
    from mappa.mappa_core.core.search.planning import cut_sweep_items

    elements = ["C", "C", "H", "H"]
    wbo = np.zeros((4, 4))
    for left, right in ((0, 1), (0, 2), (1, 3)):
        wbo[left, right] = wbo[right, left] = 1.0

    assert cut_sweep_items(
        wbo, elements=elements, heavy_only=True) == [(), ((0, 1),)]


def test_live_branch_cap_discards_only_overflowing_parent_subtree(monkeypatch):
    g_r = build_graph(["C", "O"], np.zeros((2, 2)), bond_cut=0.2)
    g_p = build_graph(
        ["C", "C", "O", "O", "O"], np.zeros((5, 5)), bond_cut=0.2)

    def fake_grow(_g_r, _g_p, seed, mapping, **_kwargs):
        if seed == 0:
            return [
                _IsoResult({0: 0}, fragment={0}),
                _IsoResult({0: 1}, fragment={0}),
            ]
        if mapping[0] == 0:
            # This parent alone would create three children under a cap of
            # two, so its subtree must be removed atomically.
            return [
                _IsoResult({1: p}, fragment={1}) for p in (2, 3, 4)
            ]
        return [_IsoResult({1: 2}, fragment={1})]

    import mappa.mappa_core.core.matching.api as fragment_mod
    monkeypatch.setattr(fragment_mod, "grow_island", fake_grow)

    branches = branch_mod.find_islands(
        g_r, g_p, [0, 1], max_branches=2)

    assert [branch.mapping for branch in branches.paths()] == [{0: 1, 1: 2}]


def test_existence_reassigns_a_correlated_witness():
    wbo_r = np.zeros((3, 3))
    wbo_r[0, 2] = wbo_r[2, 0] = 1.
    source = build_graph(['O', 'O', 'C'], wbo_r, bond_cut=.2)
    wbo_p = np.zeros((4, 4))
    wbo_p[1, 3] = wbo_p[3, 1] = 1.
    target = build_graph(['O', 'O', 'C', 'C'], wbo_p, bond_cut=.2)
    query = FragmentQuery(source, target, {0, 1, 2}, [(0, 2)], {},
                          DEFAULT_NODE_POLICY, .1)
    assert query.find(preferred={0: 0, 1: 1}) == {0: 1, 1: 0, 2: 3}


def test_branch_symmetry_record_closes_open_pool_with_mapping_owner():
    from mappa.mappa_core.core.matching.representation.records import FragmentPlacement
    branch = branch_mod._Branch()
    branch.commit(
        FragmentPlacement.from_match(_IsoResult(
            {27: 27},
            fragment={27},
            symmetry={
                "witness": {27: 27},
                "blocks": [{
                    "r_atoms": [27],
                    "p_atoms": [26, 27],
                    "extendable": False,
                    "open": True,
                    "assignments": "2",
                }],
            },
        )),
    )
    branch.mapping[22] = 26
    branch.islands_R[22] = 2

    record = _branch_symmetry_record(branch)

    assert record["blocks"][0]["r_atoms"] == [22, 27]
    assert record["blocks"][0]["p_atoms"] == [26, 27]
    assert record["blocks"][0]["open"] is False
    assert record["blocks"][0]["assignments"] == "2!"


def test_completion_retains_placements_on_different_atom_sets():
    source = build_graph(['C', 'C'], np.array([[0., 1.], [1., 0.]]), bond_cut=.2)
    weights = np.zeros((4, 4))
    weights[0, 2] = weights[2, 0] = weights[1, 3] = weights[3, 1] = 1.
    target = build_graph(['C'] * 4, weights, bond_cut=.2)
    candidates = complete(source, target, {0, 1}, [(0, 1)], .1)
    assert len(candidates) == 2
    assert {frozenset(candidate.mapping.values()) for candidate in candidates} == {
        frozenset({0, 2}), frozenset({1, 3})}
    assert represented(candidates, target, .1) == {
        ((0, a), (1, b)) for a, b in ((0, 2), (2, 0), (1, 3), (3, 1))}


def test_completion_ignores_unaccepted_source_pairs():
    weights = np.zeros((3, 3))
    weights[1, 2] = weights[2, 1] = 1.
    weights[0, 2] = weights[2, 0] = .4
    source = build_graph(['C', 'C', 'O'], weights.copy(), bond_cut=.5)
    weights[0, 2] = weights[2, 0] = 0.
    target = build_graph(['C', 'C', 'O'], weights, bond_cut=.5)
    result = complete(source, target, {0, 1, 2}, [(1, 2)], .1,
                      locked={0: 0, 1: 1})
    assert len(result) == 1 and result[0].mapping == {0: 0, 1: 1, 2: 2}


def test_completion_rejects_an_accepted_edge_mismatch():
    weights = np.zeros((3, 3))
    weights[1, 2] = weights[2, 1] = 1.
    weights[0, 2] = weights[2, 0] = .6
    source = build_graph(['C', 'C', 'O'], weights.copy(), bond_cut=.5)
    weights[0, 2] = weights[2, 0] = 0.
    target = build_graph(['C', 'C', 'O'], weights, bond_cut=.5)
    assert complete(source, target, {0, 1, 2}, [(0, 2)], .1,
                    locked={0: 0, 1: 1}) == []


def test_unique_center_has_one_compressed_hydrogen_pool():
    source = build_graph(['C', 'H'], np.array([[0., 1.], [1., 0.]]), bond_cut=.2)
    weights = np.array([[0., 1., 1.], [1., 0., 0.], [1., 0., 0.]])
    target = build_graph(['C', 'H', 'H'], weights, bond_cut=.2)
    result = complete(source, target, {0, 1}, [(0, 1)], .1)
    assert len(result) == 1
    assert [(b.r_atoms, b.p_atoms) for b in result[0].blocks] == [((1,), (1, 2))]
    assert represented(result, target, .1) == {((0, 0), (1, p)) for p in (1, 2)}


def test_planar_metal_tetramethyl_metal_carbon_branch_is_fourfold():
    elements, weights = _tetramethyl_metal_graph()
    source = build_graph(elements, weights, bond_cut=.2)
    target = build_graph(elements, weights, bond_cut=.2)
    result = complete(source, target, {0, 1}, [(0, 1)], .5)
    assert len(result) == 1
    assert represented(result, target, .5) == {((0, 0), (1, p)) for p in (1, 5, 9, 13)}


def test_core_alignment_expands_internal_branch_degeneracy():
    elements, wbo = _tetramethyl_metal_graph()

    from mappa import MolecularEndpoint, AAMSearchConfig, search_core_assignments
    endpoint = MolecularEndpoint(tuple(elements), np.zeros((len(elements), 3)), wbo)
    result = search_core_assignments(endpoint, endpoint, [1, 5],
        config=AAMSearchConfig(seed_count=1, branch_limit=100, iso_tolerance=.5))
    mappings = {assignment.pairs for assignment in result.assignments}
    expected_targets = {1, 5, 9, 13}

    assert len(mappings) == 12
    assert mappings == {
        ((1, a), (5, b))
        for a in expected_targets
        for b in expected_targets
        if a != b
    }


def test_core_projection_does_not_invent_global_orbit_swaps():
    class Branch:
        mapping = {0: 0, 1: 1}
        transitions = ()

    from mappa import AAMProblem, MolecularEndpoint
    wbo = np.zeros((2, 2))
    endpoint = MolecularEndpoint(("H", "H"), np.zeros((2, 3)), wbo)

    variants = project_path_assignments(
        Branch(), AAMProblem(endpoint, endpoint), [0, 1], assignment_limit=100)

    assert variants == [{0: 0, 1: 1}]


def test_core_projection_includes_tolerance_pool_without_exact_automorphism():
    from mappa import AAMProblem, AAMSearchConfig, MolecularEndpoint, search_core_assignments

    weights = np.array([[0., 1., 1.1], [1., 0., 0.], [1.1, 0., 0.]])
    endpoint = MolecularEndpoint(('C', 'H', 'H'), np.zeros((3, 3)), weights)
    result = search_core_assignments(endpoint, endpoint, [1, 2],
        config=AAMSearchConfig(seed_count=1, iso_tolerance=.2))
    assert {assignment.pairs for assignment in result.assignments} == {
        ((1, 1), (2, 2)), ((1, 2), (2, 1))}
    path = result.branches[0].search_path
    assert all(g.images[1] == 1 and g.images[2] == 2
               for fragment in path.hierarchy.fragments for g in fragment.target_generators)
    with pytest.raises(branch_mod.BranchLimitExceeded) as exc_info:
        project_path_assignments(path, AAMProblem(endpoint, endpoint), [1, 2],
                                 assignment_limit=1)
    assert exc_info.value.branch_count == 2


def test_core_alignment_does_not_filter_nonmatching_core_edges():
    el_r = ["N", "C", "O"]
    wbo_r = np.zeros((3, 3))
    wbo_r[0, 1] = wbo_r[1, 0] = 0.79
    wbo_r[1, 2] = wbo_r[2, 1] = 1.10

    el_t = ["N", "C", "O"]
    wbo_t = np.zeros((3, 3))
    wbo_t[0, 1] = wbo_t[1, 0] = 0.0
    wbo_t[1, 2] = wbo_t[2, 1] = 1.67

    from mappa import MolecularEndpoint, AAMSearchConfig, search_core_assignments
    result = search_core_assignments(
        MolecularEndpoint(tuple(el_r), np.zeros((3, 3)), wbo_r),
        MolecularEndpoint(tuple(el_t), np.zeros((3, 3)), wbo_t), [0, 1, 2],
        config=AAMSearchConfig(seed_count=1, branch_limit=100, iso_tolerance=.1))
    assert {assignment.pairs for assignment in result.assignments} == {((0, 0), (1, 1), (2, 2))}


def test_nauty_orbits_group_near_wbo_by_tolerance():
    pytest.importorskip("pynauty")
    wbo = np.zeros((3, 3))
    wbo[0, 1] = wbo[1, 0] = 1.0
    wbo[0, 2] = wbo[2, 0] = 1.1
    g = build_graph(["C", "H", "H"], wbo, bond_cut=0.2)

    loose = _nauty_orbits(g, wbo_tol=0.2)
    strict = _nauty_orbits(g, wbo_tol=0.05)

    assert loose[1] == loose[2]
    assert loose[0] != loose[1]
    assert strict[1] != strict[2]


def test_nauty_loose_tolerance_preserves_edges_and_methyl_orbit():
    pytest.importorskip("pynauty")
    wbo = np.zeros((4, 4))
    for h, value in enumerate((0.982056, 0.980196, 0.960315), start=1):
        wbo[0, h] = wbo[h, 0] = value
    g = build_graph(["C", "H", "H", "H"], wbo, bond_cut=0.2)

    buckets, zero_bucket = _wbo_tolerance_bucket_lookup(g, 1.0)
    orbits = _nauty_orbits(g, wbo_tol=1.0)

    assert all(buckets[(0, h)] != zero_bucket for h in (1, 2, 3))
    assert buckets[(1, 2)] == zero_bucket
    assert orbits[1] == orbits[2] == orbits[3]
    assert orbits[0] != orbits[1]


def test_nauty_masks_positive_wbo_below_graph_floor():
    pytest.importorskip("pynauty")
    wbo = np.zeros((5, 5))
    wbo[0, 1] = wbo[1, 0] = 1.0
    wbo[0, 2] = wbo[2, 0] = 1.0
    wbo[1, 3] = wbo[3, 1] = 1.0
    wbo[2, 4] = wbo[4, 2] = 1.0
    wbo[1, 4] = wbo[4, 1] = 0.16
    graph = build_graph(["B", "O", "O", "H", "H"], wbo, bond_cut=0.2)

    buckets, zero_bucket = _wbo_tolerance_bucket_lookup(graph, 1.0)
    orbits = _nauty_orbits(graph, wbo_tol=1.0)

    assert buckets[(1, 4)] == zero_bucket
    assert orbits[1] == orbits[2]
    assert orbits[3] == orbits[4]


def test_tolerance_orbits_supply_a_completed_hydrogen_pool():
    weights = np.array([[0., 1., 1.1], [1., 0., 0.], [1.1, 0., 0.]])
    source = build_graph(['C', 'H', 'H'], weights, bond_cut=.2)
    target = build_graph(['C', 'H', 'H'], weights, bond_cut=.2)
    result = complete(source, target, {0, 1}, [(0, 1)], .2)
    assert len(result) == 1
    assert [(b.r_atoms, b.p_atoms) for b in result[0].blocks] == [((1,), (1, 2))]
    assert represented(result, target, .2) == {((0, 0), (1, p)) for p in (1, 2)}


def test_planar_metal_tetramethyl_child_h_keeps_parent_degeneracy():
    elements, weights = _tetramethyl_metal_graph()
    source = build_graph(elements, weights, bond_cut=.2)
    target = build_graph(elements, weights, bond_cut=.2)
    result = complete(source, target, {0, 1, 2}, [(0, 1), (1, 2)], .5)
    assert len(result) == 1
    expected = {((0, 0), (1, c), (2, h))
                for c in (1, 5, 9, 13) for h in (c + 1, c + 2, c + 3)}
    assert len(expected) == 12 and represented(result, target, .5) == expected


def test_planar_metal_tetramethyl_two_metal_carbon_edges_are_correlated():
    elements, weights = _tetramethyl_metal_graph()
    source = build_graph(elements, weights, bond_cut=.2)
    target = build_graph(elements, weights, bond_cut=.2)
    result = complete(source, target, {0, 1, 5}, [(0, 1), (0, 5)], .5)
    assert len(result) == 1
    expected = {((0, 0), (1, a), (5, b))
                for a in (1, 5, 9, 13) for b in (1, 5, 9, 13) if a != b}
    assert len(expected) == 12 and represented(result, target, .5) == expected


def test_completion_uses_later_constraints_to_select_the_valid_assignment():
    weights = np.zeros((3, 3))
    weights[0, 2] = weights[2, 0] = 1.
    source = build_graph(['C', 'C', 'O'], weights, bond_cut=.2)
    weights = np.zeros((3, 3))
    weights[1, 2] = weights[2, 1] = 1.
    target = build_graph(['C', 'C', 'O'], weights, bond_cut=.2)
    unconstrained = complete(source, target, {0, 1}, [], .5)
    assert len(unconstrained) == 2
    assert represented(unconstrained, target, .5) == {
        ((0, 0), (1, 1)), ((0, 1), (1, 0))}
    result = complete(source, target, {0, 1, 2}, [(0, 2)], .5)
    assert len(result) == 1 and result[0].mapping == {0: 1, 1: 0, 2: 2}


def test_completion_preserves_pair_relations_within_one_vertex_orbit():
    from mappa.mappa_core.core.matching.representation.state import _SymCand
    source = build_graph(['C', 'C'], np.array([[0., 1.], [1., 0.]]), bond_cut=.2)
    weights = np.zeros((4, 4))
    for a, b in [(0, 1), (1, 2), (2, 3), (3, 0)]:
        weights[a, b] = weights[b, a] = 1.
    target = build_graph(['C'] * 4, weights, bond_cut=.2)
    adjacent = {((0, a), (1, b)) for a, b in target.edges} | {
        ((0, b), (1, a)) for a, b in target.edges}
    opposite = {((0, a), (1, (a + 2) % 4)) for a in target}
    # One vertex orbit does not make independent endpoint shuffling valid.
    assert represented([_SymCand({0: 0, 1: 1})], target, .2) == adjacent
    assert represented([_SymCand({0: 0, 1: 2})], target, .2) == opposite
    assert adjacent.isdisjoint(opposite)
    result = complete(source, target, {0, 1}, [(0, 1)], .2)
    assert len(result) == 1 and represented(result, target, .2) == adjacent


def test_completion_preserves_correlated_carbon_and_hydrogen_motion():
    elements = ['X', 'C', 'C', 'H', 'H']
    weights = np.zeros((5, 5))
    bonds = [(0, 1), (0, 2), (1, 3), (2, 4)]
    for a, b in bonds:
        weights[a, b] = weights[b, a] = 1.
    graph = build_graph(elements, weights, bond_cut=.2)
    result = complete(graph, graph, set(graph), bonds, .5)
    assert len(result) == 1
    assert represented(result, graph, .5) == {
        ((0, 0), (1, 1), (2, 2), (3, 3), (4, 4)),
        ((0, 0), (1, 2), (2, 1), (3, 4), (4, 3))}


def test_completed_weighted_constraints_preserve_every_tolerated_assignment():
    source_weights = np.zeros((4, 4))
    target_weights = np.zeros((4, 4))
    for atom, weight in enumerate((1.907, 1.869, 1.138), start=1):
        source_weights[0, atom] = source_weights[atom, 0] = weight
    for atom, weight in enumerate((1.910, 1.875, 1.127), start=1):
        target_weights[0, atom] = target_weights[atom, 0] = weight
    source = build_graph(['W', 'O', 'O', 'O'], source_weights, bond_cut=.2)
    target = build_graph(['W', 'O', 'O', 'O'], target_weights, bond_cut=.2)
    bonds = [(0, 1), (0, 2), (0, 3)]
    result = complete(source, target, set(source), bonds, 1.,
                      target_orbits={0: 0, 1: 1, 2: 1, 3: 2})
    expected = explicit(source, target, set(source), 1., bonds)
    assert len(expected) == 6 and represented(result, target, 1.) == expected


def test_locked_anchor_restricts_its_terminal_atom_to_the_same_component():
    source = build_graph(['C', 'H'], np.array([[0., 1.], [1., 0.]]), bond_cut=.2)
    weights = np.zeros((4, 4))
    weights[0, 1] = weights[1, 0] = weights[2, 3] = weights[3, 2] = 1.
    target = build_graph(['C', 'H', 'C', 'H'], weights, bond_cut=.2)
    result = complete(source, target, {0, 1}, [(0, 1)], .5, locked={0: 0})
    assert len(result) == 1 and result[0].mapping == {0: 0, 1: 1}
    assert _cand_possible_p_atoms(result[0]) == {0, 1}


def test_active_source_edge_requires_active_target_edge_even_with_loose_tolerance():
    source = build_graph(['C', 'H'], np.array([[0., .964], [.964, 0.]]), bond_cut=.2)
    weights = np.zeros((3, 3))
    weights[0, 1] = weights[1, 0] = .963
    target = build_graph(['C', 'H', 'H'], weights, bond_cut=.2)
    result = complete(source, target, {0, 1}, [(0, 1)], 1.)
    assert len(result) == 1
    assert represented(result, target, 1.) == {((0, 0), (1, 1))}


def test_subfloor_boundary_weights_remain_in_the_global_action():
    from mappa.mappa_core.core.matching.representation.state import _SymCand
    from mappa.mappa_core.core.matching.symmetry.workspace import SymmetryWorkspace
    weights = np.zeros((3, 3))
    weights[0, 1] = weights[1, 0] = 1.
    weights[0, 2] = weights[2, 0] = .4
    source = build_graph(['C', 'C', 'O'], weights, bond_cut=.5)
    target = build_graph(['C', 'C', 'O'], weights, bond_cut=.5)
    generators = SymmetryWorkspace(target, .5).conditioned_generators(_SymCand({0: 0}), ())
    assert all(np.array_equal(weights, weights[np.ix_(action, action)]) for action in generators)
    result = complete(source, target, {0, 1}, [(0, 1)], .5,
                      target_orbits={0: 0, 1: 0, 2: 1})
    relations = [represented([candidate], target, .5) for candidate in result]
    assert set.union(*relations) == {((0, 0), (1, 1)), ((0, 1), (1, 0))}
    assert sum(map(len, relations)) == 2


def test_classify_bonds_uses_lower_metal_event_threshold():
    wbo_r = np.zeros((2, 2))
    wbo_p = np.zeros((2, 2))
    wbo_p[0, 1] = wbo_p[1, 0] = 0.35
    mapping = {0: 0, 1: 1}

    organic = classify_bonds(
        mapping, wbo_r, wbo_p,
        dwbo_threshold=0.5,
        elements_R=["C", "O"], elements_P=["C", "O"],
        metal_dwbo_threshold=0.3)
    metal = classify_bonds(
        mapping, wbo_r, wbo_p,
        dwbo_threshold=0.5,
        elements_R=["Pd", "O"], elements_P=["Pd", "O"],
        metal_dwbo_threshold=0.3)

    assert organic[1] == []
    assert metal[1] == [(0, 1, 0.0, 0.35)]


def test_mechanism_certificate_distinguishes_pair_orbits_not_vertex_orbits():
    # A distance-colored K5 is vertex-transitive, but cycle edges and
    # diagonals are different pair orbits.  Endpoint-orbit IDs alone conflate
    # these two one-bond events.
    n = 5
    elements = ["C"] * n
    wbo_r = np.full((n, n), 0.75)
    np.fill_diagonal(wbo_r, 0.0)
    for atom in range(n):
        neighbor = (atom + 1) % n
        wbo_r[atom, neighbor] = wbo_r[neighbor, atom] = 1.5
    wbo_cycle = wbo_r.copy()
    wbo_cycle[0, 1] = wbo_cycle[1, 0] = 0.0
    wbo_diagonal = wbo_r.copy()
    wbo_diagonal[0, 2] = wbo_diagonal[2, 0] = 0.0
    graph_r = build_graph(elements, wbo_r, bond_cut=0.2)
    r_orbits = _nauty_orbits(graph_r, wbo_tol=0.2)
    assert len(set(r_orbits.values())) == 1
    assert (r_orbits[0], r_orbits[1]) == (r_orbits[0], r_orbits[2])
    mapping = {atom: atom for atom in range(n)}

    from mappa.mappa_core.domain import AAMProblem, MolecularEndpoint
    from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex

    def signature(product_wbo):
        coordinates = np.zeros((n, 3))
        problem = AAMProblem(MolecularEndpoint(tuple(elements), coordinates, wbo_r),
                             MolecularEndpoint(tuple(elements), coordinates, product_wbo))
        index = SignedEventIndex(problem)
        return index.describe(index.vector(mapping))

    cycle_signature = signature(wbo_cycle)
    diagonal_signature = signature(wbo_diagonal)
    assert cycle_signature['total'] == diagonal_signature['total'] == 1
    assert cycle_signature['id'] != diagonal_signature['id']


def test_generate_seed_orders_honors_trial_cap():
    elements = ["C", "C", "C", "C", "H", "H"]
    wbo = np.zeros((6, 6))
    for i, j in [(0, 1), (1, 2), (2, 3), (0, 4), (3, 5)]:
        wbo[i, j] = wbo[j, i] = 1.0
    g = build_graph(elements, wbo, bond_cut=0.2)

    orders = _generate_seed_orders(g, n_trials=2)

    assert len(orders) == 2
    assert all(g.degree[order[0]] > 0 for order in orders)


def test_generate_seed_orders_deprioritizes_common_isolated_atoms():
    elements = ["C", "C", "C", "H", "H", "H", "H"]
    wbo = np.zeros((7, 7))
    for i, j in [(0, 1), (1, 2), (2, 3)]:
        wbo[i, j] = wbo[j, i] = 1.0
    g = build_graph(elements, wbo, bond_cut=0.2)

    order = _generate_seed_orders(g, n_trials=1)[0]

    connected = {0, 1, 2, 3}
    isolated_common = {4, 5, 6}
    assert max(order.index(atom) for atom in connected) < min(
        order.index(atom) for atom in isolated_common)
    assert set(order) == set(range(7))


def test_core_atoms_do_not_reorder_seed_sequence(monkeypatch):
    elements = ["C", "C", "O"]
    wbo = np.zeros((3, 3))
    wbo[0, 1] = wbo[1, 0] = 1.0
    wbo[1, 2] = wbo[2, 1] = 1.0
    g = build_graph(elements, wbo, bond_cut=0.2)
    seen = []

    def fake_grow_island(_g_R, _g_P, seed, *_args, **_kwargs):
        seen.append(seed)
        return []

    import mappa.mappa_core.core.matching.api as fragment_mod
    monkeypatch.setattr(fragment_mod, "grow_island", fake_grow_island)

    branch_mod.find_islands(
        g, g, [0, 1, 2],
        core_R=[2],
        stop_when_core_mapped=True,
        p_orbits={0: 0, 1: 1, 2: 2},
    )

    assert seen == [0, 1, 2]


def test_find_islands_reuses_precomputed_orbits(monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', '0')
    wbo = np.zeros((2, 2))
    wbo[0, 1] = wbo[1, 0] = 1.0
    g_r = build_graph(["C", "H"], wbo, bond_cut=0.2)
    g_p = build_graph(["C", "H"], wbo, bond_cut=0.2)

    def fail_if_called(*args, **kwargs):
        raise AssertionError("orbits should be supplied by caller")

    monkeypatch.setattr(branch_mod, "_nauty_orbits", fail_if_called)
    branches = branch_mod.find_islands(
        g_r, g_p, [0, 1],
        p_orbits={0: 0, 1: 1},
    )

    assert branches
    assert next(branches.paths()).mapping == {0: 0, 1: 1}


def test_removed_overlap_scorers_are_not_public():
    import mappa
    import mappa.mappa_core.utils.io.modes
    assert not hasattr(mappa, 'bond_overlap_per_mode')
    assert not hasattr(mappa.mappa_core.utils.io.modes, 'bond_overlap_per_mode')
    assert not hasattr(mappa, 'rxn_overlap_per_mode')


def test_anchored_noop_seed_does_not_keep_pass_loop_alive(monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', '0')
    elements = ["C", "C", "C"]
    wbo_r = np.zeros((3, 3))
    wbo_r[0, 1] = wbo_r[1, 0] = 1.0
    wbo_r[1, 2] = wbo_r[2, 1] = 1.0

    wbo_p = np.zeros((3, 3))
    wbo_p[1, 2] = wbo_p[2, 1] = 1.0

    g_r = build_graph(elements, wbo_r, bond_cut=0.2)
    g_p = build_graph(elements, wbo_p, bond_cut=0.2)
    events = []

    branches = branch_mod.find_islands(
        g_r, g_p, [1, 0],
        iso_tol=1.0,
        max_branches=100,
        anchor_map={0: 0},
        events=events,
    )

    assert branches
    assert all(branch.mapping[0] == 0 for branch in branches.paths())
    assert [e["seed"] for e in events if e.get("type") == "seed_start"][:2] == [1, 0]
    assert sum(1 for e in events if e.get("type") == "pass_start") == 2


def test_approximate_orbits_are_not_an_exact_symmetry_backend():
    import mappa.mappa_core.core.matching as matcher
    import mappa.mappa_core.core.matching.symmetry.orbits as orbits
    assert not hasattr(matcher, '_color_refine_orbits')
    assert not hasattr(orbits, '_color_refine_orbits')
