"""Completed constraints must preserve the independently enumerated relation."""
from itertools import permutations

import networkx as nx
import numpy as np
import pytest

from mappa.mappa_core.core.matching.graphs import build_graph
from mappa.mappa_core.core.matching.representation.state import _SymBlock, _SymCand
from mappa.mappa_core.core.matching.policy import DEFAULT_NODE_POLICY
from mappa.mappa_core.core.matching.symmetry.orbits import _nauty_orbits
from growth_cases import represented


def graph(elements, edges):
    weights = np.zeros((len(elements), len(elements)))
    for a, b, weight in edges:
        weights[a, b] = weights[b, a] = weight
    return build_graph(elements, weights, bond_cut=0.2)


def test_historical_v2_roles_equal_direct_full_matrix_scan():
    from mappa.mappa_core.core.matching.symmetry.canonical import _CandidateAutomorphismCanonicalizer
    rng = np.random.default_rng(7321)
    nodes = (2, 5, 9, 12, 19)
    for _ in range(12):
        weights = rng.choice([0., .1, .35, 1., 1.1], size=(20, 20))
        weights = np.triu(weights, 1)
        weights += weights.T
        target = build_graph(['C'] * 20, weights, bond_cut=.2).subgraph(nodes).copy()
        canonicalizer = _CandidateAutomorphismCanonicalizer(target, wbo_tol=.5, compression_version=2)
        for images in permutations(nodes, 3):
            for pooled in (False, True):
                blocks = (_SymBlock((1,), (images[1], images[2])),) if pooled else ()
                candidate = _SymCand({0: images[0], 1: images[1]}, blocks)
                expected = canonicalizer._candidate_roles(candidate, group_domains=True)
                pool_atoms = {p for block in blocks for p in block.p_atoms}
                outside = set(nodes) - (set(candidate.mapping.values()) - pool_atoms)
                for p in nodes:
                    boundary = tuple((q, float(weights[p, q])) for q in sorted(outside)
                                     if p != q and weights[p, q] != 0)
                    role = ('outside', p) if p in outside else ('boundary', boundary)
                    expected[p] = (*expected.get(p, ()), role)
                assert canonicalizer._growth_roles(candidate, group_domains=True) == expected


def explicit(source, target, atoms, tolerance, edges=None):
    result = set()
    for images in permutations(target.nodes, len(atoms)):
        mapping = dict(zip(atoms, images))
        if any(source.nodes[r]['element'] != target.nodes[p]['element']
               for r, p in mapping.items()):
            continue
        if all(target.has_edge(mapping[a], mapping[b]) and
               abs(source.graph['wbo_matrix'][a, b] -
                   target.graph['wbo_matrix'][mapping[a], mapping[b]]) <= tolerance
               for a, b in (source.edges if edges is None else edges) if a in mapping and b in mapping):
            result.add(tuple(sorted(mapping.items())))
    return result


def complete(source, target, fragment, bonds, tolerance=.2, *, locked=None,
             target_orbits=None, bond_tolerance=None):
    """Run completion on an explicit constraint set, including impossible ones."""
    from mappa.mappa_core.core.matching.growth.completion import FragmentCompletion
    target_orbits = (_nauty_orbits(target, wbo_tol=tolerance)
                     if target_orbits is None else target_orbits)
    return FragmentCompletion(source, target, fragment, bonds, locked or {},
        target_orbits, tolerance, DEFAULT_NODE_POLICY, bond_tolerance).run()


def assert_completed_constraints(source, target, order, tolerance=.2):
    fragment = {order[0]}
    edges = []
    candidates = complete(source, target, fragment, edges, tolerance)
    assert represented(candidates, target, tolerance) == explicit(source, target, fragment, tolerance, edges)
    for atom in order[1:]:
        anchor = min(set(source.neighbors(atom)) & fragment)
        edges.append((anchor, atom))
        fragment.add(atom)
        candidates = complete(source, target, fragment, edges, tolerance)
        relations = [represented([candidate], target, tolerance) for candidate in candidates]
        expected = explicit(source, target, fragment, tolerance, edges)
        assert (set().union(*relations)) == expected, (order, atom)
        assert sum(map(len, relations)) == len(expected), (order, atom)
    return candidates


def test_growth_retains_both_identical_components():
    source = graph(['O', 'H', 'C'], [(0, 1, 1), (0, 2, 1)])
    target = graph(['O', 'H', 'C'] * 2,
                   [(0, 1, 1), (0, 2, 1), (3, 4, 1), (3, 5, 1)])
    final = assert_completed_constraints(source, target, [0, 1, 2])
    assert len(final) == 2


SMALL_GRAPHS = [g for g in nx.graph_atlas_g() if 1 <= len(g) <= 5 and nx.is_connected(g)]


@pytest.mark.parametrize('topology', SMALL_GRAPHS, ids=lambda g: f'n{len(g)}-e{g.number_of_edges()}')
def test_each_completed_constraint_set_matches_explicit_isomorphism(topology):
    target = graph(['C'] * len(topology), [(a, b, 1) for a, b in topology.edges])
    for seed in (0, len(topology) - 1):
        order = list(nx.bfs_tree(topology, seed))
        assert_completed_constraints(target, target, order)


def test_methyl_permutations_stay_compressed():
    methane = graph(['C', 'H', 'H', 'H', 'H'], [(0, p, 1) for p in range(1, 5)])
    final = assert_completed_constraints(methane, methane, [0, 1, 2, 3, 4])
    assert len(final) == 1
    assert len(final[0].blocks) == 1
    assert len(represented(final, methane, 0.2)) == 24


def test_bucket_equivalence_does_not_hide_different_future_support():
    source = graph(['C', 'C', 'O'], [(0, 1, 1), (0, 2, 1.2)])
    target = graph(['C', 'C', 'O'], [(0, 1, 1), (0, 2, 1.0), (1, 2, 1.1)])
    assert_completed_constraints(source, target, [0, 1, 2], tolerance=0.15)


def test_locked_island_merge_constrains_only_join_edge():
    source = graph(['X', 'C', 'O', 'N'],
                   [(0, 1, 1), (1, 2, 1), (2, 3, 1), (1, 3, 1)])
    target = graph(['X', 'C', 'C', 'O', 'N'],
                   [(0, 1, 1), (0, 2, 1), (1, 3, 1), (2, 3, 1),
                    (3, 4, 1), (2, 4, 1)])
    result = complete(source, target, {0, 1, 2, 3}, [(0, 1), (1, 2)],
                      locked={2: 3, 3: 4})
    assert represented(result, target, .2) == {((0, 0), (1, p), (2, 3), (3, 4)) for p in (1, 2)}


def test_unaccepted_bonds_do_not_restrict_completed_pools():
    source = graph(['X', 'C', 'C', 'C'],
                   [(0, 1, 1), (0, 2, 1), (0, 3, 1), (1, 3, 1)])
    target = graph(['X', 'C', 'C', 'C'],
                   [(0, 1, 1), (0, 2, 1), (0, 3, 1), (1, 2, 1), (2, 3, 1)])
    bonds = [(0, 1), (0, 2), (0, 3)]
    result = complete(source, target, set(source), bonds)
    assert represented(result, target, .2) == explicit(source, target, set(source), .2, bonds)


def test_weighted_constraints_preserve_exact_tolerance_checks():
    # A bucket automorphism need not preserve the source-target tolerance test.
    weighted = graph(['C'] * 4, [(0, 1, 1.2), (0, 2, .8), (0, 3, 1.2),
                                (1, 2, 1.2), (1, 3, .8), (2, 3, 1.5)])
    assert_completed_constraints(weighted, weighted, [0, 1, 2, 3], tolerance=.5)


@pytest.mark.parametrize('tolerance', [.2, .5])
def test_weighted_source_target_pairs_against_enumeration(tolerance):
    random = np.random.default_rng(1204)
    for _ in range(80):
        graphs = []
        for _side in range(2):
            weights = np.triu(random.choice([0, .8, 1, 1.2, 1.5], (4, 4)), 1)
            edges = [(a, b, weights[a, b]) for a in range(4) for b in range(a + 1, 4)
                     if weights[a, b]]
            graphs.append(graph(['C'] * 4, edges))
        if nx.is_connected(graphs[0]):
            assert_completed_constraints(*graphs, list(nx.bfs_tree(graphs[0], 0)), tolerance)




def test_standalone_fragment_matcher_returns_both_component_placements(monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', '0')
    from mappa import match_fragment, FragmentMatchConfig
    source = graph(['O', 'H', 'C'], [(0, 1, 1), (0, 2, 1)])
    target = graph(['O', 'H', 'C'] * 2,
                   [(0, 1, 1), (0, 2, 1), (3, 4, 1), (3, 5, 1)])
    result = match_fragment(source, target, seed=0,
                            config=FragmentMatchConfig(iso_tolerance=.2), events=[])
    assert not result.capped
    assert {tuple(sorted(match.items())) for match in result.matches} == {
        ((0, 0), (1, 1), (2, 2)), ((0, 3), (1, 4), (2, 5))}


def test_compression_contract_survives_typed_hierarchy_roundtrip():
    from mappa.mappa_core.domain import AAMHierarchy
    from mappa.mappa_core.core.matching.representation.state import _symmetry_state
    from mappa.mappa_core.core.matching.representation.state import candidate_from_record
    candidate = _SymCand({0: 0})
    record = {'fragments': [{'fragment': [0], 'symmetry': _symmetry_state(candidate)}]}
    hierarchy = AAMHierarchy.from_record(record)
    restored = hierarchy.to_record()['fragments'][0]['symmetry']
    assert restored['compression_version'] == 3
    assert candidate_from_record(restored).compression_version == 3
    del restored['compression_version']
    assert candidate_from_record(restored).compression_version == 1


@pytest.mark.parametrize('native', ['0', '1'])
@pytest.mark.parametrize('reverse', [False, True])
def test_ring_chain_keeps_growth_subgraph_through_public_consumers(native, reverse, monkeypatch):
    from mappa import AAMProblem, MolecularEndpoint, decode_events, query_path
    from mappa.mappa_core.core.search.branching.construction import find_islands
    from mappa.mappa_core.core.search.branching.symmetry import finalize_graph_symmetry
    from mappa.mappa_core.core.search.graph import AAMSearchGraph
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    import json

    monkeypatch.setenv('MAPPA_NATIVE', native)
    ring = graph(['C'] * 6, [(a, (a + 1) % 6, 1) for a in range(6)])
    chain = graph(['C'] * 6, [(a, a + 1, 1) for a in range(5)])
    source, target = (chain, ring) if reverse else (ring, chain)
    endpoint = lambda g: MolecularEndpoint(('C',) * 6, np.zeros((6, 3)), g.graph['wbo_matrix'])
    problem = AAMProblem(endpoint(source), endpoint(target))
    found = find_islands(source, target, list(range(6)), iso_tol=1, max_branches=100)
    found, _ = finalize_graph_symmetry(found, target, iso_tolerance=1)
    restored = AAMSearchGraph.from_record(json.loads(json.dumps(found.to_record())))
    catalogue = FinalBranchCatalogue(problem).add_graph(restored)
    assert catalogue.families and not found.capped
    for path in restored.paths():
        assert len(path.hierarchy.fragments) == 1
        fragment = path.hierarchy.fragments[0]
        assert set(fragment.r_atoms) == set(range(6))
        assert len(fragment.preserved_bonds) == 5
        tree = nx.Graph(fragment.preserved_bonds)
        assert nx.is_tree(tree)
        assert query_path(path, problem, path.mapping, source_atoms=tuple(range(6)))[0] == 'recovered'
        from mappa.mappa_core.core.search.branching.validation import RepresentativeValidationWorkspace
        RepresentativeValidationWorkspace(problem).validate_path(path)
    events = decode_events(catalogue)
    for i, family in enumerate(catalogue.families):
        family.validate_representative(problem)
        assert len(family.required_edges) == 5
        result = events.contributions(i)
        assert result['total'] == 1
        assert result['counts']['forming' if reverse else 'breaking'] == 1
        assert sum(sum(v.values()) for v in result['boundary'].values()) == 0


def test_equal_atom_sets_keep_distinct_growth_constraints():
    from dataclasses import replace
    from mappa import AAMProblem, MolecularEndpoint
    from mappa.mappa_core.core.search.branching.construction import find_islands
    from mappa.mappa_core.core.search.branching.symmetry import finalize_graph_symmetry
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    ring = graph(['C'] * 3, [(0, 1, 1), (1, 2, 1), (0, 2, 1)])
    endpoint = MolecularEndpoint(('C',) * 3, np.zeros((3, 3)), ring.graph['wbo_matrix'])
    found = find_islands(ring, ring, [0, 1, 2], iso_tol=.2, max_branches=100)
    found, _ = finalize_graph_symmetry(found, ring, iso_tolerance=.2)
    changed = replace(found, transitions=tuple(replace(e, preserved_bonds=((0, 1), (1, 2)))
                                               for e in found.transitions))
    other = replace(found, transitions=tuple(replace(e, preserved_bonds=((0, 1), (0, 2)))
                                             for e in found.transitions))
    cat = FinalBranchCatalogue(AAMProblem(endpoint, endpoint)).add_graph(changed).add_graph(other)
    assert len({f.required_edges for f in cat.families}) == 2
