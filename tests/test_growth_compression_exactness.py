"""Compression must preserve the fragment relation at every extension."""
from itertools import permutations, product
from collections import defaultdict

import networkx as nx
import numpy as np
import pytest
from sympy.combinatorics import Permutation, PermutationGroup

from mappa.frag import build_graph
from mappa.matcher import _SymBlock, _SymCand, _dedup_sym_cands, _extend_sym_cands, _nauty_orbits
from mappa.search_symmetry import SymmetryWorkspace


def graph(elements, edges):
    weights = np.zeros((len(elements), len(elements)))
    for a, b, weight in edges:
        weights[a, b] = weights[b, a] = weight
    return build_graph(elements, weights, bond_cut=0.2)


def represented(candidates, target, tolerance):
    """Enumerate tiny test groups independently of the production decoder."""
    workspace = SymmetryWorkspace(target, tolerance)
    result = set()
    for candidate in candidates:
        atoms = sorted(candidate.mapping)
        generators = workspace.conditioned_generators(candidate, ())
        group = PermutationGroup([Permutation(list(g)) for g in generators]
                                 or [Permutation(list(range(len(target))))])
        choices = [list(permutations(block.p_atoms, len(block.r_atoms)))
                   for block in candidate.blocks]
        for selection in product(*choices):
            mapping = dict(candidate.mapping)
            for block, images in zip(candidate.blocks, selection):
                mapping.update(zip(block.r_atoms, images))
            for action in group.generate_schreier_sims():
                result.add(tuple((r, action(mapping[r])) for r in atoms))
    return result


def explicit(source, target, atoms, tolerance):
    result = set()
    for images in permutations(target.nodes, len(atoms)):
        mapping = dict(zip(atoms, images))
        if any(source.nodes[r]['element'] != target.nodes[p]['element']
               for r, p in mapping.items()):
            continue
        if all(target.has_edge(mapping[a], mapping[b]) and
               abs(source.graph['wbo_matrix'][a, b] -
                   target.graph['wbo_matrix'][mapping[a], mapping[b]]) <= tolerance
               for a, b in source.edges if a in mapping and b in mapping):
            result.add(tuple(sorted(mapping.items())))
    return result


def assert_prefixes(source, target, order, tolerance=0.2):
    ro = _nauty_orbits(source, wbo_tol=tolerance)
    po = _nauty_orbits(target, wbo_tol=tolerance)
    seed = order[0]
    pools = defaultdict(list)
    for p in target:
        if source.nodes[seed]['element'] == target.nodes[p]['element']:
            pools[po[p]].append(p)
    candidates = [_SymCand({seed: values[0]},
                           (_SymBlock((seed,), tuple(values), extendable=False),)
                           if len(values) > 1 else ())
                  for values in pools.values()]
    candidates = _dedup_sym_cands(candidates, source, target, ro, po)
    fragment = {seed}
    assert represented(candidates, target, tolerance) == explicit(source, target, fragment, tolerance)
    for atom in order[1:]:
        candidates = _extend_sym_cands(candidates, fragment, atom, source, target,
                                       {}, tolerance, None, r_orbits=ro, p_orbits=po)
        fragment.add(atom)
        assert represented(candidates, target, tolerance) == explicit(source, target, fragment, tolerance), (order, atom)
    return candidates


def test_growth_retains_both_identical_components():
    source = graph(['O', 'H', 'C'], [(0, 1, 1), (0, 2, 1)])
    target = graph(['O', 'H', 'C'] * 2,
                   [(0, 1, 1), (0, 2, 1), (3, 4, 1), (3, 5, 1)])
    final = assert_prefixes(source, target, [0, 1, 2])
    assert len(final) == 2


SMALL_GRAPHS = [g for g in nx.graph_atlas_g() if 1 <= len(g) <= 5 and nx.is_connected(g)]


@pytest.mark.parametrize('topology', SMALL_GRAPHS, ids=lambda g: f'n{len(g)}-e{g.number_of_edges()}')
def test_every_prefix_matches_explicit_isomorphism(topology):
    target = graph(['C'] * len(topology), [(a, b, 1) for a, b in topology.edges])
    for seed in (0, len(topology) - 1):
        order = list(nx.bfs_tree(topology, seed))
        assert_prefixes(target, target, order)


def test_methyl_permutations_stay_compressed():
    methane = graph(['C', 'H', 'H', 'H', 'H'], [(0, p, 1) for p in range(1, 5)])
    final = assert_prefixes(methane, methane, [0, 1, 2, 3, 4])
    assert len(final) == 1
    assert len(final[0].blocks) == 1
    assert len(represented(final, methane, 0.2)) == 24


def test_bucket_equivalence_does_not_hide_different_future_support():
    source = graph(['C', 'C', 'O'], [(0, 1, 1), (0, 2, 1.2)])
    target = graph(['C', 'C', 'O'], [(0, 1, 1), (0, 2, 1.0), (1, 2, 1.1)])
    assert_prefixes(source, target, [0, 1, 2], tolerance=0.15)


def test_locked_island_merge_preserves_pool_support_for_every_added_atom():
    source = graph(['X', 'C', 'O', 'N'],
                   [(0, 1, 1), (1, 2, 1), (2, 3, 1), (1, 3, 1)])
    target = graph(['X', 'C', 'C', 'O', 'N'],
                   [(0, 1, 1), (0, 2, 1), (1, 3, 1), (2, 3, 1),
                    (3, 4, 1), (2, 4, 1)])
    candidate = _SymCand({0: 0, 1: 1}, (_SymBlock((1,), (1, 2)),))
    result = _extend_sym_cands([candidate], {0, 1}, 2, source, target,
                               {2: 3, 3: 4}, .2, {2: 1, 3: 1})
    assert represented(result, target, .2) == {((0, 0), (1, 2), (2, 3), (3, 4))}


def test_pool_join_requires_every_new_atom_choice_to_be_interchangeable():
    source = graph(['X', 'C', 'C', 'C'],
                   [(0, 1, 1), (0, 2, 1), (0, 3, 1), (1, 3, 1)])
    target = graph(['X', 'C', 'C', 'C'],
                   [(0, 1, 1), (0, 2, 1), (0, 3, 1), (1, 2, 1), (2, 3, 1)])
    candidate = _SymCand({0: 0, 1: 1, 2: 2}, (_SymBlock((1, 2), (1, 2, 3)),))
    assert represented([candidate], target, .2) == explicit(source, target, {0, 1, 2}, .2)
    result = _extend_sym_cands([candidate], {0, 1, 2}, 3, source, target, {}, .2, None)
    assert represented(result, target, .2) == explicit(source, target, {0, 1, 2, 3}, .2)


def test_weighted_prefixes_preserve_exact_tolerance_checks():
    # A bucket automorphism need not preserve the source-target tolerance test.
    weighted = graph(['C'] * 4, [(0, 1, 1.2), (0, 2, .8), (0, 3, 1.2),
                                (1, 2, 1.2), (1, 3, .8), (2, 3, 1.5)])
    assert_prefixes(weighted, weighted, [0, 1, 2, 3], tolerance=.5)


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
            assert_prefixes(*graphs, list(nx.bfs_tree(graphs[0], 0)), tolerance)


def test_support_limit_reports_incomplete_instead_of_no_match():
    from mappa.matcher.support import _support_assignments_for_value
    target = graph(['C', 'C', 'O'], [(1, 2, 1)])
    candidate = _SymCand({0: 0}, (_SymBlock((0,), (0, 1)),))
    with pytest.raises(RuntimeError, match='support state limit exceeded'):
        _support_assignments_for_value(candidate, 1, 2, [0], [(0, 1)],
                                       target, .2, max_states=1)


def test_standalone_fragment_matcher_returns_both_component_placements():
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
    from mappa.alignment.post_aam import AAMHierarchy
    from mappa.matcher import _symmetry_state
    from mappa.matcher.state import candidate_from_record
    candidate = _SymCand({0: 0})
    record = {'fragments': [{'fragment': [0], 'symmetry': _symmetry_state(candidate)}]}
    hierarchy = AAMHierarchy.from_record(record)
    restored = hierarchy.to_record()['fragments'][0]['symmetry']
    assert restored['compression_version'] == 2
    assert candidate_from_record(restored).compression_version == 2
    del restored['compression_version']
    assert candidate_from_record(restored).compression_version == 1
