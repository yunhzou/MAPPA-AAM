"""Completed constraints must preserve the independently enumerated relation."""
from itertools import permutations

import networkx as nx
import numpy as np
import pytest

from mappa.mappa_core.core.matching.graphs import build_graph
from mappa.mappa_core.core.matching.policy import DEFAULT_NODE_POLICY
from growth_cases import represented, discovery_mappings


def graph(elements, edges):
    weights = np.zeros((len(elements), len(elements)))
    for a, b, weight in edges:
        weights[a, b] = weights[b, a] = weight
    return build_graph(elements, weights, bond_cut=0.2)




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
             bond_tolerance=None, self_tolerance=1.):
    """Complete a proven explicit accepted graph; impossible queries return no result."""
    from mappa.mappa_core.core.matching.growth.completion import complete_fragment
    from mappa.mappa_core.core.matching.growth.existence import FragmentQuery
    from mappa.mappa_core.core.matching.representation.records import FragmentGrowth
    query = FragmentQuery(source, target, fragment, bonds, locked or {},
                          DEFAULT_NODE_POLICY, tolerance, bond_tolerance)
    witness = query.find()
    if witness is None:
        return []
    return complete_fragment(FragmentGrowth(set(fragment), witness, set(bonds)), source, target,
        locked or {}, tolerance=tolerance, self_tolerance=self_tolerance,
        node_policy=DEFAULT_NODE_POLICY, branch_limit=1000000, bond_tolerance=bond_tolerance)


def discovered(fragments):
    return set().union(*(discovery_mappings(fragment) for fragment in fragments))


def assert_completed_constraints(source, target, order, tolerance=.2):
    fragment = {order[0]}
    edges = []
    candidates = complete(source, target, fragment, edges, tolerance)
    assert discovered(candidates) == explicit(source, target, fragment, tolerance, edges)
    for atom in order[1:]:
        anchor = min(set(source.neighbors(atom)) & fragment)
        edges.append((anchor, atom))
        fragment.add(atom)
        candidates = complete(source, target, fragment, edges, tolerance)
        relations = [discovered([candidate]) for candidate in candidates]
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
    assert len(represented(final)) == 24


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
    assert represented(result) == {((0, 0), (1, p), (2, 3), (3, 4)) for p in (1, 2)}


def test_unaccepted_bonds_do_not_restrict_discovery():
    source = graph(['X', 'C', 'C', 'C'],
                   [(0, 1, 1), (0, 2, 1), (0, 3, 1), (1, 3, 1)])
    target = graph(['X', 'C', 'C', 'C'],
                   [(0, 1, 1), (0, 2, 1), (0, 3, 1), (1, 2, 1), (2, 3, 1)])
    bonds = [(0, 1), (0, 2), (0, 3)]
    result = complete(source, target, set(source), bonds)
    assert discovered(result) == explicit(source, target, set(source), .2, bonds)


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
    assert {match.witness for match in result.matches} == {
        ((0, 0), (1, 1), (2, 2)), ((0, 3), (1, 4), (2, 5))}


def test_compression_contract_survives_typed_hierarchy_roundtrip():
    from mappa.mappa_core.domain import AAMHierarchy
    source = graph(['C', 'C'], [(0, 1, 1)])
    fragments = complete(source, source, {0, 1}, [(0, 1)])
    hierarchy = AAMHierarchy(tuple(fragments))
    restored = AAMHierarchy.from_record(hierarchy.to_record())
    assert restored == hierarchy
    assert represented(restored.fragments) == {((0, 0), (1, 1)), ((0, 1), (1, 0))}


@pytest.mark.parametrize('native', ['0', '1'])
@pytest.mark.parametrize('reverse', [False, True])
def test_ring_chain_keeps_growth_subgraph_through_public_consumers(native, reverse, monkeypatch):
    from mappa import AAMProblem, MolecularEndpoint, decode_events, query_path
    from mappa.mappa_core.core.search.branching.construction import find_islands
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
    restored = AAMSearchGraph.from_record(json.loads(json.dumps(found.to_record())))
    catalogue = FinalBranchCatalogue(problem).add_graph(restored)
    assert catalogue.families and not found.capped
    for path in restored.paths():
        assert len(path.hierarchy.fragments) == 1
        fragment = path.hierarchy.fragments[0]
        assert set(fragment.source_atoms) == set(range(6))
        assert len(fragment.accepted_edges) == 5
        tree = nx.Graph(fragment.accepted_edges)
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
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    ring = graph(['C'] * 3, [(0, 1, 1), (1, 2, 1), (0, 2, 1)])
    endpoint = MolecularEndpoint(('C',) * 3, np.zeros((3, 3)), ring.graph['wbo_matrix'])
    found = find_islands(ring, ring, [0, 1, 2], iso_tol=.2, max_branches=100)
    from mappa.mappa_core.core.matching.representation.records import Fragment
    from mappa.mappa_core.core.symmetry.relations import DiscoveryEdge
    def with_edges(pairs):
        return replace(found, transitions=tuple(replace(edge,
            match=Fragment(replace(edge.match.relation,
                accepted_edges=tuple(DiscoveryEdge(pair, 1., .2) for pair in pairs))))
            for edge in found.transitions))
    changed = with_edges(((0, 1), (1, 2)))
    other = with_edges(((0, 1), (0, 2)))
    cat = FinalBranchCatalogue(AAMProblem(endpoint, endpoint)).add_graph(changed).add_graph(other)
    assert len(cat.families) == 2


def image_query(source, target, atoms, bonds, *, locked=None, tolerance=.2):
    from mappa.mappa_core.core.matching.growth.existence import FragmentQuery
    return FragmentQuery(source, target, atoms, bonds, locked or {},
                         DEFAULT_NODE_POLICY, tolerance)


def assert_distinct_images(source, target, atoms, bonds, *, locked=None, tolerance=.2,
                           preferred=None):
    """Compare image discovery with all valid assignments on tiny graphs."""
    locked = locked or {}
    occupied = set(locked.values())
    assignments = {
        pairs for pairs in explicit(source, target, tuple(sorted(atoms)), tolerance, bonds)
        if all(p == locked[r] if r in locked else p not in occupied for r, p in pairs)
    }
    expected = {frozenset(p for r, p in pairs) for pairs in assignments}
    query = image_query(source, target, atoms, bonds, locked=locked, tolerance=tolerance)
    witnesses = query.images(len(expected) + 1, preferred=preferred)
    images = [frozenset(witness.values()) for witness in witnesses]
    assert len(images) == len(set(images))
    assert set(images) == expected
    assert all(tuple(sorted(witness.items())) in assignments for witness in witnesses)
    return query, witnesses


def test_image_discovery_matches_independent_weighted_oracle():
    random = np.random.default_rng(6318)
    for _ in range(48):
        source = graph(['C', 'C', 'O'], [(a, b, random.choice([.8, 1., 1.2]))
            for a in range(3) for b in range(a + 1, 3) if random.random() < .7])
        target = graph(['C', 'C', 'O', 'C', 'O'], [(a, b, random.choice([.8, 1., 1.2]))
            for a in range(5) for b in range(a + 1, 5) if random.random() < .65])
        bonds = [edge for edge in source.edges if random.random() < .8]
        locked = {0: 0} if random.random() < .5 else {}
        assert_distinct_images(source, target, set(source), bonds,
                               locked=locked, tolerance=.15)


@pytest.mark.parametrize('copies', [1, 2])
def test_image_discovery_keeps_decorated_rings_as_one_choice_each(copies):
    ring = [(a, (a + 1) % 6, 1.) for a in range(6)]
    source = graph(['C'] * 6, ring)
    carbons = 6 * copies
    # Distinct external labels destroy whole-target symmetry, but cannot create
    # additional placements on the six atoms of either internal ring.
    target = graph(['C'] * carbons + [f'X{a}' for a in range(carbons)],
        [(a + 6 * i, b + 6 * i, weight) for i in range(copies) for a, b, weight in ring]
        + [(a, carbons + a, 1.) for a in range(carbons)])
    preferred = {a: carbons - 6 + a for a in range(6)}
    query = image_query(source, target, set(source), source.edges)
    witnesses = query.images(copies + 1, preferred=preferred)
    assert witnesses[0] == preferred
    assert {frozenset(witness.values()) for witness in witnesses} == {
        frozenset(range(6 * i, 6 * i + 6)) for i in range(copies)}
    assert len(witnesses) == copies
    assert query.nodes <= 4


def test_image_discovery_preserves_overlapping_sets_and_stops_at_limit():
    source = graph(['C'] * 3, [(0, 1, 1.), (1, 2, 1.)])
    target = graph(['C'] * 4, [(a, (a + 1) % 4, 1.) for a in range(4)])
    preferred = {0: 2, 1: 3, 2: 0}
    _, witnesses = assert_distinct_images(source, target, set(source), source.edges,
                                          preferred=preferred)
    assert len(witnesses) == 4
    assert witnesses[0] == preferred
    limited = image_query(source, target, set(source), source.edges)
    assert limited.images(3, preferred=preferred) == witnesses[:3]


def test_image_discovery_handles_disconnected_queries_and_external_locks():
    disconnected = graph(['C'] * 2, [])
    targets = graph(['C'] * 3, [])
    _, witnesses = assert_distinct_images(disconnected, targets, {0, 1}, [])
    assert len(witnesses) == 3
    source = graph(['C'] * 3, [(0, 1, 1.), (1, 2, 1.)])
    target = graph(['C'] * 4, [(a, a + 1, 1.) for a in range(3)])
    _, witnesses = assert_distinct_images(source, target, {1, 2}, [(1, 2)], locked={0: 0})
    assert {frozenset(witness.values()) for witness in witnesses} == {
        frozenset({1, 2}), frozenset({2, 3})}


def test_image_discovery_rejects_invalid_preferred_and_proves_absence():
    triangle = graph(['C'] * 3, [(0, 1, 1.), (1, 2, 1.), (0, 2, 1.)])
    chain = graph(['C'] * 3, [(0, 1, 1.), (1, 2, 1.)])
    query, witnesses = assert_distinct_images(triangle, chain, set(triangle), triangle.edges,
                                              preferred={0: 0, 1: 1, 2: 2})
    assert witnesses == []
    assert query.nodes == 1
    assert query.rejected == 1


def test_image_discovery_blocks_same_set_before_factorial_assignment_search():
    atoms = graph(['C'] * 12, [])
    witness = {a: a for a in atoms}
    query = image_query(atoms, atoms, set(atoms), [])
    assert query.images(2, preferred=witness) == [witness]
    assert query.nodes == 2


def repeated_arms():
    source = graph(['O', 'C', 'C'], [(0,1,1.), (1,2,1.)])
    target = graph(['O'] + ['C'] * 8,
        [(0,1+2*i,1.) for i in range(4)] + [(1+2*i,2+2*i,1.) for i in range(4)])
    return source, target


def test_global_image_orbits_share_discovery_without_losing_distinct_arms():
    from mappa.mappa_core.core.symmetry import relaxed_automorphism
    source, target = repeated_arms()
    prepared = relaxed_automorphism(['O'] + ['C'] * 8, target.graph['wbo_matrix'], tolerance=1.)
    query = image_query(source, target, tuple(source), tuple(source.edges))
    witnesses = query.images(5, global_self=prepared)
    assert {frozenset(w.values()) for w in witnesses} == {
        frozenset((0,1+2*i,2+2*i)) for i in range(4)}
    assert query.global_generators > 0
    assert query.global_images == query.global_pruned == 3
    plain = image_query(source, target, tuple(source), tuple(source.edges))
    assert {frozenset(w.values()) for w in plain.images(5)} == {frozenset(w.values()) for w in witnesses}
    assert query.nodes < plain.nodes


def test_global_image_orbits_stop_at_cap_evidence():
    from mappa.mappa_core.core.symmetry import relaxed_automorphism
    source = graph(['C'] * 3, [(0,1,1.), (1,2,1.)])
    target = graph(['C'] * 10, [(a,b,1.) for a in range(10) for b in range(a+1,10)])
    prepared = relaxed_automorphism(['C'] * 10, target.graph['wbo_matrix'], tolerance=1.)
    query = image_query(source, target, tuple(source), tuple(source.edges))
    witnesses = query.images(3, preferred={0:0,1:1,2:2}, global_self=prepared)
    assert len({frozenset(w.values()) for w in witnesses}) == len(witnesses) == 3
    assert query.nodes == 1 and query.global_images == 2


def weighted_star_self(tolerance):
    from dataclasses import replace
    from mappa.mappa_core.core.symmetry import relaxed_automorphism
    target = graph(['O','C','C','C'], [(0,1,1.), (0,2,1.4), (0,3,1.8)])
    prepared = relaxed_automorphism(['O','C','C','C'], target.graph['wbo_matrix'], tolerance=tolerance)
    # These two structural generators each exchange adjacent weights. Their
    # composition can move the endpoint weights farther than either exchange.
    return target, replace(prepared, generators=((0,2,1,3), (0,1,3,2)))


def test_global_pruning_preserves_the_actual_discovery_support_predicate():
    target, prepared = weighted_star_self(1.)
    source = graph(['O','C'], [(0,1,1.)])
    query = image_query(source, target, tuple(source), tuple(source.edges), tolerance=.5)
    witnesses = query.images(4, global_self=prepared)
    assert {frozenset(w.values()) for w in witnesses} == {frozenset((0,1)), frozenset((0,2))}
    assert query.global_generators == 1 and query.global_images == 1


def test_nontransitive_self_choices_cannot_be_used_as_an_unrestricted_group():
    target, prepared = weighted_star_self(.5)
    assert all(prepared.contains(generator) for generator in prepared.generators)
    assert not prepared.contains((0,3,2,1))
    source = graph(['O','C'], [(0,1,1.4)])
    query = image_query(source, target, tuple(source), tuple(source.edges), tolerance=1.)
    witnesses = query.images(4, global_self=prepared)
    assert {frozenset(w.values()) for w in witnesses} == {frozenset((0,i)) for i in (1,2,3)}
    assert query.global_generators == query.global_images == query.global_pruned == 0


def test_global_pruning_matches_independent_weighted_locked_image_oracle():
    import random
    from itertools import combinations
    from mappa.mappa_core.core.symmetry import relaxed_automorphism
    rng = random.Random(82351)
    transported = pruned = 0
    for trial in range(80):
        n, k = rng.choice((4,5,6)), rng.choice((2,3))
        source = graph(['C']*k, [(a,a+1,rng.choice((1.,1.3,1.6))) for a in range(k-1)])
        target = graph(['C']*n, [(a,b,rng.choice((0.,1.,1.3,1.6))) for a,b in combinations(range(n),2)])
        tolerance, self_tolerance = rng.choice((.2,.5,1.)), rng.choice((.2,.5,1.))
        locked = {0:0} if trial % 4 == 0 else ({k:n-1} if trial % 4 == 1 else {})
        valid = {pairs for pairs in explicit(source,target,range(k),tolerance) if all(
            (dict(pairs)[r] == p if r in dict(pairs) else p not in dict(pairs).values())
            for r,p in locked.items())}
        prepared = relaxed_automorphism(['C']*n,target.graph['wbo_matrix'],tolerance=self_tolerance)
        query = image_query(source,target,range(k),tuple(source.edges),locked=locked,tolerance=tolerance)
        preferred = dict(min(valid)) if valid and trial % 2 else None
        witnesses = query.images(100,preferred=preferred,global_self=prepared)
        actual = [frozenset(w.values()) for w in witnesses]
        assert len(actual) == len(set(actual))
        assert set(actual) == {frozenset(p for _,p in pairs) for pairs in valid}
        assert all(tuple(sorted(w.items())) in valid for w in witnesses)
        transported += query.global_images
        pruned += query.global_pruned
    assert transported > 0 and pruned > 0


def test_global_preparation_reuses_one_record_and_skips_complete_components(monkeypatch):
    from mappa import match_fragment, FragmentMatchConfig
    from mappa.mappa_core.core.matching.growth import completion
    monkeypatch.setenv('MAPPA_NATIVE','0')
    source, target = repeated_arms()
    original = completion._prepare_self
    preparations = []
    def prepare(graph_, atoms, *args):
        atoms = tuple(atoms)
        if graph_ is target and len(atoms) == len(target):
            preparations.append(atoms)
        return original(graph_,atoms,*args)
    monkeypatch.setattr(completion,'_prepare_self',prepare)
    for tolerance in (1.,1.,.2):
        result = match_fragment(source,target,seed=0,
            config=FragmentMatchConfig(iso_tolerance=.2,self_tolerance=tolerance))
        assert len(result.matches) == 4
    assert len(preparations) == 2
    ring = graph(['C']*6, [(a,(a+1)%6,1.) for a in range(6)])
    copies = graph(['C']*12, [(6*c+a,6*c+(a+1)%6,1.) for c in range(2) for a in range(6)])
    result = match_fragment(ring,copies,seed=0,config=FragmentMatchConfig(iso_tolerance=.2))
    assert len(result.matches) == 2 and '_relaxed_global' not in copies.graph


def test_public_completion_global_pruning_has_native_python_relation_parity(monkeypatch):
    from mappa import match_fragment, FragmentMatchConfig
    records = []
    for mode in ('0','1'):
        monkeypatch.setenv('MAPPA_NATIVE',mode)
        source,target = repeated_arms()
        profile = []
        result = match_fragment(source,target,seed=0,
            config=FragmentMatchConfig(iso_tolerance=.2),profile=profile)
        assert len(result.matches) == 4
        assert profile[0]['completion_global_images'] == 3
        records.append([fragment.to_record() for fragment in result.matches])
    assert records[0] == records[1]


def test_global_preparation_rechecks_copied_and_mutable_endpoint_weights(monkeypatch):
    from mappa import match_fragment, FragmentMatchConfig
    monkeypatch.setenv('MAPPA_NATIVE','0')
    source,target = repeated_arms()
    config = FragmentMatchConfig(iso_tolerance=.2)
    assert len(match_fragment(source,target,seed=0,config=config).matches) == 4
    original = target.graph['_relaxed_global'][3]
    copied = target.copy()
    copied.graph['wbo_matrix'] = target.graph['wbo_matrix'].copy()
    copied.graph['wbo_matrix'][0,1] = copied.graph['wbo_matrix'][1,0] = 1.9
    assert len(match_fragment(source,copied,seed=0,config=config).matches) == 3
    changed = copied.graph['_relaxed_global'][3]
    assert changed is not original and changed.weights[0,1] == 1.9
    assert target.graph['_relaxed_global'][3] is original
    copied.graph['wbo_matrix'][0,1] = copied.graph['wbo_matrix'][1,0] = 1.
    assert len(match_fragment(source,copied,seed=0,config=config).matches) == 4
    restored = copied.graph['_relaxed_global'][3]
    assert restored is not changed and restored.weights[0,1] == 1.


def test_global_preparation_supports_sparse_atom_indices_without_a_matrix():
    from mappa.mappa_core.core.matching.growth.completion import _global_self
    value = nx.Graph()
    value.add_nodes_from(((10,{'element':'O'}), (20,{'element':'C'}), (30,{'element':'C'})))
    value.add_weighted_edges_from(((10,20,1.), (10,30,1.)), weight='wbo')
    labels = {atom:value.nodes[atom]['element'] for atom in value}
    original = _global_self(value,labels,1.)
    assert original.atoms == (10,20,30)
    assert _global_self(value,labels,1.) is original
    value[10][30]['wbo'] = 1.9
    changed = _global_self(value,labels,1.)
    assert changed is not original and changed.weights[0,2] == 1.9


@pytest.mark.parametrize('backend', ['0', '1'], ids=['python', 'native'])
def test_attached_rings_keep_all_eighteen_target_sets_and_complete_discovery(backend, monkeypatch):
    from collections import defaultdict
    from mappa import FragmentMatchConfig, match_fragment

    monkeypatch.setenv('MAPPA_NATIVE', backend)
    source = graph(['C'] * 6, [(a, a + 1, 1.) for a in range(5)])
    target = graph(['C'] * 12, [(6*c + a, 6*c + (a+1) % 6, 1.)
                              for c in range(2) for a in range(6)] + [(0, 6, 1.)])
    oracle = nx.isomorphism.GraphMatcher(target, source,
        node_match=lambda left, right: left['element'] == right['element'])
    expected = defaultdict(set)
    for assignment in oracle.subgraph_monomorphisms_iter():
        expected[frozenset(assignment)].add(tuple(sorted((r, p) for p, r in assignment.items())))
    whole_rings = {frozenset(range(6)), frozenset(range(6, 12))}
    crossing = set(expected) - whole_rings
    assert len(expected) == 18 and whole_rings <= set(expected) and len(crossing) == 16
    assert all(image.intersection(range(6)) and image.intersection(range(6, 12)) for image in crossing)

    profile = []
    result = match_fragment(source, target, seed=0,
        config=FragmentMatchConfig(iso_tolerance=.2, branch_limit=18), profile=profile)
    assert profile[0]['completion_global_generators'] > 0
    assert not result.capped and result.branch_count == len(result.matches) == 18
    actual = {}
    for fragment in result.matches:
        image = frozenset(fragment.target_atoms)
        assert image not in actual
        assert fragment.source_atoms == tuple(range(6))
        assert set(fragment.accepted_edges) == set(source.edges)
        actual[image] = discovery_mappings(fragment)
    assert actual == expected


@pytest.mark.parametrize('backend', ['0', '1'], ids=['python', 'native'])
@pytest.mark.parametrize('tolerance, image_count', [(.2, 2), (1., 12)])
def test_binary_global_generators_keep_continuous_discovery_predicates(backend, tolerance, image_count, monkeypatch):
    from collections import defaultdict
    from mappa import FragmentMatchConfig, match_fragment
    from mappa.mappa_core.core.symmetry import relaxed_automorphism

    monkeypatch.setenv('MAPPA_NATIVE', backend)
    source = graph(['C'] * 5, [(a, a+1, 1.) for a in range(4)])
    target = graph(['C'] * 12, [(a, (a+1) % 12, 1. if a % 6 < 4 else 1.8)
                              for a in range(12)])
    weights = target.graph['wbo_matrix']
    global_self = relaxed_automorphism(['C'] * 12, weights, tolerance=1.)
    assert any(not np.array_equal(weights[np.ix_(g, g)], weights)
               for g in global_self.generators)
    # The same four accepted-edge predicates recur, but a broad structural
    # shuffle may change whether those predicates hold on the original WBOs.
    oracle = nx.isomorphism.GraphMatcher(target, source,
        node_match=lambda left, right: left['element'] == right['element'],
        edge_match=lambda left, right: abs(left['wbo'] - right['wbo']) <= tolerance)
    expected = defaultdict(set)
    for assignment in oracle.subgraph_monomorphisms_iter():
        expected[frozenset(assignment)].add(tuple(sorted((r, p) for p, r in assignment.items())))
    assert len(expected) == image_count
    result = match_fragment(source, target, seed=0,
        config=FragmentMatchConfig(iso_tolerance=tolerance, branch_limit=image_count))
    assert not result.capped and len(result.matches) == image_count
    actual = {}
    for fragment in result.matches:
        image = frozenset(fragment.target_atoms)
        assert image not in actual
        assert fragment.source_atoms == tuple(range(5))
        assert set(fragment.accepted_edges) == set(source.edges)
        actual[image] = discovery_mappings(fragment)
    assert actual == expected
    if backend == '1':
        monkeypatch.setenv('MAPPA_NATIVE', '0')
        reference = match_fragment(source, target, seed=0,
            config=FragmentMatchConfig(iso_tolerance=tolerance, branch_limit=image_count))
        assert [fragment.to_record() for fragment in result.matches] == [
            fragment.to_record() for fragment in reference.matches]
