"""Local relaxed relations preserve correlated branch choices and target occupancy."""
from itertools import product, permutations

import numpy as np
import pytest

from mappa import FragmentMatchConfig, match_fragment
from test_growth_compression_exactness import graph, represented, discovered, assert_completed_constraints, complete


from growth_cases import rings, ring_mappings


@pytest.mark.parametrize('order', [list(range(12)), [a for i in range(6) for a in (i, 6 + i)]])
def test_ring_constraints_preserve_all_768_matches(order):
    source, target = rings()
    fragment = {0}
    bonds = []
    for atom in order[1:]:
        bonds.append((min(set(source.neighbors(atom)) & fragment), atom))
        fragment.add(atom)
        candidates = complete(source, target, fragment, bonds, 1)
        expected = {tuple((r, p) for r, p in mapping if r in fragment) for mapping in ring_mappings()}
        assert discovered(candidates) == expected
    assert len(candidates.matches) == 1


@pytest.mark.parametrize('native', ['0', '1'])
@pytest.mark.smoke
def test_public_ring_match_compresses_64_atom_sets_with_768_assignments(native, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', native)
    source, target = rings()
    result = match_fragment(source, target, seed=0,
                            config=FragmentMatchConfig(iso_tolerance=1, branch_limit=100))
    assert not result.capped
    assert len(result.matches) == 1
    assert all(set(fragment.source_atoms) == set(source) for fragment in result.matches)
    assert represented(result) == ring_mappings()


def test_target_local_actions_fix_fragment_anchors():
    source, target = rings()
    fragments = complete(source, target, {0, 1, 6, 7}, [(0, 1), (0, 6), (1, 7)],
                         locked={0: 0})
    assert fragments.matches
    for fragment in fragments.matches:
        relation = fragment.target_self_automorph
        anchor = relation.atoms.index(0)
        assert all(action[anchor] == anchor for action in relation.generators)


def test_asymmetric_future_neighbors_keep_non_equivalent_matches():
    source = graph(['C', 'C', 'O'], [(0, 1, 1), (1, 2, 1)])
    target = graph(['C', 'C', 'C', 'O'], [(0, 1, 1), (1, 2, 1), (2, 3, 1)])
    assert_completed_constraints(source, target, [0, 1, 2])


def test_decoding_keeps_fragment_atom_sets_and_later_occupancy():
    from mappa import AAMProblem, query_path, decode_events
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    from test_relaxed_fragment_decoding import dense_events
    from mappa.mappa_core.core.search.branching.construction import find_islands
    from mappa.mappa_workflow.mapping_analytics.family import project_path_assignments
    from test_family_query import endpoint

    source = graph(['C'] * 4, [(0, 1, 1), (1, 2, 1), (2, 0, 1)])
    target = graph(['C'] * 4, [(i, (i + 1) % 4, 1) for i in range(4)])
    problem = AAMProblem(endpoint('CCCC', source.edges), endpoint('CCCC', target.edges))
    found = find_islands(source, target, [0, 1, 2, 3], iso_tol=.2, max_branches=100)
    paths = list(found.paths())
    assert len(paths) == 1
    decoded = decode_events(FinalBranchCatalogue(problem).add_graph(found))
    assert len(decoded.catalogue.families) == 1
    union = set()
    for path in paths:
        first, last = path.fragments
        expected = represented(path)
        union.update(expected)
        for images in permutations(range(4)):
            status, _ = query_path(path, problem, dict(enumerate(images)),
                                   source_atoms=tuple(range(4)), timeout_ms=2000)
            assert status == ('recovered' if tuple(enumerate(images)) in expected else 'not_recovered')
            if status == 'recovered':
                mapping = dict(enumerate(images))
                counts = decoded.tables(0).evaluate(decoded.index.vector(mapping), include_events=True)
                events = [(event['category'], event['reactant_pair'], event['product_pair'])
                          for event in counts['events']]
                assert sorted(events, key=repr) == dense_events(problem, mapping)
                assert counts['total'] == len(events)
                assert sum(sum(part.values()) for field in ('internal', 'boundary')
                           for part in counts[field].values()) == len(events)
        projected = {tuple(mapping.items()) for mapping in project_path_assignments(
            path, problem, (0, 3), assignment_limit=100)}
        assert projected == {tuple((r, p) for r, p in mapping if r in (0, 3)) for mapping in expected}
    assert union == {tuple(enumerate(images)) for images in permutations(range(4))}


def test_locked_ring_then_free_hydrogen_preserves_all_joint_choices():
    import z3
    from mappa import AAMProblem
    from mappa.mappa_core.core.search.branching.construction import find_islands
    from mappa.mappa_workflow.mapping_analytics.family import compile_path
    from test_family_query import endpoint

    cycle = [(i, (i + 1) % 3, 1) for i in range(3)]
    source = graph(['C'] * 3 + ['H'] * 4, cycle + [(i, 3 + i, 1) for i in range(3)])
    target = graph(['C'] * 3 + ['H'] * 6,
                   cycle + [(i, 3 + 2 * i + h, 1) for i in range(3) for h in range(2)])
    problem = AAMProblem(endpoint('CCC' + 'H' * 4, source.edges),
                         endpoint('CCC' + 'H' * 6, target.edges))
    found = find_islands(source, target, list(range(7)), iso_tol=1, max_branches=100)
    paths = list(found.paths())
    assert len(paths) == 3  # Remaining H positions differ once prior ring/H targets are locked.
    expected = set()
    for carbons in permutations(range(3)):
        for choices in product(range(2), repeat=3):
            images = tuple(carbons) + tuple(3 + 2 * p + h for p, h in zip(carbons, choices))
            for free in set(range(3, 9)) - set(images):
                expected.add(tuple(enumerate(images + (free,))))
    assert len(expected) == 144
    decoded = set()
    for path in paths:
        first, last = path.fragments
        assert set(first.source_atoms) == set(range(6)) and last.source_atoms == (6,)
        assert not set(first.target_atoms).intersection(last.target_atoms)
        compiled = compile_path(path, problem, {}, source_atoms=tuple(range(7)), complete_reference=False)
        compiled.solver.set(timeout=3000)
        path_mappings = set()
        while True:
            status = compiled.solver.check()
            assert status != z3.unknown, compiled.solver.reason_unknown()
            if status == z3.unsat:
                break
            mapping = tuple(map(tuple, compiled.realize(compiled.solver.model())['mapping']))
            assert mapping in expected and mapping not in path_mappings
            path_mappings.add(mapping)
            decoded.add(mapping)
            compiled.solver.add(z3.Or(*(value[0] != p
                for value, (_r, p) in zip(compiled.values, mapping))))
    assert decoded == expected


@pytest.mark.parametrize('native', ['0', '1'])
def test_hydrogen_caps_count_placement_orbits_not_individual_atoms(native, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', native)
    source = graph(['H'], [])
    equivalent = graph(['C', 'H', 'H', 'H'], [(0, p, 1.) for p in (1, 2, 3)])
    result = match_fragment(source, equivalent, seed=0,
                            config=FragmentMatchConfig(branch_limit=1))
    assert not result.capped and result.branch_count == len(result.matches) == 1
    assert represented(result) == {((0, p),) for p in (1, 2, 3)}

    different = graph(['C', 'H', 'H', 'N', 'H'],
                      [(0, 1, 1.), (0, 2, 1.), (3, 4, 1.)])
    limited = match_fragment(source, different, seed=0,
                             config=FragmentMatchConfig(branch_limit=1))
    assert limited.capped and limited.branch_count == 2 and not limited.matches
    complete_result = match_fragment(source, different, seed=0,
                                    config=FragmentMatchConfig(branch_limit=2))
    assert not complete_result.capped and len(complete_result.matches) == 2
    assert represented(complete_result) == {((0, p),) for p in (1, 2, 4)}


@pytest.mark.parametrize('native', ['0', '1'])
def test_placement_group_fixes_prior_locks_outside_current_fragment(native, monkeypatch):
    from mappa import FragmentMatchContext
    from growth_cases import group_actions
    monkeypatch.setenv('MAPPA_NATIVE', native)
    source = graph(['H', 'C'], [])
    target = graph(['C', 'H', 'C', 'H'], [(0, 1, 1.), (2, 3, 1.)])
    result = match_fragment(source, target, seed=0,
        context=FragmentMatchContext(locked_mapping={1: 0}, islands={1: 1}),
        config=FragmentMatchConfig(branch_limit=2))
    assert not result.capped and len(result.matches) == 2
    assert all(action[0] == 0 for action in group_actions(result.placement_generators, 4))
    assert represented(result) == {((0, 1),), ((0, 3),)}


@pytest.mark.parametrize('native', ['0', '1'])
@pytest.mark.parametrize('contact', [False, True], ids=['active-weight', 'subthreshold-contact'])
def test_placement_group_preserves_original_weights_not_only_binary_structure(native, contact, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', native)
    source = graph(['C', 'H'], [(0, 1, 1.)])
    elements = ['C', 'H', 'C', 'H'] + (['O'] if contact else [])
    edges = [(0, 1, 1.), (2, 3, 1. if contact else 1.3)]
    if contact:
        edges += [(0, 4, .05), (2, 4, .15)]
    target = graph(elements, edges)
    result = match_fragment(source, target, seed=0,
        config=FragmentMatchConfig(iso_tolerance=.5, self_tolerance=1., branch_limit=2))
    assert not result.capped and len(result.matches) == 2
    assert represented(result) == {((0, 0), (1, 1)), ((0, 2), (1, 3))}


def test_full_discovery_on_one_image_is_not_replaced_by_a_witness_orbit():
    source = graph(['C'] * 4, [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.), (0, 3, 1.)])
    result = complete(source, source, tuple(source), ((0, 1), (1, 2), (0, 3)),
                      tolerance=0., self_tolerance=0.)
    assert len(result.matches) == 1
    assert len(discovered(result)) == 4
    assert len(represented(result)) == 6


@pytest.mark.parametrize('native', ['0', '1'])
def test_global_placement_does_not_close_nontransitive_local_self_choices(native, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', native)
    star = graph(['O', 'C', 'C', 'C'], [(0, 1, 1.), (0, 2, 1.6), (0, 3, 2.2)])
    result = match_fragment(star, star, seed=0,
        config=FragmentMatchConfig(iso_tolerance=0., self_tolerance=.7, branch_limit=1))
    assert not result.capped and len(result.matches) == 1
    assert not result.placement_generators
    choices = represented(result)
    assert len(choices) == 5
    assert ((0, 0), (1, 3), (2, 2), (3, 1)) not in choices


def test_set_orbit_uses_correlated_motion_and_true_pointwise_stabilizers():
    from itertools import combinations
    from mappa import _group_ops
    from growth_cases import group_actions
    symmetric = ((1, 0, 2), (0, 2, 1))
    stabilizer = _group_ops.SetOrbit(3, symmetric).stabilizer((1,))
    assert set(group_actions(stabilizer, 3)) == {(0, 1, 2), (2, 1, 0)}

    correlated = ((1, 0, 3, 2),)
    orbit = _group_ops.SetOrbit(4, correlated)
    assert orbit.contains((0, 2), (1, 3))
    assert not orbit.contains((0, 2), (0, 3))
    for size in range(5):
        for first in combinations(range(4), size):
            images = {frozenset(action[p] for p in first)
                      for action in group_actions(correlated, 4)}
            for second in combinations(range(4), size):
                assert orbit.contains(first, second) == (frozenset(second) in images)


def test_placement_group_preserves_compatibility_of_future_source_atoms(monkeypatch):
    from mappa.mappa_core.core.matching.policy import CallableNodeMatchPolicy
    monkeypatch.setenv('MAPPA_NATIVE', '0')
    source = graph(['C', 'C'], [])
    target = graph(['C', 'C'], [])
    policy = CallableNodeMatchPolicy(lambda _s, r, _t, p: r == 0 or p == 1)
    result = match_fragment(source, target, seed=0,
        config=FragmentMatchConfig(branch_limit=2, node_policy=policy))
    # Both targets fit the current source atom. Only one fits the next atom,
    # so exchanging them would incorrectly identify different continuations.
    assert len(result.matches) == 2 and not result.placement_generators
    assert represented(result) == {((0, 0),), ((0, 1),)}


@pytest.mark.parametrize('generators', [
    ((1, 2, 3, 0),),
    ((1, 2, 3, 0), (0, 3, 2, 1)),
    ((1, 2, 3, 0), (1, 0, 2, 3)),
    ((1, 0, 3, 2),),
])
def test_set_orbit_and_stabilizer_equal_independent_small_group_oracle(generators):
    from itertools import combinations
    from mappa import _group_ops
    from growth_cases import group_actions
    actions = group_actions(generators, 4)
    orbit = _group_ops.SetOrbit(4, generators)
    for size in range(5):
        subsets = tuple(combinations(range(4), size))
        for first in subsets:
            images = {frozenset(action[p] for p in first) for action in actions}
            for second in subsets:
                assert orbit.contains(first, second) == (frozenset(second) in images)
            actual = set(group_actions(orbit.stabilizer(first), 4))
            assert actual == {action for action in actions if all(action[p] == p for p in first)}


@pytest.mark.parametrize('native', ['0', '1'])
def test_placement_group_preserves_active_cut_topology(native, monkeypatch):
    from growth_cases import group_actions
    monkeypatch.setenv('MAPPA_NATIVE', native)
    source = graph(['C'], [])
    target = graph(['C'] * 4, [(p, (p + 1) % 4, 1.) for p in range(4)])
    target.remove_edge(0, 1)
    result = match_fragment(source, target, seed=0,
        config=FragmentMatchConfig(branch_limit=4))
    # Original weights still describe a ring. Search topology is a path, whose
    # endpoints and middle atoms cannot share one future continuation.
    assert len(result.matches) >= 2
    actions = group_actions(result.placement_generators, 4)
    assert all({action[0], action[1]} == {0, 1} for action in actions)
    assert represented(result) == {((0, p),) for p in range(4)}


@pytest.mark.parametrize('native', ['0', '1'])
def test_external_fixed_targets_preserve_literal_offers_without_occupying_them(native, monkeypatch):
    from mappa import FragmentMatchContext
    monkeypatch.setenv('MAPPA_NATIVE', native)
    source = graph(['H'], [])
    target = graph(['C', 'H', 'H', 'H'], [(0, p, 1.) for p in (1, 2, 3)])
    config = FragmentMatchConfig(branch_limit=3)
    compressed = match_fragment(source, target, seed=0, config=config)
    fixed = match_fragment(source, target, seed=0, config=config,
        context=FragmentMatchContext(placement_fixed_targets=tuple(target)))
    assert not compressed.capped and len(compressed.matches) == 1
    assert not fixed.capped and len(fixed.matches) == 3
    assert not fixed.placement_generators
    assert represented(compressed) == represented(fixed) == {((0, p),) for p in (1, 2, 3)}


@pytest.mark.parametrize('native', ['0', '1'])
def test_placement_preparation_rechecks_in_place_active_edge_changes(native, monkeypatch):
    from growth_cases import group_actions
    monkeypatch.setenv('MAPPA_NATIVE', native)
    source = graph(['C'], [])
    target = graph(['C'] * 4, [(p, (p + 1) % 4, 1.) for p in range(4)])
    weights = target.graph['wbo_matrix'].copy()
    config = FragmentMatchConfig(branch_limit=4)
    initial = match_fragment(source, target, seed=0, config=config)
    assert len(initial.matches) == 1
    target.remove_edge(0, 1)
    changed = match_fragment(source, target, seed=0, config=config)
    assert np.array_equal(target.graph['wbo_matrix'], weights)
    assert len(changed.matches) >= 2
    assert all({action[0], action[1]} == {0, 1}
               for action in group_actions(changed.placement_generators, 4))
    assert represented(changed) == {((0, p),) for p in range(4)}
    target.add_edge(0, 1, wbo=1.)
    restored = match_fragment(source, target, seed=0, config=config)
    assert len(restored.matches) == 1
    assert represented(restored) == represented(initial)


@pytest.mark.parametrize('native', ['0', '1'])
def test_graph_view_rechecks_source_active_topology(native, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', native)
    source = graph(['C', 'O'], [(0, 1, 1.)])
    target = graph(['C', 'O'], [(0, 1, 1.)])
    config = FragmentMatchConfig(iso_tolerance=.2)
    initial = match_fragment(source, target, seed=0, config=config)
    assert initial.matches[0].source_atoms == (0, 1)
    source.remove_edge(0, 1)
    changed = match_fragment(source, target, seed=0, config=config)
    assert changed.matches[0].source_atoms == (0,)
    source.add_edge(0, 1, wbo=1.)
    restored = match_fragment(source, target, seed=0, config=config)
    assert represented(restored) == represented(initial)


@pytest.mark.parametrize('native', ['0', '1'])
@pytest.mark.parametrize('changed_field', ['weight', 'element'])
def test_graph_view_rechecks_target_weights_and_elements(native, changed_field, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', native)
    source = graph(['C', 'O'], [(0, 1, 1.)])
    target = graph(['C', 'O'], [(0, 1, 1.)])
    config = FragmentMatchConfig(iso_tolerance=.2)
    initial = match_fragment(source, target, seed=0, config=config)
    assert initial.matches[0].source_atoms == (0, 1)
    if changed_field == 'weight':
        target.graph['wbo_matrix'][0, 1] = target.graph['wbo_matrix'][1, 0] = 2.
    else:
        target.nodes[1]['element'] = 'N'
    changed = match_fragment(source, target, seed=0, config=config)
    assert changed.matches[0].source_atoms == (0,)
    target.graph['wbo_matrix'][0, 1] = target.graph['wbo_matrix'][1, 0] = 1.
    target.nodes[1]['element'] = 'O'
    restored = match_fragment(source, target, seed=0, config=config)
    assert represented(restored) == represented(initial)


def test_fixed_target_labels_skip_unused_global_preparation(monkeypatch):
    from mappa import FragmentMatchContext
    from mappa.mappa_core.core.matching.growth import completion
    monkeypatch.setenv('MAPPA_NATIVE', '0')
    source = graph(['H'], [])
    target = graph(['C', 'H', 'H', 'H'], [(0, p, 1.) for p in (1, 2, 3)])
    monkeypatch.setattr(completion, '_global_self',
                        lambda *args: pytest.fail('fixed target labels need no global preparation'))
    result = match_fragment(source, target, seed=0,
        context=FragmentMatchContext(placement_fixed_targets=tuple(target)),
        config=FragmentMatchConfig(branch_limit=3))
    assert len(result.matches) == 3 and not result.placement_generators
