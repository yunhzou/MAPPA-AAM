"""Witness improvement preserves discovery, symmetry and honest time limits."""
from dataclasses import replace
from itertools import permutations
import time

import numpy as np
import pytest

from mappa.mappa_core.core.matching.representation.records import Fragment
from mappa.mappa_core.core.symmetry import binary_symmetry, DiscoveryEdge, RelaxedIsomorphism
from mappa.mappa_core.core.symmetry.minimum import fragment_symmetry_minimum
from mappa.mappa_core.core.symmetry.objective import FragmentObjective
from mappa.mappa_core.core.symmetry.optimization import optimize_fragment
from test_binary_symmetry import weights
from test_relaxed_isomorphism import represented


def fragment(source, target, accepted, witness, *, tolerance=1., elements=None, anchors=()):
    n = len(source)
    elements = tuple(elements or 'C' * n)
    return Fragment(RelaxedIsomorphism(
        binary_symmetry(elements, source, anchors=tuple(r for r, _ in anchors)),
        binary_symmetry(elements, target, anchors=tuple(p for _, p in anchors)),
        tuple(enumerate(witness)), anchors,
        tuple((r, tuple(p for p in range(n) if elements[r] == elements[p])) for r in range(n)),
        tuple(DiscoveryEdge((a, b), source[a, b], tolerance) for a, b in accepted),
        tuple((int(a), int(b)) for a, b in zip(*np.where(np.triu(target > .2, 1)))), .2))


def assert_witness(value, answer):
    mapping = value.relation.realize_permutations(answer.discovery_witness,
        source_permutation=answer.source_operation, target_permutation=answer.target_operation)
    assert tuple(sorted(mapping.items())) == answer.witness
    objective = FragmentObjective(value)
    assert objective.count(objective.images(answer.witness)) == answer.event_count
    assert answer.lower_bound <= answer.event_count
    assert answer.complete == (answer.lower_bound == answer.event_count)


def triangle():
    matrix = weights(4, [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.), (0, 3, 1.)])
    return fragment(matrix, matrix, ((0, 1), (1, 2), (0, 3)), (1, 0, 3, 2))


def test_triangle_leaves_two_event_starting_class_and_attains_zero():
    value = triangle()
    before = value.witness
    assert not value.relation.discovery_in_witness_orbit
    baseline = fragment_symmetry_minimum(value)
    assert baseline.event_count == 2 and baseline.complete
    answer = optimize_fragment(value)
    assert answer.event_count == answer.lower_bound == 0 and answer.complete
    assert answer.work.discovery_improvements > 0
    assert value.witness == before
    assert_witness(value, answer)


@pytest.mark.parametrize('case', ['triangle_tail', 'ring_chain', 'continuous', 'anchored'])
def test_completed_atom_sets_find_an_optimum_without_growth_history(case):
    from mappa import AAMProblem, MolecularEndpoint, decode_events
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    from test_event_patterns import saved_path

    source = weights(4, [(0, 1, 1.), (1, 2, 1.), (2, 3, 1.), (3, 0, 1.)])
    target = weights(4, [(0, 1, 1.), (1, 2, 1.), (2, 3, 1.)])
    if case in ('triangle_tail', 'anchored'):
        source = target = triangle().relation.source.weights
    if case == 'continuous':
        target = weights(4, [(0, 1, 1.7), (1, 2, 1.), (2, 3, 2.2)])
    elements = ('C',) * 4
    endpoint = lambda matrix: MolecularEndpoint(elements, np.zeros((4, 3)), matrix)
    problem = AAMProblem(endpoint(source), endpoint(target))
    anchors = ((0, 1),) if case == 'anchored' else ()
    value = Fragment.from_atom_sets(problem, range(4), range(4),
        enumerate((1, 0, 3, 2)), anchors=anchors)
    objective = FragmentObjective(value)
    allowed = [image for image in permutations(range(4))
               if all(image[r] == p for r, p in anchors)]
    expected = min(objective.count(np.asarray(image)) for image in allowed)
    answer = optimize_fragment(value)
    assert answer.event_count == expected
    assert_witness(value, answer)
    # Actual decoding shares the same selected response and counts original WBOs.
    graph = saved_path(problem, (value,)).graph
    graph = replace(graph, contexts=(replace(graph.contexts[0], anchors=anchors),))
    decoded = decode_events(FinalBranchCatalogue(problem).add_graph(graph))
    assert len(decoded.fragments) == 1
    assert decoded.minimum_events == expected
    mapping = decoded.realize(0)
    assert all(mapping[r] == p for r, p in anchors)


@pytest.mark.parametrize('seed', range(32))
def test_starting_bound_matches_exhaustive_injective_fixed_bond_cost(seed):
    from mappa import AAMProblem, MolecularEndpoint
    from mappa.mappa_core.core.symmetry.optimization import _StartingMappings
    from mappa.mappa_core.core.symmetry.minimum import SymmetryWork
    rng = np.random.default_rng(seed)
    def matrix():
        upper = np.triu(rng.choice([0., .2, .5, 1., 1.4, 2.], (5, 5)), 1)
        return upper + upper.T
    elements = ('C', 'C', 'C', 'O', 'O')
    left, right = matrix(), matrix()
    endpoint = lambda x: MolecularEndpoint(elements, np.zeros((5, 3)), x)
    value = Fragment.from_atom_sets(AAMProblem(endpoint(left), endpoint(right)),
                                    range(5), range(5), enumerate(range(5)))
    objective = FragmentObjective(value)
    domains = [(i,) if rng.random() < .4 else tuple(j for j in range(5)
        if elements[i] == elements[j] and (i == j or rng.random() < .7)) for i in range(5)]
    fixed = {i for i, values in enumerate(domains) if len(values) == 1}
    expected = min(sum(abs(left[i, j]-right[image[i], image[j]]) > .5
                       for i in range(5) for j in range(i+1, 5) if i in fixed or j in fixed)
        for image in permutations(range(5)) if all(image[i] in domains[i] for i in range(5)))
    search = _StartingMappings(value, objective, time.perf_counter()+5., 100, SymmetryWork())
    assert search.search.bound(domains, search.ceiling) == expected


def test_double_triangle_continuous_weights_improve_three_to_two():
    source = weights(4, [(0, 1, 1.), (1, 2, 1.), (2, 3, 1.)])
    target = weights(4, [(0, 1, 2.), (0, 2, 1.), (0, 3, 1.), (1, 2, 1.), (1, 3, 1.)])
    value = fragment(source, target, ((0, 1), (1, 2), (2, 3)), (2, 0, 1, 3))
    assert fragment_symmetry_minimum(value).event_count == 3
    answer = optimize_fragment(value)
    assert answer.event_count == answer.lower_bound == 2 and answer.complete
    assert_witness(value, answer)


def test_complete_growth_uses_existing_symmetry_without_repeating_discovery():
    source = weights(4, [(0, 1, 1.), (1, 2, 1.), (2, 3, 2.)])
    target = weights(4, [(0, 1, 1.), (1, 2, 2.), (2, 3, 1.)])
    value = fragment(source, target, ((0, 1), (1, 2), (2, 3)), range(4))
    assert FragmentObjective(value).lower_bound() == 0
    assert value.relation.discovery_in_witness_orbit
    answer = optimize_fragment(value)
    assert answer.complete and answer.event_count == answer.lower_bound == 2
    assert answer.work.discovery_nodes == 0
    assert 'covered by solved symmetry' in answer.reason
    assert_witness(value, answer)


def test_unattainable_pair_bound_remains_explicitly_best_found():
    # All six pair weights have the same multiset, but a star is not a path.
    source = weights(4, [(0, 1, 1.), (0, 2, 1.), (0, 3, 1.)])
    target = weights(4, [(0, 1, 1.), (1, 2, 1.), (2, 3, 1.)])
    value = fragment(source, target, ((0, 1),), (0, 1, 2, 3))
    answer = optimize_fragment(value)
    assert answer.lower_bound == 0 and answer.event_count == 2
    assert not answer.complete and not answer.work.time_limit_reached
    assert 'unproved' in answer.reason
    assert_witness(value, answer)


@pytest.mark.parametrize('difference,expected', [(.5, 0), (np.nextafter(.5, 1.), 1)])
def test_weight_bound_uses_strict_threshold_and_subthreshold_contacts(difference, expected):
    source, target = weights(2, []), weights(2, [(0, 1, difference)])
    value = fragment(source, target, (), (0, 1))
    objective = FragmentObjective(value)
    assert objective.lower_bound() == expected
    answer = optimize_fragment(value)
    assert answer.complete and answer.event_count == expected


def test_connectivity_change_point_six_to_point_two_is_not_an_event():
    value = fragment(weights(2, [(0, 1, .6)]), weights(2, [(0, 1, .2)]), (), (0, 1))
    answer = optimize_fragment(value)
    assert answer.complete and answer.event_count == answer.lower_bound == 0


def test_original_event_problem_retains_sparse_indices_and_metal_threshold():
    from mappa import AAMProblem, MolecularEndpoint

    source_elements, target_elements = ('C', 'H', 'O', 'Fe'), ('Fe', 'O', 'H', 'C')
    source = binary_symmetry(source_elements, weights(4, [(1, 3, 1.)]), atoms=(1, 3))
    target = binary_symmetry(target_elements, weights(4, [(0, 2, 1.)]), atoms=(0, 2))
    value = Fragment(RelaxedIsomorphism(source, target, ((1, 2), (3, 0)), (),
        ((1, (2,)), (3, (0,))), (DiscoveryEdge((1, 3), 1., 1.),), ((0, 2),), .2))
    original = AAMProblem(
        MolecularEndpoint(source_elements, np.zeros((4, 3)), weights(4, [(1, 3, .8)])),
        MolecularEndpoint(target_elements, np.zeros((4, 3)), weights(4, [(0, 2, 1.2)])))
    assert optimize_fragment(value).event_count == 0
    answer = optimize_fragment(value, event_problem=original)
    assert answer.complete and answer.event_count == answer.lower_bound == 1
    assert answer.witness == ((1, 2), (3, 0))
    assert FragmentObjective(value, original).count((1, 0)) == 1
    assert optimize_fragment(value, event_problem=original, metal_threshold=.5).event_count == 0


def test_zero_budget_returns_original_valid_witness():
    value = triangle()
    answer = optimize_fragment(value, time_limit=0)
    assert answer.witness == value.witness and answer.event_count == 2
    assert answer.work.time_limit_reached and not answer.complete
    assert_witness(value, answer)


def test_native_preparation_delay_preserves_incumbent_and_spends_same_deadline(monkeypatch):
    from mappa.mappa_core.core.symmetry import minimum
    original = minimum.permutation_group

    def delayed(*args):
        time.sleep(.02)
        return original(*args)

    monkeypatch.setattr(minimum, 'permutation_group', delayed)
    value = triangle()
    answer = optimize_fragment(value, time_limit=.01)
    assert answer.work.time_limit_reached and not answer.complete
    assert answer.work.discovery_nodes == 0
    assert answer.event_count == 2
    assert_witness(value, answer)


def test_anchor_is_retained_while_improving_other_discoveries():
    value = triangle()
    anchored = fragment(value.source_symmetry.weights, value.target_symmetry.weights,
        ((0, 1), (1, 2), (0, 3)), (1, 0, 3, 2), anchors=((3, 2),))
    answer = optimize_fragment(anchored)
    assert dict(answer.witness)[3] == 2
    assert answer.event_count == 2 and not answer.complete
    assert_witness(anchored, answer)


@pytest.mark.parametrize('seed', range(24))
def test_tiny_weighted_relations_against_independent_complete_oracle(seed):
    rng = np.random.default_rng(seed)
    source = weights(4, [(a, b, float(rng.choice([0., .2, .8, 1., 1.6, 2.])))
                        for a in range(4) for b in range(a + 1, 4)])
    target = weights(4, [(a, b, float(rng.choice([0., .2, .8, 1., 1.6, 2.])))
                        for a in range(4) for b in range(a + 1, 4)])
    accepted = [(a, b) for a in range(4) for b in range(a + 1, 4)
                if source[a, b] > .2 and target[a, b] > .2 and abs(source[a, b] - target[a, b]) <= 1.]
    value = fragment(source, target, accepted, range(4))
    objective = FragmentObjective(value)
    allowed = represented(value.relation)
    exact = min(objective.count(images) for images in allowed)
    answer = optimize_fragment(value, time_limit=.2)
    assert tuple(dict(answer.witness)[r] for r in range(4)) in allowed
    assert objective.lower_bound() <= exact <= answer.event_count
    assert answer.event_count <= objective.count(range(4))
    if answer.complete:
        assert answer.event_count == exact
    assert_witness(value, answer)


@pytest.mark.parametrize('seed', range(12))
def test_pair_weight_lower_bound_matches_independent_pair_assignment_oracle(seed):
    rng = np.random.default_rng(seed)
    pairs = [(0, 1), (0, 2), (1, 2)]
    left, right = rng.choice([0., .2, .5, .8, 1., 1.6, 2.], size=(2, 3))
    value = fragment(weights(3, [(*pair, w) for pair, w in zip(pairs, left)]),
                     weights(3, [(*pair, w) for pair, w in zip(pairs, right)]), (), range(3))
    expected = min(sum(abs(a - b) > .5 for a, b in zip(left, perm)) for perm in permutations(right))
    assert FragmentObjective(value).lower_bound() == expected
