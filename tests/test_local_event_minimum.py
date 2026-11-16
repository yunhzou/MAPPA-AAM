"""Independent tiny mapping oracles for prepared, correlated event counts."""
from itertools import permutations, product
from collections import defaultdict

import numpy as np
import pytest
import z3

from mappa.mappa_core.core.matching.representation.records import Fragment
from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex
from mappa.mappa_core.post_processing.decoder.fragments import FragmentDecoder, PartitionPiece
from mappa.mappa_core.post_processing.decoder.local_counts import FragmentEventCounts, count_bounds
from growth_cases import saved_fragment
from test_event_interpretation import endpoint
from mappa import AAMProblem, MolecularEndpoint
from test_relaxed_fragment_decoding import problem
from test_relaxed_isomorphism import relation, represented, weights


def prepare(endpoints, actions, mapping, atoms=None):
    index = SignedEventIndex(endpoints)
    local = FragmentDecoder(index, None)
    domains, dependencies = local.programs.bind(actions, tuple(mapping.items()))
    atoms = tuple(sorted(mapping)) if atoms is None else atoms
    piece = PartitionPiece(atoms, tuple(sorted(mapping[a] for a in atoms)))
    response, _ = local.fragment(piece, mapping, actions, domains, dependencies, ())
    return local, response


def project(counts, atoms):
    # Tiny test oracle: condition every requested image; production retains a
    # queryable relation instead of materializing this Cartesian output table.
    result = {}
    for images in product(*(counts.response.atom_domains[a] for a in atoms)):
        for count in range(len(counts.response.piece.source_atoms)**2):
            answer = counts.find(count, assignments=dict(zip(atoms, images)))
            assert answer.status != 'unknown'
            if answer.status == 'found':
                result[images] = count
                break
    return result


def endpoints_for(counts, atom):
    return {images[0]: cost for images,cost in project(counts,(atom,)).items()}


def literal_count(endpoints, mapping, atoms):
    return sum(abs(endpoints.reactant.wbo[a, b] - endpoints.product.wbo[mapping[a], mapping[b]]) > .5
               for i, a in enumerate(atoms) for b in atoms[i + 1:])


def cyclic_response():
    matrix = weights(4, [(0, 1, 1.), (0, 2, 2.), (0, 3, 3.)])
    endpoints = problem(matrix, matrix)
    representative = {0: 0, 1: 2, 2: 1, 3: 3}
    cycle = (0, 2, 3, 1)
    local, response = prepare(endpoints, (('group', (cycle,)),), representative)
    choices = tuple({a: p for a, p in zip(range(4), images)} for images in (
        (0, 2, 1, 3), (0, 3, 2, 1), (0, 1, 3, 2)))
    return local, response, choices


def forbid_solver(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail('Preparing fragment event counts must not launch a solver query')
    monkeypatch.setattr(z3.Solver, 'check', forbidden)


def test_large_symmetric_zero_does_not_expand_permutations_or_launch_solver(monkeypatch):
    forbid_solver(monkeypatch)
    n = 100
    matrix = weights(n, [(0, i, 1.) for i in range(1, n)])
    cycle = (0, *range(2, n), 1)
    swap = (0, 2, 1, *range(3, n))
    local, response = prepare(problem(matrix, matrix), (('group', (cycle, swap)),), dict(enumerate(range(n))))
    assert count_bounds(local.index, response) == (0, 0)
    assert response.event_counts.minimum == response.event_counts.minimum == 0


def test_joint_cycle_minimum_is_two_when_each_pair_individually_allows_zero(monkeypatch):
    forbid_solver(monkeypatch)
    local, response, choices = cyclic_response()
    expected = {literal_count(local.index.problem, choice, response.piece.source_atoms) for choice in choices}
    assert expected == {2}
    assert count_bounds(local.index, response) == (0, 2)
    assert all(0 in term.response.categories for term in response.terms)
    assert response.event_counts.minimum == 2
    assert set(project(response.event_counts, response.piece.source_atoms).values()) == {2}


def test_prepared_response_retains_only_attainable_correlated_counts():
    local, response, choices = cyclic_response()
    assert {literal_count(local.index.problem, choice, response.piece.source_atoms)
            for choice in choices} == {2}
    expected = {tuple(choice[a] for a in response.piece.source_atoms): 2 for choice in choices}
    assert project(response.event_counts, response.piece.source_atoms) == expected
    assert endpoints_for(response.event_counts, 1) == {1: 2, 2: 2, 3: 2}


def test_repeated_branch_uses_share_bounds_and_complete_binding(monkeypatch):
    local, response, _ = cyclic_response()
    assert local.bounds(response) == (0, 2)
    monkeypatch.setattr('mappa.mappa_core.post_processing.decoder.local_counts.count_bounds',
                        lambda *_: pytest.fail('A shared response rescanned all event bounds'))
    assert local.bounds(response) == (0, 2)
    assert response.event_counts.minimum == 2
    first = local.bind(response.terms, (3,))
    assert first is local.bind(response.terms, (3,))
    assert first is not response.terms
    assert all(term.action_ids == (3,) for term in first)
    assert all(term.action_ids == (0,) for term in response.terms)
    second = local.bind(response.terms, (4,))
    assert second is not first
    assert all(term.action_ids == (4,) for term in second)
    counts = response.event_counts
    mapping = dict(response.mapping)
    repeated, _ = local.fragment(response.piece, mapping, response.actions,
                                 response.atom_domains, response.atom_actions, ())
    assert repeated is response and repeated.event_counts is counts
    assert counts.minimum == 2


@pytest.mark.parametrize('atoms', [(0, 1, 2, 3), (0, 1)])
def test_two_sided_binary_relation_matches_all_tiny_permutations(atoms):
    matrix = weights(4, [(0, 1, 1.), (0, 2, 1.6), (0, 3, 2.2)])
    value = relation(matrix, matrix, ((0, 1), (0, 2), (0, 3)))
    expected = represented(value)
    assert len(expected) == 6 and (0, 3, 2, 1) in expected
    endpoints = problem(matrix, matrix)
    local, response = prepare(endpoints, (('relation', Fragment(value)),), dict(value.witness), atoms)
    projections = {tuple(images[a] for a in atoms):
                   literal_count(endpoints, dict(enumerate(images)), atoms) for images in expected}
    assert project(response.event_counts, atoms) == projections
    assert response.event_counts.minimum == min(projections.values())


def test_endpoint_conditioned_counts_compose_without_losing_internal_correlation():
    local, response, choices = cyclic_response()
    rows = {(x, y): -1 if (x, y) == (1, 2) else 1 for x in (1, 2, 3) for y in (1, 2, 3)}
    expected = min(literal_count(local.index.problem, a, response.piece.source_atoms)
                   + literal_count(local.index.problem, b, response.piece.source_atoms)
                   + rows[a[1], b[1]] for a, b in product(choices, repeat=2))
    assert expected == 3
    endpoint = endpoints_for(response.event_counts, 1)
    actual = min(left_count + right_count + rows[left, right]
                 for left, left_count in endpoint.items() for right, right_count in endpoint.items())
    assert actual == expected


def test_local_relation_then_placement_retains_the_original_action_order():
    matrix = weights(4, [(0, 1, 1.), (0, 2, 1.), (0, 3, 1.)])
    value = relation(matrix, matrix, ((0, 1), (0, 2), (0, 3)), anchors=((0, 0), (3, 3)))
    placement = (0, 1, 3, 2)
    expected = {tuple(word[p] for p in images) for images in represented(value)
                for word in (tuple(range(4)), placement)}
    reverse_order = {tuple(images[word[p]] for p in range(4)) for images in represented(value)
                     for word in (tuple(range(4)), placement)}
    assert len(expected) == 4 and expected != reverse_order
    scoring = weights(4, [(0, 1, 1.), (0, 2, 2.), (0, 3, 3.)])
    endpoints = problem(scoring, scoring)
    local, response = prepare(endpoints,
        (('relation', Fragment(value)), ('group', (placement,))), dict(value.witness))
    expected_counts = {images: literal_count(endpoints, dict(enumerate(images)), tuple(range(4)))
                       for images in expected}
    assert project(response.event_counts, tuple(range(4))) == expected_counts


def assert_bounds_cover_mappings(local, response, mappings):
    atoms = response.piece.source_atoms
    rows = [(mapping, literal_count(local.index.problem, mapping, atoms)) for mapping in mappings]
    counts = response.event_counts
    assert counts.owner.aggregate_bound(response) <= min(count for _, count in rows)
    for atom in atoms:
        bounds = counts.endpoint_bounds(atom)
        for mapping, actual in rows:
            assert bounds[mapping[atom]] <= actual


@pytest.mark.parametrize('anchors,atoms', [
    ((), (0, 1, 2, 3)), (((0, 0),), (0, 1, 2, 3)), ((), (0, 1)),
])
def test_weight_bag_bounds_cover_every_permutation_and_partial_piece(anchors, atoms):
    source = weights(4, [(0, 1, 1.), (0, 2, 1.), (0, 3, 1.), (1, 2, 1.4)])
    target = weights(4, [(0, 1, 1.), (0, 2, 1.), (1, 2, .6), (2, 3, 1.)])
    endpoints = problem(source, target)
    fragment = saved_fragment(endpoints, anchors=anchors)
    local, response = prepare(endpoints, (('relation', fragment),), dict(fragment.witness), atoms)
    mappings = [dict(enumerate(images)) for images in permutations(range(4))
                if all(images[a] == p for a, p in anchors)]
    assert_bounds_cover_mappings(local, response, mappings)
    if len(atoms) < 4:
        assert not response.event_counts.owner.fixed_pair_bag(response)
        assert response.event_counts.owner.aggregate_bound(response) == 0


def test_weight_bag_certifies_nonzero_minimum_without_event_contraction(monkeypatch):
    source = weights(3, [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.)])
    target = weights(3, [(0, 1, 1.), (1, 2, 1.)])
    endpoints = problem(source, target)
    fragment = saved_fragment(endpoints)
    from mappa.mappa_core.core.symmetry.minimum import _Search
    monkeypatch.setattr(_Search, '_solve',
                        lambda *_: pytest.fail('The pair-weight bags already certify one event'))
    local, response = prepare(endpoints, (('relation', fragment),), dict(fragment.witness))
    assert count_bounds(local.index, response) == (0, 1)
    assert response.event_counts.owner.aggregate_bound(response) == 1
    assert response.event_counts.minimum == 1
    expected = {literal_count(endpoints, dict(enumerate(images)), (0, 1, 2))
                for images in permutations(range(3))}
    assert expected == {1}


def test_incident_weight_bags_strengthen_only_the_conditioned_endpoint():
    matrix = weights(3, [(0, 1, 1.)])
    endpoints = problem(matrix, matrix)
    fragment = saved_fragment(endpoints)
    local, response = prepare(endpoints, (('relation', fragment),), dict(fragment.witness))
    assert response.event_counts.lower == response.event_counts.minimum == 0
    assert response.event_counts.endpoint_bounds(0) == {0: 0, 1: 0, 2: 2}
    # Mapping a bonded endpoint onto the isolated atom forces two changes.
    mappings = [dict(enumerate(images)) for images in permutations(range(3))]
    assert_bounds_cover_mappings(local, response, mappings)
    # The selected binary orbit cannot move a bonded atom onto an isolated atom.
    assert endpoints_for(response.event_counts, 0) == {0: 0, 1: 0}


@pytest.mark.parametrize('seed', range(12))
def test_counted_pair_bags_preserve_all_weights_and_assignment_bound(seed):
    from scipy.optimize import linear_sum_assignment
    from mappa.mappa_core.utils.bonds import bond_event_threshold

    rng = np.random.default_rng(seed)
    elements = ('C', 'C', 'C', 'O', 'O', 'Fe')
    matrices = []
    for _ in range(2):
        values = rng.choice([0., .2, .3, .5, .5000001, .8, 1., 2.], size=(6, 6))
        matrix = np.triu(values, 1)
        matrices.append(matrix + matrix.T)
    endpoints = AAMProblem(*(MolecularEndpoint(elements, np.zeros((6, 3)), matrix)
                             for matrix in matrices))
    local = FragmentDecoder(SignedEventIndex(endpoints), None)
    counts, atoms = local.counts, tuple(range(6))
    for anchor in (None, *atoms):
        bags = []
        for side, matrix in zip(('reactant', 'product'), matrices):
            expected = defaultdict(list)
            for a in atoms:
                for b in atoms[a + 1:]:
                    expected[tuple(sorted((elements[a], elements[b]))), anchor in (a, b)].append(matrix[a, b])
            actual = counts.bag(side, atoms, anchor)
            for key, (values, multiplicities) in actual.items():
                assert np.array_equal(np.repeat(values, multiplicities), sorted(expected[key]))
            bags.append(expected)
        cost = 0
        for (labels, incident), left in bags[0].items():
            a, b = (elements.index(label) for label in labels)
            limit = bond_event_threshold(elements, a, b, default_threshold=.5, metal_threshold=.3)
            matrix = np.abs(np.subtract.outer(left, bags[1][labels, incident])) > limit
            rows, columns = linear_sum_assignment(matrix)
            cost += matrix[rows, columns].sum()
        assert counts.distance(counts.bag('reactant', atoms, anchor),
                               counts.bag('product', atoms, anchor)) == cost


def test_sweep_endpoint_cache_does_not_retain_all_nonincident_pairs():
    n = 120
    matrix = weights(n, [(i, (i + 1) % n, 1.) for i in range(n)])
    local = FragmentDecoder(SignedEventIndex(problem(matrix, matrix)), None)
    atoms = tuple(range(n))
    for atom in atoms:
        bag = local.counts.bag('reactant', atoms, atom)
        assert sum(int(counts.sum()) for _, counts in bag.values()) == n * (n - 1) // 2
    arrays = [array for collection in (local.counts.bags, local.counts.incident_bags)
              for bag in collection.values() for pair in bag.values() for array in pair]
    assert sum(array.nbytes for array in arrays) < 64 * (n + 1)


def test_weight_preserving_placement_can_move_the_entire_fragment_atom_set():
    source = weights(6, [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.)])
    target = weights(6, [(0, 1, 1.), (1, 2, 1.), (3, 4, 1.), (4, 5, 1.)])
    endpoints = problem(source, target)
    fragment = saved_fragment(endpoints, source_atoms=(0, 1, 2))
    placement = (3, 4, 5, 0, 1, 2)
    local, response = prepare(endpoints,
        (('relation', fragment), ('group', (placement,))), dict(fragment.witness))
    assert response.event_counts.owner.fixed_pair_bag(response)
    assert response.event_counts.owner.aggregate_bound(response) == 1
    mappings = [dict(enumerate(images)) for atoms in ((0, 1, 2), (3, 4, 5))
                for images in permutations(atoms)]
    assert_bounds_cover_mappings(local, response, mappings)
    assert response.event_counts.endpoint_bounds(0) == dict.fromkeys(range(6), 1)


def test_count_only_contraction_excludes_weight_preserving_placement_choices(monkeypatch):
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints
    source = weights(6, [(0, 1, 1.), (1, 2, 1.)])
    target = weights(6, [(0, 1, 1.), (1, 2, 1.), (3, 4, 1.), (4, 5, 1.)])
    endpoints = problem(source, target)
    fragment = saved_fragment(endpoints, source_atoms=(0, 1, 2),
                              witness=((0, 1), (1, 0), (2, 2)))
    placement = ('group', ((3, 4, 5, 0, 1, 2), (2, 1, 0, 3, 4, 5)))
    _, response = prepare(endpoints, (('relation', fragment), placement), dict(fragment.witness))
    assert response.event_counts.bounds == (0, 0)
    monkeypatch.setattr(FiniteConstraints, 'act',
        lambda *_: pytest.fail('Weight-preserving placement entered the count contraction'))
    counts = response.event_counts
    assert counts.minimum == 0
    assert endpoints_for(counts, 1) == {1: 0, 4: 0}


def test_weight_preserving_placements_transport_correlated_images_in_order():
    matrix = weights(3, [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.)])
    endpoints = problem(matrix, matrix)
    actions = (('group', ((1, 0, 2),)), ('group', ((0, 2, 1),)))
    _, response = prepare(endpoints, actions, {0: 0, 1: 1})
    assert project(response.event_counts, (0, 1)) == {
        (0, 1): 0, (1, 0): 0, (0, 2): 0, (2, 0): 0}
    assert endpoints_for(response.event_counts, 0) == {0: 0, 1: 0, 2: 0}


def test_placement_certificate_uses_event_weights_instead_of_matching_weights():
    matching = problem(weights(3, [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.)]),
                       weights(3, [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.)]))
    fragment = saved_fragment(matching, source_atoms=(0, 1), anchors=((0, 0), (1, 1)))
    events = problem(weights(3, [(0, 1, 1.)]), weights(3, [(0, 1, .4), (0, 2, 1.)]))
    placement = ('group', ((0, 2, 1),))
    local, response = prepare(events, (('relation', fragment), placement), dict(fragment.witness))
    # The matching graph permits the placement; original event weights do not
    # preserve its selected pair bag. Reusing that bag would wrongly assert one.
    assert not response.event_counts.owner.preserves_weights(placement)
    assert response.event_counts.owner.aggregate_bound(response) == 0
    assert response.event_counts.minimum == 0
    mappings = [{0: 0, 1: 1}, {0: 0, 1: 2}]
    assert {literal_count(events, mapping, (0, 1)) for mapping in mappings} == {0, 1}
    assert_bounds_cover_mappings(local, response, mappings)


def test_weight_bag_certificate_respects_metal_and_ordinary_thresholds(monkeypatch):
    elements = ('Fe', 'C', 'C')
    endpoints = AAMProblem(endpoint(elements, [(0, 1, 1.), (0, 2, 1.), (1, 2, .4)]),
                           endpoint(elements, [(0, 1, 1.), (0, 2, 1.4)]))
    fragment = saved_fragment(endpoints)
    from mappa.mappa_core.core.symmetry.minimum import _Search
    monkeypatch.setattr(_Search, '_solve',
                        lambda *_: pytest.fail('The element-specific weight bags certify one event'))
    local, response = prepare(endpoints, (('relation', fragment),), dict(fragment.witness))
    assert count_bounds(local.index, response) == (0, 1)
    assert response.event_counts.owner.aggregate_bound(response) == response.event_counts.minimum == 1
    for images in ((0, 1, 2), (0, 2, 1)):
        actual = sum(abs(endpoints.reactant.wbo[a, b] - endpoints.product.wbo[images[a], images[b]])
                     > (.3 if a == 0 or b == 0 else .5)
                     for a in range(3) for b in range(a + 1, 3))
        assert actual == 1
        for atom in range(3):
            assert response.event_counts.endpoint_bounds(atom)[images[atom]] <= actual
