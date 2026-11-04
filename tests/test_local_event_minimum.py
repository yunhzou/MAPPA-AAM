"""Independent tiny mapping oracles for lazy, correlated fragment minima."""
from itertools import permutations, product

import pytest
import z3

from mappa.mappa_core.core.matching.representation.records import Fragment
from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex
from mappa.mappa_core.post_processing.decoder.fragments import FragmentDecoder, PartitionPiece
from mappa.mappa_core.post_processing.decoder.local_counts import compile_local, count_bounds
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


def test_large_symmetric_zero_uses_witness_certificate_without_compilation(monkeypatch):
    n = 100
    matrix = weights(n, [(0, i, 1.) for i in range(1, n)])
    cycle = (0, *range(2, n), 1)
    swap = (0, 2, 1, *range(3, n))
    local, response = prepare(problem(matrix, matrix), (('group', (cycle, swap)),), dict(enumerate(range(n))))
    monkeypatch.setattr('mappa.mappa_core.post_processing.decoder.local_counts.compile_local',
                        lambda *_: pytest.fail('A proved zero minimum compiled its enormous symmetry'))
    assert count_bounds(local.index, response) == (0, 0)
    assert local.minimum(response) == local.minimum(response) == 0
    assert len(local.minima) == 1


def test_joint_cycle_minimum_is_two_when_each_pair_individually_allows_zero(monkeypatch):
    local, response, choices = cyclic_response()
    expected = {literal_count(local.index.problem, choice, response.piece.source_atoms) for choice in choices}
    assert expected == {2}
    assert count_bounds(local.index, response) == (0, 2)
    assert all(0 in term.response.categories for term in response.terms)
    assert local.minimum(response) == 2


def test_bounded_minimum_requires_an_attainable_count_and_includes_upper_endpoint():
    local, response, choices = cyclic_response()
    assert {literal_count(local.index.problem, choice, response.piece.source_atoms)
            for choice in choices} == {2}
    query = compile_local(local.index, (response,))
    assert query.minimum(query.counts[0], 0, 1) is None
    assert query.minimum(query.counts[0], 0, 2) == 2
    assert query.minimum(query.counts[0], 2, 2) == 2
    assert query.minimum(query.counts[0], 0, 0) is None


def test_repeated_branch_uses_share_bounds_and_complete_binding(monkeypatch):
    local, response, _ = cyclic_response()
    assert local.bounds(response) == (0, 2)
    monkeypatch.setattr('mappa.mappa_core.post_processing.decoder.local_counts.count_bounds',
                        lambda *_: pytest.fail('A shared response rescanned all event bounds'))
    assert local.bounds(response) == (0, 2)
    assert local.minimum(response) == 2
    first = local.bind(response.terms, (3,))
    assert first is local.bind(response.terms, (3,))
    assert first is not response.terms
    assert all(term.action_ids == (3,) for term in first)
    assert all(term.action_ids == (0,) for term in response.terms)
    second = local.bind(response.terms, (4,))
    assert second is not first
    assert all(term.action_ids == (4,) for term in second)
    monkeypatch.setattr('mappa.mappa_core.post_processing.decoder.local_counts.compile_local',
                        lambda *_: pytest.fail('A shared numeric minimum was recomputed'))
    assert local.minimum(response) == 2


@pytest.mark.parametrize('atoms', [(0, 1, 2, 3), (0, 1)])
def test_two_sided_nontransitive_relation_matches_all_tiny_permutations(atoms):
    matrix = weights(4, [(0, 1, 1.), (0, 2, 1.6), (0, 3, 2.2)])
    value = relation(matrix, matrix, ((0, 1), (0, 2), (0, 3)), self_tolerance=.7)
    expected = represented(value)
    assert len(expected) == 5 and (0, 3, 2, 1) not in expected
    endpoints = problem(matrix, matrix)
    local, response = prepare(endpoints, (('relation', Fragment(value)),), dict(value.witness), atoms)
    query = compile_local(local.index, (response,))
    projections = {tuple(images[a] for a in atoms) for images in expected}
    for images in permutations(range(4), len(atoms)):
        query.solver.push()
        query.solver.add(*(query.images[0][a] == p for a, p in zip(atoms, images)))
        status = query.solver.check()
        assert (status == z3.sat) == (images in projections)
        if status == z3.sat:
            actual = query.solver.model().eval(query.counts[0]).as_long()
            assert actual == literal_count(endpoints, dict(zip(atoms, images)), atoms)
        query.solver.pop()


def test_independent_copies_preserve_internal_correlation_and_joint_pair_cost():
    local, response, choices = cyclic_response()
    query = compile_local(local.index, (response, response))
    rows = {(x, y): -1 if (x, y) == (1, 2) else 1 for x in (1, 2, 3) for y in (1, 2, 3)}
    correction = query.pair_cost(query.images[0][1], query.images[1][1], rows)
    objective = sum(query.counts) + correction
    expected = min(literal_count(local.index.problem, a, response.piece.source_atoms)
                   + literal_count(local.index.problem, b, response.piece.source_atoms)
                   + rows[a[1], b[1]] for a, b in product(choices, repeat=2))
    assert expected == 3
    assert query.minimum(objective, -1, 5) == expected
    with pytest.raises(ValueError, match='cover both endpoint domains'):
        query.pair_cost(query.images[0][1], query.images[1][1], {(1, 2): -1})


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
    query = compile_local(local.index, (response,))
    for images in permutations(range(4)):
        query.solver.push()
        query.solver.add(*(query.images[0][a] == p for a, p in enumerate(images)))
        status = query.solver.check()
        assert (status == z3.sat) == (images in expected)
        if status == z3.sat:
            assert query.solver.model().eval(query.counts[0]).as_long() == literal_count(
                endpoints, dict(enumerate(images)), tuple(range(4)))
        query.solver.pop()


def test_inconclusive_local_proof_is_not_a_minimum(monkeypatch):
    local, response, _ = cyclic_response()
    query = compile_local(local.index, (response,))
    monkeypatch.setattr(query.solver, 'check', lambda: z3.unknown)
    monkeypatch.setattr(query.solver, 'reason_unknown', lambda: 'test interruption')
    with pytest.raises(RuntimeError, match='inconclusive: test interruption'):
        query.minimum(query.counts[0], 0, 2)
