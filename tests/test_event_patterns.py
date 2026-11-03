"""Bounded checks: no mapping-permutation or group-element enumeration."""
from collections import defaultdict
import math
from pathlib import Path
import random
import sys

import numpy as np
import pynauty
import pytest
from sympy.combinatorics import Permutation, PermutationGroup

from mappa import AAMProblem, MolecularEndpoint
from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex
from mappa.mappa_core.core.symmetry.responses import response_bounds
from reference_event_decoder import extract_path_events
from mappa.mappa_core.core.search.graph import AAMSearchGraph, SearchContext, SearchState, FragmentTransition, SearchStop
from growth_cases import saved_fragment

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bench'))
from reference_dense_events import DeltaPatterns


def problem(elements, r_edges, p_edges):
    def endpoint(edges):
        w = np.zeros((len(elements), len(elements)))
        for a, b, value in edges:
            w[a, b] = w[b, a] = value
        return MolecularEndpoint(tuple(elements), np.zeros((len(elements), 3)), w)
    return AAMProblem(endpoint(r_edges), endpoint(p_edges))


def raw_record(value):
    return {side: dict(elements=list(ep.elements), coordinates=ep.coordinates.tolist(), wbo=ep.wbo.tolist())
            for side, ep in (('reactant', value.reactant), ('product', value.product))}


def saved_path(value, fragments=None):
    """Construct recorded path evidence directly; never run AAM search here."""
    n = value.source_atom_count
    identity = tuple((i, i) for i in range(n))
    parts = fragments if fragments is not None else (saved_fragment(value, anchors=identity),)
    states = [SearchState(0, 0, (), (), ())]
    transitions = []
    assigned, owners = {}, {}
    for step, fragment in enumerate(parts):
        assigned.update(fragment.witness)
        owners.update((i, step) for i in fragment.source_atoms)
        states.append(SearchState(step + 1, 0, tuple(sorted(assigned.items())), tuple(sorted(owners.items())), ()))
        transitions.append(FragmentTransition(step, step, step + 1,
            fragment.source_atoms[0] if fragment.source_atoms else None, (step, 0), fragment, step))
    graph = AAMSearchGraph((SearchContext(tuple(range(n)), tuple(range(value.target_atom_count)), tuple(range(n)),
                                        iso_tolerance=1.0),), (0,), tuple(states), tuple(transitions),
                          (SearchStop(len(parts), 'objective_met'),))
    path = next(graph.paths())
    assert path.mapping == assigned
    return path


@pytest.mark.parametrize('tau', [0.3, 0.5, 1.0])
def test_response_bounds_match_dense_float_comparisons(tau):
    values = np.array(sorted({x for v in (0., .1, .3, .5, .7, 1., 1.4, 1.8)
                             for x in (math.nextafter(v, -math.inf), v, math.nextafter(v, math.inf))}))
    for weight in values:
        a, b = response_bounds(values, weight, tau)
        dense = np.where(values - weight < -tau, -1, np.where(values - weight > tau, 1, 0))
        assert np.array_equal(dense, [-1] * a + [0] * (b - a) + [1] * (len(values) - b))


def test_sparse_and_dense_groups_and_event_partitions_agree():
    rng = random.Random(117)
    for sample in range(8):
        elements = ('C', 'O', 'O', 'H', 'H') if sample % 2 else ('V', 'O', 'O', 'H', 'H')
        edges = [(a, b, rng.choice((0., .1, .3, .7, 1., 1.4, 1.8)))
                 for a in range(5) for b in range(a + 1, 5)]
        product = [(a, b, rng.choice((0., .1, .3, .7, 1., 1.4, 1.8))) for a, b, _ in edges]
        value = problem(elements, edges, product)
        legacy, index = DeltaPatterns(raw_record(value)), SignedEventIndex(value)
        dense_generators = [tuple(g[:5]) for g in pynauty.autgrp(legacy.base)[0]]
        dense_group = PermutationGroup([Permutation(g, size=5) for g in dense_generators] or [Permutation(list(range(5)))])
        sparse_group = PermutationGroup([Permutation(g, size=5) for g in index.source_generators] or [Permutation(list(range(5)))])
        # Generator containment proves equality; no closure/group enumeration.
        assert all(sparse_group.contains(g) for g in dense_group.generators)
        assert all(dense_group.contains(g) for g in sparse_group.generators)
        maps = [(0, 1, 2, 3, 4), (0, 2, 1, 3, 4), (0, 1, 2, 4, 3)]
        old, new = [legacy.describe(m) for m in maps], [index.describe(m) for m in maps]
        for a, b in zip(old, new):
            assert a['events'] == b['events'] and a['total'] == b['total']
        assert [[a['id'] == b['id'] for b in old] for a in old] == [[a['id'] == b['id'] for b in new] for a in new]
        assert index.base.number_of_vertices <= legacy.base.number_of_vertices


def test_equivalent_hydrogen_events_merge_but_signs_do_not():
    value = problem('CHH', [(0, 1, 1), (0, 2, 1)], [(0, 1, 1), (0, 2, .1)])
    index = SignedEventIndex(value)
    first, second = index.describe([0, 1, 2]), index.describe([0, 2, 1])
    assert first['events'] != second['events'] and first['id'] == second['id']
    value = problem('COO', [(0, 1, 1), (0, 2, 1.4)], [(0, 1, 1.4), (0, 2, 1.6)])
    assert SignedEventIndex(value).describe([0, 1, 2])['id'] != SignedEventIndex(value).describe([0, 2, 1])['id']


def test_window_retains_nonminimum_events_after_zero_event_seed():
    value = problem('COO', [(0, 1, 1), (0, 2, 1.4)], [(0, 1, 1.4), (0, 2, 1.6)])
    path = saved_path(value, (saved_fragment(value),))
    index = SignedEventIndex(value)
    result = extract_path_events(path, value, index, max_events=1)
    assert result['complete'], result
    assert {p['total'] for p in result['patterns']} == {0, 1}
    result = extract_path_events(path, value, index, max_events=0)
    assert result['complete'] and [p['total'] for p in result['patterns']] == [0]


def test_threshold_change_reuses_original_search_path():
    value = problem('COO', [(0, 1, 1), (0, 2, 1.4)], [(0, 1, 1.4), (0, 2, 1.6)])
    path = saved_path(value, (saved_fragment(value),))
    strict = extract_path_events(path, value, SignedEventIndex(value, threshold=.5))
    loose = extract_path_events(path, value, SignedEventIndex(value, threshold=.7))
    assert strict['complete'] and loose['complete']
    assert {p['total'] for p in strict['patterns']} == {0, 1}
    assert {p['total'] for p in loose['patterns']} == {0}
    assert path.context.iso_tolerance == 1.0


def test_raw_metal_events_remain_distinct_under_binary_search_constraints():
    raw = problem(('V', 'O', 'O'), [(0, 1, 2.2), (0, 2, .7)], [(0, 1, .8), (0, 2, 2.1)])
    search = problem(('V', 'O', 'O'), [(0, 1, 1), (0, 2, 1)], [(0, 1, 1), (0, 2, 1)])
    path = saved_path(search, (saved_fragment(search),))
    result = extract_path_events(path, search, SignedEventIndex(raw), max_events=2)
    assert result['complete'] and {p['total'] for p in result['patterns']} == {0, 2}
    assert np.array_equal(search.reactant.wbo, [[0, 1, 1], [1, 0, 0], [1, 0, 0]])


def test_cross_fragment_events_are_not_dropped():
    edges = [(0, 1, 1), (0, 2, 1.8), (1, 2, 1)]
    value = problem('COO', edges, edges)
    path = saved_path(value, tuple(saved_fragment(value, source_atoms=atoms)
                                  for atoms in ((1, 2), (0,))))
    result = extract_path_events(path, value, SignedEventIndex(value), max_events=2)
    assert result['complete'] and {p['total'] for p in result['patterns']} == {0, 2}


def test_one_zero_event_witness_does_not_certify_the_whole_relation():
    edges = [(0, 1, 1), (0, 2, 1.4), (0, 3, 1.8)]
    value = problem('COOO', edges, edges)
    path = saved_path(value, (saved_fragment(value),))
    index = SignedEventIndex(value)
    assert index.describe([0, 2, 1, 3])['total'] == index.describe([0, 1, 3, 2])['total'] == 0
    assert not index.invariant_path(path)
    result = extract_path_events(path, value, index)
    assert result['complete'] and any(p['total'] > 0 for p in result['patterns'])


def test_large_invariant_fragment_never_enters_solver(monkeypatch):
    from mappa import decode_events
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    # The completed fragment represents 12! hydrogen assignments.
    edges = [(0, i, 1) for i in range(1, 13)]
    value = problem(('C',) + ('H',) * 12, edges, edges)
    path = saved_path(value, (saved_fragment(value),))
    def forbidden(*args, **kwargs):
        raise AssertionError('invariant pool must remain compressed')
    monkeypatch.setattr('mappa.mappa_workflow.mapping_analytics.family.compile_path', forbidden)
    decoded = decode_events(FinalBranchCatalogue(value).add_graph(path.graph))
    assert decoded.contributions(0)['total'] == 0
    assert decoded.tables(0).statistics['projected_pair_states'] == 0


def test_reject_partial_invalid_and_asymmetric_inputs():
    value = problem('COO', [], [])
    index = SignedEventIndex(value)
    assert index.counts([]).shape == (0, 2)
    for vector in ([], [0, 1], [0, 1, 1], [1, 0, 2], [0., 1., 2.]):
        with pytest.raises(ValueError):
            index.describe(vector)
    for tau in (0, -1, float('nan')):
        with pytest.raises(ValueError):
            SignedEventIndex(value, threshold=tau)
    unbalanced = AAMProblem(value.reactant, problem('CO', [], []).product)
    assert SignedEventIndex(unbalanced).counts([[0, 1, -1]]).shape == (1, 2)


def test_correlated_fragment_relation_does_not_allow_independent_atom_swaps():
    edges = [(0, 1, 1.), (1, 2, 1.8), (2, 3, 1.8)]
    value = problem('CCCC', edges, edges)
    path = saved_path(value, (saved_fragment(value,
        accepted_edges=((0, 1), (1, 2), (2, 3))),))
    result = extract_path_events(path, value, SignedEventIndex(value))
    assert result['complete'], result
    assert {p['total'] for p in result['patterns']} == {0, 2}
    # Reversal moves both pairs together. Swapping only the outer pair breaks
    # the path and lies outside the complete saved relationship.
    assert SignedEventIndex(value).describe([3, 1, 2, 0])['total'] == 4


def test_symbolic_orbit_equality_matches_canonical_identity():
    import z3
    from mappa.mappa_workflow.mapping_analytics.event_query import _same_event_orbit
    value = problem('CHHH', [(0, a, 1.) for a in (1, 2, 3)], [(0, 1, .1), (0, 2, 1.), (0, 3, 1.)])
    index = SignedEventIndex(value)
    pattern = index.describe([0, 1, 2, 3])
    for vector in ([0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 2, 1]):
        candidate = index.describe(vector)
        terms = {tuple((int(a), int(b))): 0 for a, b in zip(index.a, index.b)}
        for kind, sign in [('broken', -1), ('formed', 1)]:
            terms.update({tuple(pair): sign for pair in candidate['events'][kind]})
        solver = z3.Solver()
        solver.set(timeout=1000)
        solver.add(_same_event_orbit(terms, z3.IntVal(candidate['total']), pattern, index))
        assert (solver.check() == z3.sat) == (pattern['id'] == candidate['id'])
        # Same event location with reversed sign is a distinct output class.
        reversed_terms = {k: -v for k, v in terms.items()}
        solver = z3.Solver()
        solver.set(timeout=1000)
        solver.add(_same_event_orbit(reversed_terms, z3.IntVal(candidate['total']), pattern, index))
        assert solver.check() == z3.unsat
