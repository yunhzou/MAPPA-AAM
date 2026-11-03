"""Conditional thresholds preserve exact matches and native/Python parity."""
from itertools import permutations

import numpy as np
import pytest

from mappa import BondTolerance, FragmentMatchConfig, FragmentMatchContext, match_fragment
from mappa.mappa_core.utils.bonds import bond_event_threshold
from mappa.mappa_core.core.matching.growth import native
from test_growth_compression_exactness import graph, discovered, complete


def event_threshold(source, a, b):
    elements = {n: source.nodes[n]['element'] for n in (a, b)}
    return bond_event_threshold(elements, a, b, default_threshold=.5, metal_threshold=.3)


@pytest.fixture(params=['0', '1'])
def backend(request, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', request.param)
    calls = []
    original = native.grow_island
    def observe(*args, **kwargs):
        calls.append(True)
        return original(*args, **kwargs)
    monkeypatch.setattr(native, 'grow_island', observe)
    if request.param == '1':
        assert native.built(), 'native parity tests require a built current engine'
    yield request.param
    assert bool(calls) == (request.param == '1')


def match(source, target, policy, anchor=0, **kwargs):
    result = match_fragment(source, target, seed=anchor,
        context=FragmentMatchContext(locked_mapping={anchor: anchor}, islands={anchor: 0}),
        config=FragmentMatchConfig(iso_tolerance=.5, allow_mapped_seed=True,
                                   bond_tolerance=policy), **kwargs)
    assert not result.capped
    return result


def test_mixed_metal_organic_fragment(backend):
    source = graph(['Co', 'C', 'O'], [(0, 1, 1.), (1, 2, 1.)])
    target = graph(['Co', 'C', 'O'], [(0, 1, 1.4), (1, 2, 1.4)])
    policy = BondTolerance(event_threshold, inclusive=False)
    result = match(source, target, policy, anchor=1)
    assert len(result.matches) == 1
    assert result.matches[0].source_atoms == (1, 2)
    assert dict(result.matches[0].witness) == {1: 1, 2: 2}
    assert (0, 1) in result.deferred_edges


@pytest.mark.parametrize('element,weight,threshold', [('C', 1., .5), ('Co', .3, .3)])
@pytest.mark.parametrize('position', [-1, 0, 1])
@pytest.mark.parametrize('inclusive', [False, True])
def test_threshold_boundary(backend, element, weight, threshold, position, inclusive):
    target_weight = weight + threshold
    if position:
        target_weight = np.nextafter(target_weight, -np.inf if position < 0 else np.inf)
    source = graph([element, 'O'], [(0, 1, weight)])
    target = graph([element, 'O'], [(0, 1, target_weight)])
    result = match(source, target, BondTolerance(event_threshold, inclusive=inclusive))
    difference = abs(target_weight - weight)
    expected = difference <= threshold if inclusive else difference < threshold
    assert any(m.source_atoms == (0, 1) for m in result.matches) == expected


def test_rule_is_evaluated_once_per_bond_and_uniform_policy_is_identical(backend):
    source = graph(['C'] * 3, [(0, 1, 1.), (0, 2, 1.4)])
    target = graph(['C'] * 3, [(0, 1, 1.), (0, 2, 1.4)])
    seen = []
    def tolerance(g, a, b):
        seen.append((a, b))
        return .5
    legacy = match(source, target, None)
    conditional = match(source, target, BondTolerance(tolerance))
    assert seen == [(0, 1), (0, 2)]
    assert legacy.matches == conditional.matches


def test_source_bond_attributes_keep_distinct_assignment_constraints(backend):
    source = graph(['O', 'C', 'C', 'C'], [(0, i, 1.) for i in range(1, 4)])
    target = graph(['O', 'C', 'C', 'C'], [(0, 1, .875), (0, 2, 1.125), (0, 3, 1.375)])
    for atom, tolerance in [(1, .2), (2, .4), (3, .2)]:
        source[0][atom]['precision'] = tolerance
    policy = BondTolerance(lambda g, a, b: g[a][b]['precision'], inclusive=False)
    result = match(source, target, policy)
    assert all(m.source_atoms == tuple(range(4)) for m in result.matches)
    actual = discovered(result.matches)
    expected = {((0, 0), *zip(range(1, 4), images)) for images in permutations(range(1, 4))
                if all(abs(1. - target.graph['wbo_matrix'][0, p]) < source[0][r]['precision']
                       for r, p in zip(range(1, 4), images))}
    assert len(expected) == 2 and actual == expected


def test_completed_constraints_match_explicit_conditional_relation():
    source = graph(['O', 'C', 'C', 'C'], [(0, i, 1.) for i in range(1, 4)])
    target = graph(['O', 'C', 'C', 'C'], [(0, 1, .875), (0, 2, 1.125), (0, 3, 1.375)])
    limits = {1: .2, 2: .4, 3: .2}
    policy = BondTolerance(lambda g, a, b: limits[max(a, b)], inclusive=False).resolve(source)
    fragment = {0}
    bonds = []
    for atom in (2, 1, 3):
        fragment.add(atom)
        bonds.append((0, atom))
        candidates = complete(source, target, fragment, bonds, .5,
                              locked={0: 0}, bond_tolerance=policy)
        atoms = sorted(fragment - {0})
        expected = {((0, 0), *zip(atoms, images)) for images in permutations(range(1, 4), len(atoms))
                    if all(abs(1. - target.graph['wbo_matrix'][0, p]) < limits[r]
                           for r, p in zip(atoms, images))}
        assert discovered(candidates) == expected


def test_noncontiguous_atom_indices(backend):
    source = graph(['C'] * 7, [(2, 4, 1.), (4, 6, 1.)]).subgraph([2, 4, 6]).copy()
    target = graph(['C'] * 7, [(2, 4, 1.4), (4, 6, 1.4)]).subgraph([2, 4, 6]).copy()
    policy = BondTolerance(lambda g, a, b: .3 if a == 2 else .5, inclusive=False)
    result = match(source, target, policy, anchor=4)
    assert {m.source_atoms for m in result.matches} == {(4, 6)}


def test_fragment_cache_keeps_its_conditional_policy(monkeypatch):
    from mappa.mappa_core.core.search.branching.reuse import FragmentRepair
    monkeypatch.setenv('MAPPA_NATIVE', '1')
    assert native.built()
    source = graph(['Co', 'C', 'O'], [(0, 1, 1.), (1, 2, 1.)])
    target = graph(['Co', 'C', 'O'], [(0, 1, 1.4), (1, 2, 1.4)])
    policy = BondTolerance(event_threshold, inclusive=False)
    session = FragmentRepair(source, target, bond_tolerance=policy)
    view = session.for_cut(())
    config = FragmentMatchConfig(allow_mapped_seed=True, bond_tolerance=policy)
    context = FragmentMatchContext(locked_mapping={1: 1}, islands={1: 0},
                                   fragment_cache=view)
    replayed = match_fragment(view.source, target, seed=1, context=context, config=config)
    assert replayed.matches == match(source, target, policy, anchor=1).matches
    with pytest.raises(ValueError, match='bond tolerances differ'):
        match_fragment(view.source, target, seed=1, context=context,
                       config=FragmentMatchConfig(allow_mapped_seed=True))


def test_trace_uses_conditional_thresholds(monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', '0')
    source = graph(['Co', 'O'], [(0, 1, 1.)])
    target = graph(['Co', 'O'], [(0, 1, 1.4)])
    events = []
    result = match(source, target, BondTolerance(event_threshold, inclusive=False), events=events)
    assert result.matches[0].source_atoms == (0,)
    assert events
