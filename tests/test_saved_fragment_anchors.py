"""Saved anchors constrain rematching; a stored witness does not."""
from dataclasses import replace

import pytest

from mappa import BondTolerance, FragmentMatchConfig, FragmentMatchContext, match_fragment
from mappa.mappa_core.domain import AAMHierarchy, AtomPermutation, FragmentMatch, SymmetryDomain
from mappa.mappa_core.core.matching.growth import native
from mappa.mappa_core.post_processing.final_branches import FinalFamily
from mappa.mappa_core.core.search.branching.construction import find_islands
from mappa.mappa_core.core.matching.representation.state import candidate_from_record
from test_bond_tolerance import event_threshold
from test_growth_compression_exactness import graph, represented


def saved_fragment(mapping, *, pools=(), fixed=(), generators=(), cuts=()):
    return FragmentMatch(0, 0, tuple(mapping),
        representative_assignments=tuple(mapping.items()),
        symmetry_domains=tuple(SymmetryDomain(rs, ps, 'sym_block') for rs, ps in pools),
        exact_fixed=fixed, target_generators=tuple(AtomPermutation(g) for g in generators),
        deferred_edges=cuts)


def test_unpooled_assignments_are_fixed_without_explicit_flags():
    fragment = saved_fragment({2: 8, 4: 6}, cuts=((2, 4),))
    context = FragmentMatchContext.from_saved_fragment(fragment)
    assert context.locked_mapping == {2: 8, 4: 6}
    assert context.islands == {2: 2, 4: 4}
    assert context.deferred_edges == ((2, 4),)


def test_explicit_pool_lock_leaves_other_pool_members_free():
    fragment = saved_fragment({0: 0, 1: 1, 2: 2},
        pools=(((0, 1, 2), (0, 1, 2)),), fixed=(0,))
    assert fragment.fixed_assignments == ((0, 0),)


def test_only_remaining_pool_assignment_is_fixed():
    fragment = saved_fragment({0: 0, 1: 1},
        pools=(((0, 1), (0, 1)),), fixed=(0,))
    assert fragment.fixed_assignments == ((0, 0), (1, 1))


def test_open_pool_with_one_source_atom_is_still_free():
    fragment = saved_fragment({0: 0, 1: 1}, pools=(((1,), (1, 2)),))
    assert fragment.fixed_assignments == ((0, 0),)


def test_jointly_moving_atoms_are_not_anchors_even_with_local_pool_lock():
    fragment = saved_fragment(dict(enumerate(range(5))), fixed=(0, 4),
        generators=((1, 0, 3, 2, 4),))
    assert fragment.fixed_assignments == ((4, 4),)
    assert FragmentMatchContext.from_saved_fragment(fragment).locked_mapping == {4: 4}
    assert fragment.target_generators == (AtomPermutation((1, 0, 3, 2, 4)),)


def test_ring_rotation_has_no_fixed_atom():
    # All six C atoms and their six attached H atoms rotate together.
    rotation = tuple((i + 1) % 6 for i in range(6)) + tuple(6 + (i + 1) % 6 for i in range(6))
    fragment = saved_fragment(dict(enumerate(range(12))), generators=(rotation,))
    assert fragment.fixed_assignments == ()


def test_anchors_survive_serialization_and_target_transport():
    fragment = saved_fragment({0: 0, 1: 1, 2: 2},
        pools=(((1, 2), (1, 2)),), generators=((0, 2, 1),))
    hierarchy = AAMHierarchy((fragment,))
    restored = AAMHierarchy.from_record(hierarchy.to_record())
    assert restored.fragments[0].fixed_assignments == ((0, 0),)
    transported = hierarchy.relabel_target({0: 2, 1: 0, 2: 1})
    assert transported.fragments[0].fixed_assignments == ((0, 2),)


def test_unfinalized_symmetry_is_not_assumed_fixed():
    fragment = replace(saved_fragment({0: 0}), target_generators=None)
    with pytest.raises(ValueError, match='Finalize saved fragment symmetry'):
        FragmentMatchContext.from_saved_fragment(fragment)


def test_final_family_uses_saved_actions_without_freezing_the_witness():
    family = FinalFamily(tuple(enumerate(range(6))), (),
        (('pool', (1, 2)), ('group', ((0, 1, 2, 4, 3, 5),))), (.2, 1.))
    assert family.fixed_assignments == ((0, 0), (5, 5))


@pytest.mark.parametrize('backend', ['0', '1'])
def test_scheduler_splits_mixed_threshold_fragment(backend, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', backend)
    if backend == '1':
        assert native.built()
    source = graph(['Co', 'C', 'O'], [(0, 1, 1.), (1, 2, 1.)])
    target = graph(['Co', 'C', 'O'], [(0, 1, 1.4), (1, 2, 1.4)])
    result = find_islands(source, target, [1, 0, 2], iso_tol=.5,
        anchor_map={1: 1}, bond_tolerance=BondTolerance(event_threshold, inclusive=False))
    paths = list(result.paths())
    assert paths and not any(stop.reason == 'capped' for stop in result.stops)
    for path in paths:
        assert path.mapping == {0: 0, 1: 1, 2: 2}
        owners = dict(result.states[path.terminal].islands)
        assert owners[1] == owners[2] != owners[0]


@pytest.mark.parametrize('backend', ['0', '1'])
def test_tight_rematch_rediscover_mapping_instead_of_pinning_saved_witness(backend, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', backend)
    if backend == '1':
        assert native.built()
    source = graph(['O', 'C', 'C'], [(0, 1, 1.), (0, 2, 2.)])
    target = graph(['O', 'C', 'C'], [(0, 1, 1.), (0, 2, 2.)])
    fragment = saved_fragment({0: 0, 1: 2, 2: 1},
        pools=(((1, 2), (1, 2)),), fixed=(0,))
    context = FragmentMatchContext.from_saved_fragment(fragment)
    assert context.locked_mapping == {0: 0}
    result = match_fragment(source, target, seed=0, context=context,
        config=FragmentMatchConfig(allow_mapped_seed=True,
                                  bond_tolerance=BondTolerance(event_threshold, inclusive=False)))
    assert not result.capped
    assert len(result.matches) == 1
    assert result.matches[0].mapping == {0: 0, 1: 1, 2: 2}
    assert result.matches[0].fragment == frozenset((0, 1, 2))
    assert fragment.representative_assignments == ((0, 0), (1, 2), (2, 1))


@pytest.mark.parametrize('backend', ['0', '1'])
def test_rematch_rediscovers_compressed_symmetry(backend, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', backend)
    if backend == '1':
        assert native.built()
    source = graph(['C', 'H', 'H', 'H'], [(0, h, 1.) for h in (1, 2, 3)])
    target = graph(['C', 'H', 'H', 'H'], [(0, h, 1.) for h in (1, 2, 3)])
    fragment = saved_fragment({0: 0, 1: 2, 2: 3, 3: 1},
        pools=(((1, 2, 3), (1, 2, 3)),), fixed=(0,))
    result = match_fragment(source, target, seed=0,
        context=FragmentMatchContext.from_saved_fragment(fragment),
        config=FragmentMatchConfig(allow_mapped_seed=True,
                                  bond_tolerance=BondTolerance(event_threshold, inclusive=False)))
    assert not result.capped and len(result.matches) == 1
    assert len(represented([candidate_from_record(result.matches[0].symmetry)], target, .5)) == 6
