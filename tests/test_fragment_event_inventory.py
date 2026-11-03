"""Shared local decoding must preserve every original branch's joint choices."""
from dataclasses import replace
from collections import Counter
from itertools import product
import gc
import pickle
from weakref import WeakSet, ref

import pytest

from mappa import decode_events
from mappa.mappa_core.post_processing.decoder.events import (
    FragmentEventTables, SignedEventIndex,
)
from mappa.mappa_core.post_processing.decoder.fragments import FragmentInventory
from mappa.mappa_core.core.symmetry.actions import _ActionSequence
from mappa.mappa_core.core.scheduler.storage import BundleStore
from mappa.mappa_core.utils.io.records import decoded_events_record
from test_event_interpretation import dense_categories
from test_fragment_conversion import endpoint, saved


def two_branches(*, shift_actions=False):
    ep = endpoint('CCCCCC', [(0, 1, 1.), (0, 2, 1.8), (3, 4, 1.), (3, 5, 1.8), (2, 3, 1.)])
    actions = (('pool', (1, 2)), ('pool', (4, 5)))
    source = saved(ep, ep, dict(enumerate(range(6))), actions=actions,
                   pairs=(((0, 1, 2), (0, 1, 2)), ((3, 4, 5), (3, 4, 5))))
    second = replace(source.families[0], mapping=((0, 0), (1, 1), (2, 2), (3, 3), (4, 5), (5, 4)),
                     actions=actions[::-1] if shift_actions else actions, provenance=[])
    source.families.append(second)
    next(iter(source.branches.values())).add(1)
    return source


def correlated_branches():
    ep = endpoint('CCCC', [(0, 1, 1.), (2, 3, 1.), (0, 2, 1.8), (1, 3, .2)])
    actions = (('group', ((1, 0, 3, 2),)), ('pool', (1, 2)))
    source = saved(ep, ep, dict(enumerate(range(4))), actions=actions,
                   pairs=(((0, 1), (0, 1)), ((2, 3), (2, 3))))
    first = source.families[0]
    source.families.extend((
        replace(first, required_edges=((0, 1),), policy=(.2, .1),
                provenance=[dict(cuts=((0, 1),))]),
        replace(first, mapping=((0, 1), (1, 0), (2, 2), (3, 3)),
                provenance=[dict(cuts=((0, 2),))]),
        replace(first, actions=actions[::-1], provenance=[]),
        replace(first, actions=(actions[0], ('pool', (0, 1))), provenance=[]),
    ))
    next(iter(source.branches.values())).update((1, 2, 3, 4))
    source.coupled_branches[0] = {tuple(((a,), (a,)) for a in range(4))}
    return source


def assert_all_choices(decoded):
    for fid, family in enumerate(decoded.catalogue.families):
        choices = [(None, tuple(reversed(action[1]))) for action in family.actions]
        for choice in product(*choices):
            mapping = decoded.realize(fid, choice)
            expected = dense_categories(decoded.catalogue.problem, mapping)
            result = decoded.contributions(fid, choice, include_events=True)
            actual = [(e['category'], e['reactant_pair'], e['product_pair']) for e in result['events']]
            assert sorted(actual, key=repr) == expected
            assert result['total'] == len(expected)


def test_inventory_is_complete_before_any_event_decoding(monkeypatch):
    source = two_branches()
    collected = []
    collect = FragmentInventory.__init__
    check = SignedEventIndex.target_action_responses

    def inventory(self, catalogue):
        collect(self, catalogue)
        collected.append(self)

    def responses(self, *args, **kwargs):
        assert len(collected) == 1 and len(collected[0].fragments) == 2
        assert len(collected[0].uses) == 2
        return check(self, *args, **kwargs)

    monkeypatch.setattr(FragmentInventory, '__init__', inventory)
    monkeypatch.setattr(SignedEventIndex, 'target_action_responses', responses)
    decoded = decode_events(source)
    assert decoded.inventory is collected[0]
    assert not decoded.local.fragments and not decoded.local.bonds and not decoded._tables
    assert decoded.inventory.uses[0][0][0] is decoded.inventory.uses[1][0][0]


@pytest.mark.parametrize('shift_actions', [False, True])
def test_shared_fragment_results_bind_original_branch_action_positions(shift_actions):
    source = two_branches(shift_actions=shift_actions)
    before = source.to_record()
    decoded = decode_events(source)
    first, second = decoded.tables(0), decoded.tables(1)
    # A has the same local program even when unrelated B changes its witness
    # and the positions of A's actions in the complete program move.
    assert first.fragment_responses[0][0] is second.fragment_responses[0][0]
    assert first.fragment_responses[1][0] is not second.fragment_responses[1][0]
    assert len(decoded.local.fragments) == 3
    assert first.fragment_responses[0][1] == (0,)
    assert second.fragment_responses[0][1] == ((1,) if shift_actions else (0,))
    record = decoded_events_record(decoded)
    first_group, second_group = [record['tables'][row['tables'][0]]['internal'][0]
                                 for row in record['families']]
    assert (first_group is second_group) == (not shift_actions)
    assert first_group['terms'][0]['action_ids'] == (0,)
    assert second_group['terms'][0]['action_ids'] == ((1,) if shift_actions else (0,))
    assert_all_choices(decoded)
    assert source.to_record() == before


def test_cut_annotations_and_constraints_are_not_lost_through_sharing():
    source = two_branches()
    source.families[1] = replace(source.families[1], mapping=source.families[0].mapping,
        required_edges=((0, 1),), policy=(.2, .1), provenance=[dict(cuts=((1, 2), (2, 3)))])
    decoded = decode_events(source)
    first, second = decoded.tables(0), decoded.tables(1)
    assert first.cut_bonds == frozenset()
    assert second.cut_bonds == frozenset(((1, 2), (2, 3)))
    assert any(t.atoms == (1, 2) and t.manual_cut for t in second.internal[0])
    assert not any(t.atoms == (1, 2) for t in first.internal[0])
    record = decoded_events_record(decoded)
    first_group, second_group = [record['tables'][row['tables'][0]]['internal'][0]
                                 for row in record['families']]
    assert first_group is not second_group
    assert any(t['atoms'] == (1, 2) and t['manual_cut'] for t in second_group['terms'])
    assert decoded.contributions(0, [(2, 1), None])['total'] > 0
    with pytest.raises(ValueError, match='original fragment constraints'):
        decoded.contributions(1, [(2, 1), None])


def test_shared_preparation_keeps_joint_choices_cuts_and_constraints():
    source = correlated_branches()
    before = source.to_record()
    decoded = decode_events(source)
    decoded_events_record(decoded)
    allowed = []
    for fid, family in enumerate(source.families):
        permitted = 0
        alternatives = [(None, (1,)) if action[0] == 'group'
                        else (None, tuple(reversed(action[1])))
                        for action in family.actions]
        for choices in product(*alternatives):
            # Enumerate only this four-atom oracle's two binary operations.
            # The group swaps both fragments together, never independently.
            mapping = dict(family.mapping)
            for action, choice in zip(family.actions, choices):
                if choice is not None:
                    if action[0] == 'group':
                        permutation = (1, 0, 3, 2)
                    else:
                        permutation = list(range(4))
                        a, b = action[1]
                        permutation[a], permutation[b] = b, a
                    mapping = {a: permutation[p] for a, p in mapping.items()}
            valid = all(abs(source.problem.reactant.wbo[a, b] -
                            source.problem.product.wbo[mapping[a], mapping[b]]) <= .1
                        for a, b in family.required_edges)
            if not valid:
                with pytest.raises(ValueError, match='original fragment constraints'):
                    decoded.contributions(fid, choices)
                continue
            permitted += 1
            assert decoded.realize(fid, choices) == mapping
            expected = dense_categories(source.problem, mapping)
            for partition in range(len(decoded.fragment_partitions[fid])):
                result = decoded.contributions(fid, choices, partition=partition, include_events=True)
                events = [(e['category'], e['reactant_pair'], e['product_pair'])
                          for e in result['events']]
                assert sorted(events, key=repr) == expected
                assert result['total'] == len(expected)
        allowed.append(permitted)
    assert allowed == [4, 2, 4, 4, 4]
    internal = decoded.tables(1)
    boundary = decoded.tables(2)
    assert internal.cut_bonds == frozenset(((0, 1),))
    assert any(term.atoms == (0, 1) and term.manual_cut for term in internal.internal[0])
    assert boundary.cut_bonds == frozenset(((0, 2),))
    assert any(term.atoms == (0, 2) and term.manual_cut for term in boundary.boundary[0, 1])
    assert decoded.tables(0).atom_domains[0] == (0, 1, 2)
    assert decoded.tables(3).atom_domains[0] == (0, 1)
    # The second action still controls the result when it does not enlarge a domain.
    assert decoded.tables(4).atom_domains[0] == (0, 1)
    assert decoded.tables(4).atom_actions[0] == (0, 1)
    assert decoded.tables(4).atom_actions[2] == (0,)
    assert source.to_record() == before


def test_preparation_scales_with_unique_programs_and_images(monkeypatch):
    source = correlated_branches()
    prepared = Counter()
    prepare = _ActionSequence._prepare_atom

    def observed(self, image):
        prepared[tuple(self), image] += 1
        return prepare(self, image)

    monkeypatch.setattr(_ActionSequence, '_prepare_atom', observed)
    decoded = decode_events(source)
    record = decoded_events_record(decoded)
    expected = {(family.actions, image) for family in source.families
                for _, image in family.mapping}
    assert len(expected) == 12  # Three ordered programs, four initial images each.
    assert prepared == Counter({key: 1 for key in expected})
    assert len(record['tables']) == 6  # Distinct cuts, constraints and partitions remain.
    assert not decoded._tables
    for fid, partitions in enumerate(decoded.fragment_partitions):
        for partition in range(len(partitions)):
            decoded.tables(fid, partition)
    assert prepared == Counter({key: 1 for key in expected})


def test_preparation_eviction_bounds_retention_without_changing_rows(monkeypatch):
    ep = endpoint('C' * 128, [])
    source = saved(ep, ep, dict(enumerate(range(128))), actions=(('pool', (0, 1)),))
    source.families.extend(replace(source.families[0], actions=(('pool', (0, atom)),))
                           for atom in range(2, 34))
    next(iter(source.branches.values())).update(range(33))
    prepared = Counter()
    prepare = _ActionSequence._prepare_atom

    def observed(self, image):
        prepared[tuple(self), image] += 1
        return prepare(self, image)

    monkeypatch.setattr(_ActionSequence, '_prepare_atom', observed)
    decoded = decode_events(source)
    first = decoded.tables(0)
    expected = first.atom_domains[0], first.atom_actions[0]
    assert expected == ((0, 1), (0,))
    decoded_events_record(decoded)
    # 33 genuine programs over 128 images exceed the retained-row allowance.
    assert sum(prepared.values()) == 4224
    rows = decoded.index._prepared_atoms
    assert len(rows.rows) <= 4096 and 0 < rows.bytes <= 16 * 2**20
    program = decoded.local.programs.program(source.families[0].actions)
    assert program.atom(0) == expected
    assert prepared[source.families[0].actions, 0] == 2
    assert len(rows.rows) <= 4096 and 0 < rows.bytes <= 16 * 2**20


def test_session_does_not_retain_released_preparation():
    from mappa.mappa_core.post_processing.decoder.api import EventDecodeSession

    source = correlated_branches()
    session = EventDecodeSession(source.problem)
    prepared = []
    for _ in range(3):
        decoded = session.decode(source)
        decoded_events_record(decoded)
        prepared.append(ref(decoded.local.programs))
        del decoded
    gc.collect()
    assert all(reference() is None for reference in prepared)
    # Releasing branch preparation must leave the session usable.
    assert decoded_events_record(session.decode(source))['catalogue'] == source.to_record()


def test_export_and_bundle_share_records_without_mutating_runtime_tables(tmp_path):
    decoded = decode_events(two_branches())
    runtime = decoded.tables(0)
    domains, actions = dict(runtime.atom_domains), dict(runtime.atom_actions)
    record = decoded_events_record(decoded)
    first, second = [record['tables'][row['tables'][0]] for row in record['families']]
    a = next(region['terms'] for region in first['internal'] if region['fragment'] == 0)
    b = next(region['terms'] for region in second['internal'] if region['fragment'] == 0)
    assert a is b and all(x is y for x, y in zip(a, b))
    restored = pickle.loads(pickle.dumps(record, protocol=5))
    assert restored == record
    first, second = [restored['tables'][row['tables'][0]] for row in restored['families']]
    assert first['internal'][0] is second['internal'][0]
    store = BundleStore(tmp_path / 'bundles', scratch=tmp_path)
    restored = store.load(store.put(record))
    assert restored == record
    first, second = [restored['tables'][row['tables'][0]] for row in restored['families']]
    assert first['internal'][0] is second['internal'][0]
    assert first['components'] is second['components']
    exported = record['tables'][0]
    assert all('atom_domains' not in table and 'atom_actions' not in table
               for table in record['tables'])
    exported['internal'][0]['fragment'] = -1
    assert runtime.atom_domains == domains and runtime.atom_actions == actions
    assert -1 not in runtime.internal
    assert_all_choices(decoded)


def test_export_distinguishes_group_owners_with_shared_terms():
    source = two_branches()
    source.coupled_branches[0] = {(((0,), (0,)), ((1, 2), (1, 2)), ((3, 4, 5), (3, 4, 5)))}
    decoded = decode_events(source)
    record = decoded_events_record(decoded)
    groups = []
    for number in record['families'][0]['tables']:
        table = record['tables'][number]
        owner = table['fragments'].index((3, 4, 5))
        groups.append(next(row for row in table['internal'] if row['fragment'] == owner))
    assert {row['fragment'] for row in groups} == {1, 2}
    assert groups[0]['terms'] is groups[1]['terms']
    assert groups[0] is not groups[1]
    assert_all_choices(decoded)


def test_export_releases_working_tables_and_does_not_compile_evaluation(monkeypatch):
    source = two_branches()
    source.families.extend((replace(source.families[0]),
                            replace(source.families[1], required_edges=((0, 1),))))
    next(iter(source.branches.values())).update((2, 3))
    decoded = decode_events(source)
    alive = WeakSet()
    construct = FragmentEventTables.__init__

    def observed(self, *args, **kwargs):
        assert not alive, 'export retained a previous working table'
        construct(self, *args, **kwargs)
        alive.add(self)

    def forbidden(self):
        pytest.fail('export compiled an evaluation kernel')

    with monkeypatch.context() as patch:
        patch.setattr(FragmentEventTables, '__init__', observed)
        patch.setattr(FragmentEventTables, '_evaluation', property(forbidden))
        record = decoded_events_record(decoded)
    assert not alive and not decoded._tables
    assert [row['tables'] for row in record['families']] == [[0], [1], [0], [2]]
    assert len(record['tables']) == 3
    assert record['catalogue'] == source.to_record()
    assert_all_choices(decoded)
    table = decoded.tables(0)
    assert table is decoded.tables(2)
    assert '_evaluation' in vars(table)


def test_inventory_keeps_alternative_and_coupled_partitions():
    source = two_branches()
    source.coupled_branches[0] = {tuple(((a,), (a,)) for a in range(6))}
    decoded = decode_events(source)
    assert len(decoded.inventory.fragments) == 8
    assert len(decoded.fragment_partitions[0]) == 2
    for partition in range(2):
        actual = decoded.contributions(0, [(2, 1), (5, 4)], partition=partition, include_events=True)
        expected = dense_categories(source.problem, decoded.realize(0, [(2, 1), (5, 4)]))
        assert actual['total'] == len(expected)
