"""Shared local decoding must preserve every original branch's joint choices."""
from dataclasses import replace
from itertools import product
import pickle
from weakref import WeakSet

import pytest

from mappa import decode_events
from mappa.mappa_core.post_processing.decoder.events import FragmentEventTables, SignedEventIndex
from mappa.mappa_core.post_processing.decoder.fragments import FragmentInventory
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
    assert decoded.contributions(0, [(2, 1), None])['total'] > 0
    with pytest.raises(ValueError, match='original fragment constraints'):
        decoded.contributions(1, [(2, 1), None])


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
    assert first['internal'][0]['terms'] is second['internal'][0]['terms']
    store = BundleStore(tmp_path / 'bundles', scratch=tmp_path)
    restored = store.load(store.put(record))
    assert restored == record
    first, second = [restored['tables'][row['tables'][0]] for row in restored['families']]
    assert first['atom_domains'] is second['atom_domains']
    assert first['atom_actions'] is second['atom_actions']
    assert first['atom_domains'][1] is first['atom_domains'][2]
    assert first['atom_actions'][1] is first['atom_actions'][2]
    assert first['internal'][0]['terms'] is second['internal'][0]['terms']
    exported = record['tables'][0]
    assert exported['atom_domains'] is not runtime.atom_domains
    assert exported['atom_actions'] is not runtime.atom_actions
    exported['atom_domains'][0] = (-1,)
    exported['atom_actions'][0] = (-1,)
    assert runtime.atom_domains == domains and runtime.atom_actions == actions
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
