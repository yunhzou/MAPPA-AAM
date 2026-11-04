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
from mappa.mappa_core.core.symmetry.actions import PreparedActionPrograms
from mappa.mappa_core.core.scheduler.storage import BundleStore
from mappa.mappa_core.utils.io.records import decoded_events_record
from mappa import AAMProblem
from growth_cases import saved_fragment
from test_event_interpretation import dense_categories, endpoint, saved, swap


def two_branches(*, shift_actions=False):
    ep = endpoint('CCCCCC', [(0, 1, 1.), (0, 2, 1.8), (3, 4, 1.), (3, 5, 1.8), (2, 3, 1.)])
    fragments = tuple(saved_fragment(AAMProblem(ep, ep), source_atoms=atoms) for atoms in ((1, 2), (4, 5)))
    actions = tuple(('relation', fragment) for fragment in fragments)
    source = saved(ep, ep, dict(enumerate(range(6))), fragments=fragments,
                   pairs=(((0, 1, 2), (0, 1, 2)), ((3, 4, 5), (3, 4, 5))))
    second = replace(source.families[0], mapping=((0, 0), (1, 1), (2, 2), (3, 3), (4, 5), (5, 4)),
                     actions=actions[::-1] if shift_actions else actions, provenance=[])
    source.families.append(second)
    next(iter(source.branches.values())).add(1)
    return source


def correlated_branches():
    ep = endpoint('CCCC', [(0, 1, 1.), (2, 3, 1.), (0, 2, 1.8), (1, 3, .2)])
    problem = AAMProblem(ep, ep)
    joint = saved_fragment(problem, accepted_edges=((0, 1), (0, 2), (2, 3)),
                           iso_tolerance=0., self_tolerance=0.)
    pair = saved_fragment(problem, source_atoms=(1, 2))
    fixed_pair = saved_fragment(problem, source_atoms=(1, 2), anchors=((1, 1), (2, 2)))
    actions = (('relation', joint), ('relation', pair))
    source = saved(ep, ep, dict(enumerate(range(4))), fragments=(joint, pair),
                   pairs=(((0, 1), (0, 1)), ((2, 3), (2, 3))))
    first = source.families[0]
    source.families.extend((
        replace(first, actions=(actions[0], ('relation', fixed_pair)),
                provenance=[dict(cuts=((0, 1),))]),
        replace(first, mapping=((0, 2), (1, 3), (2, 0), (3, 1)),
                provenance=[dict(cuts=((0, 2),))]),
        replace(first, actions=actions[::-1], provenance=[]),
        replace(first, actions=(actions[0], ('relation', saved_fragment(problem, source_atoms=(0, 1)))),
                provenance=[]),
    ))
    next(iter(source.branches.values())).update((1, 2, 3, 4))
    source.coupled_branches[0] = {tuple(((a,), (a,)) for a in range(4))}
    return source


def assert_all_choices(decoded):
    for fid, family in enumerate(decoded.catalogue.families):
        choices = [(None, swap(fragment)) for _, fragment in family.actions]
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
    from mappa.mappa_core.core.symmetry.actions import PreparedActionPrograms
    check = PreparedActionPrograms.program

    def inventory(self, catalogue):
        collect(self, catalogue)
        collected.append(self)

    def responses(self, *args, **kwargs):
        assert len(collected) == 1 and len(collected[0].fragments) == 2
        assert len(collected[0].uses) == 2
        return check(self, *args, **kwargs)

    monkeypatch.setattr(FragmentInventory, '__init__', inventory)
    monkeypatch.setattr(PreparedActionPrograms, 'program', responses)
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
    constrained = saved_fragment(source.problem, source_atoms=(0, 1, 2),
        anchors=((0, 0),), accepted_edges=((0, 1),), iso_tolerance=.1, self_tolerance=.1)
    source.families[1] = replace(source.families[1], mapping=source.families[0].mapping,
        actions=(('relation', constrained), source.families[1].actions[1]),
        provenance=[dict(cuts=((1, 2), (2, 3)))])
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
    assert decoded.contributions(0, [swap(source.families[0].actions[0][1]), None])['total'] > 0
    with pytest.raises(ValueError, match='accepted edge'):
        decoded.contributions(1, [dict(discovery={0: 0, 1: 2, 2: 1}), None])


def test_shared_preparation_keeps_joint_choices_cuts_and_constraints():
    source = correlated_branches()
    before = source.to_record()
    decoded = decode_events(source)
    decoded_events_record(decoded)
    allowed = []
    for fid, family in enumerate(source.families):
        permitted = 0
        alternatives = [(None, dict(discovery={0: 2, 1: 3, 2: 0, 3: 1}))
                        if len(fragment.source_atoms) == 4 else (None, swap(fragment))
                        for _, fragment in family.actions]
        for choices in product(*alternatives):
            # Only the two explicit four-atom test choices are expanded here.
            mapping = dict(family.mapping)
            valid = True
            for action, choice in zip(family.actions, choices):
                if choice is not None:
                    fragment = action[1]
                    discovery = choice['discovery']
                    valid &= all(discovery[a] == p for a, p in fragment.fixed_assignments)
                    transport = {p: discovery[a] for a, p in fragment.witness}
                    mapping = {a: transport.get(p, p) for a, p in mapping.items()}
            if not valid:
                with pytest.raises(ValueError, match='anchor'):
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
    assert decoded.tables(0).atom_domains[0] != decoded.tables(3).atom_domains[0] or \
        decoded.tables(0).atom_actions[0] != decoded.tables(3).atom_actions[0]
    # The second action still controls the result when it does not enlarge a domain.
    assert decoded.tables(4).atom_actions[0] == (0, 1)
    assert decoded.tables(4).atom_actions[2] == (0, 1)
    assert source.to_record() == before


def test_independent_fragments_reuse_preparation_across_distinct_family_orders(monkeypatch):
    count = 12
    pairs = tuple((2*i, 2*i + 1) for i in range(count))
    ep = endpoint('C' * (2*count), [(a, b, 1.) for a, b in pairs] +
                  [(2*i + 1, 2*i + 2, 1.) for i in range(count - 1)])
    problem = AAMProblem(ep, ep)
    fragments = tuple(saved_fragment(problem, source_atoms=atoms, accepted_edges=(atoms,))
                      for atoms in pairs)
    mapping = dict(enumerate(range(2*count)))
    source = saved(ep, ep, mapping, fragments=fragments,
                   pairs=tuple((atoms, atoms) for atoms in pairs))
    first = source.families[0]
    source.families.extend(replace(first, actions=first.actions[i:] + first.actions[:i],
                                   provenance=[]) for i in range(1, count))
    next(iter(source.branches.values())).update(range(count))
    before = source.to_record()
    prepared = Counter()
    action_domains = SignedEventIndex.action_domains

    def observed(index, action):
        prepared[action] += 1
        return action_domains(index, action)

    def unnecessary_union(*args):
        pytest.fail('independent singleton rows repeated domain composition')

    monkeypatch.setattr(SignedEventIndex, 'action_domains', observed)
    monkeypatch.setattr(PreparedActionPrograms, '_combine_domain', unnecessary_union)
    decoded = decode_events(source)
    record = decoded_events_record(decoded)
    assert prepared == Counter({action: 1 for action in first.actions})
    assert len(decoded.local.fragments) == count
    assert len(record['tables']) == count and not decoded._tables
    shared = decoded.tables(0).fragment_responses
    for fid, family in enumerate(source.families):
        table = decoded.tables(fid)
        for piece, (response, binding) in enumerate(table.fragment_responses):
            assert response is shared[piece][0]
            assert binding == (family.actions.index(('relation', fragments[piece])),)
        for selected in (None, *range(count)):
            choices = [None] * count
            expected_mapping = dict(mapping)
            if selected is not None:
                a, b = pairs[selected]
                expected_mapping[a], expected_mapping[b] = b, a
                position = family.actions.index(('relation', fragments[selected]))
                choices[position] = swap(fragments[selected])
            assert decoded.realize(fid, choices) == expected_mapping
            expected = dense_categories(problem, expected_mapping)
            actual = decoded.contributions(fid, choices, include_events=True)
            events = [(e['category'], e['reactant_pair'], e['product_pair'])
                      for e in actual['events']]
            assert sorted(events, key=repr) == expected
            assert actual['total'] == len(expected)
    assert prepared == Counter({action: 1 for action in first.actions})
    assert source.to_record() == before


def test_ordered_joint_placements_preserve_cut_and_unmatched_boundary_events():
    matching = endpoint('CCCCCC', [])
    fragment = saved_fragment(AAMProblem(matching, matching), source_atoms=(0, 1))
    source = saved(matching, matching, dict(enumerate(range(4))), fragments=(fragment,),
                   pairs=(((0, 1), (0, 1)), ((2, 3), (2, 3))))
    exchange_first = (2, 3, 0, 1, 4, 5)
    exchange_last = (0, 1, 4, 5, 2, 3)
    orders = ((exchange_first, exchange_last), (exchange_last, exchange_first))
    original = source.families[0]
    source.families = [replace(original, actions=(('relation', fragment),
        ('group', (first,)), ('group', (second,))), provenance=[dict(cuts=((0, 2),))])
                      for first, second in orders]
    next(iter(source.branches.values())).add(1)
    before = source.to_record()
    # Both placements preserve matching weights. Original event weights make
    # their order and changing target occupancy observable.
    event_problem = AAMProblem(
        endpoint('CCCCCC', [(0, 2, 1.), (1, 3, 1.8), (0, 4, 1.)]),
        endpoint('CCCCCC', [(0, 2, 1.), (1, 3, 1.8), (2, 4, .9),
                            (3, 5, 1.6), (4, 5, .7)]))
    decoded = decode_events(source, event_problem=event_problem)
    totals = []
    for fid, (first, second) in enumerate(orders):
        outcomes = {}
        for local, move_first, move_second in product((False, True), repeat=3):
            choices = (swap(fragment) if local else None,
                       (1,) if move_first else (), (1,) if move_second else ())
            images = [1, 0, 2, 3] if local else [0, 1, 2, 3]
            if move_first:
                images = [first[p] for p in images]
            if move_second:
                images = [second[p] for p in images]
            mapping = dict(enumerate(images))
            assert decoded.realize(fid, choices) == mapping
            expected = dense_categories(event_problem, mapping)
            actual = decoded.contributions(fid, choices, include_events=True)
            events = [(e['category'], e['reactant_pair'], e['product_pair'])
                      for e in actual['events']]
            assert sorted(events, key=repr) == expected
            assert actual['total'] == len(expected)
            outcomes[local, move_first, move_second] = actual['total']
        table = decoded.tables(fid)
        cut = next(term for term in table.boundary[0, 1] if term.atoms == (0, 2))
        assert cut.manual_cut and cut.action_ids == (0, 1, 2)
        assert table.components == (((0, 1), (0, 1, 2)),)
        assert len(table.product_boundary_pairs) == 5
        totals.append(outcomes[False, True, True])
    assert totals == [8, 4]
    record = decoded_events_record(decoded)
    assert record['catalogue'] == before == source.to_record()
    assert pickle.loads(pickle.dumps(record, protocol=5)) == record


def test_sparse_binding_does_not_close_a_relaxed_action_with_itself():
    ep = endpoint('CCCC', [(0, 1, 1.), (1, 2, 1.), (2, 3, 1.)])
    fragment = saved_fragment(AAMProblem(ep, ep), anchors=((0, 0),),
        domains=((0, (0,)), (1, (1, 2)), (2, (2, 3)), (3, (3,))))
    programs = PreparedActionPrograms(SignedEventIndex(AAMProblem(ep, ep)))
    action = ('relation', fragment)
    mapping = ((1, 1),)
    domains, positions = programs.bind((action,), mapping)
    assert domains == {1: (1, 2)} and positions == {1: (0,)}
    domains, positions = programs.bind((action, action), mapping)
    assert domains == {1: (1, 2, 3)} and positions == {1: (0, 1)}


def test_branch_assembly_requests_only_boundaries_after_unique_internal_preparation(monkeypatch):
    source = correlated_branches()
    observed = Counter()
    neighboring_pairs = SignedEventIndex.neighboring_pairs
    def inspect(index, domains, *, owners=None):
        observed['internal' if owners is None else 'boundary'] += 1
        pairs = neighboring_pairs(index, domains, owners=owners)
        if owners is not None:
            assert all(owners[a] != owners[b] for a, b in pairs)
        return pairs
    monkeypatch.setattr(SignedEventIndex, 'neighboring_pairs', inspect)
    decoded = decode_events(source)
    record = decoded_events_record(decoded)
    assert observed['internal'] == len(decoded.local.fragments)
    assert observed['boundary'] == len(record['tables'])
    for fid, partitions in enumerate(decoded.fragment_partitions):
        for partition in range(len(partitions)):
            table = decoded.tables(fid, partition)
            pairs = neighboring_pairs(decoded.index, table.atom_domains)
            pairs.update((a, b) for a, b in table.cut_bonds
                         if a in table.atom_domains or b in table.atom_domains)
            assert table.statistics['neighboring_pairs'] == len(pairs)


def test_domain_preparation_eviction_preserves_results_and_bounds_retention(monkeypatch):
    ep = endpoint('C' * 15, [])
    problem = AAMProblem(ep, ep)
    programs = PreparedActionPrograms(SignedEventIndex(problem))
    fragment = saved_fragment(problem, source_atoms=(0, 1))
    action = programs.program((('relation', fragment),))[0]
    prepared = Counter()
    combine = PreparedActionPrograms._combine_domain

    def observed(action, domain):
        prepared[domain] += 1
        return combine(action, domain)

    monkeypatch.setattr(PreparedActionPrograms, '_combine_domain', staticmethod(observed))
    first = (0, 2)
    assert programs._advance(action, first) == (0, 1, 2)
    assert programs._advance(action, first) == (0, 1, 2)
    assert prepared[first] == 1
    # The same local swap sees distinct upstream bounds; their identity must
    # not depend on complete branch histories, and retention remains bounded.
    for mask in range(2, 4098):
        domain = (0,) + tuple(atom + 2 for atom in range(13) if mask & (1 << atom))
        assert programs._advance(action, domain) == (0, 1) + domain[1:]
    assert (action, first) not in programs._domains
    assert programs._advance(action, first) == (0, 1, 2)
    assert prepared[first] == 2
    assert len(programs._domains) <= 4096
    assert 0 < programs._domain_bytes <= 16 * 2**20
    assert programs._domain_bytes == sum(size for _, size in programs._domains.values())


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
                            replace(source.families[1], provenance=[dict(cuts=((0, 1),))])))
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
    assert len(decoded.inventory.fragments) == 2
    assert len(decoded.inventory.pieces) == 8
    assert len(decoded.fragment_partitions[0]) == 2
    for partition in range(2):
        choice = [swap(fragment) for _, fragment in source.families[0].actions]
        actual = decoded.contributions(0, choice, partition=partition, include_events=True)
        expected = dense_categories(source.problem, decoded.realize(0, choice))
        assert actual['total'] == len(expected)
