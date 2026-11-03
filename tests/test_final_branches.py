"""Reduction preserves complete local relations, order, provenance and sharing."""
from dataclasses import replace
from itertools import permutations
import json

import pytest

from growth_cases import saved_fragment
from test_event_patterns import problem, saved_path
from mappa import decode_events
from mappa.mappa_core.core.search.graph import AAMSearchGraph, SearchContext, SearchState, FragmentTransition, SearchStop
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue, canonical_actions
from mappa.mappa_workflow.mapping_analytics.family import query_path


def independent_path(order, mapping=None, labels=(7, 19)):
    value = problem('HHHH', [], [])
    mapping = mapping or dict(enumerate(range(4)))
    parts = ((0, 1), (2, 3))
    states, transitions, assigned, owners = [SearchState(0, 0, (), (), ())], [], {}, {}
    for step, part_id in enumerate(order):
        atoms = parts[part_id]
        assigned.update((r, mapping[r]) for r in atoms)
        owners.update((r, labels[part_id]) for r in atoms)
        fragment = saved_fragment(value, source_atoms=atoms,
            target_atoms=tuple(sorted(mapping[r] for r in atoms)),
            witness={r: mapping[r] for r in atoms})
        states.append(SearchState(step + 1, 0, tuple(sorted(assigned.items())),
                                  tuple(sorted(owners.items())), ()))
        transitions.append(FragmentTransition(step, step, step + 1, atoms[0], (step, 0), fragment, labels[part_id]))
    graph = AAMSearchGraph((SearchContext(tuple(range(4)), tuple(range(4)), tuple(order), iso_tolerance=1.),),
        (0,), tuple(states), tuple(transitions), (SearchStop(2, 'objective_met'),))
    return next(graph.paths())


def recovered(path, value):
    return {image for image in permutations(range(value.target_atom_count), len(path.mapping))
        if query_path(path, value, dict(zip(sorted(path.mapping), image)),
                      source_atoms=tuple(sorted(path.mapping)))[0] == 'recovered'}


def test_order_and_owner_labels_do_not_create_branches_or_families():
    value = problem('HHHH', [], [])
    first, second = independent_path((0, 1)), independent_path((1, 0), labels=(123, 2))
    saved = FinalBranchCatalogue(value).add_graph(first.graph, 'first').add_graph(second.graph, 'second')
    assert saved.path_count == 2 and saved.branch_count == len(saved.families) == 1
    assert len(saved.families[0].provenance) == 2
    assert recovered(first, value) == recovered(saved.families[0].as_path(value), value)


def test_different_target_pairs_remain_separate_and_overlap_is_not_lost():
    value = problem('HHHH', [], [])
    first = independent_path((0, 1))
    second = independent_path((1, 0), mapping={0: 2, 1: 3, 2: 0, 3: 1})
    graph = AAMSearchGraph.combine((first.graph, second.graph))
    assert len(graph.branches()) == 2
    saved = FinalBranchCatalogue(value).add_graph(graph)
    assert saved.branch_count == len(saved.families) == 2
    assert not recovered(first, value).intersection(recovered(second, value))


def test_overlapping_relations_keep_order_and_are_not_absorbed():
    value = problem('CCC', [], [])
    a = saved_fragment(value, source_atoms=(0, 1), target_atoms=(0, 1))
    b = saved_fragment(value, source_atoms=(1, 2), target_atoms=(1, 2))
    left, right = saved_path(value, (a, b)), saved_path(value, (b, a))
    saved = FinalBranchCatalogue(value).add_graph(left.graph).add_graph(right.graph)
    assert len(saved.families) == 2
    assert canonical_actions((('relation', a), ('relation', a))) == (('relation', a), ('relation', a))
    assert recovered(left, value) != recovered(right, value)
    for path, family in zip((left, right), saved.families):
        assert recovered(path, value) == recovered(family.as_path(value), value)


def test_discovery_and_self_constraints_are_part_of_family_identity():
    value = problem('CCC', [(0, 1, 1.)], [(0, 1, 1.)])
    free = saved_fragment(value)
    constrained = saved_fragment(value, accepted_edges=((0, 1),), iso_tolerance=0., self_tolerance=0.)
    saved = FinalBranchCatalogue(value).add_graph(saved_path(value, (free,)).graph)
    saved.add_graph(saved_path(value, (constrained,)).graph)
    assert saved.branch_count == 1 and len(saved.families) == 2
    assert len(recovered(saved.families[0].as_path(value), value)) == 6
    assert len(recovered(saved.families[1].as_path(value), value)) == 2


def test_anchors_remain_attached_to_the_local_relation():
    value = problem('HHHH', [], [])
    fragment = saved_fragment(value, anchors=((0, 0), (1, 1)))
    saved = FinalBranchCatalogue(value).add_graph(saved_path(value, (fragment,)).graph)
    assert saved.families[0].fixed_assignments == ((0, 0), (1, 1))
    assert recovered(saved.families[0].as_path(value), value) == {(0, 1, 2, 3), (0, 1, 3, 2)}


def test_discarded_history_pruning_preserves_catalogue_and_boundary_counts():
    from mappa.mappa_core.core.search.branching.construction import _Branch
    from mappa.mappa_core.core.search.graph import SearchGraphBuilder
    value = problem('CCC', [(0, 1, 1.)], [(0, 2, 1.)])
    builder = SearchGraphBuilder(SearchContext((0, 1, 2), (0, 1, 2), (0,), cuts=((0, 1),)))
    root = _Branch(builder)
    discarded = root.fork()
    discarded.commit(saved_fragment(value, source_atoms=(0,), target_atoms=(1,), witness={0: 1}))
    builder.stop(discarded, 'capped', stage='combined_live_leaves', count=2, limit=1)
    live = root.fork()
    fragment = saved_fragment(value, anchors=((0, 0),))
    live.commit(fragment)
    builder.stop(live, 'objective_met')
    before = FinalBranchCatalogue(value).add_graph(builder.finish())
    builder.retain([live])
    after = FinalBranchCatalogue(value).add_graph(builder.finish())
    assert len(builder.finish().states) == 2
    assert before.branches == after.branches and before.path_count == after.path_count == 1
    assert before.families[0].actions == after.families[0].actions
    for saved in (before, after):
        decoded = decode_events(saved)
        assert [decoded.contributions(0, [dict(discovery=dict(enumerate(image)))])['total']
                for image in ((0, 1, 2), (0, 2, 1))] == [2, 0]
        assert recovered(saved.families[0].as_path(value), value) == {(0, 1, 2), (0, 2, 1)}


@pytest.mark.parametrize('transport', ['json', 'bundle'])
def test_catalogue_restore_shares_fragments_without_changing_choices(tmp_path, transport):
    from mappa.execution import BundleStore
    from mappa.mappa_core.core.scheduler.results import CatalogueResult
    from mappa import AAMSearchConfig
    value = problem('HHHH', [], [])
    first, second = independent_path((0, 1)), independent_path((1, 0))
    saved = FinalBranchCatalogue(value).add_graph(first.graph).add_graph(second.graph)
    family = saved.families[0]
    saved.families.append(replace(family, provenance=[dict(cuts=((0, 1),))]))
    next(iter(saved.branches.values())).add(1)
    expected = json.loads(json.dumps(saved.to_record()))
    assert len(expected['fragments']) == 2
    if transport == 'json':
        restored = FinalBranchCatalogue.from_record(value, expected)
    else:
        store = BundleStore(tmp_path / 'bundles', scratch=tmp_path)
        restored = store.load(store.put(CatalogueResult('test', saved, (), (), AAMSearchConfig(), 'native'))).catalogue
    assert json.loads(json.dumps(restored.to_record())) == expected
    assert restored.families[0].mapping is restored.families[1].mapping
    assert restored.families[0].actions is restored.families[1].actions
    assert all(a is b for a, b in zip(restored.families[0].fragments, restored.families[1].fragments))
    assert recovered(restored.families[0].as_path(value), value) == recovered(first, value)
    assert all(decode_events(restored).contributions(i)['total'] == 0 for i in range(2))
    other = problem('HHHH', [(0, 1, 1.)], [])
    with pytest.raises(ValueError, match='different endpoint'):
        FinalBranchCatalogue.from_record(other, expected)


def test_fragment_record_restoration_shares_source_preparation():
    from mappa.mappa_core.core.matching.representation.records import Fragment
    value = problem('CCCC', [], [])
    a = saved_fragment(value, source_atoms=(0, 1), target_atoms=(0, 1))
    b = saved_fragment(value, source_atoms=(0, 1), target_atoms=(2, 3), witness={0: 2, 1: 3})
    pool = {}
    left = Fragment.from_record(json.loads(json.dumps(a.to_record())), self_pool=pool)
    right = Fragment.from_record(json.loads(json.dumps(b.to_record())), self_pool=pool)
    assert left.source_self_automorph is right.source_self_automorph
    assert left.target_atoms != right.target_atoms


def shared_source_catalogue():
    from mappa.mappa_core.core.matching.representation.records import Fragment
    value = problem('CCCCCC', [(0, 1, 1.), (0, 2, 2.), (3, 4, 1.), (3, 5, 2.)],
                    [(0, 1, 1.), (0, 2, 2.), (3, 4, 1.), (3, 5, 2.)])
    first = saved_fragment(value, source_atoms=(0, 1, 2), accepted_edges=((0, 1), (0, 2)),
                           self_tolerance=.5)
    second = saved_fragment(value, source_atoms=(0, 1, 2), target_atoms=(3, 4, 5),
        witness={0: 3, 1: 4, 2: 5}, accepted_edges=((0, 1), (0, 2)), self_tolerance=.5)
    second = Fragment(replace(second.relation, source=first.source_self_automorph))
    return FinalBranchCatalogue(value).add_graph(saved_path(value, (first,)).graph).add_graph(
        saved_path(value, (second,)).graph)


def test_binary_catalogue_preserves_shared_preparation_without_record_expansion(tmp_path, monkeypatch):
    import gc
    from weakref import ref
    from mappa import AAMSearchConfig
    from mappa.execution import BundleStore
    from mappa.mappa_core.core.scheduler.results import CatalogueResult
    import mappa.mappa_core.core.symmetry.relaxed as relaxed

    saved = shared_source_catalogue()
    expected = json.loads(json.dumps(saved.to_record()))
    original = ref(saved.families[0].fragments[0].source_self_automorph)
    assert original().constraints
    assert original() is saved.families[1].fragments[0].source_self_automorph
    def forbidden(*args, **kwargs):
        pytest.fail('binary transport must preserve prepared objects without JSON or re-preparation')
    store = BundleStore(tmp_path / 'bundles', scratch=tmp_path)
    with monkeypatch.context() as guard:
        guard.setattr(FinalBranchCatalogue, 'to_record', forbidden)
        guard.setattr(relaxed, '_constraints', forbidden)
        reference = store.put(CatalogueResult('shared', saved, (), (), AAMSearchConfig(), 'native'))
        del saved
        gc.collect()
        assert original() is None
        restored = store.load(reference).catalogue
    first, second = (family.fragments[0] for family in restored.families)
    assert first.source_self_automorph is second.source_self_automorph
    assert first.target_atoms != second.target_atoms
    assert not first.source_self_automorph.weights.flags.writeable
    assert not first.target_self_automorph.weights.flags.writeable
    assert json.loads(json.dumps(restored.to_record())) == expected
    assert '_evaluation' not in vars(first.source_self_automorph)
    assert len(restored._family_ids) == len(restored.families) == 2
    decoded = decode_events(restored)
    for fid, fragment in enumerate((first, second)):
        assert decoded.contributions(fid)['total'] == 0
        mapping = dict(fragment.witness)
        mapping[1], mapping[2] = mapping[2], mapping[1]
        assert decoded.contributions(fid, [dict(discovery=mapping)])['total'] == 2
    del first, second, fragment, decoded, restored
    gc.collect()


def test_decoded_record_and_bundle_share_self_payloads_across_distinct_fragments(tmp_path):
    from mappa.execution import BundleStore
    from mappa.mappa_core.utils.io.records import decoded_events_record
    record = decoded_events_record(decode_events(shared_source_catalogue()))
    first, second = record['catalogue']['fragments']
    assert first['target']['atoms'] != second['target']['atoms']
    assert first['source'] is second['source']
    store = BundleStore(tmp_path / 'bundles', scratch=tmp_path)
    restored = store.load(store.put(record))
    a, b = restored['catalogue']['fragments']
    assert a['source'] is b['source']
    assert a['source']['weights'] is b['source']['weights']
    assert restored == record


@pytest.mark.smoke
def test_public_grouping_is_eager_and_preserves_distinct_witnesses(monkeypatch):
    from mappa import AAMResult, AAMSearchConfig, AAMSearchMetrics
    value = problem('COO', [], [])
    a = saved_fragment(value, witness={0: 0, 1: 1, 2: 2})
    b = saved_fragment(value, witness={0: 0, 1: 2, 2: 1})
    first, second = saved_path(value, (a,)), saved_path(value, (b,))
    graph = AAMSearchGraph.combine((first.graph, second.graph))
    result = AAMResult(value, AAMSearchConfig(), graph, AAMSearchMetrics.from_record({}, 0.))
    assert len(result.branches) == result.metrics.retained_branch_count == 1
    assert len(result.branches[0].paths) == 2
    assert len(graph.literal_branches()) == 2
    monkeypatch.setattr(AAMSearchGraph, 'branches', lambda *_: pytest.fail('branch groups were rebuilt'))
    assert len(result.branches) == 1
