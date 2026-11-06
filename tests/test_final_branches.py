"""Atom-set deduplication retains one original fragment and its shared response."""
from dataclasses import replace
from itertools import permutations
import json

import pytest

from growth_cases import saved_fragment
from test_event_patterns import problem, saved_path
from mappa import decode_events
from mappa.mappa_core.core.search.graph import AAMSearchGraph, SearchContext, SearchState, FragmentTransition, SearchStop
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
from mappa.mappa_core.core.symmetry.actions import canonical_actions
from mappa.mappa_core.core.search.queries.family import query_path


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


def test_assignment_pairs_are_shared_before_publication_and_after_restore():
    import pickle

    value = problem('HHHH', [], [])
    graphs = [independent_path((0, 1), mapping=mapping).graph for mapping in
              ({0: 0, 1: 1, 2: 2, 3: 3}, {0: 0, 1: 2, 2: 1, 3: 3})]
    catalogue = FinalBranchCatalogue(value)
    for graph in graphs:
        catalogue.add_graph(graph)
    record = catalogue.to_record()
    for saved in (catalogue, pickle.loads(pickle.dumps(catalogue)),
                  FinalBranchCatalogue.from_record(value, record)):
        assert saved.to_record() == record
        assert len(saved.families) == 2
        pairs = {}
        for family in saved.families:
            for pair in family.mapping:
                assert pairs.setdefault(pair, pair) is pair
        assert saved.families[0].mapping != saved.families[1].mapping


def test_raw_overlapping_relations_keep_their_order():
    value = problem('CCC', [], [])
    a = saved_fragment(value, source_atoms=(0, 1), target_atoms=(0, 1))
    b = saved_fragment(value, source_atoms=(1, 2), target_atoms=(1, 2))
    left, right = saved_path(value, (a, b)), saved_path(value, (b, a))
    assert canonical_actions((('relation', a), ('relation', a))) == (('relation', a), ('relation', a))
    assert recovered(left, value) != recovered(right, value)
    for original in (left, right):
        saved = FinalBranchCatalogue(value).add_graph(original.graph)
        assert recovered(saved.families[0].as_path(value), value) == recovered(original, value)


def test_growth_constraints_do_not_split_a_completed_atom_set_fragment():
    value = problem('CCC', [(0, 1, 1.)], [(0, 1, 1.)])
    free = saved_fragment(value)
    constrained = saved_fragment(value, accepted_edges=((0, 1),), iso_tolerance=0.)
    saved = FinalBranchCatalogue(value).add_graph(saved_path(value, (free,)).graph)
    saved.add_graph(saved_path(value, (constrained,)).graph)
    assert saved.branch_count == len(saved.families) == 1
    assert saved.families[0].fragments == (free,)
    assert saved.families[0].fragments[0] is free
    assert sum(len(f.provenance) for f in saved.families) == 2
    assert recovered(saved_path(value, (free,)), value) == recovered(saved.families[0].as_path(value), value)


def test_duplicate_starting_witness_is_discarded_and_later_search_finds_alignment():
    from test_fragment_optimization import triangle
    original = triangle()
    value = problem('CCCC', [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.), (0, 3, 1.)],
                    [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.), (0, 3, 1.)])
    improved = replace(original, relation=replace(original.relation,
        witness=tuple(enumerate(range(4)))))
    record = original.to_record()
    saved = FinalBranchCatalogue(value)
    saved.add_graph(saved_path(value, (original,)).graph)
    saved.add_graph(saved_path(value, (improved,)).graph)
    assert saved.branch_count == len(saved.families) == 1
    assert saved.families[0].fragments == (original,)
    assert saved.families[0].fragments[0] is original
    assert len(saved.to_record()['fragments']) == 1
    decoded = decode_events(saved, event_slack=5)
    assert len(decoded.inventory.fragments) == 1
    assert decoded.minimum_events == 0
    assert original.to_record() == record
    answer = decoded.local.prepare(original).optimization
    assert answer.work.discovery_improvements > 0


@pytest.mark.parametrize('relabeled_frame', [False, True])
def test_branches_share_one_fragment_and_use_its_retained_witness(relabeled_frame):
    value = problem('HHHH', [], [])
    shared = saved_fragment(value, source_atoms=(0, 1), target_atoms=(0, 1))
    duplicate = (shared.relabel_target({0: 1, 1: 0}) if relabeled_frame else
                 saved_fragment(value, source_atoms=(0, 1), target_atoms=(0, 1),
                                witness={0: 1, 1: 0}))
    together = saved_fragment(value, source_atoms=(2, 3), target_atoms=(2, 3))
    split = tuple(saved_fragment(value, source_atoms=(i,), target_atoms=(i,)) for i in (2, 3))
    saved = FinalBranchCatalogue(value).add_graph(saved_path(value, (shared, together)).graph)
    saved.add_graph(saved_path(value, (duplicate, *split)).graph)
    assert saved.branch_count == len(saved.families) == 2
    assert len(saved.to_record()['fragments']) == 4
    for family in saved.families:
        assert next(f for f in family.fragments if f.source_atoms == (0, 1)) is shared
        assert dict(family.mapping)[0] == 0 and dict(family.mapping)[1] == 1
        family.validate_representative(value)
    assert len(decode_events(saved, event_slack=5).inventory.fragments) == 4


def test_selected_scope_supplies_the_retained_fragment():
    value = problem('CC', [], [])
    first = saved_fragment(value)
    second = saved_fragment(value, witness={0: 1, 1: 0})
    graph = AAMSearchGraph.combine((saved_path(value, (first,)).graph,
                                    saved_path(value, (second,)).graph))
    paths = tuple(graph.paths())
    saved = FinalBranchCatalogue(value).add_graph(graph, paths=paths[1:])
    assert saved.families[0].fragments[0] is second
    assert saved.families[0].mapping == second.witness
    assert saved.path_count == 1


def test_explicit_endpoint_anchors_remain_fixed_after_reduction():
    value = problem('HHHH', [], [])
    fragment = saved_fragment(value, anchors=((0, 0), (1, 1)))
    graph = saved_path(value, (fragment,)).graph
    graph = replace(graph, contexts=(replace(graph.contexts[0], anchors=((0, 0), (1, 1))),))
    saved = FinalBranchCatalogue(value).add_graph(graph)
    assert saved.families[0].fixed_assignments == ((0, 0), (1, 1))
    assert recovered(saved.families[0].as_path(value), value) == {(0, 1, 2, 3), (0, 1, 3, 2)}


def test_reduction_shares_final_pairs_without_symmetry_or_starting_search(monkeypatch):
    from mappa.mappa_core.core.symmetry import structure, optimization
    from mappa.mappa_core.core.search.branching.construction import _Branch
    from mappa.mappa_core.core.search.graph import SearchGraphBuilder

    value = problem('CCC', [(0, 1, 1.), (1, 2, 1.)], [(0, 1, 1.), (1, 2, 1.)])
    builder = SearchGraphBuilder(SearchContext((0, 1, 2), (0, 1, 2), (0,), iso_tolerance=1.))
    branch = _Branch(builder)
    branch.commit(saved_fragment(value, source_atoms=(0, 1), accepted_edges=((0, 1),)))
    branch.commit(saved_fragment(value, source_atoms=(1, 2), anchors=((1, 1),), accepted_edges=((1, 2),)))
    builder.stop(branch, 'objective_met')
    graph = builder.finish()
    other = saved_path(value, (saved_fragment(value, witness={0: 2, 1: 1, 2: 0}),)).graph
    def forbidden(*args, **kwargs):
        pytest.fail('Reduction must only aggregate records, not solve a fragment')
    with monkeypatch.context() as guard:
        guard.setattr(structure, 'structural_generators', forbidden)
        guard.setattr(optimization, 'optimize_fragment', forbidden)
        saved = FinalBranchCatalogue(value).add_graph(graph).add_graph(other)
    assert saved.path_count == 2 and saved.branch_count == 1
    assert len(saved.families) == 2
    assert {id(f) for f in saved.families[0].fragments} == {id(edge.match) for edge in graph.transitions}
    assert saved.families[1].fragments == (other.transitions[0].match,)
    for original, family in zip((graph, other), saved.families):
        from growth_cases import represented
        assert represented(original) == represented(family.as_path(value))


@pytest.mark.parametrize('recorded_decisions', [False, True])
def test_million_growth_histories_reduce_and_decode_without_unfolding(monkeypatch, recorded_decisions):
    import pickle
    from mappa.mappa_core.core.search.graph import SearchHistories
    bonds = [(0, 1, 1.), (1, 2, 1.), (0, 2, 1.)]
    value = problem('CCC', bonds, bonds)
    decisions = (tuple(saved_fragment(value, accepted_edges=edges)
                       for edges in (((0, 1), (1, 2)), ((0, 2), (1, 2))))
                 if recorded_decisions else (None, None))
    mapping, owners = tuple(enumerate(range(3))), tuple((a, 1) for a in range(3))
    states = tuple(SearchState(i, 0, mapping, owners, ()) for i in range(21))
    edges = tuple(FragmentTransition(2*i+j, i, i+1, None, (i, j), decisions[j])
                  for i in range(20) for j in range(2))
    graph = AAMSearchGraph((SearchContext((0, 1, 2), (0, 1, 2), (0,)),), (0,), states,
                           edges, (SearchStop(20, 'objective_met'),))
    def forbidden(*args, **kwargs):
        pytest.fail('Reduction, publication and sweep setup must not unfold growth paths')
    monkeypatch.setattr(AAMSearchGraph, 'paths', forbidden)
    monkeypatch.setattr(SearchHistories, 'paths', forbidden)
    catalogue = FinalBranchCatalogue(value).add_graph(graph)
    assert catalogue.path_count == 2**20
    assert catalogue.branch_count == len(catalogue.families) == 1
    if recorded_decisions:
        assert set(catalogue.families[0].fragments) == {decisions[0]}
    origins = catalogue.families[0].provenance
    assert len(origins) == 2**20 and len(origins.references) == 1
    assert len(origins.references[0][0].parents) == 21
    encoded = pickle.dumps(catalogue)
    assert len(encoded) < 10_000
    restored = pickle.loads(encoded)
    assert restored.path_count == 2**20
    encoded = json.dumps(restored.to_record())
    assert len(encoded) < 10_000
    restored = FinalBranchCatalogue.from_record(value, json.loads(encoded))
    assert len(restored.families[0].provenance) == 2**20
    decoded = decode_events(restored)
    assert decoded.minimum_events == 0
    from mappa.mappa_core.utils.io.records import decoded_events_record
    assert decoded_events_record(decoded)['catalogue']['input_paths'] == 2**20


def test_native_history_counts_exceed_machine_integer_without_expansion():
    from mappa.mappa_core.core.search.graph import SearchHistories
    states = tuple(SearchState(i, 0, (), (), ()) for i in range(81))
    edges = tuple(FragmentTransition(2*i+j, i, i+1, None, (), None)
                  for i in range(80) for j in range(2))
    graph = AAMSearchGraph((SearchContext((), (), ()),), (0,), states,
                          edges, (SearchStop(80, 'objective_met'),))
    histories, endings = SearchHistories.collect(graph)
    assert len(histories.parents) == 81
    assert histories.counts[endings[0][2]] == 2**80


def test_compact_origins_preserve_each_path_and_selected_scope():
    from mappa.mappa_core.core.search.graph import SearchHistories
    value = problem('CC', [], [])
    # The first two decisions rejoin before a placement-only final decision.
    states = tuple(SearchState(i, 0, ((0, 0), (1, 1)), ((0, 1), (1, 2)), ()) for i in range(4))
    edges = (FragmentTransition(0, 0, 1, None, (), None),
             FragmentTransition(1, 0, 2, None, (), None),
             FragmentTransition(2, 1, 3, None, (), None, placement_generators=((1, 0),)),
             FragmentTransition(3, 2, 3, None, (), None))
    graph = AAMSearchGraph((SearchContext((0, 1), (0, 1), (0,)),), (0,), states,
                           edges, (SearchStop(3, 'objective_met'),))
    expected = [path.transitions for path in graph.paths()]
    histories, endings = SearchHistories.collect(graph)
    assert [path for _, _, node in endings for path in histories.paths(node)] == expected
    chosen = tuple(graph.paths())[-1:]
    selected = FinalBranchCatalogue(value).add_graph(graph, paths=chosen)
    assert selected.path_count == 1
    assert selected.families[0].provenance[0]['transitions'] == expected[-1]
    saved = FinalBranchCatalogue(value).add_graph(graph)
    restored = FinalBranchCatalogue.from_record(value, json.loads(json.dumps(saved.to_record())))
    restored.add_graph(graph)
    assert restored.path_count == 4 and len(restored.families) == len(saved.families)
    assert sorted(tuple(origin['transitions']) for f in restored.families for origin in f.provenance) == sorted(expected*2)


def test_placement_alternatives_keep_original_histories():
    value = problem('CC', [], [])
    states = tuple(SearchState(i, 0, ((0, 0), (1, 1)), ((0, 1), (1, 1)), ()) for i in range(5))
    pairs = ((0, 1), (0, 2), (1, 3), (2, 3), (0, 3), (3, 4))
    edges = tuple(FragmentTransition(i, a, b, None, (), None,
        placement_generators=((1, 0),) if i in (2, 4) else ()) for i, (a, b) in enumerate(pairs))
    graph = AAMSearchGraph((SearchContext((0, 1), (0, 1), (0,)),), (0,), states,
                           edges, (SearchStop(4, 'objective_met'),))
    saved = FinalBranchCatalogue(value).add_graph(graph)
    assert saved.branch_count == 1 and len(saved.families) == 2
    assert sorted(origin['transitions'] for family in saved.families for origin in family.provenance) == sorted(
        path.transitions for path in graph.paths())


@pytest.mark.parametrize('seed', range(16))
def test_bulk_atom_sets_match_independent_partition_grouping(seed):
    import random
    from mappa._group_ops import fragment_compositions
    rng = random.Random(seed)
    states, expected = [], []
    for number in range(12):
        atoms = rng.sample(range(10), rng.randrange(1, 11))
        mapping = dict(zip(atoms, rng.sample(range(10), len(atoms))))
        owners = {atom: rng.randrange(4) for atom in atoms if rng.random() < .8}
        groups = {}
        for atom in atoms:
            groups.setdefault(('island', owners[atom]) if atom in owners else ('single', atom), set()).add(atom)
        expected.append(tuple(sorted((tuple(sorted(group)), tuple(sorted(mapping[a] for a in group)))
                                     for group in groups.values())))
        states.append(SearchState(number, number % 2, tuple(mapping.items()), tuple(owners.items()), ()))
    records, compositions = fragment_compositions(states, range(len(states)), (0, 1))
    actual = [tuple(records[index][1:3] for index in ids) for _, ids in compositions]
    assert actual == expected
    for policy, source, target, witness in records:
        assert tuple(sorted(dict(witness))) == source
        assert tuple(sorted(dict(witness).values())) == target


def test_native_projection_reads_shared_blocks_without_expanding_states(monkeypatch):
    from mappa._group_ops import fragment_compositions
    from mappa.mappa_core.core.search.pairs import SharedPairs

    prefix = ((8, 3), (0, 7))
    tails = (((2, 5), (6, 1)), ((2, 1), (6, 5)), ((4, 2),))
    owners = SharedPairs((((8, 0), (0, 1)), ((2, 0), (6, 1))), 4)
    states = [SearchState(i, i % 2, SharedPairs((prefix, tail), 2 + len(tail)), owners, ())
              for i, tail in enumerate(tails * 2)]
    plain = [replace(state, mapping=tuple(state.mapping), islands=tuple(state.islands))
             for state in states]
    expected = fragment_compositions(plain, range(len(plain)), (0, 1))
    def forbidden(*args):
        pytest.fail('Native projection must reuse blocks without expanding each shared sequence')
    monkeypatch.setattr(SharedPairs, '__iter__', forbidden)
    assert fragment_compositions(states, range(len(states)), (0, 1)) == expected


@pytest.mark.parametrize('seed', range(32))
@pytest.mark.parametrize('mixed', [False, True])
def test_representative_compositions_keep_valid_mappings_and_atom_set_scope(seed, mixed):
    import numpy as np
    from mappa import search_aam, AAMSearchConfig
    from growth_cases import represented
    rng = np.random.default_rng(seed)
    pairs = [(a, b) for a in range(4) for b in range(a + 1, 4)]
    source = [(a, b, float(rng.choice([0., 1., 2.]))) for a, b in pairs]
    target = [(a, b, float(rng.choice([0., 1., 2.]))) for a, b in pairs]
    value = problem('CCOH' if mixed else 'CCCC', source, target)
    anchors = ((2, 2),) if mixed and seed % 2 else ()
    original = search_aam(value, AAMSearchConfig(seed_count=2, iso_tolerance=2., branch_limit=100,
                                              anchors=anchors))
    saved = FinalBranchCatalogue(value).add_aam(original)
    before = represented(original.graph)
    after = set().union(*(represented(family.as_path(value)) for family in saved.families))
    # One original relation is retained per atom-set pair. The discarded
    # discovery histories are deliberately not an exhaustive mapping union.
    assert after and after <= before
    assert saved.branch_count == len(original.graph.branches())
    fragments = [fragment for family in saved.families for fragment in family.fragments]
    assert len({id(f) for f in fragments}) == len({f.atom_sets for f in fragments})
    for family in saved.families:
        family.validate_representative(value)


@pytest.mark.parametrize('offset', [0, 6])
def test_reduced_ring_retains_its_symmetry_coverage_and_skips_extra_search(offset):
    from mappa.mappa_core.core.symmetry.optimization import optimize_fragment
    atoms = tuple(range(offset, offset + 6))
    bonds = tuple((i + offset, (i + 1) % 6 + offset) for i in range(6))
    value = problem('C' * (offset + 6), [(a, b, w) for (a, b), w in zip(bonds, (2, 2, 1, 1, 1, 1))],
        [(a, b, w) for (a, b), w in zip(bonds, (2, 1, 2, 1, 1, 1))])
    fragment = saved_fragment(value, source_atoms=atoms, target_atoms=atoms, accepted_edges=bonds)
    before = fragment.to_record()
    saved = FinalBranchCatalogue(value).add_graph(saved_path(value, (fragment,)).graph)
    original, = saved.families[0].fragments
    assert original is fragment and original.to_record() == before
    assert original.relation.discovery_in_witness_orbit
    answer = optimize_fragment(original)
    assert answer.event_count == answer.lower_bound == 2 and answer.complete
    assert answer.work.discovery_nodes == 0
    decoded = decode_events(saved, event_slack=5)
    assert decoded.minimum_events == 2


def test_discarded_history_pruning_preserves_catalogue_and_boundary_counts():
    from mappa.mappa_core.core.search.branching.construction import _Branch
    from mappa.mappa_core.core.search.graph import SearchGraphBuilder
    value = problem('CCC', [(0, 1, 1.)], [(0, 2, 1.)])
    builder = SearchGraphBuilder(SearchContext((0, 1, 2), (0, 1, 2), (0,),
        cuts=((0, 1),), anchors=((0, 0),)))
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
    assert left.source_symmetry is right.source_symmetry
    assert left.target_atoms != right.target_atoms


def shared_source_catalogue():
    from mappa.mappa_core.core.matching.representation.records import Fragment
    value = problem('CCCCCC', [(0, 1, 1.), (0, 2, 2.), (3, 4, 1.), (3, 5, 2.)],
                    [(0, 1, 1.), (0, 2, 2.), (3, 4, 1.), (3, 5, 2.)])
    first = saved_fragment(value, source_atoms=(0, 1, 2), accepted_edges=((0, 1), (0, 2)))
    second = saved_fragment(value, source_atoms=(0, 1, 2), target_atoms=(3, 4, 5),
        witness={0: 3, 1: 4, 2: 5}, accepted_edges=((0, 1), (0, 2)))
    second = Fragment(replace(second.relation, source=first.source_symmetry))
    return FinalBranchCatalogue(value).add_graph(saved_path(value, (first,)).graph).add_graph(
        saved_path(value, (second,)).graph)


def test_binary_catalogue_preserves_shared_preparation_without_record_expansion(tmp_path, monkeypatch):
    import gc
    from weakref import ref
    from mappa import AAMSearchConfig
    from mappa.execution import BundleStore
    from mappa.mappa_core.core.scheduler.results import CatalogueResult
    from mappa.mappa_core.core.symmetry.relations import RelaxedIsomorphism
    import mappa.mappa_core.core.symmetry.structure as structure

    saved = shared_source_catalogue()
    expected = json.loads(json.dumps(saved.to_record()))
    original = ref(saved.families[0].fragments[0].source_symmetry)
    assert original().generators
    assert original() is saved.families[1].fragments[0].source_symmetry
    def forbidden(*args, **kwargs):
        pytest.fail('binary transport must preserve prepared objects without JSON or re-preparation')
    store = BundleStore(tmp_path / 'bundles', scratch=tmp_path)
    with monkeypatch.context() as guard:
        guard.setattr(FinalBranchCatalogue, 'to_record', forbidden)
        guard.setattr(structure, 'structural_generators', forbidden)
        guard.setattr(RelaxedIsomorphism, '__post_init__', forbidden)
        reference = store.put(CatalogueResult('shared', saved, (), (), AAMSearchConfig(), 'native'))
        del saved
        gc.collect()
        assert original() is None
        restored = store.load(reference).catalogue
    first, second = (family.fragments[0] for family in restored.families)
    assert first.source_symmetry is second.source_symmetry
    assert first.target_atoms != second.target_atoms
    assert not first.source_symmetry.weights.flags.writeable
    assert not first.target_symmetry.weights.flags.writeable
    assert json.loads(json.dumps(restored.to_record())) == expected
    assert '_evaluation' not in vars(first.source_symmetry)
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


def placed_fragment_path():
    value = problem('CCCC', [(0, 1, 1.)], [(0, 1, 1.), (2, 3, 1.)])
    fragment = saved_fragment(value, source_atoms=(0, 1), accepted_edges=((0, 1),),
                              iso_tolerance=0.)
    path = saved_path(value, (fragment,))
    graph = replace(path.graph, transitions=(replace(path.graph.transitions[0],
                    placement_generators=((2, 3, 0, 1),)),))
    return value, next(graph.paths())


@pytest.mark.parametrize('transport', ['json', 'bundle'])
def test_placement_groups_preserve_distinct_images_and_shared_transport(tmp_path, transport):
    from mappa.execution import BundleStore
    value, path = placed_fragment_path()
    saved = FinalBranchCatalogue(value).add_graph(path.graph)
    family = saved.families[0]
    assert saved.branch_count == 1 and len(family.fragments) == 1
    assert tuple(kind for kind, _ in family.actions) == ('relation', 'group')
    assert family.fixed_assignments == ()
    expected = {(0, 1), (1, 0), (2, 3), (3, 2)}
    assert recovered(path, value) == expected
    assert recovered(family.as_path(value), value) == expected
    saved.families.append(replace(family, provenance=[dict(cuts=((0, 2),))]))
    next(iter(saved.branches.values())).add(1)
    record = json.loads(json.dumps(saved.to_record()))
    assert len(record['fragments']) == len(record['placements']) == 1
    if transport == 'json':
        restored = FinalBranchCatalogue.from_record(value, record)
    else:
        store = BundleStore(tmp_path / 'placement-bundle', scratch=tmp_path)
        restored = store.load(store.put(saved))
    assert json.loads(json.dumps(restored.to_record())) == record
    assert restored.families[0].actions is restored.families[1].actions
    assert restored.families[0].actions[1][1] is restored.families[1].actions[1][1]
    assert restored.families[0].actions[1].support == frozenset(range(4))
    assert recovered(restored.families[0].as_path(value), value) == expected


def test_placement_only_saved_transition_preserves_order_without_a_fragment():
    value, path = placed_fragment_path()
    saved = FinalBranchCatalogue(value).add_graph(path.graph)
    family = saved.families[0]
    second_group = ((0, 1, 3, 2),)
    family.actions += (('group', second_group),)
    restored_path = family.as_path(value)
    assert any(edge.match is None and edge.placement_generators
               for edge in restored_path.graph.transitions)
    roundtrip = FinalBranchCatalogue(value).add_graph(restored_path.graph)
    assert roundtrip.families[0].actions == family.actions
    assert recovered(restored_path, value) == {(0, 1), (1, 0), (2, 3), (3, 2)}
    assert dict(restored_path.realize(placement_choices={0: (1,), 1: (1,)}).mapping) == {0: 3, 1: 2}


def test_placement_graph_transport_combine_and_realization_preserve_the_program():
    value, path = placed_fragment_path()
    restored = AAMSearchGraph.from_record(json.loads(json.dumps(path.graph.to_record())))
    combined = AAMSearchGraph.combine((path.graph, restored))
    for current in combined.paths():
        edge_id, = current.transitions
        for image in ((0, 1), (1, 0)):
            for word in ((), (1,)):
                mapping = dict(current.realize({edge_id: dict(discovery=dict(zip((0, 1), image)))},
                                               {edge_id: word}).mapping)
                assert mapping == dict(zip((0, 1), tuple(p + (2 if word else 0) for p in image)))
        assert recovered(current, value) == {(0, 1), (1, 0), (2, 3), (3, 2)}


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
