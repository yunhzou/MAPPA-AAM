"""Search-graph contracts independent of optional event/geometry selection."""
import itertools
import json
import random
from dataclasses import replace

import numpy as np
import pytest

from mappa import AAMProblem, AAMSearchConfig, MolecularEndpoint, search_aam
from mappa.mappa_core.utils.io.records import aam_from_record, aam_record
from mappa.mappa_core.core.search.branching.construction import _Branch, find_islands
from mappa.mappa_core.core.matching.graphs import build_graph
from mappa.mappa_core.core.matching.representation.records import Fragment, FragmentMatchConfig, FragmentMatchContext, FragmentMatchResult
from mappa.mappa_core.core.matching.api import match_fragment
from mappa.mappa_core.core.matching.representation.records import IslandBranchLimitExceeded
from mappa.mappa_core.core.search.graph import AAMSearchGraph, SearchContext, SearchGraphBuilder


def network(elements, bonds=()):
    weights = np.zeros((len(elements), len(elements)))
    for left, right in bonds:
        weights[left, right] = weights[right, left] = 1
    return build_graph(elements, weights, bond_cut=0.2)


def iso(mapping):
    from mappa.mappa_core.core.symmetry.structure import binary_symmetry
    from mappa.mappa_core.core.symmetry.relations import RelaxedIsomorphism
    source, target = tuple(sorted(mapping)), tuple(sorted(mapping.values()))
    base = binary_symmetry(('C',) * len(source), np.zeros((len(source), len(source))))
    return Fragment(RelaxedIsomorphism(replace(base, atoms=source), replace(base, atoms=target),
        tuple(mapping.items()), (), tuple((atom, target) for atom in source), (), (), .2))


def growth_result(*fragments):
    return FragmentMatchResult(fragments, False, len(fragments), 100, 0.)


def problem():
    weights = np.array([[0., 1.], [1., 0.]])
    endpoint = MolecularEndpoint(('H', 'H'), [[0, 0, 0], [1, 0, 0]], weights)
    return AAMProblem(endpoint, endpoint, 'hydrogen')


def path_signatures(graph):
    return [(path.context, tuple(path.mapping.items()), tuple(
        (edge.seed, edge.step, edge.match, edge.island)
        for edge in (graph.transitions[index] for index in path.transitions)))
        for path in graph.paths()]


def test_prefix_fork_reconvergence_does_not_copy_or_cross_contexts():
    graph = network(['C'] * 4)
    recorder = SearchGraphBuilder(SearchContext(tuple(graph), tuple(graph), (0, 1, 2, 3)))
    prefix = _Branch(recorder)
    prefix.commit(iso({0: 0}))
    discarded = prefix.fork()
    discarded.commit(iso({2: 3}))
    recorder.stop(discarded, 'capped', stage='combined_live_leaves', count=3, limit=2)
    left, right = prefix.fork(), prefix.fork()
    left.commit(iso({1: 1}))
    left.commit(iso({2: 2}))
    left.commit(iso({1: 1, 2: 2}))
    right.commit(iso({1: 1, 2: 2}))
    assert left.islands_R == right.islands_R
    left.merge_exact_paths(right)
    left.commit(iso({3: 3}))
    recorder.stop(left, 'objective_met')
    expected = path_signatures(recorder.finish())
    recorder.retain([left])
    result = recorder.finish()
    assert path_signatures(result) == expected
    assert all(edge.match is None or dict(edge.match.witness) != {2: 3}
               for edge in result.transitions)
    assert result.states[left.node].mapping == tuple(sorted(left.mapping.items()))
    cap, = [stop for stop in result.stops if stop.reason == 'capped']
    assert (cap.state, result.stop_context(cap), cap.occurrences) == (None, 0, 1)
    paths = tuple(result.paths())
    assert len(paths) == 2
    assert all(p.transitions[0] == 0 for p in paths)
    assert paths[0].transitions[-1] == paths[1].transitions[-1]
    assert len(result.branches()) == 1  # same final fragment pairs
    assert len(result.branches()[0].paths) == 2
    assert len(result.literal_branches()) == 2
    assert all(edge.source < edge.target for edge in result.transitions)
    combined = AAMSearchGraph.combine((result, result))
    assert len(tuple(combined.paths())) == 4  # not a Cartesian product
    assert len(combined.branches()) == 1
    assert len(combined.branches()[0].paths) == 4
    assert len(combined.literal_branches()) == 2


def test_rejected_subtree_removes_its_history_without_losing_shared_survivor(monkeypatch):
    import mappa.mappa_core.core.search.branching.construction as construction
    from mappa.mappa_core.core.matching.representation.records import FragmentMatchResult

    def grow(source, target, *, seed, context, config, **kwargs):
        placements = ([iso({0: 0, 1: 1}), iso({0: 0, 1: 2}), iso({0: 0})]
                      if seed == 0 else
                      [iso({0: 0, 1: value}) for value in (1, 3, 4)])
        return FragmentMatchResult(tuple(placements), False, len(placements),
                                   config.branch_limit, 0.)

    monkeypatch.setattr(construction, 'match_fragment', grow)
    source = network(['C', 'O'])
    target = network(['C', 'O', 'O', 'O', 'O'])
    result = find_islands(source, target, [0, 1], max_branches=3)
    paths = list(result.paths())
    assert [path.mapping for path in paths] == [{0: 0, 1: 1}, {0: 0, 1: 2}]
    assert [[list(edge.match.source_atoms) for edge in
             (result.transitions[index] for index in path.transitions) if edge.match]
            for path in paths] == [[[0, 1]], [[0, 1]]]
    assert all(state.mapping in ((), ((0, 0), (1, 1)), ((0, 0), (1, 2)))
               for state in result.states)
    cap, = [stop for stop in result.stops if stop.reason == 'capped']
    assert (cap.state, result.stop_context(cap), cap.seed, cap.step, cap.stage,
            cap.count, cap.limit, cap.occurrences) == (
                None, 0, 1, (1, 1), 'combined_live_leaves', 4, 3, 1)


def test_retention_releases_discarded_payloads_and_aggregates_stop_evidence():
    import gc
    import weakref

    recorder = SearchGraphBuilder(SearchContext((0, 1), tuple(range(102)), (0, 1)))
    live = _Branch(recorder)
    live.commit(iso({0: 0}))
    payloads = []
    for target in range(1, 101):
        discarded = live.fork()
        match = iso({1: target})
        payloads.append(weakref.ref(match))
        discarded.commit(match)
        recorder.stop(discarded, 'capped', stage='combined_live_leaves', count=2, limit=1)
        recorder.retain([live])
        assert len(recorder.states) == 2 and len(recorder.transitions) == 1
        del discarded, match
    gc.collect()
    assert all(reference() is None for reference in payloads)
    result = recorder.finish()
    cap, = result.stops
    assert (cap.state, cap.count, cap.limit, cap.occurrences) == (None, 2, 1, 100)
    assert not result.terminals and not list(result.paths())
    recorder.stop(live, 'stalled')
    recorder.retain([])
    assert [path.mapping for path in recorder.finish().paths()] == [{0: 0}]


def test_online_admission_reuses_equal_continuation(monkeypatch):
    import mappa.mappa_core.core.matching.api as fragment_module
    calls = []
    def grow(source, target, seed, mapping, **kwargs):
        calls.append((seed, tuple(mapping)))
        if seed == 0:
            return growth_result(iso({0: 0}), iso({0: 0}))
        return growth_result(iso({1: 1}))
    monkeypatch.setattr(fragment_module, 'grow_island', grow)
    graph = network(['C', 'O'])
    result = find_islands(graph, graph, [0, 1], max_branches=1)
    assert calls == [(0, ()), (1, (0,))]
    assert not result.capped
    assert len(tuple(result.paths())) == 2
    assert len(result.branches()) == 1


def test_sibling_commits_share_snapshots_without_mutating_the_parent():
    parent = _Branch(anchor_map={0: 10, 2: 12})
    before = parent.state_key()
    left, right = parent.fork(), parent.fork()
    assert left.islands_R is right.islands_R is parent.islands_R
    left.commit(iso({1: 11}))
    right.commit(iso({1: 13}))
    assert parent.state_key() == before
    assert parent.mapping == {0: 10, 2: 12}
    assert parent.islands_R == {0: 1, 2: 2}
    # Only the source partition is shared: different target witnesses survive.
    assert left.state_key()[1] is right.state_key()[1]
    assert left.state_key() != right.state_key()
    assert left.state_key()[0][0] is before[0][0]
    old_node = left.node
    left.commit(iso({1: 11}))
    assert left.node == old_node  # a no-op doesn't invent a new decision
    assert len(parent.graph.states) == 3


def test_deferred_bonds_belong_to_continuation_states():
    parent = _Branch()
    for mapping, cuts in [({0: 1}, ((0, 2), (1, 3))),
                          ({0: 2}, ((0, 2), (1, 3))),
                          ({0: 3}, ((0, 2),))]:
        branch = parent.fork()
        branch.commit(iso(mapping), deferred_edges=cuts)
        parent.graph.stop(branch, 'stalled')
    graph = parent.graph.finish()
    a, b, c = [graph.states[e.target].deferred_edges for e in graph.transitions]
    assert a == b == ((0, 2), (1, 3))
    assert c == ((0, 2),)
    assert all(not hasattr(e.match, 'deferred_edges') for e in graph.transitions)
    assert len(graph.stops) == 3  # retained partial results keep their own records
    record = graph.to_record()
    restored = AAMSearchGraph.from_record(json.loads(json.dumps(record)), copy=False)
    assert restored.to_record() == record
    x, y, z = [restored.states[e.target].deferred_edges for e in restored.transitions]
    assert x is y and x[0] is z[0]


@pytest.mark.parametrize('seed', range(20))
def test_cached_partitions_and_events_match_literal_commit(seed):
    """Compare the optimized commit against the pre-optimization definition."""
    rng = random.Random(seed)
    atoms = list(range(0, 36, 3))
    assignment = dict(zip(atoms, rng.sample(range(50, 150), len(atoms))))
    anchors = dict((a, assignment[a]) for a in rng.sample(atoms, 3))
    branch = _Branch(anchor_map=anchors)
    mapping, islands = dict(anchors), dict(branch.islands_R)
    cuts, next_iid = set(), branch.next_iid
    for _ in range(30):
        selected = rng.sample(atoms, rng.randrange(1, 6))
        pairs = {a: assignment[a] for a in selected}
        new_cuts = frozenset({tuple(sorted(rng.sample(atoms, 2)))})
        match = iso(pairs)
        pairs = dict(match.witness)
        touched = {islands[a] for a in pairs if a in islands}
        iid = min(touched) if touched else next_iid
        added, relabeled = [], []
        for a, b in pairs.items():
            if a not in mapping:
                mapping[a] = b
                added.append((a, b))
            elif islands[a] != iid:
                relabeled.append((a, islands[a]))
            islands[a] = iid
        for a, label in list(islands.items()):
            if label in touched and label != iid:
                relabeled.append((a, label))
                islands[a] = iid
        groups = {}
        for a, label in islands.items():
            groups.setdefault(label, []).append(a)
        groups = sorted(tuple(sorted(group)) for group in groups.values())
        islands = {a: label for label, group in enumerate(groups, 1) for a in group}
        next_iid = len(groups) + 1
        cuts.update(new_cuts)
        expected_event = dict(type='island_locked', island_idx=iid, pairs=added,
                              merged_with=sorted(touched - {iid}), relabeled=sorted(relabeled),
                              mapped_total=len(mapping))
        events = []
        branch.commit(match, events=events, deferred_edges=new_cuts)
        assert [{**event, 'relabeled': sorted(event['relabeled'])}
                for event in events] == [expected_event]
        assert branch.mapping == mapping
        assert list(branch.islands_R.items()) == list(islands.items())
        assert branch.islands_P == {mapping[a]: label for a, label in islands.items()}
        assert branch.next_iid == next_iid
        assert branch.state_key() == (tuple(sorted(mapping.items())),
                                      tuple(sorted(islands.items())), tuple(sorted(cuts)))



def test_graph_hierarchy_references_prepared_fragments_without_reconstruction():
    graph = network(['C'] * 3)
    recorder = SearchGraphBuilder(SearchContext(tuple(graph), tuple(graph), (0, 1, 2)))
    prefix = _Branch(recorder)
    first = iso({0: 0})
    prefix.commit(first)
    live, capped = prefix.fork(), prefix.fork()
    second = iso({1: 1})
    live.commit(second)
    capped.commit(iso({2: 2}))
    recorder.stop(live, 'objective_met')
    recorder.stop(capped, 'capped')
    recorder.retain([live])
    result = recorder.finish()
    assert len(result.transitions) == 2
    assert result.capped
    assert result.ancestor_transitions(result.terminals) == frozenset({0, 1})
    path = next(result.paths())
    assert path.hierarchy.fragments[0] is first
    assert path.hierarchy.fragments[1] is second
    combined = AAMSearchGraph.combine((result, result))
    assert combined.transitions[0].match is combined.transitions[2].match is first
    restored = AAMSearchGraph.from_record(json.loads(json.dumps(result.to_record())))
    assert next(restored.paths()).hierarchy == path.hierarchy


def test_completed_fragment_cap_is_recorded_without_claiming_a_match(monkeypatch):
    import mappa.mappa_core.core.matching.api as fragment_module
    def grow(source, target, seed, mapping, **kwargs):
        raise IslandBranchLimitExceeded(2, 1, seed=seed)
    monkeypatch.setattr(fragment_module, 'grow_island', grow)
    graph = network(['C', 'C'])
    result = find_islands(graph, graph, [0], max_branches=1)
    assert not result.terminals
    assert result.capped
    assert all(s.stage == 'fragment_growth' and s.count == 2 and s.limit == 1
               for s in result.stops)
    assert len(result.states) == 1 and not result.transitions
    assert all(s.state is None and result.stop_context(s) == 0 for s in result.stops)
    assert sum(s.occurrences for s in result.stops) == 1


def test_combined_cap_preserves_partial_matches_and_stall_status():
    source = network(['C', 'O', 'H'])
    target = network(['C', 'C', 'O', 'O', 'N', 'F', 'Cl', 'Br'],
                     [(0, 4), (1, 5), (2, 6), (3, 7)])
    result = find_islands(source, target, [0, 1, 2], iso_tol=.1, max_branches=3)
    assert [path.mapping for path in result.paths()] == [
        {0: 0, 1: 2}, {0: 0, 1: 3}]
    cap, = [stop for stop in result.stops if stop.reason == 'capped']
    assert (cap.stage, cap.count, cap.limit, cap.occurrences) == (
        'combined_live_leaves', 4, 3, 1)
    assert [stop.reason for stop in result.stops if stop.state is not None] == [
        'stalled', 'stalled']


def test_partial_fragment_api_and_sparse_original_atom_indices():
    source = network(['C', 'O', 'H'], [(0, 1), (0, 2)]).subgraph([0, 2]).copy()
    target = network(['C', 'O', 'H', 'C', 'H'], [(0, 2), (3, 4)]).subgraph([0, 2, 3, 4]).copy()
    config = FragmentMatchConfig()
    context = FragmentMatchContext()
    result = match_fragment(source, target, seed=0, context=context, config=config)
    assert result.matches and not result.capped
    assert all(set(m.source_atoms) == {0, 2} for m in result.matches)
    graph = AAMSearchGraph.initial_fragment_search(source, target, 0, result, config)
    restored = AAMSearchGraph.from_record(json.loads(json.dumps(graph.to_record())))
    for path in restored.paths():
        realization = path.realize()
        mapping = dict(realization.mapping)
        assert set(mapping) == {0, 2}
        assert target.has_edge(mapping[0], mapping[2])


def test_tiny_exact_symmetry_matches_exhaustive_ring_automorphisms():
    bonds = [(0, 1), (1, 2), (2, 3), (3, 0)]
    graph = network(['C'] * 4, bonds)
    result = find_islands(graph, graph, [0], iso_tol=.5)
    path = next(result.paths())
    edge = next(e for e in result.transitions if e.match is not None)
    generators = edge.match.target_symmetry.generators
    orbit = {tuple(path.mapping[a] for a in range(4))}
    queue = list(orbit)
    for mapping in queue:
        for generator in generators:
            moved = tuple(generator[a] for a in mapping)
            if moved not in orbit:
                orbit.add(moved)
                queue.append(moved)
    exact = {p for p in itertools.permutations(range(4))
             if all(graph.has_edge(p[a], p[b]) for a, b in bonds)}
    assert orbit == exact
    assert len(exact) == 8
    for seed in range(10):
        word = tuple(1 + random.Random(seed + i).randrange(len(generators)) for i in range(4))
        sample = path.realize({edge.id: {'target_word': word}})
        assert tuple(dict(sample.mapping)[a] for a in range(4)) in exact


def test_saved_raw_aam_is_independent_of_postprocessing(tmp_path, monkeypatch):
    import mappa.mappa_core.post_processing.decoder.api as postprocessing
    def forbidden(*args, **kwargs):
        raise AssertionError('raw search must not score mechanisms')
    monkeypatch.setattr(postprocessing, 'decode_events', forbidden)
    base = problem()
    enriched = AAMProblem(replace(base.reactant, energy=-1.2, metadata={'origin': 'saved-test'}),
                          base.product)
    result = search_aam(enriched, AAMSearchConfig(seed_count=2), intermediate_dir=tmp_path)
    restored = aam_from_record(json.loads((tmp_path / 'aam.json').read_text()))
    assert restored.problem.reactant.energy == -1.2
    assert restored.problem.reactant.metadata['origin'] == 'saved-test'
    assert restored.graph.to_record() == aam_from_record(
        json.loads(json.dumps(aam_record(result)))).graph.to_record()
    assert list(tmp_path.glob('*.bundle'))
    assert not list(tmp_path.glob('cut_*'))
    assert len(restored.branches) == len(result.branches)


def test_serial_and_workers_produce_identical_search_graphs():
    config = AAMSearchConfig(seed_count=2)
    serial = search_aam(problem(), config, workers=1)
    parallel = search_aam(problem(), config, workers=2)
    assert serial.graph.to_record() == parallel.graph.to_record()


def test_local_relation_choices_transport_witnesses_and_preserve_anchors():
    from mappa.mappa_core.core.symmetry.structure import binary_symmetry
    from mappa.mappa_core.core.symmetry.relations import RelaxedIsomorphism
    graph = network(['C'] * 4)
    recorder = SearchGraphBuilder(SearchContext(tuple(graph), tuple(graph), (0, 2)))
    branch = _Branch(recorder)
    branch.commit(iso({0: 0, 1: 1}))
    anchored = binary_symmetry(('C',) * 4, np.zeros((4, 4)),  anchors=(0, 1))
    later = Fragment(RelaxedIsomorphism(anchored, anchored, tuple(enumerate(range(4))),
        ((0, 0), (1, 1)), ((0, (0,)), (1, (1,)), (2, (2, 3)), (3, (2, 3))), (), (), .2))
    branch.commit(later)
    recorder.stop(branch, 'objective_met')
    path = next(recorder.finish().paths())
    sample = path.realize({0: {'discovery': {0: 1, 1: 0}},
                           1: {'discovery': {0: 0, 1: 1, 2: 3, 3: 2}}})
    assert dict(sample.mapping) == {0: 1, 1: 0, 2: 3, 3: 2}
    assert dict(sample.hierarchy.fragments[0].witness) == {0: 1, 1: 0}
    with pytest.raises(ValueError, match='domain|anchor'):
        path.realize({1: {'discovery': {0: 1, 1: 0, 2: 2, 3: 3}}})


def test_capped_sibling_does_not_discard_a_successful_path(monkeypatch):
    import mappa.mappa_core.core.matching.api as fragment_module
    def grow(source, target, seed, mapping, **kwargs):
        if seed == 0:
            return growth_result(iso({0: 0}), iso({0: 1}))
        if mapping[0] == 0:
            raise IslandBranchLimitExceeded(3, 2, seed=seed)
        return growth_result(iso({1: 0}))
    monkeypatch.setattr(fragment_module, 'grow_island', grow)
    graph = network(['C', 'C'])
    result = find_islands(graph, graph, [0, 1], max_branches=2)
    assert result.capped
    assert [p.mapping for p in result.paths()] == [{0: 1, 1: 0}]


def test_offline_viewer_uses_saved_graph_without_new_matching(tmp_path, monkeypatch):
    import re
    import shutil
    import subprocess
    import mappa.mappa_core.core.matching.api as fragment_module
    from mappa.presentation.viewers import write_aam_bundle
    result = search_aam(problem(), AAMSearchConfig(seed_count=1))
    def forbidden(*args, **kwargs):
        raise AssertionError('viewer must not invoke AAM')
    monkeypatch.setattr(fragment_module, 'grow_island', forbidden)
    write_aam_bundle(result, tmp_path)
    page = (tmp_path / 'search.html').read_text()
    assert '__DATA__' not in page and '__LIBRARY__' not in page
    assert '<script src=' not in page
    assert 'Next recorded history' in page
    if shutil.which('node'):
        script = tmp_path / 'viewer.js'
        script.write_text(re.findall(r'<script>(.*?)</script>', page, re.S)[-1])
        subprocess.run(['node', '--check', str(script)], check=True, capture_output=True)
