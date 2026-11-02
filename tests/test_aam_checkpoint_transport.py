import json
from dataclasses import replace
import pytest

from mappa.mappa_core.core.scheduler import search as aam
from mappa.mappa_core.core.search.session import _GrowthCounts
from mappa.mappa_core.core.search.graph import (
    AAMSearchGraph, SearchContext, SearchGraphBuilder, SearchState, SearchStop)
from mappa.mappa_core.utils.io.checkpoints import raw_cut_paths, read_raw_cut, write_raw_cut


def test_cut_worker_persists_graph_and_returns_only_reference(tmp_path, monkeypatch):
    from mappa.mappa_core.core.scheduler.bundles import read_blob
    graph = SearchGraphBuilder(SearchContext((), (), ())).finish()
    monkeypatch.setattr(aam, '_CHECKPOINT_STORE', None)
    monkeypatch.setattr(aam, '_search_cut', lambda cut: (graph, {
        'search_seconds': 1., 'max_live_branches': 2, 'max_growth_candidates': 3}))
    for index in range(20):
        returned, reference, counts = aam._search_cut_task((index, (), str(tmp_path)))
        assert returned == index and read_blob(reference, tmp_path) == graph
        assert counts['checkpoint'] == reference
    assert len(list(tmp_path.glob('*.bundle'))) == 1
    _, owned, _ = aam._search_cut_task((21, (), str(tmp_path)), in_process=True)
    assert owned is graph


def test_in_memory_worker_and_online_profile_reducer(monkeypatch):
    graph = SearchGraphBuilder(SearchContext((), (), ())).finish()
    monkeypatch.setattr(aam, '_search_cut', lambda cut: (graph, {'search_seconds': 0.}))
    assert aam._search_cut_task((4, (), None)) == (4, graph, {
        'search_seconds': 0., 'checkpoint_seconds': 0.})
    reducer = _GrowthCounts()
    rows = [{}, {'max_cands_before': 3}, {'result': 'subtree_branch_cap'}, {'max_cands_before': 2}]
    for row in rows: reducer.append(row)
    assert reducer.maximum == max(row.get('max_cands_before', 0) for row in rows)


def test_raw_cut_formats_and_order_are_explicit(tmp_path):
    graph = SearchGraphBuilder(SearchContext((), (), ())).finish()
    binary=tmp_path/'cut_00002.raw.pkl.gz';text=tmp_path/'cut_00001.json'
    write_raw_cut(graph,binary);write_raw_cut(graph,text)
    assert raw_cut_paths(tmp_path)==[text,binary]
    assert read_raw_cut(binary)==read_raw_cut(text)==graph
    with pytest.raises(ValueError,match='Unsupported raw cut format'):
        write_raw_cut(graph,tmp_path/'unknown.bin')
    write_raw_cut(graph,tmp_path/'cut_00002.json')
    with pytest.raises(ValueError,match='Duplicate raw cut'):
        raw_cut_paths(tmp_path)


def test_weighted_cap_summaries_roundtrip_without_inventing_paths(tmp_path):
    context = SearchContext((0,), (0, 1), (0,))
    graph = AAMSearchGraph((context,), (0,), (SearchState(0, 0, (), (), ()),), (), (
        SearchStop(None, 'capped', 0, 'combined_live_leaves', 3, 2, (1, 0),
                   context=0, occurrences=7),))
    for suffix in ('raw.pkl.gz', 'json'):
        path = tmp_path / f'cut_00000.{suffix}'
        write_raw_cut(graph, path)
        restored = read_raw_cut(path)
        assert restored == graph
        assert not restored.terminals and not list(restored.paths())
        stop, = restored.stops
        assert (restored.stop_context(stop), stop.count, stop.occurrences) == (0, 3, 7)
    combined = AAMSearchGraph.combine((graph, graph))
    assert [combined.stop_context(stop) for stop in combined.stops] == [0, 1]
    assert sum(stop.occurrences for stop in combined.stops) == 14
    assert not list(combined.paths())


def test_legacy_attached_stops_keep_nonzero_context_through_combine():
    contexts = (SearchContext((0,), (0,), (0,)),
                SearchContext((0,), (0,), (0,), cuts=((0, 1),)))
    legacy = AAMSearchGraph(contexts, (0, 1), (
        SearchState(0, 0, (), (), ()), SearchState(1, 1, (), (), ())), (), (
            SearchStop(1, 'capped', 0, 'fragment_growth', 3, 2, (1, 0)),))
    record = legacy.to_record()
    for stop in record['stops']:
        stop.pop('context')
        stop.pop('occurrences')
    restored = AAMSearchGraph.from_record(json.loads(json.dumps(record)))
    assert restored == legacy and restored.stop_context(restored.stops[0]) == 1
    detached = replace(legacy, stops=(replace(legacy.stops[0], state=None,
                                             context=1, occurrences=4),))
    combined = AAMSearchGraph.combine((restored, detached, restored))
    assert [stop.state for stop in combined.stops] == [1, None, 5]
    assert [combined.stop_context(stop) for stop in combined.stops] == [1, 3, 5]
    assert [stop.occurrences for stop in combined.stops] == [1, 4, 1]
    assert AAMSearchGraph.from_record(combined.to_record()) == combined
