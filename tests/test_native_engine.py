"""Differential tests: the compiled growth engine against the Python engine.

Every ``grow_island`` call made by a full ``search_aam`` run on the
tetraphenylmethane benchmark is executed by both engines and the results
(mappings, deferred edges, fragments, symmetry states, cap behaviour) must be
identical; the typed search results of the two modes must agree as well.
"""
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "bench"))

from mappa.mappa_core.core.matching.growth import native  # noqa: E402

pytestmark = pytest.mark.skipif(
    not native.built(), reason="native growth engine not built")


def _result_key(result):
    matches = [fragment.to_record() for fragment in result.matches]
    return (matches, tuple(sorted(result.deferred_edges)), result.capped,
            result.branch_count, result.branch_limit)


def _search(case, workers=1):
    from native_cases import CASES
    from mappa.mappa_core.core.search.api import search_aam
    from mappa.mappa_core.domain import AAMProblem, AAMSearchConfig

    R, P = CASES[case]()
    return search_aam(AAMProblem(R, P, name=case), AAMSearchConfig(),
                      workers=workers, execution='reference')


def test_every_growth_call_agrees_with_python(monkeypatch):
    import mappa.mappa_core.core.matching.api as fragment_mod
    import mappa.mappa_core.core.matching.growth.engine as island_mod
    from mappa.mappa_core.core.matching.representation.records import IslandBranchLimitExceeded
    from mappa.mappa_core.core.matching.policy import as_node_match_policy

    original = island_mod.grow_island
    stats = {"calls": 0, "native": 0}

    def both(g_R, g_P, seed, mapping, **kw):
        stats["calls"] += 1
        monkeypatch.setenv("MAPPA_NATIVE", "0")
        try:
            try:
                expected = ("ok", _result_key(original(g_R, g_P, seed, dict(mapping), **kw)))
            except IslandBranchLimitExceeded as exc:
                expected = ("raised", (exc.count, exc.limit))
        finally:
            monkeypatch.setenv("MAPPA_NATIVE", "1")
        policy = as_node_match_policy(kw.get("node_policy"))
        if native.available():
            native.validate_inputs(policy, kw.get("events"))
            stats["native"] += 1
            try:
                out = native.grow_island(
                    g_R, g_P, seed, dict(mapping),
                    graph_floor=kw.get("graph_floor", 0.2),
                    iso_tol=kw.get("iso_tol", 0.5),
                    min_lock_size=kw.get("min_lock_size", 1),
                    max_branches=kw.get("max_branches", 1_000_000),
                    islands_R=kw.get("islands_R"),
                    prior_deferred_edges=kw.get("prior_deferred_edges"),
                    allow_mapped_seed=kw.get("allow_mapped_seed", False),
                    profile=None, profile_context=None,
                    self_tolerance=kw.get("self_tolerance", 1.0))
                got = ("ok", _result_key(out))
            except IslandBranchLimitExceeded as exc:
                got = ("raised", (exc.count, exc.limit))
            assert got == expected, (seed, len(mapping))
        return original(g_R, g_P, seed, dict(mapping), **kw)

    monkeypatch.setattr(fragment_mod, "grow_island", both)
    result = _search("tetraphenyl")
    assert result.branches
    # Search policy changes can change the call count. Require differential
    # coverage of every actual growth request, not a historical count floor.
    assert stats["calls"] > 0
    assert stats["native"] == stats["calls"]


def test_search_results_agree_between_engines(monkeypatch):
    monkeypatch.setenv("MAPPA_NATIVE", "1")
    with_native = _search("tetraphenyl").graph.to_record()
    monkeypatch.setenv("MAPPA_NATIVE", "0")
    with_python = _search("tetraphenyl").graph.to_record()
    assert with_native == with_python


def test_reused_target_rejects_source_elements_introduced_later():
    """A cached target must safely handle new catalog element codes."""
    from mappa import _engine

    target = _engine.TargetGraph(
        ["NativeTargetC", "NativeTargetH"],
        [[0.0, 1.0], [1.0, 0.0]],
        0.2,
        [(0, 1)],
    )
    for index in range(8):
        source = _engine.SourceGraph(
            ["NativeTargetC", f"PreviouslyUnseenElement{index}"],
            [[0.0, 1.0], [1.0, 0.0]],
            0.2,
            [(0, 1)],
        )
        result = _engine.grow_island(
            source, target, 0, [-1, -1], 0.2, 0.5, 1, 100,
            None, [], False,
        )
        assert result["capped"] is False
        assert result["fragments"]


def test_sparse_atom_ids_use_native_engine_and_restore_original_ids():
    """Augmented component graphs retain parent atom IDs across native calls."""
    import networkx as nx
    import numpy as np


    source = nx.Graph()
    source.add_nodes_from(((10, {"element": "C"}),
                           (20, {"element": "O"})))
    source.add_edge(10, 20, wbo=1.0)
    source_wbo = np.zeros((21, 21))
    source_wbo[10, 20] = source_wbo[20, 10] = 1.0
    source.graph.update(wbo_matrix=source_wbo, bond_cut=0.2)

    target = nx.Graph()
    target.add_nodes_from(((100, {"element": "C"}),
                           (110, {"element": "O"}),
                           (120, {"element": "C"})))
    target.add_edge(100, 110, wbo=1.0)
    target_wbo = np.zeros((121, 121))
    target_wbo[100, 110] = target_wbo[110, 100] = 1.0
    target.graph.update(wbo_matrix=target_wbo, bond_cut=0.2)

    result = native.grow_island(
        source, target, 10, {}, graph_floor=0.2, iso_tol=0.5,
        min_lock_size=1, max_branches=100, islands_R=None,
        prior_deferred_edges=None, allow_mapped_seed=False,
        profile=None, profile_context=None,
    )

    assert [dict(match.witness) for match in result.matches] == [{10: 100, 20: 110}]
    assert result.matches[0].source_atoms == (10, 20)


@pytest.mark.parametrize('selection', [None, '1'])
def test_missing_selected_engine_is_an_error(selection, monkeypatch):
    from mappa.mappa_core.core.matching.growth import native
    monkeypatch.setattr(native, '_engine', None)
    if selection is None:
        monkeypatch.delenv('MAPPA_NATIVE', raising=False)
    else:
        monkeypatch.setenv('MAPPA_NATIVE', selection)
    with pytest.raises(ImportError, match='not built'):
        native.available()
    monkeypatch.setenv('MAPPA_NATIVE', '0')
    assert not native.available()


def test_invalid_backend_is_not_an_implicit_selection(monkeypatch):
    from mappa.mappa_core.core.matching.growth import native
    monkeypatch.setenv('MAPPA_NATIVE', 'auto')
    with pytest.raises(ValueError, match='must be 0'):
        native.available()


@pytest.mark.parametrize('request_kind', ['trace', 'custom_policy'])
def test_unsupported_native_request_never_dispatches_to_python(request_kind, monkeypatch):
    from mappa import FragmentMatchConfig, match_fragment
    from mappa.mappa_core.core.matching.graphs import build_graph
    from mappa.mappa_core.core.matching.growth import engine
    import numpy as np
    graph = build_graph(('C', 'C'), np.array([[0., 1.], [1., 0.]]))
    config = FragmentMatchConfig(
        node_policy=('element',) if request_kind == 'custom_policy' else None)
    events = [] if request_kind == 'trace' else None
    original = engine.grow_island_python
    monkeypatch.setattr(engine, 'grow_island_python',
        lambda *args, **kwargs: pytest.fail('native request silently selected Python'))
    monkeypatch.setenv('MAPPA_NATIVE', '1')
    with pytest.raises(ValueError, match='explicitly select MAPPA_NATIVE=0'):
        match_fragment(graph, graph, seed=0, config=config, events=events)
    monkeypatch.setattr(engine, 'grow_island_python', original)
    monkeypatch.setenv('MAPPA_NATIVE', '0')
    assert match_fragment(graph, graph, seed=0, config=config, events=events).matches


def test_retired_symmetry_repair_kernel_is_not_exported():
    from mappa import _engine
    assert not hasattr(_engine, 'repair_group')
    assert not hasattr(_engine, 'repair_scores')
    assert not hasattr(_engine, 'round12')


def test_native_execution_error_is_propagated_without_retry(monkeypatch):
    from growth_cases import graph
    from mappa import match_fragment
    from mappa.mappa_core.core.matching.growth import engine, native

    failure = RuntimeError('native matching failure')

    def fail_native(*args, **kwargs):
        raise failure

    def forbid_python(*args, **kwargs):
        pytest.fail('native execution failure was retried through Python')

    monkeypatch.setenv('MAPPA_NATIVE', '1')
    monkeypatch.setattr(native._engine, 'grow_island', fail_native)
    monkeypatch.setattr(engine, 'grow_island_python', forbid_python)
    ring = graph(6, [(i, (i + 1) % 6) for i in range(6)])
    with pytest.raises(RuntimeError) as caught:
        match_fragment(ring, ring, seed=0)
    assert caught.value is failure
