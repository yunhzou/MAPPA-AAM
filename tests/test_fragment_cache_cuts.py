"""Exact differential tests, including correlated symmetry and cap stops."""
import numpy as np
import pytest

from mappa.mappa_core.core.matching.reuse import FragmentPool
from mappa.mappa_core.core.matching.graphs import build_graph
from mappa.mappa_core.core.matching.growth import native
from mappa.mappa_core.core.search.branching.construction import find_islands

pytestmark = pytest.mark.skipif(not native.built(), reason="native engine not built")


from growth_cases import graph, raw

@pytest.mark.parametrize('tol',[.5,1.])
def test_every_cut_matches_fresh_native(tol):
    rng = np.random.default_rng(819)
    for trial in range(4):
        edges = [(i,i+1,float(rng.choice([.8,1.2,1.8]))) for i in range(7)]
        edges += [(0,4,1.),(2,6,1.)]
        source=graph(8,edges)
        target=graph(8,edges if trial==0 else [(a,b,w) for a,b,w in edges[:-1]])
        pool=FragmentPool(source,target)
        cuts=[()]+[(e[:2],) for e in edges]+[((1,2),(4,5))]
        for seed in (0,4):
            for cap in (1,3,100):
                for cut in cuts:
                    expected=raw(source,target,seed,[-1]*8,cut,cap=cap,tol=tol)
                    actual=raw(source,target,seed,[-1]*8,cut,pool=pool,cap=cap,tol=tol)
                    assert actual==expected,(trial,seed,cap,cut)


def test_unread_component_can_change_source_symmetry_and_still_reuse():
    source=graph(6,[(0,1),(1,2),(3,4),(4,5)])
    target=graph(6,[(0,1),(1,2),(3,4),(4,5)])
    pool=FragmentPool(source,target)
    for cut in ((),((4,5),)):
        assert raw(source,target,0,[-1]*6,cut,pool=pool)==raw(source,target,0,[-1]*6,cut)
    assert pool.engine.exists([0, 1], [(0, 1)], [-1]*6, 1.)
    assert pool.stats()['fragment_hits'] > 0


def test_locks_islands_and_boundaries_are_respected_by_fresh_growth():
    source=graph(8,[(i,i+1) for i in range(7)]+[(0,3)])
    target=graph(8,[(i,i+1) for i in range(7)]+[(0,3)])
    pool=FragmentPool(source,target)
    for mapped in ((0,1),(1,0)):
        mapping=[*mapped]+[-1]*6
        for islands in ([(0,1),(1,1)],[(0,1),(1,2)]):
            for deferred in ((),((1,2),)):
                for cut in ((),((5,6),)):
                    kw=dict(islands=islands,deferred=deferred)
                    assert raw(source,target,7,mapping,cut,pool=pool,**kw)==raw(source,target,7,mapping,cut,**kw)


def test_cache_budget_does_not_limit_search():
    source=graph(4,[(0,1),(1,2),(2,3)])
    pool=FragmentPool(source,source,cache_bytes=1)
    for cut in ((),((1,2),)):
        assert raw(source,source,0,[-1]*4,cut,pool=pool)==raw(source,source,0,[-1]*4,cut)
    assert pool.stats()['entries']==0


def test_snapshot_preserves_completion_measurements():
    ring = graph(6, [(r, (r + 1) % 6) for r in range(6)])
    original = FragmentPool(ring, ring)
    expected = raw(ring, ring, 0, [-1] * 6, (), pool=original, tol=.5)
    assert expected['profile']['completion_generated'] > 0
    restored = FragmentPool(ring, ring)
    restored.loads(original.dumps())
    actual = raw(ring, ring, 0, [-1] * 6, (), pool=restored, tol=.5)
    assert actual == expected
    assert restored.engine.exists([0, 1], [(0, 1)], [-1]*6, .5)
    assert restored.stats()['fragment_hits'] > 0


def test_invalid_cut_does_not_change_session_topology():
    source=graph(4,[(0,1),(1,2),(2,3)])
    pool=FragmentPool(source,source)
    for _ in range(2):
        with pytest.raises(ValueError,match='not an active edge'):
            raw(source,source,0,[-1]*4,((0,1),(0,3)),pool=pool)
    for cut in ((),((1,2),)):
        assert raw(source,source,0,[-1]*4,cut,pool=pool)==raw(source,source,0,[-1]*4,cut)


def test_full_search_graph_is_identical():
    source=graph(8,[(i,i+1) for i in range(7)]+[(0,3)])
    target=graph(8,[(i,i+1) for i in range(7)]+[(2,5)])
    pool=FragmentPool(source,target)
    for cut in ((),((3,4),),((0,1),),((0,3),)):
        view=pool.for_cut(cut)
        options=dict(cuts=cut,max_branches=100)
        expected=find_islands(view.source,target,list(range(8)),**options)
        actual=find_islands(view.source,target,list(range(8)),fragment_cache=view,**options)
        assert actual==expected
