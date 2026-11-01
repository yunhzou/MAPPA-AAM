"""Differential tests for dependency-directed conditional fragment repair."""
import numpy as np
import pytest

from mappa.mappa_core.core.search.branching.construction import find_islands
from mappa.mappa_core.core.search.branching.reuse import FragmentRepair
from mappa.mappa_core.core.matching.growth import native
from mappa.mappa_core.core.matching.symmetry.orbits import _nauty_orbits
from mappa.mappa_core.core.search.branching.symmetry import finalize_graph_symmetry
from test_cut_replay import graph, raw as native_raw

pytestmark=pytest.mark.skipif(not native.built(),reason='native engine not built')


@pytest.fixture(params=['incremental', 'existence'])
def growth_policy(request):
    return request.param


def raw(*args, **kwargs):
    result = native_raw(*args, **kwargs)
    result.pop('profile')  # Timing and operation counts are not fragment results.
    return result


def test_reuse_later_fragment_after_independent_locked_mapping_changes(growth_policy):
    source=graph(10,[(0,1),(1,2),(2,3),(4,5),(6,7),(7,8),(8,9)])
    orbits=_nauty_orbits(source,wbo_tol=1.)
    repair=FragmentRepair(source,source,orbits,growth_policy=growth_policy)
    for mapping,islands in (([-1]*4+[4,5]+[-1]*4,[(4,1),(5,1)]),
                            ([-1]*4+[5,4]+[-1]*4,[(4,1),(5,2)])):
        assert raw(source,source,orbits,0,mapping,(),islands=islands,replay=repair)==raw(
            source,source,orbits,0,mapping,(),islands=islands)
    assert repair.stats()['history_hits']==1


def test_rebase_unrelated_previous_boundaries(growth_policy):
    source=graph(10,[(0,1),(1,2),(2,3),(4,5),(6,7),(7,8),(8,9)])
    orbits=_nauty_orbits(source,wbo_tol=1.)
    repair=FragmentRepair(source,source,orbits,growth_policy=growth_policy)
    for deferred in (((4,5),),((6,7),),()):
        assert raw(source,source,orbits,0,[-1]*10,(),deferred=deferred,replay=repair)==raw(
            source,source,orbits,0,[-1]*10,(),deferred=deferred)
    assert repair.stats()['boundary_rebases']==2


def test_changed_occupancy_cannot_hide_new_seed_alternatives(growth_policy):
    source=graph(6,[(0,1),(1,2),(3,4),(4,5)])
    orbits=_nauty_orbits(source,wbo_tol=1.)
    repair=FragmentRepair(source,source,orbits,growth_policy=growth_policy)
    expected=[]
    for mapping in ([-1]*3+[3,4,5],[-1]*6):
        fresh=raw(source,source,orbits,0,mapping,())
        assert raw(source,source,orbits,0,mapping,(),replay=repair)==fresh
        expected.append(fresh)
    assert expected[0]!=expected[1]
    assert repair.stats()['full_hits']==0


@pytest.mark.parametrize('tol',[.5,1.])
@pytest.mark.parametrize('budget',[1,64*1024*1024])
def test_changed_contexts_and_cuts_match_fresh_results(tol,budget,growth_policy):
    rng=np.random.default_rng(93852)
    for trial in range(5):
        elements=['C']*6+['H']*3
        edges=[(i,i+1,float(rng.choice([.8,1.2,1.8]))) for i in range(5)]
        edges += [(0,6,1.),(2,7,1.),(4,8,1.)]
        if trial%2:edges.append((0,4,1.))
        source=graph(9,edges,elements)
        target=graph(9,edges[:-1]+[(1,5,.8)],elements)
        orbits=_nauty_orbits(target,wbo_tol=tol)
        repair=FragmentRepair(source,target,orbits,cache_bytes=budget,growth_policy=growth_policy)
        for _ in range(40):
            mapping=[-1]*9
            locked=list(map(int,rng.choice(6,3,replace=False)))
            for a,b in zip(locked,rng.permutation([1,3,5])):mapping[a]=int(b)
            islands=[(r,1+int(rng.integers(2))) for r in locked]
            cut=() if rng.integers(3)==0 else (edges[int(rng.integers(len(edges)))][:2],)
            deferred=() if rng.integers(2)==0 else (edges[int(rng.integers(len(edges)))][:2],)
            seed=int(rng.integers(9));cap=int(rng.choice([1,3,100]))
            kw=dict(cap=cap,tol=tol,islands=islands,deferred=deferred)
            assert raw(source,target,orbits,seed,mapping,cut,replay=repair,**kw)==raw(
                source,target,orbits,seed,mapping,cut,**kw),(trial,mapping,islands,cut)


def test_complete_compressed_search_graphs_and_generators(growth_policy):
    source=graph(9,[(0,1),(1,2),(2,3),(3,4),(4,5),(0,6),(2,7),(4,8)],['C']*6+['H']*3)
    target=graph(9,[(0,1),(1,2),(2,3),(3,4),(4,5),(1,6),(3,7),(5,8)],['C']*6+['H']*3)
    orbits=_nauty_orbits(target,wbo_tol=1.)
    repair=FragmentRepair(source,target,orbits,growth_policy=growth_policy)
    for order in (list(range(9)),list(reversed(range(9))),[4,2,6,0,1,3,5,7,8]):
        for cut in ((),((2,3),),((0,6),),((0,1),)):
            view=repair.for_cut(cut)
            options=dict(cuts=cut,max_branches=100,p_orbits=orbits,iso_tol=1.)
            expected=find_islands(view.source,target,order,**options)
            actual=find_islands(view.source,target,order,growth_replay=view,growth_policy=growth_policy,**options)
            assert actual==expected
            expected,_=finalize_graph_symmetry(expected,target,iso_tolerance=1.)
            actual,_=finalize_graph_symmetry(actual,target,iso_tolerance=1.)
            assert actual==expected


def test_eviction_bounds_retention_without_changing_results(growth_policy):
    source=graph(8,[(i,(i+1)%8) for i in range(8)])
    orbits=_nauty_orbits(source,wbo_tol=1.)
    probe=FragmentRepair(source,source,orbits,growth_policy=growth_policy)
    raw(source,source,orbits,0,[-1]*8,(),replay=probe)
    budget=2*probe.stats()['resident_bytes']
    repair=FragmentRepair(source,source,orbits,cache_bytes=budget,growth_policy=growth_policy)
    for _ in range(3):
        for seed in range(8):
            assert raw(source,source,orbits,seed,[-1]*8,(),replay=repair)==raw(
                source,source,orbits,seed,[-1]*8,())
            stats=repair.stats()
            assert 0 < stats['resident_bytes'] <= stats['peak_resident_bytes'] <= budget
    assert repair.stats()['evictions'] > 0
    assert repair.stats()['entries'] <= 2
    before=repair.stats()
    raw(source,source,orbits,7,[-1]*8,(),replay=repair)
    after=repair.stats()
    assert after['full_hits']==before['full_hits']+1
    assert after['resident_bytes']==before['resident_bytes']


def test_oversized_results_are_returned_without_retention(growth_policy):
    source=graph(8,[(i,i+1) for i in range(7)])
    orbits=_nauty_orbits(source,wbo_tol=1.)
    repair=FragmentRepair(source,source,orbits,cache_bytes=1,growth_policy=growth_policy)
    for _ in range(3):
        assert raw(source,source,orbits,0,[-1]*8,(),replay=repair)==raw(
            source,source,orbits,0,[-1]*8,())
    stats=repair.stats()
    assert stats['entries']==stats['resident_bytes']==stats['peak_resident_bytes']==0
    assert stats['oversized_results']==3


def test_public_fragment_cache_reuses_results_and_rejects_policy_mismatch():
    from mappa import match_fragment, FragmentMatchConfig, FragmentMatchContext
    from test_existence_growth import record
    source=graph(6,[(i,(i+1)%6) for i in range(6)])
    orbits=_nauty_orbits(source,wbo_tol=1.)
    cache=FragmentRepair(source,source,orbits,growth_policy='existence')
    view=cache.for_cut(())
    context=FragmentMatchContext(growth_replay=view)
    config=FragmentMatchConfig(iso_tolerance=1.,growth_policy='existence')
    expected=record(match_fragment(view.source,source,seed=0,config=config))
    for _ in range(3):
        assert record(match_fragment(view.source,source,seed=0,context=context,config=config))==expected
    assert cache.stats()['full_hits']==2
    with pytest.raises(ValueError,match='growth_policy'):
        match_fragment(view.source,source,seed=0,context=context,
                       config=FragmentMatchConfig(growth_policy='incremental'))


def test_snapshot_preserves_results_dependencies_and_caps(growth_policy):
    source=graph(9,[(0,1),(1,2),(2,3),(3,4),(4,5),(0,6),(2,7),(4,8)],['C']*6+['H']*3)
    orbits=_nauty_orbits(source,wbo_tol=1.)
    cache=FragmentRepair(source,source,orbits,growth_policy=growth_policy)
    requests=[(0,[-1]*9,(),dict(cap=cap)) for cap in (1,3,100)]
    requests += [(0,[-1]*6+[6,7,8],((2,3),),dict(cap=3,islands=[(6,1),(7,1),(8,2)])),
                 (2,[-1]*9,((0,1),),dict(cap=100,deferred=((2,3),)))]
    for seed,mapping,cut,options in requests:
        raw(source,source,orbits,seed,mapping,cut,replay=cache,**options)
    restored=FragmentRepair(source,source,orbits,growth_policy=growth_policy)
    restored.loads(cache.dumps())
    for seed,mapping,cut,options in requests:
        assert raw(source,source,orbits,seed,mapping,cut,replay=restored,**options)==raw(
            source,source,orbits,seed,mapping,cut,**options)
    assert restored.stats()['full_hits']==len(requests)
    # A changed lock must still invalidate a saved entry after transport.
    mapping=[-1]*8+[8]
    assert raw(source,source,orbits,0,mapping,(),replay=restored)==raw(source,source,orbits,0,mapping,())
    assert restored.stats()['full_hits']==len(requests)
    before=restored.stats()['entries']
    with pytest.raises(ValueError,match='Truncated'):
        restored.loads(cache.dumps()[:-1])
    assert restored.stats()['entries']==before
    changed=FragmentRepair(source,source,orbits,growth_policy='existence' if growth_policy=='incremental' else 'incremental')
    with pytest.raises(ValueError,match='context differs'):
        changed.loads(cache.dumps())
    tiny=FragmentRepair(source,source,orbits,cache_bytes=1,growth_policy=growth_policy)
    tiny.loads(cache.dumps())
    assert tiny.stats()['resident_bytes']==tiny.stats()['entries']==0


def test_snapshot_is_portable_across_process_element_tables(tmp_path):
    import subprocess, sys
    from mappa import AAMSearchConfig
    from mappa.mappa_core.core.search.session import SearchSession
    from test_reused_aam_api import problem
    search=SearchSession(problem(), AAMSearchConfig(growth_policy='existence'))
    search.search_cut(())
    snapshot=tmp_path/'fragment.cache'
    snapshot.write_bytes(search.repair.dumps())
    script='''
import sys
from pathlib import Path
from mappa import AAMSearchConfig
from mappa.mappa_core.core.search.session import SearchSession
from mappa.mappa_core.core.matching.growth.native import source_graph
from test_cut_replay import graph
from test_reused_aam_api import problem
source_graph(graph(2, [(0,1)], ['Xe','Ne']))
search=SearchSession(problem(), AAMSearchConfig(growth_policy='existence'))
search.repair.loads(Path(sys.argv[1]).read_bytes())
actual,_=search.search_cut(())
fresh,_=SearchSession(problem(), AAMSearchConfig(growth_policy='existence')).search_cut(())
assert actual==fresh
stats=search.repair.stats()
assert stats['full_hits']==stats['calls'] and stats['calls']>0
'''
    subprocess.run([sys.executable,'-c',script,str(snapshot)],check=True)


def test_snapshot_import_retains_local_answers_without_duplicate_entries(growth_policy):
    source=graph(8,[(i,i+1) for i in range(7)])
    orbits=_nauty_orbits(source,wbo_tol=1.)
    remote=FragmentRepair(source,source,orbits,growth_policy=growth_policy)
    local=FragmentRepair(source,source,orbits,growth_policy=growth_policy)
    raw(source,source,orbits,0,[-1]*8,(),replay=remote)
    raw(source,source,orbits,7,[-1]*8,(),replay=local)
    assert local.loads(remote.dumps())  # This worker has an unpublished answer.
    entries=local.stats()['entries']
    assert local.loads(remote.dumps())
    assert local.stats()['entries']==entries
    before=local.stats()['full_hits']
    for seed in (0,7):
        assert raw(source,source,orbits,seed,[-1]*8,(),replay=local)==raw(source,source,orbits,seed,[-1]*8,())
    assert local.stats()['full_hits']==before+2
    assert not remote.loads(local.dumps())  # The union covers all its local answers.


def test_recent_cut_hit_does_not_scan_older_cut_results(growth_policy):
    source=graph(20,[(i,i+1) for i in range(19)])
    orbits=_nauty_orbits(source,wbo_tol=1.)
    cache=FragmentRepair(source,source,orbits,growth_policy=growth_policy)
    for edge in range(18):
        raw(source,source,orbits,0,[-1]*20,((edge,edge+1),),replay=cache)
    before=cache.stats()
    for _ in range(10):
        assert raw(source,source,orbits,0,[-1]*20,((17,18),),replay=cache)==raw(
            source,source,orbits,0,[-1]*20,((17,18),))
    after=cache.stats()
    assert after['full_hits']-before['full_hits']==10
    assert after['dependency_checks']-before['dependency_checks']==10


def test_reused_older_entry_becomes_first_choice_for_next_request(growth_policy):
    source=graph(10,[(i,i+1) for i in (0,1,2,3,5,6,7,8)])
    orbits=_nauty_orbits(source,wbo_tol=1.)
    cache=FragmentRepair(source,source,orbits,growth_policy=growth_policy)
    for cut in (((8,9),),((2,3),),((7,8),)):
        assert raw(source,source,orbits,0,[-1]*10,cut,replay=cache)==raw(source,source,orbits,0,[-1]*10,cut)
    before=cache.stats()
    for _ in range(10):
        assert raw(source,source,orbits,0,[-1]*10,((6,7),),replay=cache)==raw(
            source,source,orbits,0,[-1]*10,((6,7),))
    after=cache.stats()
    assert after['full_hits']-before['full_hits']==10
    assert after['dependency_checks']-before['dependency_checks']==10
