"""Differential tests for dependency-directed conditional fragment repair."""
import numpy as np
import pytest
from mappa.mappa_core.core.matching.policy import BondTolerance

from mappa.mappa_core.core.search.branching.construction import find_islands
from mappa.mappa_core.core.search.branching.reuse import FragmentRepair
from mappa.mappa_core.core.matching.growth import native
from growth_cases import graph, raw as native_raw

pytestmark=pytest.mark.skipif(not native.built(),reason='native engine not built')


def raw(*args, **kwargs):
    result = native_raw(*args, **kwargs)
    result.pop('profile')  # Timing and operation counts are not fragment results.
    return result


def test_reuse_later_fragment_after_independent_locked_mapping_changes():
    source=graph(10,[(0,1),(1,2),(2,3),(4,5),(6,7),(7,8),(8,9)])
    repair=FragmentRepair(source,source)
    for mapping,islands in (([-1]*4+[4,5]+[-1]*4,[(4,1),(5,1)]),
                            ([-1]*4+[5,4]+[-1]*4,[(4,1),(5,2)])):
        assert raw(source,source,0,mapping,(),islands=islands,replay=repair)==raw(
            source,source,0,mapping,(),islands=islands)
    assert repair.stats()['history_hits']==1


def test_rebase_unrelated_previous_boundaries():
    # O cannot extend, while C-N has two distinct placements. A remote prior
    # boundary changes continuation without changing either completed fragment.
    source=graph(6,[(0,1),(0,2),(3,4),(4,5)],['C','O','N','F','Cl','Br'])
    target=graph(7,[(0,1),(2,3),(0,4),(2,5)],['C','N','C','N','F','Cl','Br'])
    repair=FragmentRepair(source,target)
    initial=raw(source,target,0,[-1]*6,(),deferred=((3,4),),replay=repair)
    assert [item['target']['atoms'] for item in initial['fragments']]==[[0,1],[2,3]]
    assert all(item['source']['atoms']==[0,2] and item['accepted_edges']==[(0,1)]
               for item in initial['fragments'])
    assert set(initial['deferred_edges'])=={(0,1),(3,4)}
    assert initial['fragments'][0]['source'] is initial['fragments'][1]['source']
    snapshot=repair.dumps()
    repair=FragmentRepair(source,target)
    repair.loads(snapshot)
    for deferred in (((4,5),),()):
        actual=raw(source,target,0,[-1]*6,(),deferred=deferred,replay=repair)
        assert actual==raw(source,target,0,[-1]*6,(),deferred=deferred)
        assert actual['fragments']==initial['fragments']
        assert actual['fragments'][0]['source'] is actual['fragments'][1]['source']
        assert set(actual['deferred_edges'])=={(0,1),*deferred}
    assert set(initial['deferred_edges'])=={(0,1),(3,4)}
    assert repair.stats()['full_hits']==2
    assert repair.stats()['boundary_rebases']==2


def test_changed_occupancy_cannot_hide_new_seed_alternatives():
    source=graph(6,[(0,1),(1,2),(3,4),(4,5)])
    repair=FragmentRepair(source,source)
    expected=[]
    for mapping in ([-1]*3+[3,4,5],[-1]*6):
        fresh=raw(source,source,0,mapping,())
        assert raw(source,source,0,mapping,(),replay=repair)==fresh
        expected.append(fresh)
    assert expected[0]!=expected[1]
    assert repair.stats()['full_hits']==0


@pytest.mark.parametrize('tol',[.5,1.])
@pytest.mark.parametrize('budget',[1,64*1024*1024])
def test_changed_contexts_and_cuts_match_fresh_results(tol,budget):
    rng=np.random.default_rng(93852)
    for trial in range(5):
        elements=['C']*6+['H']*3
        edges=[(i,i+1,float(rng.choice([.8,1.2,1.8]))) for i in range(5)]
        edges += [(0,6,1.),(2,7,1.),(4,8,1.)]
        if trial%2:edges.append((0,4,1.))
        source=graph(9,edges,elements)
        target=graph(9,edges[:-1]+[(1,5,.8)],elements)
        repair=FragmentRepair(source,target,cache_bytes=budget)
        for _ in range(40):
            mapping=[-1]*9
            locked=list(map(int,rng.choice(6,3,replace=False)))
            for a,b in zip(locked,rng.permutation([1,3,5])):mapping[a]=int(b)
            islands=[(r,1+int(rng.integers(2))) for r in locked]
            cut=() if rng.integers(3)==0 else (edges[int(rng.integers(len(edges)))][:2],)
            deferred=() if rng.integers(2)==0 else (edges[int(rng.integers(len(edges)))][:2],)
            seed=int(rng.integers(9));cap=int(rng.choice([1,3,100]))
            kw=dict(cap=cap,tol=tol,islands=islands,deferred=deferred)
            assert raw(source,target,seed,mapping,cut,replay=repair,**kw)==raw(
                source,target,seed,mapping,cut,**kw),(trial,mapping,islands,cut)


def test_complete_compressed_search_graphs_and_generators():
    source=graph(9,[(0,1),(1,2),(2,3),(3,4),(4,5),(0,6),(2,7),(4,8)],['C']*6+['H']*3)
    target=graph(9,[(0,1),(1,2),(2,3),(3,4),(4,5),(1,6),(3,7),(5,8)],['C']*6+['H']*3)
    repair=FragmentRepair(source,target)
    for order in (list(range(9)),list(reversed(range(9))),[4,2,6,0,1,3,5,7,8]):
        for cut in ((),((2,3),),((0,6),),((0,1),)):
            view=repair.for_cut(cut)
            options=dict(cuts=cut,max_branches=100,iso_tol=1.)
            expected=find_islands(view.source,target,order,**options)
            actual=find_islands(view.source,target,order,fragment_cache=view,**options)
            assert actual==expected


def test_eviction_bounds_retention_without_changing_results():
    source=graph(8,[(i,(i+1)%8) for i in range(8)])
    probe=FragmentRepair(source,source)
    raw(source,source,0,[-1]*8,(),replay=probe)
    budget=2*probe.stats()['resident_bytes']
    repair=FragmentRepair(source,source,cache_bytes=budget)
    for _ in range(3):
        for seed in range(8):
            assert raw(source,source,seed,[-1]*8,(),replay=repair)==raw(
                source,source,seed,[-1]*8,())
            stats=repair.stats()
            assert 0 < stats['resident_bytes'] <= stats['peak_resident_bytes'] <= budget
    assert repair.stats()['evictions'] > 0
    before=repair.stats()
    raw(source,source,7,[-1]*8,(),replay=repair)
    after=repair.stats()
    assert after['full_hits']==before['full_hits']+1
    assert after['resident_bytes']==before['resident_bytes']


def test_shared_self_preparation_survives_transport_and_releases_on_eviction():
    from mappa import match_fragment, FragmentMatchConfig, FragmentMatchContext
    source=graph(21,[(i,(i+1)%20) for i in range(20)],['C']*20+['O'])
    cache=FragmentRepair(source,source)
    raw(source,source,0,[-1]*21,(),replay=cache)
    first=cache.stats()
    for seed in range(1,6):
        assert raw(source,source,seed,[-1]*21,(),replay=cache)==raw(source,source,seed,[-1]*21,())
    shared=cache.stats()
    assert shared['entries']==6 and shared['shared_self_records']==2
    assert shared['shared_self_bytes']==first['shared_self_bytes']
    assert shared['resident_bytes'] < 6*first['resident_bytes']
    restored=FragmentRepair(source,source)
    restored.loads(cache.dumps())
    assert restored.stats()['shared_self_records']==2
    for seed in range(6):
        assert raw(source,source,seed,[-1]*21,(),replay=restored)==raw(source,source,seed,[-1]*21,())
    assert restored.stats()['full_hits']==6

    # Small oxygen-only entries eventually evict every ring entry. The ring's
    # shared matrices and registry records must leave with their last owner.
    bounded=FragmentRepair(source,source,cache_bytes=first['resident_bytes'])
    raw(source,source,0,[-1]*21,(),replay=bounded)
    view=bounded.for_cut(())
    context=FragmentMatchContext(fragment_cache=view)
    expected=match_fragment(source,source,seed=20).matches[0].to_record()
    for random_seed in range(30):
        result=match_fragment(view.source,source,seed=20,context=context,
                             config=FragmentMatchConfig(random_seed=random_seed))
        assert result.matches[0].to_record()==expected
        assert bounded.stats()['resident_bytes']<=first['resident_bytes']
    released=bounded.stats()
    assert released['evictions']>0 and released['shared_self_records']==2
    assert released['shared_self_bytes']<first['shared_self_bytes']


def test_oversized_results_are_returned_without_retention():
    source=graph(8,[(i,i+1) for i in range(7)])
    repair=FragmentRepair(source,source,cache_bytes=1)
    for _ in range(3):
        assert raw(source,source,0,[-1]*8,(),replay=repair)==raw(
            source,source,0,[-1]*8,())
    stats=repair.stats()
    assert stats['entries']==stats['resident_bytes']==stats['peak_resident_bytes']==0
    assert stats['oversized_results']==3


def test_public_fragment_cache_reuses_complete_results():
    from mappa import match_fragment, FragmentMatchConfig, FragmentMatchContext
    from test_fragment_growth import record
    source=graph(6,[(i,(i+1)%6) for i in range(6)])
    cache=FragmentRepair(source,source)
    view=cache.for_cut(())
    context=FragmentMatchContext(fragment_cache=view)
    config=FragmentMatchConfig(iso_tolerance=1.)
    expected=record(match_fragment(view.source,source,seed=0,config=config))
    for _ in range(3):
        assert record(match_fragment(view.source,source,seed=0,context=context,config=config))==expected
    assert cache.stats()['full_hits']==2


def test_snapshot_preserves_results_dependencies_and_caps():
    source=graph(9,[(0,1),(1,2),(2,3),(3,4),(4,5),(0,6),(2,7),(4,8)],['C']*6+['H']*3)
    cache=FragmentRepair(source,source)
    requests=[(0,[-1]*9,(),dict(cap=cap)) for cap in (1,3,100)]
    requests += [(0,[-1]*6+[6,7,8],((2,3),),dict(cap=3,islands=[(6,1),(7,1),(8,2)])),
                 (2,[-1]*9,((0,1),),dict(cap=100,deferred=((2,3),)))]
    for seed,mapping,cut,options in requests:
        raw(source,source,seed,mapping,cut,replay=cache,**options)
    restored=FragmentRepair(source,source)
    restored.loads(cache.dumps())
    for seed,mapping,cut,options in requests:
        assert raw(source,source,seed,mapping,cut,replay=restored,**options)==raw(
            source,source,seed,mapping,cut,**options)
    assert restored.stats()['full_hits']==len(requests)
    # A changed lock must still invalidate a saved entry after transport.
    mapping=[-1]*8+[8]
    assert raw(source,source,0,mapping,(),replay=restored)==raw(source,source,0,mapping,())
    assert restored.stats()['full_hits']==len(requests)
    before=restored.stats()['entries']
    with pytest.raises(ValueError,match='Truncated'):
        restored.loads(cache.dumps()[:-1])
    assert restored.stats()['entries']==before
    changed=FragmentRepair(source,source,bond_tolerance=BondTolerance(.25))
    with pytest.raises(ValueError,match='context differs'):
        changed.loads(cache.dumps())
    tiny=FragmentRepair(source,source,cache_bytes=1)
    tiny.loads(cache.dumps())
    assert tiny.stats()['resident_bytes']==tiny.stats()['entries']==0


def test_snapshot_is_portable_across_process_element_tables(tmp_path):
    import os, subprocess, sys
    from pathlib import Path
    from mappa import AAMSearchConfig
    from mappa.mappa_core.core.search.session import SearchSession
    from test_reused_aam_api import problem
    search=SearchSession(problem(), AAMSearchConfig())
    search.search_cut(())
    snapshot=tmp_path/'fragment.cache'
    snapshot.write_bytes(search.repair.dumps())
    script='''
import sys
from pathlib import Path
from mappa import AAMSearchConfig
from mappa.mappa_core.core.search.session import SearchSession
from mappa.mappa_core.core.matching.growth.native import source_graph
from growth_cases import graph
from test_reused_aam_api import problem
source_graph(graph(2, [(0,1)], ['Xe','Ne']))
search=SearchSession(problem(), AAMSearchConfig())
search.repair.loads(Path(sys.argv[1]).read_bytes())
actual,_=search.search_cut(())
fresh,_=SearchSession(problem(), AAMSearchConfig()).search_cut(())
assert actual==fresh
stats=search.repair.stats()
assert stats['full_hits']==stats['calls'] and stats['calls']>0
'''
    subprocess.run([sys.executable,'-c',script,str(snapshot)],check=True,
        env={**os.environ, 'PYTHONPATH': str(Path(__file__).parent) + os.pathsep + os.environ.get('PYTHONPATH', '')})


def test_incompatible_snapshot_layout_does_not_replace_resident_results():
    import sys
    source=graph(4,[(0,1),(1,2),(2,3)])
    repair=FragmentRepair(source,source)
    expected=raw(source,source,0,[-1]*4,(),replay=repair)
    snapshot=repair.dumps()
    header_end=8+int.from_bytes(snapshot[:8],sys.byteorder)
    # Prior snapshots began directly with context identity. Keeping the same
    # identity must not let an incompatible entry layout reach the parser.
    previous_layout=snapshot[header_end:]
    entries=repair.stats()['entries']
    with pytest.raises(ValueError,match='Unsupported fragment cache schema'):
        repair.loads(previous_layout)
    assert repair.stats()['entries']==entries
    assert raw(source,source,0,[-1]*4,(),replay=repair)==expected
    assert repair.stats()['full_hits']==1


def test_snapshot_import_retains_local_answers_without_duplicate_entries():
    source=graph(8,[(i,i+1) for i in range(7)])
    remote=FragmentRepair(source,source)
    local=FragmentRepair(source,source)
    raw(source,source,0,[-1]*8,(),replay=remote)
    raw(source,source,7,[-1]*8,(),replay=local)
    assert local.loads(remote.dumps())  # This worker has an unpublished answer.
    entries=local.stats()['entries']
    assert local.loads(remote.dumps())
    assert local.stats()['entries']==entries
    before=local.stats()['full_hits']
    for seed in (0,7):
        assert raw(source,source,seed,[-1]*8,(),replay=local)==raw(source,source,seed,[-1]*8,())
    assert local.stats()['full_hits']==before+2
    assert not remote.loads(local.dumps())  # The union covers all its local answers.


def test_recent_cut_hit_does_not_scan_older_cut_results():
    source=graph(20,[(i,i+1) for i in range(19)])
    cache=FragmentRepair(source,source)
    for edge in range(18):
        raw(source,source,0,[-1]*20,((edge,edge+1),),replay=cache)
    before=cache.stats()
    for _ in range(10):
        assert raw(source,source,0,[-1]*20,((17,18),),replay=cache)==raw(
            source,source,0,[-1]*20,((17,18),))
    after=cache.stats()
    assert after['full_hits']-before['full_hits']==10
    assert after['dependency_checks']-before['dependency_checks']==10


def test_reused_older_entry_becomes_first_choice_for_next_request():
    source=graph(10,[(i,i+1) for i in (0,1,2,3,5,6,7,8)])
    cache=FragmentRepair(source,source)
    for cut in (((8,9),),((2,3),),((7,8),)):
        assert raw(source,source,0,[-1]*10,cut,replay=cache)==raw(source,source,0,[-1]*10,cut)
    before=cache.stats()
    for _ in range(10):
        assert raw(source,source,0,[-1]*10,((6,7),),replay=cache)==raw(
            source,source,0,[-1]*10,((6,7),))
    after=cache.stats()
    assert after['full_hits']-before['full_hits']==10
    assert after['dependency_checks']-before['dependency_checks']==10


def test_cache_reuses_binary_symmetry_without_restricting_discovery():
    from mappa import match_fragment, FragmentMatchConfig, FragmentMatchContext
    from growth_cases import discovery_mappings, represented
    source = graph(3, [(0, 1, 1.), (0, 2, 2.)], ['O', 'C', 'C'])
    cache = FragmentRepair(source, source)
    view = cache.for_cut(())
    context = FragmentMatchContext(fragment_cache=view)
    strict = FragmentMatchConfig(iso_tolerance=.1)
    relaxed = FragmentMatchConfig(iso_tolerance=.1)
    initial = match_fragment(view.source, source, seed=0, context=context, config=strict)
    expanded = match_fragment(view.source, source, seed=0, context=context, config=relaxed)
    assert len(initial.matches) == len(expanded.matches) == 1
    assert discovery_mappings(initial.matches[0]) == discovery_mappings(expanded.matches[0]) == {
        ((0, 0), (1, 1), (2, 2))}
    assert represented(initial) == represented(expanded) == {
        ((0, 0), (1, 1), (2, 2)), ((0, 0), (1, 2), (2, 1))}
    assert cache.stats()['full_hits'] == 1
    restored = FragmentRepair(source, source)
    restored.loads(cache.dumps())
    restored_view = restored.for_cut(())
    for config, expected in ((strict, initial), (relaxed, expanded)):
        actual = match_fragment(restored_view.source, source, seed=0,
            context=FragmentMatchContext(fragment_cache=restored_view), config=config)
        assert [fragment.key for fragment in actual.matches] == [fragment.key for fragment in expected.matches]
    assert restored.stats()['full_hits'] == 2


def test_external_placement_context_survives_cache_reuse_and_transport():
    from mappa import FragmentMatchConfig, FragmentMatchContext, match_fragment
    from growth_cases import represented
    from test_fragment_growth import record
    source = graph(1, [], ['H'])
    target = graph(4, [(0, p) for p in (1, 2, 3)], ['C', 'H', 'H', 'H'])
    config = FragmentMatchConfig(branch_limit=3)
    expected = {}
    cache = FragmentRepair(source, target)
    view = cache.for_cut(())
    for fixed, count in (((), 1), (tuple(target), 3)):
        expected[fixed] = match_fragment(source, target, seed=0, config=config,
            context=FragmentMatchContext(placement_fixed_targets=fixed))
        actual = match_fragment(view.source, target, seed=0, config=config,
            context=FragmentMatchContext(fragment_cache=view, placement_fixed_targets=fixed))
        assert len(actual.matches) == count
        assert record(actual) == record(expected[fixed])
        assert represented(actual) == {((0, p),) for p in (1, 2, 3)}
    restored = FragmentRepair(source, target)
    restored.loads(cache.dumps())
    restored_view = restored.for_cut(())
    for fixed in ((), tuple(target)):
        actual = match_fragment(restored_view.source, target, seed=0, config=config,
            context=FragmentMatchContext(fragment_cache=restored_view, placement_fixed_targets=fixed))
        assert record(actual) == record(expected[fixed])
    assert restored.stats()['full_hits'] == 2
