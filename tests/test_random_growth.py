"""Random selection has no bond-order, hydrogen, or boundary-degree preference."""
from collections import Counter
import json

import numpy as np
import pytest

from mappa.mappa_core.core.matching.frontier import GrowthRandom, _pop_random_atom
from mappa.mappa_core.core.matching.growth import grow_island
from mappa.mappa_core.core.matching import native
from mappa.frag import build_graph
from mappa.mappa_core.core.matching.orbits import _nauty_orbits
from mappa import AAMProblem, AAMSearchConfig, MolecularEndpoint, search_aam
from mappa.mappa_core.core.search.planning import checkpoint_manifest


@pytest.mark.smoke
def test_selection_ignores_weights_and_samples_atoms_not_edges():
    # Atom 3 has three incident edges, atom 4 only one. Equal atom probability.
    pool=[(-3.,0,3),(-2.,1,3),(-1.,2,3),(-.5,0,4)]
    counts=Counter()
    for seed in range(2000):
        choice=_pop_random_atom(list(pool),set(),{0,1,2},GrowthRandom(seed))
        changed=[(-(.6 if i%2 else 8.),u,n) for i,(_,u,n) in enumerate(pool)]
        other=_pop_random_atom(changed,set(),{0,1,2},GrowthRandom(seed))
        assert choice[1:]==other[1:]
        counts[choice[2]]+=1
    assert 850<counts[3]<1150 and 850<counts[4]<1150


def molecule():
    # Unique element correspondence keeps matching feasible regardless of order.
    elements=('N','C','O','H','S')
    w=np.zeros((5,5))
    for a,b,v in [(0,1,3.),(0,2,2.),(0,3,1.),(1,4,1.)]:w[a,b]=w[b,a]=v
    return build_graph(elements,w,bond_cut=.2), MolecularEndpoint(elements,np.zeros((5,3)),w)


def trace(seed):
    source,_=molecule();events=[]
    grow_island(source,source,0,{},events=events,iso_tol=.1,random_seed=seed)
    return [e['edge']['ext_atom'] for e in events if e['type']=='pop']


def test_actual_growth_can_choose_hydrogen_before_triple_bond():
    orders=[trace(seed) for seed in range(30)]
    assert {order[0] for order in orders}=={1,2,3}
    assert len({tuple(order) for order in orders})>1
    assert orders[7]==trace(7)


@pytest.mark.skipif(not native.built(),reason='native engine required')
def test_nondefault_seeds_agree_in_python_native_and_repair(monkeypatch):
    from mappa.mappa_core.core.search.reuse import CutReplay, FragmentRepair
    source,_=molecule();orbits=_nauty_orbits(source,wbo_tol=.1)
    r=native.source_graph(source);p=native.target_graph(source,orbits)
    sessions=[CutReplay(source,source,orbits),FragmentRepair(source,source,orbits)]
    for seed in [1,2,7,42,2**64-1,1]:
        expected=native._engine.grow_island(r.graph,p.graph,0,[-1]*5,.2,.1,1,100,None,(),False,seed)
        for session in sessions:
            actual=session.engine.grow((),0,[-1]*5,.2,.1,1,100,None,(),False,seed)
            actual.pop('reused_extensions')
            assert actual==expected
        def result(mode):
            monkeypatch.setenv('MAPPA_NATIVE',mode)
            return [(dict(x),x.fragment,x.deferred_edges,x.symmetry) for x in grow_island(
                source,source,0,{},p_orbits=orbits,iso_tol=.1,random_seed=seed)]
        assert result('0')==result('1')


def test_new_search_records_seed_and_rejects_old_policy(tmp_path):
    _,ep=molecule();problem=AAMProblem(ep,ep)
    config=AAMSearchConfig(random_seed=137,cut_floor=10,iso_tolerance=.1)
    result=search_aam(problem,config)
    assert all(c.random_seed==137 and c.growth_policy=='uniform_boundary_atom_splitmix64_v1'
               for c in result.graph.contexts)
    assert result.graph==search_aam(problem,config).graph
    old=checkpoint_manifest(problem,config)
    old.pop('growth_policy');old['schema']='rxn_core.aam_checkpoints/v3'
    (tmp_path/'manifest.json').write_text(json.dumps(old))
    with pytest.raises(ValueError,match='configuration differs'):
        search_aam(problem,config,intermediate_dir=tmp_path,resume=True)
    old=checkpoint_manifest(problem,config)
    old.pop('matching_policy');old['schema']='rxn_core.aam_checkpoints/v4'
    (tmp_path/'manifest.json').write_text(json.dumps(old))
    with pytest.raises(ValueError,match='configuration differs'):
        search_aam(problem,config,intermediate_dir=tmp_path,resume=True)


def test_old_pickle_policy_is_labeled_and_cannot_replay_as_new():
    import pickle
    from mappa.mappa_core.core.search.graph import SearchContext
    from mappa.search_trajectory import build_trajectory
    context=SearchContext((0,),(0,),(0,))
    object.__delattr__(context,'growth_policy')
    object.__delattr__(context,'random_seed')
    restored=pickle.loads(pickle.dumps(context))
    assert restored.growth_policy=='legacy_wbo_priority'
    _,ep=molecule();result=search_aam(AAMProblem(ep,ep),AAMSearchConfig(cut_floor=10))
    from dataclasses import replace
    object.__setattr__(result.graph,'contexts',tuple(replace(c,growth_policy='legacy_wbo_priority')
                       for c in result.graph.contexts))
    with pytest.raises(ValueError,match='Growth policy differs'):
        build_trajectory([dict(aam=result,context=0)])


@pytest.mark.skipif(not native.built(),reason='native engine required')
def test_randomized_partial_matching_parity_across_streams(monkeypatch):
    from mappa.mappa_core.core.matching.records import IslandBranchLimitExceeded
    rng=np.random.default_rng(2071)
    for trial in range(8):
        edges=[(i,i+1) for i in range(7)]+[(0,3),(2,6)]
        matrices=[]
        for side in range(2):
            w=np.zeros((8,8))
            for a,b in edges:w[a,b]=w[b,a]=rng.choice([.6,1.,1.8,3.])
            matrices.append(w)
        source,target=[build_graph(['C']*5+['H']*3,w,bond_cut=.2) for w in matrices]
        orbits=_nauty_orbits(target,wbo_tol=.5)
        for random_seed in (0,1,2,7,137):
            def run(mode):
                monkeypatch.setenv('MAPPA_NATIVE',mode)
                try:
                    return [(dict(x),x.fragment,x.deferred_edges,x.symmetry) for x in grow_island(
                        source,target,trial,{},p_orbits=orbits,iso_tol=.5,max_branches=100,
                        random_seed=random_seed)]
                except IslandBranchLimitExceeded as exc:
                    return ('capped',exc.count,exc.limit)
            assert run('0')==run('1'),(trial,random_seed)


def test_selected_native_backend_failure_does_not_fall_through(monkeypatch):
    source,_=molecule()
    monkeypatch.setattr(native,'applicable',lambda *a,**k:True)
    monkeypatch.setattr(native,'grow_island',lambda *a,**k:None)
    with pytest.raises(ValueError,match='Selected native growth cannot represent'):
        grow_island(source,source,0,{})
