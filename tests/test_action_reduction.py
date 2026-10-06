from dataclasses import replace
from mappa.action_reduction import absorb_subgroups
from mappa.final_branches import FinalFamily
from mappa.family_query import query_path
from mappa.event_patterns import SignedEventIndex,extract_path_events
from test_event_patterns import problem

S3=('group',((1,0,2),(1,2,0)))
C3=('group',((1,2,0),))
T01=('group',((1,0,2),))
T12=('group',((0,2,1),))

def test_nested_groups_absorb_in_either_order():
    assert absorb_subgroups((C3,T01,S3),3)==(S3,)
    assert absorb_subgroups((S3,C3,T01),3)==(S3,)
    assert absorb_subgroups((('pool',(0,1)),S3,('pool',(1,2))),3)==(S3,)

def test_equal_orbits_do_not_prove_containment():
    assert absorb_subgroups((C3,T01),3)==(C3,T01)
    assert absorb_subgroups((T01,C3),3)==(T01,C3)

def test_noncommuting_subgroups_are_not_reordered_or_closed():
    assert absorb_subgroups((T01,T12),3)==(T01,T12)
    assert absorb_subgroups((T12,T01),3)==(T12,T01)

def test_disjoint_groups_remain_a_product():
    a=('pool',(0,1));b=('pool',(2,3))
    assert absorb_subgroups((a,b),4)==(a,b)

def test_free_pool_support_does_not_expand_locked_atoms():
    assert absorb_subgroups((('pool',(1,2)),),3)==(('pool',(1,2)),)
    assert absorb_subgroups((('pool',(1,2)),T01),3)==(('pool',(1,2)),T01)

def test_constraints_and_event_alternatives_survive_reduction(monkeypatch):
    p=problem('OOO',[(0,1,1)],[(0,1,1),(0,2,.4)])
    f=FinalFamily(((0,0),(1,1),(2,2)),((0,1),),(C3,T01,S3),(.2,1.))
    simplified=replace(f,actions=(S3,))
    with monkeypatch.context() as patch:
        patch.setattr("mappa.action_reduction.absorb_subgroups",lambda actions,degree:actions)
        original=f.as_path(p)
    for mapping,expected in [({0:0,1:1,2:2},'recovered'),({0:1,1:0,2:2},'recovered'),({0:1,1:2,2:0},'not_recovered')]:
        assert query_path(original,p,mapping,source_atoms=(0,1,2))[0]==expected
        assert query_path(f.as_path(p),p,mapping,source_atoms=(0,1,2))[0]==expected
    a=extract_path_events(original,p,SignedEventIndex(p),seconds=5,max_patterns=None)
    b=extract_path_events(simplified.as_path(p),p,SignedEventIndex(p),seconds=5,max_patterns=None)
    assert a['complete'] and b['complete']
    assert {v['id'] for v in a['patterns']}=={v['id'] for v in b['patterns']}
    assert f.required_edges==simplified.required_edges
