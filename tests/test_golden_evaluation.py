import importlib.util
from itertools import product
from pathlib import Path
import pynauty
import pytest

spec=importlib.util.spec_from_file_location('golden_evaluation',Path(__file__).parents[1]/'bench/golden_evaluation.py')
E=importlib.util.module_from_spec(spec);spec.loader.exec_module(E)


def feature(n, distinct=False):
    return dict(heavy=list(range(n)),colors=[('C',i if distinct else 0) for i in range(n)],bonds=[])


def test_joint_chemical_equivalence_is_not_independent_atom_membership():
    f=feature(4);f['bonds']=[(0,1,(1.,'')),(1,2,(1.,'')),(2,3,(1.,''))]
    original=pynauty.certificate(E.colored_graph([f,f],dict(enumerate(range(4)))))
    reversed_=pynauty.certificate(E.colored_graph([f,f],dict(enumerate([3,2,1,0]))))
    wrong=pynauty.certificate(E.colored_graph([f,f],dict(enumerate([0,2,1,3]))))
    assert original==reversed_ and original!=wrong


@pytest.mark.parametrize('finite_domain',[False,True])
def test_symbolic_conditioned_group_product_matches_tiny_exhaustive_oracle(finite_domain):
    f=feature(3,True)
    groups=[((1,0,2),),((0,2,1),)]
    possible=set()
    for first,second in product([(0,1,2),(1,0,2)],[(0,1,2),(0,2,1)]):
        possible.add(tuple(first[second[i]] for i in range(3)))
    from itertools import permutations
    for values in permutations(range(3)):
        status,_=E.symbolic_path_query(dict(enumerate(range(3))),groups,[f,f],dict(enumerate(values)),5000,
                                     finite_domain=finite_domain)
        assert (status=='recovered') == (values in possible)


@pytest.mark.parametrize('finite_domain',[False,True])
def test_symbolic_query_preserves_unmatched_domain_and_free_hydrogen_projection(finite_domain):
    f=feature(3,True)
    status,_=E.symbolic_path_query({0:0,1:1},[((1,0,2),)],[f,f],{0:1,1:0},5000,finite_domain=finite_domain)
    assert status=='recovered'
    status,_=E.symbolic_path_query({0:0,1:1},[((1,0,2),)],[f,f],{0:1,2:0},5000,finite_domain=finite_domain)
    assert status=='not_recovered'


def test_inputs_are_unmapped_explicit_H_and_reference_indices_are_separate():
    problem,features,reference=E.prepare('[CH3:7][Br:2].[OH2:5]>>[CH3:7][OH:5].[BrH:2]')
    assert problem.balanced
    assert problem.source_atom_count==8
    assert len(E.project(reference,features))==3
    assert all(':' not in e.metadata['unmapped_smiles'] for e in (problem.reactant,problem.product))


@pytest.mark.parametrize('extra_pair',[False,True])
def test_one_sided_reference_requires_exact_relation_not_subset(extra_pair):
    from types import SimpleNamespace as S
    import numpy as np
    from mappa import AAMProblem
    from mappa.mappa_core.domain import MolecularEndpoint
    endpoint=MolecularEndpoint(('C','C'),np.zeros((2,3)),np.zeros((2,2)))
    mapping={0:0,1:1} if extra_pair else {0:0}
    graph=S(terminals=[0],states=[S(mapping=tuple(mapping.items()))],capped=False,
            paths=lambda:iter(()))
    f=feature(2,True)
    result=E.evaluate(S(graph=graph,problem=AAMProblem(endpoint,endpoint)),[f,f],{0:0})
    assert not result['reference_annotation_complete']
    assert result['top1_correct'] is (not extra_pair)
    assert (result['reference_recovery']=='recovered') is (not extra_pair)


def test_reference_recovery_queries_real_compressed_aam_not_only_representative():
    import numpy as np
    from mappa import AAMProblem,AAMSearchConfig,search_aam
    from mappa.mappa_core.domain import MolecularEndpoint
    from dataclasses import replace
    from itertools import permutations
    def endpoint(n):return MolecularEndpoint(('O',)*n,np.zeros((n,3)),np.zeros((n,n)))
    result=search_aam(AAMProblem(endpoint(2),endpoint(4)),AAMSearchConfig(seed_count=1))
    terminal=result.graph.terminals[0]
    result=replace(result,graph=replace(result.graph,stops=tuple(
        stop for stop in result.graph.stops if stop.state==terminal)))
    mapping=dict(result.graph.states[terminal].mapping)
    reference=next(dict(enumerate(p)) for p in permutations(range(4),2) if dict(enumerate(p))!=mapping)
    features=[feature(2,True),feature(4,True)]
    evaluated=E.evaluate(result,features,reference,seconds=10)
    assert not evaluated['representative_recovery']
    assert evaluated['reference_recovery']=='recovered'
    assert evaluated['symbolic_queries']>0


def test_H_witness_differences_retain_full_feasibility_queries():
    from types import SimpleNamespace as S
    import numpy as np
    from mappa import AAMProblem
    from mappa.mappa_core.domain import MolecularEndpoint
    endpoint=MolecularEndpoint(('C','C','H','H'),np.zeros((4,3)),np.zeros((4,4)))
    mappings=[{0:0,1:1,2:2,3:3},{0:0,1:1,2:3,3:2}]
    graph=S(terminals=[0,1],states=[S(mapping=tuple(m.items())) for m in mappings],capped=False,
        transitions=[S(source=0),S(source=0)],
        paths=lambda:iter([S(mapping=m,terminal=i,transitions=[i],context=S(cuts=()),fragments=()) for i,m in enumerate(mappings)]),
        fragment_placement=lambda edge:S(r_atoms=(),representative_assignments=(),exact_fixed=(),deferred_edges=(),
            symmetry_domains=(),target_generators=[S(images=(1,0,2,3)),S(images=(0,1,3,2))]))
    calls=[]
    def query(*args,**kwargs):calls.append(args);return 'unknown',None
    f=feature(2,True)
    result=E.evaluate(S(graph=graph,problem=AAMProblem(endpoint,endpoint)),[f,f],{0:1,1:0},query=query)
    assert result['reference_recovery']=='unknown'
    assert len(calls)==result['symbolic_queries']==2


def test_explicit_H_donor_competition_depends_on_seed_order():
    """Completing a seed order does not enumerate every reactant subset."""
    import numpy as np
    from mappa.mappa_core.core.matching.graphs import build_graph
    from mappa.mappa_core.core.search.branching.construction import find_islands
    source=np.zeros((9,9));target=np.zeros((6,6))
    for a,b in [(0,1),(0,3),(0,4),(0,5),(1,8),(2,6),(2,7)]:source[a,b]=source[b,a]=1
    for a,b in [(0,1),(0,2),(0,3),(0,4),(1,5)]:target[a,b]=target[b,a]=1
    r=build_graph(('C','O','O')+('H',)*6,source)
    p=build_graph(('C','O')+('H',)*4,target)
    def recovered(order):
        graph=find_islands(r,p,order,iso_tol=1.,max_branches=100)
        return any(dict(graph.states[t].mapping).get(0)==0 and
                   dict(graph.states[t].mapping).get(2)==1 for t in graph.terminals)
    assert not recovered(list(range(9)))
    assert recovered([2,0,1,3,4,5,6,7,8])


def test_compressed_membership_matches_independent_cyclic_oracle():
    from itertools import permutations
    from mappa import AAMProblem
    from mappa.mappa_core.post_processing.final_branches import FinalFamily
    from test_family_query import endpoint
    problem=AAMProblem(endpoint('CCCC',[]),endpoint('CCCC',[]))
    family=FinalFamily(tuple(enumerate(range(4))),(),(('group',((1,2,3,0),)),),(.2,1.))
    f=[feature(4,True),feature(4,True)]
    allowed={tuple((i+shift)%4 for i in range(4)) for shift in range(4)}
    for values in permutations(range(4)):
        result=E.reference_membership(problem,[(0,family.mapping)],
            lambda:iter([(0,family.as_path(problem))]),f,dict(enumerate(values)))
        assert (result['status']=='recovered') == (values in allowed)
        if values in allowed:
            assert dict(result['witness']['mapping'])==dict(enumerate(values))
            if values!=tuple(range(4)):
                assert result['witness']['method']=='compressed_choices'
                assert 'encoding' in result['queries'][0]['timing']


def test_membership_representative_stops_before_path_construction():
    from mappa import AAMProblem
    from test_family_query import endpoint
    problem=AAMProblem(endpoint('CC',[]),endpoint('CCC',[]))
    def forbidden():raise AssertionError('No query needed for a verified witness')
    result=E.reference_membership(problem,[(0,{0:1,1:2})],forbidden,
        [feature(2,True),feature(3,True)],{0:1,1:2})
    assert result['status']=='recovered' and result['queries']==[]


def test_membership_preserves_exact_partial_relation():
    from mappa import AAMProblem
    from mappa.mappa_core.post_processing.final_branches import FinalFamily
    from test_family_query import endpoint
    problem=AAMProblem(endpoint('CC',[]),endpoint('CC',[]))
    f=[feature(2,True),feature(2,True)]
    family=FinalFamily(((0,0),(1,1)),(),(),(.2,1.))
    result=E.reference_membership(problem,[(0,family.mapping)],
        lambda:iter([(0,family.as_path(problem))]),f,{0:0})
    assert result['status']=='not_recovered'


def test_event_filter_does_not_reject_other_choices_in_compressed_branch():
    from mappa import AAMProblem
    from mappa.mappa_core.post_processing.final_branches import FinalFamily
    from test_family_query import endpoint
    problem=AAMProblem(endpoint('CCC',[(0,1)]),endpoint('CCC',[(0,1)]))
    f=feature(3,True);f['bonds']=[(0,1,(1.,''))]
    family=FinalFamily(((0,1),(1,2),(2,0)),(),(('group',((1,2,0),)),),(.2,1.))
    result=E.reference_membership(problem,[(0,family.mapping)],
        lambda:iter([(0,family.as_path(problem))]),[f,f],{0:0,1:1,2:2})
    assert result['status']=='recovered'
    assert result['event_filter']['rejected_representatives']==1
    assert result['event_filter']['reference_count']==0
    assert dict(result['witness']['mapping'])=={0:0,1:1,2:2}


def test_event_filter_partial_boundary_and_metal_threshold():
    f=dict(heavy=[0,1,2],colors=[('Fe',0),('C',0),('C',0)],
           bonds=[(0,1,(.4,'')),(1,2,(1.,''))])
    p=dict(f,bonds=[])
    index=E.reference_event_index([f,p])
    # Metal boundary counts at .3; wholly unmatched C-C material is ignored.
    assert tuple(index.counts([index.vector({0:0})])[0])==(1,0)


def test_event_filter_ignores_unannotated_hydrogen_correspondence():
    problem,features,reference=E.prepare('[CH3:1][OH:2]>>[CH2:1]=[O:2]')
    index=E.reference_event_index(features)
    assert index.n==index.m==2
    assert index.counts([index.vector(E.project(reference,features))])[0].sum()==1
