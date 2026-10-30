from dataclasses import replace
import pytest
import sys
from pathlib import Path
from test_event_patterns import problem,saved_path
from mappa.mappa_core.core.search.graph import AAMSearchGraph, SearchContext, SearchState, FragmentTransition, SearchStop
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue, canonical_actions, final_fragment_pairs
from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex
from reference_event_decoder import extract_path_events
from mappa.mappa_workflow.mapping_analytics.family import query_path


def independent_path(order, mapping=None, labels=(7,19)):
    mapping=mapping or {i:i for i in range(4)}
    parts=((0,1),(2,3));states=[SearchState(0,0,(),(),())];transitions=[];assigned={};owners={}
    for j,part_id in enumerate(order):
        atoms=parts[part_id];assigned.update({r:mapping[r] for r in atoms});owners.update({r:labels[part_id] for r in atoms})
        states.append(SearchState(j+1,0,tuple(sorted(assigned.items())),tuple(sorted(owners.items())),()))
        g=list(range(4));p,q=[mapping[r] for r in atoms];g[p],g[q]=q,p
        match=dict(fragment=atoms,deferred_edges=(),symmetry=dict(witness={r:mapping[r] for r in atoms},blocks=[],exact_fixed=[],automorph_generators=[g]))
        transitions.append(FragmentTransition(j,j,j+1,atoms[0],(j,0),match,()))
    graph=AAMSearchGraph((SearchContext(tuple(range(4)),tuple(range(4)),tuple(order),iso_tolerance=1),),(0,),tuple(states),tuple(transitions),(SearchStop(2,'objective_met'),))
    return next(graph.paths())


def test_order_and_owner_labels_do_not_create_branches_or_families():
    value=problem('HHHH',[],[]);a=independent_path((0,1));b=independent_path((1,0),labels=(123,2))
    index=FinalBranchCatalogue(value).add_graph(a.graph,'a').add_graph(b.graph,'b')
    assert index.path_count==2 and len(index.branches)==1 and len(index.families)==1
    assert len(index.families[0].provenance)==2
    assert final_fragment_pairs(a.graph.states[a.terminal])==final_fragment_pairs(b.graph.states[b.terminal])


def test_different_product_fragment_pairs_remain_separate():
    value=problem('HHHH',[],[]);a=independent_path((0,1));b=independent_path((0,1),{0:2,1:3,2:0,3:1})
    index=FinalBranchCatalogue(value).add_graph(a.graph).add_graph(b.graph)
    assert len(index.branches)==2


def test_noncommuting_choices_remain_alternatives_inside_one_branch():
    value=problem('HHH',[],[]);a=saved_path(value,[(1,0,2)]);b=saved_path(value,[(0,2,1)])
    def combine(first,second):
        e=first.graph.transitions[0];f=second.graph.transitions[0]
        states=(SearchState(0,0,(),(),()),SearchState(1,0,(),(),()),replace(first.graph.states[-1],id=2))
        return replace(first.graph,states=states,transitions=(replace(e,source=0,target=1),replace(f,id=1,source=1,target=2)),stops=(SearchStop(2,'objective_met'),))
    left=combine(a,b);right=combine(b,a);index=FinalBranchCatalogue(value).add_graph(left).add_graph(right)
    assert len(index.branches)==1 and len(index.families)==2
    assert index.families[0].actions!=index.families[1].actions
    for old,family in zip((next(left.paths()),next(right.paths())),index.families):
        new=family.as_path(value)
        for mapping in ({0:0,1:1,2:2},{0:1,1:2,2:0},{0:2,1:0,2:1}):
            assert query_path(old,value,mapping,source_atoms=(0,1,2))[0]==query_path(new,value,mapping,source_atoms=(0,1,2))[0]


def test_pool_anchor_locks_survive_flattening():
    value=problem('HHHH',[],[]);old=saved_path(value,blocks=[dict(r_atoms=[1,2,3],p_atoms=[1,2,3])])
    root=replace(old.graph.states[0],mapping=((1,1),),islands=((1,0),))
    graph=replace(old.graph,states=(root,old.graph.states[1]));old=next(graph.paths())
    index=FinalBranchCatalogue(value).add_graph(graph);family=index.families[0]
    assert family.actions==( ('pool',(2,3)), )
    new=family.as_path(value)
    for mapping in ({0:0,1:1,2:3,3:2},{0:0,1:2,2:1,3:3}):
        assert query_path(old,value,mapping,source_atoms=(0,1,2,3))[0]==query_path(new,value,mapping,source_atoms=(0,1,2,3))[0]


def test_constraint_variants_are_retained_within_one_fragment_combination():
    value=problem('COO',[(0,1,1),(0,2,1)],[(0,1,1),(0,2,.4)])
    old=saved_path(value,[(0,2,1)]);other=replace(old.graph,
        transitions=tuple(replace(e, preserved_bonds=((0,2),)) for e in old.graph.transitions))
    index=FinalBranchCatalogue(value).add_graph(old.graph).add_graph(other)
    assert len(index.branches)==1 and len(index.families)==2
    for graph,family in zip((old.graph,other),index.families):
        a=extract_path_events(next(graph.paths()),value,SignedEventIndex(value),max_patterns=None)
        b=extract_path_events(family.as_path(value),value,SignedEventIndex(value),max_patterns=None)
        assert a['complete'] and b['complete']
        assert {p['id'] for p in a['patterns']}=={p['id'] for p in b['patterns']}


def test_intrinsic_group_has_no_prefix_dependence_and_preserves_match_tests():
    value=problem('COO',[(0,1,1.3),(0,2,1.3)],[(0,1,.2),(0,2,1.2)])
    index=FinalBranchCatalogue(value);pair=((0,1,2),(0,1,2));group=index.rebuild_symmetry((.2,1.),pair)
    assert group['generators']==() # swapping O atoms would change match admissibility
    assert index.rebuild_symmetry((.2,1.),pair) is group
    value=problem('COO',[(0,1,1),(0,2,1)],[(0,1,.5),(0,2,1)])
    index=FinalBranchCatalogue(value);group=index.rebuild_symmetry((.2,1.),pair)
    assert (0,2,1) in group['generators']
    for g in group['generators']:
        for a in pair[0]:
            for b in pair[0]:
                if value.reactant.wbo[a,b]<.2:continue
                w=value.reactant.wbo[a,b]
                for p in pair[1]:
                    for q in pair[1]:
                        before=value.product.wbo[p,q];after=value.product.wbo[g[p],g[q]]
                        assert (before>=.2 and abs(w-before)<=1+1e-9)==(after>=.2 and abs(w-after)<=1+1e-9)


def test_intrinsic_group_is_not_mistaken_for_observed_family_coverage():
    value=problem('COO',[(0,1,1),(0,2,1)],[(0,1,.5),(0,2,1)])
    old=saved_path(value);index=FinalBranchCatalogue(value).add_graph(old.graph)
    index.rebuild_all_symmetries()
    assert any(g['generators'] for g in index.symmetries.values())
    # Rebuilding a group must not grant a shuffle absent from all saved families.
    status,_=query_path(index.families[0].as_path(value),value,{0:0,1:2,2:1},source_atoms=(0,1,2))
    assert status=='not_recovered'










def test_public_branch_grouping_preserves_distinct_mapping_alternatives():
    value=problem('COO',[],[]);a=saved_path(value)
    b=replace(a.graph,states=(a.graph.states[0],replace(a.graph.states[1],mapping=((0,0),(1,2),(2,1)))))
    combined=AAMSearchGraph.combine((a.graph,b))
    assert len(combined.branches())==1
    assert len(combined.literal_branches())==2
    assert {tuple(sorted(p.mapping.items())) for p in combined.branches()[0].paths}=={((0,0),(1,1),(2,2)),((0,0),(1,2),(2,1))}


def test_flat_archive_roundtrip_needs_no_original_search_graph():
    import json
    value=problem('COO',[(0,1,1),(0,2,.9)],[(0,1,.51),(0,2,.49)])
    catalogue=FinalBranchCatalogue(value).add_graph(saved_path(value,[(0,2,1)]).graph)
    record=json.loads(json.dumps(catalogue.to_record()))
    restored=FinalBranchCatalogue.from_record(value,record)
    assert restored.branch_count==catalogue.branch_count
    assert restored.families[0].mapping==catalogue.families[0].mapping
    result=extract_path_events(restored.families[0].as_path(value),value,SignedEventIndex(value),max_patterns=None)
    assert result['complete'] and {p['total'] for p in result['patterns']}=={0,1}
    other=problem('COO',[(0,1,1)],[(0,1,.51),(0,2,.49)])
    with pytest.raises(ValueError,match='different endpoint'):
        FinalBranchCatalogue.from_record(other,record)




@pytest.mark.smoke
def test_aam_result_groups_atom_sets_before_return_without_anchor_key(monkeypatch):
    from mappa.mappa_core.domain import AAMResult, AAMSearchConfig, AAMSearchMetrics
    value = problem('COO', [], [])
    a = saved_path(value)
    # Different fixed mappings and policies belong to one atom-set branch.
    b = replace(a.graph,
        contexts=(replace(a.graph.contexts[0], iso_tolerance=.5),),
        states=(a.graph.states[0], replace(a.graph.states[1], mapping=((0,0),(1,2),(2,1)))))
    graph = AAMSearchGraph.combine((a.graph, b))
    result = AAMResult(value, AAMSearchConfig(), graph,
                       AAMSearchMetrics.from_record({'retained_branch_count': 2}, 0.))
    assert 'branches' in vars(result)  # constructed eagerly, before access/decoding
    assert len(result.branches) == result.metrics.retained_branch_count == 1
    assert result.branches[0].fragment_pairs == (((0,1,2),(0,1,2)),)
    assert {tuple(sorted(p.mapping.items())) for p in result.branches[0].paths} == {
        ((0,0),(1,1),(2,2)), ((0,0),(1,2),(2,1))}
    assert {p.context.iso_tolerance for p in result.branches[0].paths} == {1., .5}
    def forbidden(*args):
        raise AssertionError('returned branch groups must not be rebuilt on access')
    monkeypatch.setattr(AAMSearchGraph, 'branches', forbidden)
    assert len(result.branches) == 1


@pytest.mark.smoke
def test_aam_grouping_distinguishes_target_atom_sets():
    a = independent_path((0,1))
    b = independent_path((1,0), mapping={0:2,1:3,2:0,3:1})
    graph = AAMSearchGraph.combine((a.graph,b.graph))
    assert len(graph.branches()) == 2
    assert len({branch.fragment_pairs for branch in graph.branches()}) == 2
