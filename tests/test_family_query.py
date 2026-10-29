from itertools import permutations
import numpy as np

from mappa import AAMProblem,AAMSearchConfig,search_aam,query_path
from mappa.mappa_core.domain import MolecularEndpoint
from mappa.mappa_core.core.search.graph import AAMSearchGraph, SearchContext, SearchState, FragmentTransition, SearchStop


def endpoint(elements,edges):
    w=np.zeros((len(elements),len(elements)))
    for a,b in edges:w[a,b]=w[b,a]=1
    return MolecularEndpoint(tuple(elements),np.zeros((len(elements),3)),w)


def synthetic_path(blocks,automorph_blocks=()):
    graph=AAMSearchGraph((SearchContext((0,1),(0,1,2),(0,1)),),(0,),
        (SearchState(0,0,(),(),()),SearchState(1,0,((0,0),(1,1)),((0,1),(1,1)),())),
        (FragmentTransition(0,0,1,0,(1,0),dict(fragment=[0,1],symmetry=dict(
            witness={0:0,1:1},blocks=blocks,exact_fixed=[],automorph_blocks=automorph_blocks,
            automorph_generators=[]))),),(SearchStop(1,'objective_met'),))
    return next(graph.paths())


def test_open_pools_keep_same_pair_and_different_pair_assignments():
    problem=AAMProblem(endpoint('OO',[]),endpoint('OOOO',[(0,1),(2,3)]))
    result=search_aam(problem,AAMSearchConfig(seed_count=1))
    path=next(result.graph.paths())
    # Exhaustive oracle only for this tiny unit test, never production search.
    for image in permutations(range(4),2):
        status,witness=query_path(path,problem,dict(enumerate(image)),source_atoms=(0,1),timeout_ms=2000)
        assert status=='recovered'
        assert dict(witness['mapping'])==dict(enumerate(image))


def test_correlated_fragment_is_not_independent_atom_membership():
    molecule=endpoint('CCCC',[(0,1),(1,2),(2,3)])
    problem=AAMProblem(molecule,molecule)
    result=search_aam(problem,AAMSearchConfig(seed_count=1,cut_floor=10,iso_tolerance=.5))
    paths=list(result.graph.paths())
    for image in [(0,1,2,3),(3,2,1,0)]:
        assert any(query_path(path,problem,dict(enumerate(image)),source_atoms=tuple(range(4)),
            target_generators=((3,2,1,0),),timeout_ms=2000)[0]=='recovered' for path in paths)
    assert all(query_path(path,problem,{0:0,1:2,2:1,3:3},source_atoms=tuple(range(4)),
        target_generators=((3,2,1,0),),timeout_ms=2000)[0]=='not_recovered' for path in paths)


def test_full_hydrogen_mapping_is_returned_and_one_to_one():
    molecule=endpoint('OHOH',[(0,1),(2,3)])
    problem=AAMProblem(molecule,molecule)
    result=search_aam(problem,AAMSearchConfig(seed_count=1,cut_floor=10))
    # Different target molecules are separate branches. Check their union,
    # without assuming the first branch is the swapped-component placement.
    answers=[query_path(path,problem,{0:2,2:0},source_atoms=(0,2),timeout_ms=3000)
             for path in result.graph.paths()]
    assert all(status in ('recovered','not_recovered') for status,_ in answers)
    witnesses=[witness for status,witness in answers if status=='recovered']
    assert witnesses
    witness=witnesses[0]
    mapping=dict(witness['mapping'])
    assert len(mapping)==len(set(mapping.values()))==4
    assert molecule.wbo[mapping[0],mapping[1]]==molecule.wbo[mapping[2],mapping[3]]==1


def test_one_target_cannot_supply_two_reference_atoms():
    molecule=endpoint('OO',[]);problem=AAMProblem(molecule,molecule)
    path=next(search_aam(problem,AAMSearchConfig(seed_count=1)).graph.paths())
    assert query_path(path,problem,{0:0,1:0},source_atoms=(0,1),timeout_ms=2000)[0]=='not_recovered'


def test_pool_choices_do_not_move_fixed_fragment_slots():
    problem=AAMProblem(endpoint('OO',[]),endpoint('OOO',[]))
    path=synthetic_path([dict(r_atoms=[1],p_atoms=[0,1,2])])
    assert query_path(path,problem,{0:0,1:2},source_atoms=(0,1))[0]=='recovered'
    assert query_path(path,problem,{0:2,1:0},source_atoms=(0,1))[0]=='not_recovered'


def test_automorphism_summary_is_not_an_independent_domain():
    problem=AAMProblem(endpoint('OO',[]),endpoint('OOO',[]))
    path=synthetic_path([],automorph_blocks=[dict(r_atoms=[0,1],p_atoms=[0,1,2])])
    assert query_path(path,problem,{0:0,1:2},source_atoms=(0,1))[0]=='not_recovered'


def test_heavy_projection_positive_cannot_replace_full_H_check():
    problem=AAMProblem(endpoint('OH',[(0,1)]),endpoint('OHO',[(0,1)]))
    path=synthetic_path([dict(r_atoms=[0],p_atoms=[0,2])])
    assert query_path(path,problem,{0:2},source_atoms=(0,),projected_atoms=(0,))[0]=='recovered'
    assert query_path(path,problem,{0:2},source_atoms=(0,))[0]=='not_recovered'


def test_reference_propagation_keeps_correlated_permutation_choices():
    import z3
    from mappa.mappa_workflow.mapping_analytics.constraints import ReferenceConstraints
    for target in permutations(range(3)):
        constraints=ReferenceConstraints()
        values=constraints.act([constraints.constant(i) for i in range(3)],[(1,2,0)],'ring')
        for value,wanted in zip(values,target):constraints.restrict(value[0],{wanted})
        solver=z3.Solver();constraints.emit(solver)
        expected=target in ((0,1,2),(1,2,0),(2,0,1))
        assert (solver.check()==z3.sat)==expected


def test_reference_propagation_dynamic_lookup_preserves_all_solutions():
    import z3
    from mappa.mappa_workflow.mapping_analytics.constraints import ReferenceConstraints
    # The reference fixes the output, not a freely chosen input slot.
    for selector in (0,1):
        for left in (2,3):
            for right in (2,3):
                constraints=ReferenceConstraints()
                choice=constraints.symbol((0,1))
                a=constraints.symbol((2,3));b=constraints.symbol((2,3))
                output=constraints.lookup(choice,{0:a,1:b})
                constraints.restrict(output[0],{2})
                constraints.restrict(choice[0],{selector})
                constraints.restrict(a[0],{left});constraints.restrict(b[0],{right})
                solver=z3.Solver();constraints.emit(solver)
                assert (solver.check()==z3.sat)==((left,right)[selector]==2)


def test_reference_propagation_overlapping_pools_matches_exhaustive_oracle():
    import z3
    from mappa.mappa_workflow.mapping_analytics.constraints import ReferenceConstraints
    # Two correlated permutations on overlapping target pools, with one slot
    # fixed in the second operation. Compare every output, including negatives.
    expected=set()
    for first in permutations((0,1,2)):
        for second in permutations((1,2,3)):
            table=dict(zip((1,2,3),second))
            if table[3]!=3:continue
            expected.add(tuple(table.get(v,v) for v in (*first,3)))
    for target in permutations(range(4)):
        constraints=ReferenceConstraints()
        values=[constraints.constant(i) for i in range(4)]
        for pool,fixed in (((0,1,2),()),((1,2,3),(3,))):
            table={p:constraints.symbol(pool) for p in pool}
            constraints.distinct(v[0] for v in table.values())
            for p in fixed:constraints.restrict(table[p][0],{p})
            values=[constraints.lookup(v,table) for v in values]
        for value,wanted in zip(values,target):constraints.restrict(value[0],{wanted})
        solver=z3.Solver();constraints.emit(solver)
        assert (solver.check()==z3.sat)==(target in expected)


def test_tolerant_symmetry_retains_all_small_query_mappings():
    mol=endpoint('CCCHHH',[(0,1),(1,2),(0,2),(0,3),(1,4),(2,5)])
    weights=mol.wbo.copy()
    weights[0,1]=weights[1,0]=2
    mol=MolecularEndpoint(mol.elements,mol.coordinates,weights)
    problem=AAMProblem(mol,mol)
    aam=search_aam(problem,AAMSearchConfig(seed_count=1,cut_floor=10,iso_tolerance=1))
    paths=list(aam.graph.paths())
    # Every tolerance-compatible placement must remain in the union.
    from mappa import query_path
    for p in permutations(range(3)):
        mapping=dict(enumerate((*p,*(i+3 for i in p))))
        assert any(query_path(path,problem,mapping,source_atoms=tuple(range(6)),
                              timeout_ms=2000)[0]=='recovered' for path in paths)
