from types import SimpleNamespace as NS
import numpy as np
import pytest
from mappa import MolecularEndpoint, AtomAssignment, AAMSearchConfig
from mappa.bond_diagnostics import endpoint_bond_changes

def fixture(mapping=(2,0,1)):
    xyz=np.array([[0.,0,0],[1.,0,0],[3.,0,0]])
    w=np.zeros((3,3));w[0,1]=w[1,0]=1
    r=MolecularEndpoint(('C',)*3,xyz,w)
    tw=w.copy();tw[0,2]=tw[2,0]=.8;tw[0,1]=tw[1,0]=.3
    order=np.argsort(mapping)
    t=MolecularEndpoint(('C',)*3,xyz[order],tw[np.ix_(order,order)])
    reaction=NS(problem=NS(reactant=r,product=r),config=AAMSearchConfig())
    mechanism=NS(mapping=AtomAssignment.from_mapping({0:0,1:1,2:2}),broken_bonds=((0,1),),formed_bonds=(),core_atoms=(0,1))
    target=NS(molecule=t,reactant_atom_mapping=mapping)
    assignment=AtomAssignment.from_mapping({0:mapping[0],1:mapping[1]})
    return reaction,mechanism,target,assignment

def test_full_mapping_finds_outside_core_change_in_both_directions():
    d=endpoint_bond_changes(*fixture())
    assert d['outside_core_assessed']
    for side in ('reactant_to_ts','product_to_ts'):
        changes=d[side]
        assert [c['kind'] for c in changes]==['breaking','forming']
        extra=changes[1]
        assert extra['reactant_atoms']==[0,2] and extra['ts_atoms']==[2,1]
        assert extra['touches_outside_reaction_core'] and not extra['intended_reaction_pair']
        assert extra['delta_wbo']==pytest.approx(.8)

def test_partial_mapping_does_not_claim_outside_core_coverage():
    args=fixture();args[2].reactant_atom_mapping=None
    d=endpoint_bond_changes(*args)
    assert not d['outside_core_assessed']
    assert len(d['reactant_to_ts'])==1

def test_invalid_full_mapping_is_rejected():
    args=fixture();args[2].reactant_atom_mapping=(0,0,1)
    with pytest.raises(ValueError,match='bijection'):endpoint_bond_changes(*args)

def test_mapping_disagreement_is_explicit():
    args=fixture();args[2].reactant_atom_mapping=(0,1,2)
    assert not endpoint_bond_changes(*args)['core_assignment_agrees_with_atom_identity']


def test_metal_pairs_use_existing_lower_threshold():
    args=fixture((0,1,2))
    r=args[0].problem.reactant; t=args[2].molecule
    elements=('Au','C','C')
    r2=MolecularEndpoint(elements,r.coordinates,r.wbo)
    tw=t.wbo.copy();tw[0,2]=tw[2,0]=.35
    args[0].problem.reactant=args[0].problem.product=r2
    args[2].molecule=MolecularEndpoint(elements,t.coordinates,tw)
    changes=endpoint_bond_changes(*args)['reactant_to_ts']
    assert any(c['reactant_atoms']==[0,2] and c['threshold']==.3 for c in changes)
