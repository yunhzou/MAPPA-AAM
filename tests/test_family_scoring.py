from itertools import permutations
import numpy as np

from mappa import AAMProblem,AAMSearchConfig,search_aam
from mappa.domain import MolecularEndpoint
from mappa.family_scoring import invariant_score,validate_representative,bond_events
from golden_evaluation import prepare
from publication_analysis import ranker


def endpoint(elements,edges):
    w=np.zeros((len(elements),len(elements)))
    for a,b,v in edges:w[a,b]=w[b,a]=v
    return MolecularEndpoint(tuple(elements),np.zeros((len(elements),3)),w)


def test_tolerant_symmetry_retains_all_small_query_mappings():
    mol=endpoint('CCCHHH',[(0,1,2),(1,2,1),(0,2,1),(0,3,1),(1,4,1),(2,5,1)])
    problem=AAMProblem(mol,mol)
    aam=search_aam(problem,AAMSearchConfig(seed_count=1,cut_floor=10,iso_tolerance=1))
    paths=list(aam.graph.paths())
    scores=[bond_events(problem,dict(enumerate((*p,*(i+3 for i in p)))))['total']
            for p in permutations(range(3))]
    assert set(scores)=={0,2}
    # Tolerance-compatible placements need not share an exact target orbit.
    # Every placement must remain in the union, with a proven family minimum.
    from mappa import query_path
    for p in permutations(range(3)):
        mapping=dict(enumerate((*p,*(i+3 for i in p))))
        assert any(query_path(path,problem,mapping,source_atoms=tuple(range(6)),
                              timeout_ms=2000)[0]=='recovered' for path in paths)
    for path in paths:
        validate_representative(path,problem)



def test_explicit_h_symmetry_certified_without_solver():
    problem,_,_=prepare('CC(C)(C)C>>CC(C)(C)C')
    path=next(search_aam(problem,AAMSearchConfig(seed_count=1,cut_floor=10)).graph.paths())
    invariant, tested = invariant_score(path,problem)
    assert invariant and tested > 0
    assert bond_events(problem,path.mapping)['total']==0 and len(path.mapping)==17


def test_event_score_matches_existing_ranker_including_unmatched_atoms():
    problem,_,_=prepare('CCO.O>>CC=O.O')
    mapping={i:i for i in range(min(problem.source_atom_count,problem.target_atom_count)-1)}
    for reverse in (False,True):
        actual=problem if not reverse else AAMProblem(problem.product,problem.reactant)
        oriented=mapping if not reverse else {v:k for k,v in mapping.items()}
        expected=ranker(actual)(oriented)[0][2]
        assert bond_events(problem,mapping,reverse=reverse)['total']==expected
