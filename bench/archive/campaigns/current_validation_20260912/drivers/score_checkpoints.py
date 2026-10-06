"""Verify reference membership cut by cut, without a full merge.

Cut graphs form a disjoint union. Final symmetry depends on the same fixed
target and each transition's locked mapping and candidate record, so it can
be finalized per cut. A positive witness belongs to the complete union. This
adapter certifies absence only when every scheduled cut has a negative verdict.
"""
import sys,time,json,gc
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'engine/src'),str(ROOT/'engine/bench')]
from run import read,save,sha

def inputs(seed,case,direction):
    from mappa import AAMProblem,MolecularEndpoint,AAMSearchConfig
    from mappa.search_orientation import AAMSearchPlan
    raw=read(ROOT/'golden-inputs'/str(case)/'input.json')
    problem=AAMProblem(*(MolecularEndpoint(**raw[s]) for s in ('reactant','product')),raw.get('name',''))
    config=AAMSearchConfig(seed_count=seed,branch_limit=100)
    reverse=direction=='P_to_R'
    plan=AAMSearchPlan(problem,AAMProblem(problem.product,problem.reactant,problem.name) if reverse else problem,config,reverse)
    return plan,config,ROOT/'runs/golden'/f'seed{seed}'/f'case{case}'/direction

def finalized(graph,problem,config):
    from mappa.frag import build_graph
    from mappa.search_symmetry import finalize_graph_symmetry
    from mappa.conditioned_symmetry import ConditionedSymmetryWorkspace
    target=build_graph(problem.product.elements,problem.product.wbo,bond_cut=config.graph_floor)
    workspace=ConditionedSymmetryWorkspace(target,config.iso_tolerance)
    return finalize_graph_symmetry(graph,target,iso_tolerance=config.iso_tolerance,workspace=workspace)[0]

def selftest():
    from mappa.artifacts import read_raw_cut,raw_cut_paths
    from mappa.search_graph import AAMSearchGraph
    from mappa.aam import checkpoint_manifest
    rows=[]
    for case in (0,1):
        plan,config,folder=inputs(2,case,'R_to_P')
        assert read(folder/'cuts/manifest.json')==checkpoint_manifest(plan.problem,config)
        files=raw_cut_paths(folder/'cuts')[:2]
        raw=[read_raw_cut(p) for p in files]
        together=finalized(AAMSearchGraph.combine(raw),plan.problem,config)
        separate=AAMSearchGraph.combine([finalized(read_raw_cut(p),plan.problem,config) for p in files])
        assert together==separate
        rows.append(dict(case=case,cuts={str(p):sha(p) for p in files},states=len(together.states),
                         transitions=len(together.transitions),exact_graph_equality=True))
    save(ROOT/'checkpoint-verification-tests.json',dict(status='passed',driver_sha256=sha(Path(__file__)),
         rows=rows,scope='Exact equality of combined-then-finalized and separately-finalized-then-combined graphs, including states, contexts, transition actions, constraints and stops. No atom-mapping or permutation enumeration.'))

def child(seed,case,direction):
    from golden_checkpoint_evaluation import evaluate_checkpoints
    plan,config,folder=inputs(seed,case,direction)
    reference=read(ROOT/'golden-inputs'/str(case)/'reference.json')
    result=evaluate_checkpoints(folder/'cuts',plan,reference['features'],reference['mapping'],
        seconds=275,query_timeout_ms=5000,
        progress=lambda row:save(folder/'checkpoint-verification.json',row))
    save(folder/'checkpoint-verification.json',result)

if __name__=='__main__':
    if sys.argv[1]=='test':selftest()
    else:child(int(sys.argv[1]),int(sys.argv[2]),sys.argv[3])
