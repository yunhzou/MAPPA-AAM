"""Verify a positive reference witness in a completed cut, without a full merge.

Cut graphs form a disjoint union. Final symmetry depends on the same fixed
target and each transition's locked mapping and candidate record, so it can
be finalized per cut. A positive witness belongs to the complete union. This
adapter never infers absence from an interrupted or partially checked search.
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
    from mappa.artifacts import read_raw_cut,raw_cut_paths
    from mappa.aam import checkpoint_manifest
    from mappa.domain import AAMResult,AAMSearchMetrics
    from golden_evaluation import evaluate_planned
    proof=read(ROOT/'checkpoint-verification-tests.json')
    assert proof['status']=='passed' and proof['driver_sha256']==sha(Path(__file__))
    started=time.perf_counter();cpu=time.process_time();deadline=started+275
    plan,config,folder=inputs(seed,case,direction)
    assert read(folder/'cuts/manifest.json')==checkpoint_manifest(plan.problem,config)
    reference=read(ROOT/'golden-inputs'/str(case)/'reference.json')
    files=raw_cut_paths(folder/'cuts');rows=[];hit=None
    for path in files:
        if time.perf_counter()>deadline:break
        before=time.perf_counter();graph=finalized(read_raw_cut(path),plan.problem,config)
        aam=AAMResult(plan.problem,config,graph,AAMSearchMetrics.from_record({},0))
        result=evaluate_planned(aam,plan,reference['features'],reference['mapping'],
                                seconds=min(20,max(.01,deadline-time.perf_counter())),query_timeout_ms=5000)
        row=dict(cut=str(path),sha256=sha(path),reference_recovery=result['reference_recovery'],
                 wall_seconds=time.perf_counter()-before)
        rows.append(row)
        if result['reference_recovery']=='recovered':
            hit=result;break
        del aam,graph,result;gc.collect()
        save(folder/'checkpoint-verification.json',dict(reference_recovery='unknown',checked=rows,
             available_cuts=len(files),cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-started))
    result=dict(hit or {},reference_recovery='recovered' if hit else 'unknown',
                evaluation_scope='Positive family membership in completed cuts; ranking and candidate-count fields, if present, are cut-local.',
                checkpoint_verification=dict(checked=rows,available_cuts=len(files),
                    driver_sha256=sha(Path(__file__)),proof_sha256=sha(ROOT/'checkpoint-verification-tests.json'),
                    scope='Positive witness in a completed, independently finalized raw cut; no claim that an interrupted combined search completed.'),
                cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-started)
    save(folder/'checkpoint-verification.json',result)

if __name__=='__main__':
    if sys.argv[1]=='test':selftest()
    else:child(int(sys.argv[1]),int(sys.argv[2]),sys.argv[3])
