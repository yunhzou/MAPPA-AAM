"""Current-engine trajectory diagnostics; reference-directed rows are feasibility only."""
import sys,json,time,hashlib
from pathlib import Path
ROOT=Path(sys.argv[1]);OUT=Path(sys.argv[2]);CASE=int(sys.argv[3])
sys.path[:0]=[str(ROOT/'engine/src'),str(ROOT/'engine/bench')]
import pynauty
from mappa import AAMProblem,MolecularEndpoint,AAMSearchConfig,AAMSearchPlan
from mappa.alignment.branch import find_islands,_generate_seed_orders
from mappa.matcher.policy import AttributeNodeMatchPolicy
from mappa.frag import build_graph
from mappa.aam import _initialize_search,_search_cut
from mappa.search_symmetry import finalize_graph_symmetry
from mappa.domain import AAMResult,AAMSearchMetrics
from golden_evaluation import colored_graph,project,evaluate_planned
raw=json.loads((ROOT/f'golden-inputs/{CASE}/input.json').read_text());ref=json.loads((ROOT/f'golden-inputs/{CASE}/reference.json').read_text())
p=AAMProblem(*(MolecularEndpoint(**raw[s]) for s in ('reactant','product')),raw.get('name',''));reverse=CASE==1285
config=AAMSearchConfig(seed_count=1,branch_limit=2000)
plan=AAMSearchPlan(p,AAMProblem(p.product,p.reactant,p.name) if reverse else p,config,reverse)
reference=plan.to_search_mapping(dict(ref['mapping']))
graphs=[build_graph(e.elements,e.wbo,config.graph_floor) for e in (plan.problem.reactant,plan.problem.product)]
for side,g in enumerate(graphs):
 for atom in g:
  labels=reference if side==0 else {p:r for r,p in reference.items()}
  label=(reference[atom] if side==0 else atom) if atom in labels else f'{side}:{atom}'
  g.nodes[atom]['reference_label']='H' if g.nodes[atom]['element']=='H' else str(label)
order=_generate_seed_orders(graphs[0],1,rng_seed=42)[0]
records=[];saved=[]
def certificate(mapping):return pynauty.certificate(colored_graph(ref['features'],project(plan.to_input_mapping(mapping),ref['features'])))
def changes(mapping):
 a,b=plan.problem.reactant.wbo,plan.problem.product.wbo
 return sum(a[r,s]!=b[mapping[r],mapping[s]] for r in reference for s in reference if r<s)
for policy in ['ordinary','ordinary_native','reference_constrained']:
 events=None if policy=='ordinary_native' else [];started=time.perf_counter()
 g=find_islands(*graphs,order,graph_floor=config.graph_floor,iso_tol=config.iso_tolerance,max_branches=2000,events=events,node_policy=AttributeNodeMatchPolicy(('element','reference_label')) if policy=='reference_constrained' else None)
 conflicts=[]
 for event in events or []:
  if event['type']!='commit':continue
  a,b=event['edge']['frag_atom'],event['edge']['ext_atom']
  if a in reference and b in reference:
   x=float(plan.problem.reactant.wbo[a,b]);y=float(plan.problem.product.wbo[reference[a],reference[b]])
   if abs(x-y)>config.iso_tolerance:conflicts.append(dict(source_edge=[a,b],target_edge=[reference[a],reference[b]],weights=[x,y],fragment=event['fragment']))
 mappings=[dict(g.states[t].mapping) for t in g.terminals]
 records.append(dict(policy=policy,reference_directed=policy=='reference_constrained',capped=g.capped,terminals=len(mappings),exact_reference_witnesses=sum(all(m.get(r)==p for r,p in reference.items()) for m in mappings),reference_equivalent_witnesses=sum(certificate(m)==certificate(reference) for m in mappings),heavy_pair_changes=[changes(m) for m in mappings if set(reference)<=set(m)],trace_conflicts=conflicts,seconds=time.perf_counter()-started,mappings=[sorted(m.items()) for m in mappings]))
 saved.append(g)
result=dict(case=CASE,direction=plan.direction,source='current frozen benchmark engine',native_python_graph_equal=saved[0]==saved[1],reference_heavy_pair_changes=changes(reference),rows=records,reference_mapping=sorted(reference.items()),seed_order=order,benchmark_score_unchanged=True)
if CASE==986:
 altered=dict(reference);altered[11],altered[12]=altered[12],altered[11]
 result['diagnostic_label_swap']=dict(source_indices=[11,12],original_rdf_labels=[15,16],heavy_pair_changes=changes(altered),ordinary_equivalent=any(certificate(m)==certificate(altered) for m in [dict(saved[0].states[t].mapping) for t in saved[0].terminals]))
if CASE==1285:
 config=AAMSearchConfig(seed_count=10,branch_limit=2000);plan=AAMSearchPlan(plan.input_problem,plan.problem,config,True)
 cuts=((7,8),(9,10),(10,20));_initialize_search(plan.problem,config,'reused_native');g,counts=_search_cut(cuts)
 g,_=finalize_graph_symmetry(g,build_graph(plan.problem.product.elements,plan.problem.product.wbo,config.graph_floor),iso_tolerance=config.iso_tolerance)
 ev=evaluate_planned(AAMResult(plan.problem,config,g,AAMSearchMetrics.from_record(counts,0)),plan,ref['features'],ref['mapping'],seconds=60,query_timeout_ms=5000)
 result['reference_selected_multicut']=dict(cuts=cuts,seed_count=10,branch_cap=2000,capped=g.capped,reference_recovery=ev['reference_recovery'],reference_directed_cut_selection=True,benchmark_score_unchanged=True)
OUT.mkdir(parents=True,exist_ok=True);(OUT/f'cause{CASE}.json').write_text(json.dumps(result,indent=2,default=int)+'\n');print(json.dumps(result,default=int),flush=True)
