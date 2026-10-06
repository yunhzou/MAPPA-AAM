from pathlib import Path
import sys,json,time,gzip
from dataclasses import replace
from run import S,PACKAGE
sys.path.insert(0,str(S/'engine/src'))
from run import save,read
from mappa.artifacts import read_aam_checkpoint,write_aam_checkpoint
from mappa.frag import build_graph
from mappa.matcher import _nauty_orbits
from mappa.native_search import find_islands_native
from mappa.search_graph import AAMSearchGraph
from mappa.search_symmetry import finalize_graph_symmetry
from mappa.domain import AAMResult,AAMSearchMetrics
from mappa.final_branches import FinalBranchCatalogue
from mappa.event_patterns import SignedEventIndex,extract_path_events
from mappa.competition import compete_fragments,CompetitionConfig

def run(case):
 c=time.process_time();out=S/f'seed-control/{case}';out.mkdir(parents=True,exist_ok=True)
 a=read_aam_checkpoint(S/f'runs/{case}/1/baseline.pkl.gz');p=a.problem;orders=list(dict.fromkeys(tuple(c.seed_order) for c in a.graph.contexts))
 source=build_graph(p.reactant.elements,p.reactant.wbo,bond_cut=.2);target=build_graph(p.product.elements,p.product.wbo,bond_cut=.2);po=_nauty_orbits(target,wbo_tol=1)
 graphs=[find_islands_native(source,target,order,graph_floor=.2,iso_tol=1,max_branches=2000,p_orbits=po,cuts=()) for order in orders]
 graph=AAMSearchGraph.combine(graphs);graph,_=finalize_graph_symmetry(graph,target,iso_tolerance=1)
 cfg=replace(a.config,sweep_cuts=False,seed_count=len(orders));b=AAMResult(p,cfg,graph,AAMSearchMetrics.from_record({},0));write_aam_checkpoint(b,out/'baseline.pkl.gz')
 search_cpu=time.process_time()-c;cat=FinalBranchCatalogue(p).add_aam(b,'baseline');base_ids=set(range(len(cat.families)))
 t=time.process_time();r=compete_fragments(b,CompetitionConfig(proposal_mode='prefix',operation_budget=128));comp_cpu=time.process_time()-t
 for i,repair in enumerate(r.repairs):
  assert all(not c.cuts for c in repair.graph.contexts)
  cat.add_aam(repair,f'prefix:{i}')
 idx=SignedEventIndex(p);K=read(PACKAGE/'controls.json')['per_case'][case]['window'];patterns={};baseline_patterns={};cert=[]
 for i,f in enumerate(cat.families):
  f.validate_representative(p);z=extract_path_events(f.as_path(p),p,idx,max_events=K,seconds=20,max_patterns=None);cert.append(dict(family=i,complete=z['complete'],reason=z['reason']))
  for item in z['patterns']:
   assert idx.describe(item['mapping'])['id']==item['id'];patterns[item['id']]=dict(item,family=i)
   if i in base_ids:baseline_patterns[item['id']]=item
 result=dict(case=case,orders=orders,order_count=len(orders),source_cut_contexts=len(a.graph.contexts),search_cpu=search_cpu,competition_cpu=comp_cpu,total_cpu=time.process_time()-c,complete=all(x['complete'] for x in cert),baseline_patterns=baseline_patterns,patterns=patterns,certificates=cert,competition_counts=r.counts,offers=r.offers,one_fragment_terminals=sum(len(set(dict(graph.states[t].islands).values()))==1 for t in graph.terminals),full_terminals=sum(len(graph.states[t].mapping)==p.source_atom_count for t in graph.terminals))
 save(out/'result.json',result);(out/'catalogue.json.gz').write_bytes(gzip.compress(json.dumps(cat.to_record()).encode(),mtime=0));print(case,'orders',len(orders),'min',min([x['total'] for x in patterns.values()]or[None]),'complete',result['complete'],'cpu',result['total_cpu'])
if __name__=='__main__':run(int(sys.argv[1]))
