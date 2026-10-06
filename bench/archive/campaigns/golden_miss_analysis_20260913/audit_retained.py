"""Inspect saved ten-seed cuts using necessary invariants; no new AAM search."""
import sys,json,time,gc,hashlib,resource
from pathlib import Path
from collections import Counter
ROOT=Path(sys.argv[1]);OUT=Path(sys.argv[2]);CASE=int(sys.argv[3]);DIRECTION=sys.argv[4]
sys.path[:0]=[str(ROOT/'engine/src'),str(ROOT/'engine/bench')]
import numpy as np
from golden_evaluation import endpoint_generators,project,exact_action
from mappa import AAMProblem,MolecularEndpoint,AAMSearchConfig
from mappa.aam import checkpoint_manifest
from mappa.artifacts import raw_cut_paths,read_raw_cut
from mappa.alignment.sweep import cut_sweep_items
from mappa.frag import build_graph
from mappa.search_symmetry import finalize_graph_symmetry
from mappa.conditioned_symmetry import ConditionedSymmetryWorkspace

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
class Orbits:
 def __init__(self,n):self.parent=list(range(n))
 def find(self,i):
  while self.parent[i]!=i:self.parent[i]=self.parent[self.parent[i]];i=self.parent[i]
  return i
 def join(self,a,b):
  a,b=self.find(a),self.find(b)
  self.parent[max(a,b)]=min(a,b)
 def add(self,g):
  for i,j in enumerate(g):self.join(i,j)
 def labels(self):return tuple(self.find(i) for i in range(len(self.parent)))
raw=read(ROOT/f'golden-inputs/{CASE}/input.json');reference=read(ROOT/f'golden-inputs/{CASE}/reference.json')
problem=AAMProblem(*(MolecularEndpoint(**raw[s]) for s in ('reactant','product')),raw.get('name',''))
features=reference['features'];mapping=dict(reference['mapping'])
if DIRECTION=='P_to_R':
 problem=AAMProblem(problem.product,problem.reactant,problem.name);features=features[::-1];mapping={p:r for r,p in mapping.items()}
config=AAMSearchConfig(seed_count=10,branch_limit=100)
folder=ROOT/f'runs/golden/seed10/case{CASE}/{DIRECTION}/cuts'
assert read(folder/'manifest.json')==checkpoint_manifest(problem,config)
expected=cut_sweep_items(problem.reactant.wbo,config.cut_floor)
paths=raw_cut_paths(folder)
assert [int(p.name.split('_')[1].split('.')[0]) for p in paths]==list(range(len(expected)))
ref=project(mapping,features);nref=len(ref)
ro,po=Orbits(len(features[0]['heavy'])),Orbits(len(features[1]['heavy']))
for g in endpoint_generators(features[0]):ro.add(g)
for g in endpoint_generators(features[1]):po.add(g)
rlabels,plabels=ro.labels(),po.labels()
source_signature=Counter(rlabels[r] for r in ref)
reference_pairs=Counter((rlabels[r],plabels[p]) for r,p in ref.items())
target=build_graph(problem.product.elements,problem.product.wbo,bond_cut=config.graph_floor)
target_index={a:i for i,a in enumerate(features[1]['heavy'])}
rows=[];started=time.perf_counter();cpu=time.process_time()
for i,path in enumerate(paths):
 graph=read_raw_cut(path)
 assert all(tuple(c.cuts)==tuple(expected[i]) for c in graph.contexts)
 terminals=graph.terminals
 projected=[(t,project(graph.states[t].mapping,features)) for t in terminals]
 cardinality=Counter(len(m) for _,m in projected)
 same_card=[(t,m) for t,m in projected if len(m)==nref]
 source_ok=[(t,m) for t,m in same_card if Counter(rlabels[r] for r in m)==source_signature]
 row=dict(index=i,sha256=sha(path),terminals=len(terminals),cardinality_histogram=dict(cardinality),same_cardinality=len(same_card),source_orbit_selection_matches=len(source_ok),cap_stops=Counter(s.stage for s in graph.stops if s.reason=='capped'),states=len(graph.states),transitions=len(graph.transitions))
 row['representative_exact_orbit_pair_matches']=sum(Counter((rlabels[r],plabels[p]) for r,p in m.items())==reference_pairs for _,m in source_ok)
 if source_ok:
  workspace=ConditionedSymmetryWorkspace(target,config.iso_tolerance)
  graph,_=finalize_graph_symmetry(graph,target,iso_tolerance=config.iso_tolerance,states=[t for t,m in source_ok],workspace=workspace)
  del workspace
  allowed=Orbits(len(plabels))
  for g in endpoint_generators(features[1]):allowed.add(g)
  non_endpoint=0;actions=0;free_blocks=0
  selected=graph.ancestor_transitions([t for t,m in source_ok])
  for e in graph.transitions:
   if e.id not in selected or e.match is None:continue
   symmetry=e.match['symmetry']
   assert symmetry.get('automorph_group_source')=='conditioned_search_transition'
   for generator in symmetry.get('automorph_generators',()):
    projected_g=tuple(target_index[generator[a]] for a in features[1]['heavy'])
    actions+=1;non_endpoint+=not exact_action(projected_g,features[1]);allowed.add(projected_g)
   for b in symmetry.get('blocks',()):
    if b.get('source')=='exact_automorph_group':continue
    atoms=[target_index[a] for a in b['p_atoms'] if a in target_index]
    if len(atoms)>1:
     free_blocks+=1
     for a in atoms[1:]:allowed.join(atoms[0],a)
  labels=allowed.labels();sig=Counter((rlabels[r],labels[p]) for r,p in ref.items())
  row.update(heavy_action_records=actions,non_endpoint_actions=non_endpoint,free_heavy_blocks=free_blocks,relaxed_target_orbit_pair_matches=sum(Counter((rlabels[r],labels[p]) for r,p in m.items())==sig for _,m in source_ok))
 else:row['relaxed_target_orbit_pair_matches']=0
 rows.append(row)
 del graph,projected,same_card,source_ok;gc.collect()
OUT.mkdir(parents=True,exist_ok=True)
result=dict(case=CASE,direction=DIRECTION,seed_count=10,branch_cap=100,expected_cuts=len(expected),checked_cuts=len(rows),reference_pairs=nref,heavy_counts=[len(f['heavy']) for f in features],cuts=rows,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-started,peak_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024 if sys.platform=='linux' else 1024**2),input_sha256=sha(ROOT/f'golden-inputs/{CASE}/input.json'),reference_sha256=sha(ROOT/f'golden-inputs/{CASE}/reference.json'),driver_sha256=sha(__file__))
for key in ['terminals','same_cardinality','source_orbit_selection_matches','representative_exact_orbit_pair_matches','relaxed_target_orbit_pair_matches','heavy_action_records','non_endpoint_actions','free_heavy_blocks']:
 result[key]=sum(r.get(key,0) for r in rows)
result['cardinality_histogram']=dict(sum((Counter({int(k):v for k,v in r['cardinality_histogram'].items()}) for r in rows),Counter()))
result['cap_stops']=dict(sum((Counter(r['cap_stops']) for r in rows),Counter()))
(OUT/f'case{CASE}-{DIRECTION}.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='cuts'}),flush=True)
