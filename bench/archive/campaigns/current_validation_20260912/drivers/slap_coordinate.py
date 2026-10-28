"""Fresh SLAP mapping from the archived, verified original XYZ adjacencies."""
import sys,os,json,time,random,collections
from pathlib import Path
from types import SimpleNamespace
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parent;RUN=ROOT/'slap-coordinate';REPO=Path('/Users/yunhengz/Desktop/AAM Writing')
sys.path[:0]=[str(ROOT/'slap-upstream/src'),str(ROOT/'engine/src'),str(ROOT/'engine/bench')]
from resources import read,save,sha,MemoryGuard,bounded_process

def prepare():
 from rdkit import Chem
 rows=read(REPO/'reports/holdout_cap1000_seed1_20260910/slap_xyz.json.gz');calls=read(REPO/'reports/holdout_slap_xyz_sweep_20260910/native_call_records.json.gz');sources={}
 for case in range(140):
  raw=read(ROOT.parent/'full140_inputs'/str(case)/'input.json');cs=rows[str(case)]['candidates'];graphs=cs[0]['graphs'];original=[]
  for side,g in enumerate(graphs):
   edges=sorted(tuple(e) for e in g['edges']);elements=g['elements']
   assert elements==[Chem.GetPeriodicTable().GetAtomicNumber(e) for e in raw[('reactant','product')[side]]['elements']]
   assert all(sorted(tuple(e) for e in c['graphs'][side]['edges'])==edges and c['graphs'][side]['elements']==elements for c in cs)
   direction=('R_to_P','P_to_R')[side];edge_list=[tuple(row['cut']) for row in calls[f'{case}/{direction}'] if row['cut'] is not None]
   assert edge_list==[(a,b) for a,b,w in edges]
   assert all(w==1 for a,b,w in edges)
   original.append(dict(elements=elements,edges=edges))
  save(RUN/'inputs'/str(case)/'graphs.json',original);save(RUN/'inputs'/str(case)/'input.json',raw)
  sources[str(case)]=sha(RUN/'inputs'/str(case)/'graphs.json')
 save(RUN/'manifest.json',dict(cases=140,upstream_commit='ea248fd9494f52f4865193e87a98cc92c62b5f9e',adapter_sha256=sha(Path(__file__)),input_sha256=sources,archived_graph_sha256=sha(REPO/'reports/holdout_cap1000_seed1_20260910/slap_xyz.json.gz'),scope='Fresh uncut and both-direction single-edge SLAP sweeps. Original XYZ-derived adjacency and atom order restored from archived native graphs; atom-number labels reinitialized. Verified edges against every archived cut. XYZ preparation is not rerun or included in timing.'))

def child(case):
 import numpy as np
 from slapmapper.core import LabeledGraph
 from slapmapper.aam import SlapAAM
 from slap_edge_sweep import cut_graphs
 from holdout_slap_sweep import label_key
 from holdout_cap_seed_benchmark import score_slap
 random.seed(20260910);np.random.seed(20260910);graphs=read(RUN/'inputs'/str(case)/'graphs.json')
 native=[];sweep=[];seen=set();rows=[];cpu=0.;start=time.perf_counter()
 for rev in (False,True):
  random.seed(20260910);np.random.seed(20260910);gs=graphs[::-1] if rev else graphs;base=[]
  for g in gs:
   adjacency={i:{} for i in range(len(g['elements']))}
   for a,b,w in g['edges']:adjacency[a][b]=adjacency[b][a]=w
   base.append(LabeledGraph(adjacency,list(g['elements'])))
  targets=[i for i,z in enumerate(gs[0]['elements']) if z>1];mapper=SlapAAM(binary=True)
  for ordinal,edge in enumerate([None,*[(a,b) for a,b,w in gs[0]['edges']]]):
   t=time.process_time();mapper.get_maps(cut_graphs(base,edge),break_sym_targets=targets,base=0);cpu+=time.process_time()-t
   for result in mapper.results:
    labels=[list(map(int,g.labels)) for g in result['lgp']]
    if rev:labels.reverse()
    c=dict(graphs=[dict(labels=l) for l in labels]);key=label_key(c)
    if not rev and ordinal==0:native.append(c)
    if key not in seen:seen.add(key);sweep.append(c)
   rows.append(dict(direction='P_to_R' if rev else 'R_to_P',cut=edge,candidates=len(mapper.results)))
 # Use the same symbolic H-label refinement as the archived comparator.
 for method,candidates in [('native_slap',native),('slap_sweep',sweep)]:
  d=RUN/method
  if not (d/'inputs').exists():d.mkdir(parents=True,exist_ok=True);(d/'inputs').symlink_to(RUN/'inputs',target_is_directory=True)
  save(d/'slap_xyz'/f'{case}.json',dict(status='mapped',candidates=candidates))
  (d/'slap_scoring').mkdir(exist_ok=True)
  score_slap(SimpleNamespace(run=d,index=case))
 save(RUN/'results'/f'{case}.json',dict(case=case,complete=True,mapping_cpu_seconds=cpu,wall_seconds=time.perf_counter()-start,native_families=len(native),sweep_families=len(sweep),calls=rows))

def main():
 # Create shared directories before worker launch.
 for method in ('native_slap','slap_sweep'):
  d=RUN/method;d.mkdir(parents=True,exist_ok=True)
  if not (d/'inputs').exists():(d/'inputs').symlink_to(RUN/'inputs',target_is_directory=True)
  (d/'slap_scoring').mkdir(exist_ok=True)
 guard=MemoryGuard(3072,8192,6144);env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1');start=time.perf_counter()
 def run(c):
  d=RUN/'execution'/str(c);d.mkdir(parents=True,exist_ok=True)
  with (d/'run.log').open('w') as log:state,peak=bounded_process([sys.executable,__file__,str(c)],env,log,300,guard)
  row=dict(case=c,status=state,peak_mib=peak);save(d/'status.json',row);print(json.dumps(row),flush=True);return row
 with ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(run,range(140)))
 save(RUN/'execution.json',dict(rows=rows,wall_seconds=time.perf_counter()-start,peak_mib=guard.peak_total_kib/1024))
if __name__=='__main__':
 if sys.argv[1]=='prepare':prepare()
 elif sys.argv[1]=='all':main()
 else:child(int(sys.argv[1]))
