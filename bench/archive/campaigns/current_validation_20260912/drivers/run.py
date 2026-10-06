"""Fixed-source fresh validation; search inputs exclude reference map labels."""
import os,sys,json,time,hashlib,platform,random,gzip
from pathlib import Path
from dataclasses import asdict
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parent; ENGINE=ROOT/'engine'; WORK=ROOT.parent
sys.path[:0]=[str(ENGINE/'src'),str(ENGINE/'bench')]
from run_event_campaign import MemoryGuard,bounded_process,save,read,sha

def prepare_inputs():
 from golden_evaluation import prepare
 rows=[]
 for line in (ROOT/'golden-data/audit.jsonl').read_text().splitlines():
  row=json.loads(line);p,f,r=prepare(row['mapped_reaction']);d=ROOT/'golden-inputs'/str(row['index'])
  def ep(e):return dict(elements=e.elements,coordinates=e.coordinates.tolist(),wbo=e.wbo.tolist(),label=e.label,metadata=dict(e.metadata))
  save(d/'input.json',dict(reactant=ep(p.reactant),product=ep(p.product),name=f'golden_{row["index"]}'))
  save(d/'reference.json',dict(features=f,mapping=sorted(r.items())))
  rows.append(dict(case=row['index'],atoms=max(p.source_atom_count,p.target_atom_count),input_sha256=sha(d/'input.json')))
 ordered=sorted(rows,key=lambda r:(r['atoms'],r['case']));pilot=[ordered[round(i*(len(rows)-1)/19)]['case'] for i in range(20)]
 save(ROOT/'inputs.json',dict(records=rows,pilot=pilot,selection='20 evenly spaced explicit-atom-size order statistics; no outcomes used'))
 save(ROOT/'engine.json',dict(commit='eb1a7d150e9c3c80da79639ea59c300d8f834038',python=sys.version,platform=platform.platform(),sources={str(p.relative_to(ENGINE)):sha(p) for p in ENGINE.rglob('*') if p.is_file() and p.suffix in ('.py','.cpp','.h','.so') and 'build' not in p.parts},dataset=read(ROOT/'golden-data/manifest.json')))
 print(json.dumps(dict(records=len(rows),pilot=pilot)),flush=True)

def child(dataset,seed,case,direction,phase):
 import numpy as np
 from mappa import AAMProblem,MolecularEndpoint,AAMSearchConfig,search_aam,search_aam_checkpoints
 from mappa.search_orientation import AAMSearchPlan
 from golden_checkpoint_evaluation import evaluate_checkpoints
 random.seed(42);np.random.seed(42)
 raw=read((ROOT/'golden-inputs'/str(case) if dataset=='golden' else WORK/'full140_inputs'/str(case))/'input.json')
 p=AAMProblem(*(MolecularEndpoint(**raw[s]) for s in ('reactant','product')),raw.get('name',''))
 cfg=AAMSearchConfig(seed_count=seed,branch_limit=100 if dataset=='golden' else 2000)
 rev=direction=='P_to_R';plan=AAMSearchPlan(p,AAMProblem(p.product,p.reactant,p.name) if rev else p,cfg,rev)
 folder=ROOT/'runs'/dataset/f'seed{seed}'/f'case{case}'/direction
 start=time.perf_counter();cpu=time.process_time()
 if phase=='search':
  if dataset=='golden':
   a=search_aam_checkpoints(plan.problem,cfg,execution='reused_native',workers=1,intermediate_dir=folder/'cuts',resume=(folder/'cuts/manifest.json').exists())
   archive=dict(archive_kind='raw_cut_checkpoints',manifest_sha256=sha(folder/'cuts/manifest.json'))
   capped=a.capped
  else:
   a=search_aam(plan.problem,cfg,execution='reused_native',workers=1,intermediate_dir=folder/'cuts',archive_format='checkpoint',resume=(folder/'cuts/manifest.json').exists())
   archive=dict(archive_sha256=sha(folder/'cuts/aam.pkl.gz'));capped=a.graph.capped
  save(folder/'search.json',dict(case=case,seed=seed,direction=direction,config=asdict(cfg),metrics=asdict(a.metrics),wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,capped=capped,**archive))
 else:
  ref=read(ROOT/'golden-inputs'/str(case)/'reference.json')
  record=read(folder/'search.json')
  if record.get('archive_kind')=='raw_cut_checkpoints':
   assert record['manifest_sha256']==sha(folder/'cuts/manifest.json')
  result=evaluate_checkpoints(folder/'cuts',plan,ref['features'],ref['mapping'],seconds=220,query_timeout_ms=5000)
  result.update(cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start)
  save(folder/'evaluation.json',result)

def campaign(dataset,seeds,cases,workers,label):
 workers=min(workers,8,max(1,(os.cpu_count() or 2)-1));guard=MemoryGuard(3072,8192,6144)
 env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1')
 tasks=[(s,c,d) for s in seeds for c in cases for d in (('R_to_P','P_to_R') if dataset=='golden' else ('R_to_P',))]
 out=ROOT/'campaigns'/label;out.mkdir(parents=True,exist_ok=True)
 identity=dict(dataset=dataset,seeds=seeds,cases=cases,engine_sha256=sha(ROOT/'engine.json'),driver_sha256=sha(Path(__file__)),workers=workers,watchdog=300,memory_mib=dict(worker=3072,total=8192,reserve=6144))
 if (out/'manifest.json').exists():assert read(out/'manifest.json')==identity
 else:save(out/'manifest.json',identity)
 start=time.perf_counter()
 def run(task):
  s,c,d=task;folder=ROOT/'runs'/dataset/f'seed{s}'/f'case{c}'/d;folder.mkdir(parents=True,exist_ok=True)
  reports=[]
  for phase,target in [('search','search.json')]+([('score','evaluation.json')] if dataset=='golden' else []):
   path=folder/f'{phase}_execution.json'
   if (folder/target).exists() and path.exists() and read(path)['status']=='passed':reports.append(read(path));continue
   with (folder/f'{phase}.log').open('w') as log:
    t=time.perf_counter();state,peak=bounded_process([sys.executable,str(Path(__file__)),'child',dataset,str(s),str(c),d,phase],env,log,300,guard)
   row=dict(phase=phase,status=state,peak_mib=peak,elapsed_seconds=time.perf_counter()-t);save(path,row);reports.append(row)
   if state!='passed':break
  row=dict(seed=s,case=c,direction=d,phases=reports)
  print(json.dumps(row),flush=True);return row
 rows=[]
 with ThreadPoolExecutor(max_workers=workers) as pool:
  for f in as_completed([pool.submit(run,t) for t in tasks]):
   rows.append(f.result());save(out/'progress.json',dict(completed=len(rows),total=len(tasks),wall_seconds=time.perf_counter()-start,peak_mib=guard.peak_total_kib/1024,rows=rows))
 save(out/'execution.json',dict(wall_seconds=time.perf_counter()-start,peak_mib=guard.peak_total_kib/1024,rows=rows))
if __name__=='__main__':
 if sys.argv[1]=='prepare':prepare_inputs()
 elif sys.argv[1]=='child':child(sys.argv[2],int(sys.argv[3]),int(sys.argv[4]),sys.argv[5],sys.argv[6])
 else:
  dataset=sys.argv[2];seeds=list(map(int,sys.argv[3].split(',')));arg=sys.argv[4]
  cases=read(ROOT/'inputs.json')['pilot'] if arg=='pilot' else list(range(1851 if dataset=='golden' else 140)) if arg=='all' else list(map(int,arg.split(',')))
  campaign(dataset,seeds,cases,int(sys.argv[5]),sys.argv[6])
