"""Paired sweep x competition ablation with complete per-family window certificates."""
from pathlib import Path
import sys,os,json,gzip,time,subprocess,signal,hashlib,resource
from dataclasses import asdict
from concurrent.futures import ThreadPoolExecutor,as_completed
PACKAGE=Path(__file__).resolve().parent
S=Path(os.environ.get("MAPPA_EXPERIMENT_WORK",PACKAGE/"work")).resolve()
sys.path.insert(0,str(S/'engine/src'))
def read(p):return json.loads(Path(p).read_text())
def save(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_suffix('.tmp');q.write_text(json.dumps(x)+'\n');q.replace(p)
def worker(case,sweep,stage):
 from mappa import AAMProblem,MolecularEndpoint,AAMSearchConfig,search_aam
 from mappa.artifacts import read_aam_checkpoint,write_aam_checkpoint
 out=S/f'runs/{case}/{int(sweep)}';out.mkdir(parents=True,exist_ok=True)
 c=time.process_time();w=time.perf_counter()
 if stage=='search':
  inp=read(PACKAGE/f'inputs/{case}.json');p=AAMProblem(*(MolecularEndpoint(**inp[k]) for k in ['reactant','product']),inp.get('name',''))
  cfg=AAMSearchConfig(seed_count=1,branch_limit=2000,iso_tolerance=1,graph_floor=.2,cut_floor=.2,sweep_cuts=sweep,random_seed=42,event_threshold=.5,metal_event_threshold=.3)
  a=search_aam(p,cfg,execution='reused_native',workers=1)
  write_aam_checkpoint(a,out/'baseline.pkl.gz')
  result=dict(config=asdict(cfg),metrics=asdict(a.metrics),full_terminals=sum(len(a.graph.states[t].mapping)==p.source_atom_count for t in a.graph.terminals),contexts=len(a.graph.contexts),capped=a.graph.capped)
 elif stage in ['local','prefix']:
  from mappa.competition import compete_fragments,CompetitionConfig
  a=read_aam_checkpoint(out/'baseline.pkl.gz');cfg=CompetitionConfig(proposal_mode=stage,operation_budget=128)
  (out/stage).mkdir(exist_ok=True)
  r=compete_fragments(a,cfg)
  for i,a in enumerate(r.repairs):
   if not sweep: assert all(not context.cuts for context in a.graph.contexts)
   write_aam_checkpoint(a,out/f'{stage}/repair{i}.pkl.gz')
  result=dict(config=asdict(cfg),counts=r.counts,repairs=len(r.repairs),pending=r.pending,offers=r.offers)
 elif stage=='decode':
  from mappa.final_branches import FinalBranchCatalogue
  from mappa.event_patterns import SignedEventIndex,extract_path_events
  a=read_aam_checkpoint(out/'baseline.pkl.gz');p=a.problem;cat=FinalBranchCatalogue(p).add_aam(a,'baseline');modes=['none','local','prefix'];family_ids={};counts={}
  def group_ids(mode):
   return [i for i,f in enumerate(cat.families) if any(x['archive']=='baseline' or x['archive'].startswith(mode+':') for x in f.provenance)]
  family_ids['none']=list(range(len(cat.families)));counts['none']=dict(branches=cat.branch_count,families=len(cat.families))
  for mode in ['local','prefix']:
   for i in range(read(out/f'{mode}.json')['repairs']):cat.add_aam(read_aam_checkpoint(out/f'{mode}/repair{i}.pkl.gz'),f'{mode}:{i}')
  for mode in ['local','prefix']:family_ids[mode]=group_ids(mode)
  # Families common to arms are decoded once; each arm retains its own membership.
  K=read(PACKAGE/'controls.json')['per_case'][case]['window'];idx=SignedEventIndex(p);done={};journal=out/'families.jsonl'
  if journal.exists():
   for line in journal.read_text().splitlines():
    z=json.loads(line)
    if z['complete']:done[z['family']]=z
  reverse={m:set(v) for m,v in family_ids.items()}
  roles={i:[m for m in modes if i in reverse[m]] for i in range(len(cat.families))}
  deadline=w+275
  with journal.open('a',buffering=1) as stream:
   for i,f in enumerate(cat.families):
    if i in done:continue
    if time.perf_counter()>=deadline:break
    f.validate_representative(p)
    fc=time.process_time()
    z=extract_path_events(f.as_path(p),p,idx,max_events=K,seconds=min(20,deadline-time.perf_counter()),max_patterns=None)
    for item in z['patterns']:assert idx.describe(item['mapping'])['id']==item['id']
    z=dict(family=i,roles=roles[i],extract_cpu=time.process_time()-fc,**z);stream.write(json.dumps(z)+'\n')
    if z['complete']:done[i]=z
    if i%100==0:save(out/'decode-progress.json',dict(completed=len(done),families=len(cat.families)))
  arms={}
  for mode in modes:
   pats={}
   for i in family_ids[mode]:
    if i in done:
     for item in done[i]['patterns']:pats.setdefault(item['id'],dict(item,family=i))
   complete=reverse[mode]<=done.keys()
   arms[mode]=dict(complete=complete,families=len(family_ids[mode]),completed_families=len(reverse[mode]&done.keys()),patterns=pats,minimum=min([v['total'] for v in pats.values()] or [None]),no_full_mapping=not family_ids[mode],minimum_lower_bound=K+1 if complete and family_ids[mode] and not pats else None)
  (out/'catalogue.json.gz').write_bytes(gzip.compress(json.dumps(cat.to_record()).encode(),mtime=0))
  result=dict(window=K,arms=arms,complete=all(a['complete'] for a in arms.values()),families=len(cat.families),family_ids=family_ids,certificates=[{k:v for k,v in z.items() if k!='patterns'} for z in done.values()])
 else:raise ValueError(stage)
 result.update(case=case,sweep=sweep,stage=stage,cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w,peak_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20)
 save(out/f'{stage}.json',result)
 print(json.dumps(dict(case=case,sweep=sweep,stage=stage,cpu=result['cpu_seconds'],complete=result.get('complete'))),flush=True)
def launch(case,sweep,stage):
 import psutil
 out=S/f'runs/{case}/{int(sweep)}';out.mkdir(parents=True,exist_ok=True)
 env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
 begin=time.perf_counter();peak=0;status=None
 with (out/f'{stage}.log').open('a') as log:
  p=subprocess.Popen([sys.executable,__file__,'worker',str(case),str(int(sweep)),stage],env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
  while p.poll() is None:
   try:peak=max(peak,psutil.Process(p.pid).memory_info().rss/2**20)
   except psutil.NoSuchProcess:break
   if time.perf_counter()-begin>300 or peak>3000:
    status='watchdog' if peak<=3000 else 'memory_limit';os.killpg(p.pid,signal.SIGKILL);break
   time.sleep(.2)
  code=p.wait()
 status=status or ('passed' if code==0 else 'failed')
 result=dict(case=case,sweep=sweep,stage=stage,status=status,exit_code=code,peak_mib=peak,wall_seconds=time.perf_counter()-begin)
 save(out/f'{stage}-execution.json',result);print(json.dumps(result),flush=True);return result
def pipeline(args):
 case,sweep=args;rows=[]
 for stage in ['search','local','prefix','decode']:
  out=S/f'runs/{case}/{int(sweep)}'
  if (out/f'{stage}-execution.json').exists() and read(out/f'{stage}-execution.json')['status']=='passed':continue
  row=launch(case,sweep,stage);rows.append(row)
  if row['status']!='passed':break
 return rows
if __name__=='__main__':
 if len(sys.argv)>1 and sys.argv[1]=='worker':worker(int(sys.argv[2]),bool(int(sys.argv[3])),sys.argv[4])
 else:
  jobs=[(c,s) for c in read(PACKAGE/'manifest.json')['cases'] for s in [False,True]]
  with ThreadPoolExecutor(max_workers=2) as pool:
   rows=[r for group in pool.map(pipeline,jobs) for r in group]
  save(S/'execution.json',rows)
