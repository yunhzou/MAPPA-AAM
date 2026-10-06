"""Controlled, paired Golden search/verification. Each reaction owns its outputs."""
import os,sys,json,time,hashlib,socket,platform,subprocess,signal,resource,random,gzip,traceback
from pathlib import Path
from dataclasses import asdict
ROOT=Path(__file__).resolve().parent
BASE=Path(os.environ.get('AAM_VALIDATION_BASE',str(ROOT.parent/'validation_20260912')))
sys.path[:0]=[str(ROOT/'engine/src'),str(ROOT/'engine/bench'),str(BASE/'slap-upstream/src')]
ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1')

def read(p):return json.loads(Path(p).read_text())
def save(p,v):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(v)+'\n');q.replace(p)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def configurations():
 return [dict(key=f'mappa_c{cap}_s{seed}_{"sweep" if sweep else "uncut"}',method='mappa',cap=cap,seed=seed,sweep=sweep) for cap in (100,2000) for seed,sweep in [(1,False),(1,True),(2,True),(3,True),(10,True)]]+[dict(key='slap_'+('sweep' if sweep else 'uncut'),method='slap',sweep=sweep) for sweep in (False,True)]
def setting(key):return next(c for c in configurations() if c['key']==key)
def folder(case,key,direction):return ROOT/'results'/str(case)/key/direction

def plan_for(case,cfg,direction):
 from mappa import AAMProblem,MolecularEndpoint,AAMSearchConfig
 from mappa.search_orientation import AAMSearchPlan
 raw=read(BASE/'golden-inputs'/str(case)/'input.json')
 problem=AAMProblem(*(MolecularEndpoint(**raw[s]) for s in ('reactant','product')),raw.get('name',''))
 config=AAMSearchConfig(seed_count=cfg['seed'],branch_limit=cfg['cap'],sweep_cuts=cfg['sweep'])
 reverse=direction=='P_to_R'
 return AAMSearchPlan(problem,AAMProblem(problem.product,problem.reactant,problem.name) if reverse else problem,config,reverse)

def begin(dest,**info):
 # Saved before the call so an interrupted run still has a measured CPU baseline.
 baseline=time.process_time();start=time.perf_counter()
 save(dest/'search-start.json',dict(cpu_baseline=baseline,wall_start=start,**info))
 return baseline,start

def mappa_search(case,cfg,direction,dest):
 from mappa import search_aam_checkpoints
 import mappa.aam as core
 plan=plan_for(case,cfg,direction)
 compute={'cpu_seconds':0.,'wall_seconds':0.}
 # Read-only timers around unchanged native search calls; output handling remains
 # inside the separately reported full API-call timer.
 def timed(fn):
  def wrapper(*a,**kw):
   c=time.process_time();w=time.perf_counter()
   try:return fn(*a,**kw)
   finally:compute['cpu_seconds']+=time.process_time()-c;compute['wall_seconds']+=time.perf_counter()-w
  return wrapper
 core._initialize_search=timed(core._initialize_search)
 core._search_cut=timed(core._search_cut)
 c,w=begin(dest,config=asdict(plan.config))
 result=search_aam_checkpoints(plan.problem,plan.config,execution='reused_native',workers=1,intermediate_dir=dest/'cuts',resume=False)
 api_cpu=time.process_time()-c;api_wall=time.perf_counter()-w
 save(dest/'search.json',dict(cpu_seconds=api_cpu,wall_seconds=api_wall,compute=compute,config=asdict(plan.config),metrics=asdict(result.metrics),capped=result.capped,manifest_sha256=sha(dest/'cuts/manifest.json'),search_complete=True))

def mappa_score(case,cfg,direction,dest):
 from golden_checkpoint_evaluation import evaluate_checkpoints
 plan=plan_for(case,cfg,direction)
 if not (dest/'cuts/manifest.json').exists():save(dest/'evaluation.json',dict(reference_recovery='unknown',reason='no search checkpoint'));return
 ref=read(BASE/'golden-inputs'/str(case)/'reference.json')
 result=evaluate_checkpoints(dest/'cuts',plan,ref['features'],ref['mapping'],seconds=275,query_timeout_ms=5000,progress=lambda x:save(dest/'verification-progress.json',x))
 save(dest/'evaluation.json',result)

def slap_search(case,cfg,direction,dest):
 import numpy as np
 from slapmapper.core import SlapMapper
 from slapmapper.aam._smiles import smiles2lgp,get_numbered_rxn_smiles
 import slap_adapter as m
 if not hasattr(os,'sched_getaffinity'):os.sched_getaffinity=lambda pid:{0}
 from types import SimpleNamespace
 # Inputs and imports excluded from the API-call timer for both mappers.
 run=dest/'slap';run.mkdir(parents=True,exist_ok=True)
 for name,src in [('tasks.json','slap-tasks.json'),('inputs.json','slap-inputs.json')]:
  p=run/name
  if not p.exists():p.symlink_to(ROOT/src)
 slots=[4*case+i for i in ([0,1] if direction=='R_to_P' else [2,3])]
 c,w=begin(dest,config=cfg)
 for slot in slots:m.case(SimpleNamespace(run=run,slot=slot,sweep_cuts=cfg['sweep']))
 api_cpu=time.process_time()-c;api_wall=time.perf_counter()-w
 totals=dict(graph_cpu=0.,mapping_cpu=0.,export_cpu=0.,encoding_cpu=0.);count=errors=0
 for slot in slots:
  task=read(ROOT/'slap-tasks.json')[slot]
  for row in m.records(m.folder(run,task)/'records.jsonl'):
   count+=1;errors+=row['status']!='mapped'
   for k in totals:totals[k]+=row[k]
 save(dest/'search.json',dict(cpu_seconds=api_cpu,wall_seconds=api_wall,stages=totals,cuts=count,mapping_errors=errors,search_complete=errors==0,config=cfg))

def slap_score(case,cfg,direction,dest):
 from golden_competitors import signatures
 import slap_adapter as m
 ref=json.loads((BASE/'golden-data/audit.jsonl').read_text().splitlines()[case]);ep0,sig0=signatures(ref['mapped_reaction'])
 slots=[4*case+i for i in ([0,1] if direction=='R_to_P' else [2,3])];tasks=read(ROOT/'slap-tasks.json')
 hits=[];complete=True;invalid=0;errors=0;unique=set()
 for slot in slots:
  task=tasks[slot];rows=m.records(m.folder(dest/'slap',task)/'records.jsonl')
  expected=task['expected_cuts'] if cfg['sweep'] else 1
  complete &= [r['ordinal'] for r in rows]==list(range(expected))
  for row in rows:
   errors+=row['status']!='mapped'
   for ordinal,candidate in enumerate(row['candidates']):
    try:
     ep,sig=signatures(candidate['mapped_rxn']);assert ep==ep0
    except Exception:invalid+=1;continue
    unique.add(sig)
    if sig==sig0:hits.append(dict(mode=task['mode'],cut=row['ordinal'],candidate=ordinal,mapped_rxn=candidate['mapped_rxn']))
 verdict='recovered' if hits else 'not_recovered' if complete and not errors and not invalid else 'unknown'
 save(dest/'evaluation.json',dict(reference_recovery=verdict,hits=hits,complete=bool(complete and not errors and not invalid),invalid=invalid,errors=errors,unique_classes=len(unique)))

def child(case,key,direction,phase):
 import numpy as np
 random.seed(42);np.random.seed(42)
 cfg=setting(key);dest=folder(case,key,direction);dest.mkdir(parents=True,exist_ok=True)
 c=time.process_time();w=time.perf_counter()
 globals()[cfg['method']+'_'+phase](case,cfg,direction,dest)
 save(dest/(phase+'-child.json'),dict(cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w))

def execute(case,key,direction,phase):
 import psutil
 dest=folder(case,key,direction);dest.mkdir(parents=True,exist_ok=True);record=dest/(phase+'-execution.json')
 if record.exists():return read(record) # Complete configuration records only; never resume search timing.
 start=time.perf_counter();usage=resource.getrusage(resource.RUSAGE_CHILDREN);peak=0.;status='passed'
 with (dest/(phase+'.log')).open('w') as log:
  proc=subprocess.Popen([sys.executable,__file__,'child',str(case),key,direction,phase],env=ENV,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
  while proc.poll() is None:
   try:
    process=psutil.Process(proc.pid);rss=process.memory_info().rss+sum(x.memory_info().rss for x in process.children(recursive=True));peak=max(peak,rss/2**20)
   except (psutil.NoSuchProcess,psutil.AccessDenied):pass
   if time.perf_counter()-start>300:status='timeout'
   elif peak>6144:status='memory_limit'
   if status!='passed':
    try:os.killpg(proc.pid,signal.SIGKILL)
    except ProcessLookupError:pass
    proc.wait();break
   time.sleep(.05)
  code=proc.wait()
 if status=='passed' and code!=0:status='failed'
 used=resource.getrusage(resource.RUSAGE_CHILDREN);cpu=(used.ru_utime+used.ru_stime)-(usage.ru_utime+usage.ru_stime)
 row=dict(status=status,exit_code=code,worker_cpu_seconds=cpu,elapsed_seconds=time.perf_counter()-start,peak_mib=peak,host=socket.gethostname(),affinity=sorted(os.sched_getaffinity(0)),job_id=os.environ.get('SLURM_JOB_ID'),watchdog_seconds=300,memory_limit_mib=6144)
 if phase=='search':
  if (dest/'search.json').exists():row['api_cpu_seconds']=read(dest/'search.json')['cpu_seconds'];row['api_complete']=True
  elif (dest/'search-start.json').exists():row['api_cpu_seconds']=max(0.,cpu-read(dest/'search-start.json')['cpu_baseline']);row['api_complete']=False
  else:row['api_cpu_seconds']=None;row['api_complete']=False
 save(record,row);return row

def run_case(case):
 assert os.environ.get('SLURM_JOB_PARTITION')=='cpu_short' or os.environ.get('AAM_LOCAL_SMOKE')=='1'
 if hasattr(os,'sched_getaffinity'):
  cores=sorted(os.sched_getaffinity(0));os.sched_setaffinity(0,{cores[0]})
 machine=dict(host=socket.gethostname(),platform=platform.platform(),python=sys.version,affinity=sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else [],threads={k:os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']},cpu_info=subprocess.check_output(['lscpu','-J'],text=True) if sys.platform=='linux' else platform.processor())
 save(ROOT/'results'/str(case)/'machine.json',machine)
 configs=configurations();random.Random(20260915+case).shuffle(configs)
 for cfg in configs:
  for direction in ['R_to_P','P_to_R']:
   for phase in ['search','score']:
    row=execute(case,cfg['key'],direction,phase)
    print(json.dumps(dict(case=case,config=cfg['key'],direction=direction,phase=phase,status=row['status'])),flush=True)
    if row['status']=='failed':raise RuntimeError('Code failure; inspect logs before proceeding')
 save(ROOT/'results'/str(case)/'done.json',dict(case=case,status='processed',configurations=[c['key'] for c in configs]))

def verify_sources():
 m=read(ROOT/'manifest.json')
 for rel,digest in m['source_sha256'].items():assert sha(ROOT/rel)==digest,rel
 assert sha(ROOT/'slap_adapter.py')==m['slap_adapter_sha256']
 assert sha(ROOT/'inputs.json')==m['inputs_sha256']
 print('Source hashes verified',flush=True)

if __name__=='__main__':
 if sys.argv[1]=='child':child(int(sys.argv[2]),sys.argv[3],sys.argv[4],sys.argv[5])
 elif sys.argv[1]=='case':verify_sources();run_case(int(sys.argv[2]))
 elif sys.argv[1]=='check':verify_sources()
