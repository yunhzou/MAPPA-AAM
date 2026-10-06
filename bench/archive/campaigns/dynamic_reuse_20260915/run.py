"""Bounded failure-case competition ablation; references enter only aggregation."""
from pathlib import Path
import sys,os,json,gzip,time,resource,subprocess
from concurrent.futures import ThreadPoolExecutor,as_completed
PACKAGE=Path(__file__).resolve().parent
S=Path(os.environ.get('MAPPA_EXPERIMENT_WORK',PACKAGE/'work')).resolve()
INPUTS=PACKAGE.parent/'unrestricted_competition_20260915/inputs'
sys.path.insert(0,str(S/'engine/src'))

def save(p,v):
 p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(v,indent=2)+'\n');tmp.replace(p)

def worker(case,mode,stage):
 from mappa.artifacts import read_aam_checkpoint,write_aam_checkpoint
 from mappa.final_branches import FinalBranchCatalogue
 from mappa.event_patterns import SignedEventIndex,extract_path_events
 from mappa.competition import compete_fragments,CompetitionConfig
 out=S/'runs'/str(case)/mode;out.mkdir(parents=True,exist_ok=True)
 start=time.perf_counter();cpu=time.process_time()
 src=INPUTS/f'case{case}.pkl.gz'
 a=read_aam_checkpoint(src);assert a.config.seed_count==1 and a.config.branch_limit==2000
 if stage=='search':
  c=compete_fragments(a,CompetitionConfig(proposal_mode='prefix' if mode=='prefix512' else mode, operation_budget=512 if mode=='prefix512' else 128))
  search_cpu=time.process_time()-cpu;search_wall=time.perf_counter()-start
  for i,r in enumerate(c.repairs):write_aam_checkpoint(r,out/f'takeover{i}.pkl.gz')
  save(out/'search.json',dict(case=case,mode=mode,cpu_seconds=search_cpu,wall_seconds=search_wall,counts=c.counts,offers=c.offers,pending=c.pending,repairs=len(c.repairs),peak_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2))
  print(json.dumps(dict(stage=stage,case=case,mode=mode,cpu=search_cpu,repairs=len(c.repairs),max_eaten=max([len(o['eaten']) for o in c.offers] or [0]))),flush=True);return
 cat=FinalBranchCatalogue(a.problem).add_aam(a,'baseline')
 if mode!='none':
  info=json.loads((out/'search.json').read_text())
  for i in range(info['repairs']):cat.add_aam(read_aam_checkpoint(out/f'takeover{i}.pkl.gz'),f'takeover{i}')
 cat.rebuild_all_symmetries();catalogue_cpu=time.process_time()-cpu
 (out/'catalogue.json.gz').write_bytes(gzip.compress(json.dumps(cat.to_record()).encode(),mtime=0))
 # Fixed pre-existing windows, not derived from outcomes in this experiment.
 window=8 if case==64 else 5
 index=SignedEventIndex(a.problem);patterns={};certificates=[];bad=[]
 for i,f in enumerate(cat.families):
  if time.perf_counter()-start>270:break
  try:f.validate_representative(a.problem)
  except AssertionError:bad.append(i);continue
  result=extract_path_events(f.as_path(a.problem),a.problem,index,max_events=window,seconds=min(10,270-(time.perf_counter()-start)),max_patterns=None)
  certificates.append(dict(family=i,complete=result['complete'],reason=result['reason'],solver_queries=result['solver_queries']))
  for p in result['patterns']:
   assert index.describe(p['mapping'])['id']==p['id']
   patterns.setdefault(p['id'],dict(p,family=i))
  if i%50==0:save(out/'decode-progress.json',dict(families_done=i+1,families=len(cat.families),patterns=len(patterns)))
 complete=len(certificates)==len(cat.families) and all(c['complete'] for c in certificates) and not bad
 save(out/'decode.json',dict(case=case,mode=mode,window=window,complete=complete,paths=cat.path_count,branches=cat.branch_count,families=len(cat.families),patterns=patterns,certificates=certificates,invalid_families=bad,catalogue_cpu=catalogue_cpu,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start,peak_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2))
 print(json.dumps(dict(stage=stage,case=case,mode=mode,complete=complete,branches=cat.branch_count,families=len(cat.families),patterns=len(patterns),cpu=time.process_time()-cpu)),flush=True)

def launch(case,mode,stage):
 out=S/'runs'/str(case)/mode;out.mkdir(parents=True,exist_ok=True)
 env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
 start=time.perf_counter();peak=0;status=None
 import psutil
 with (out/f'{stage}.log').open('w') as log:
  p=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'worker',str(case),mode,stage],env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
  proc=psutil.Process(p.pid)
  while p.poll() is None:
   try:rss=proc.memory_info().rss/1024**2;peak=max(peak,rss)
   except psutil.NoSuchProcess:break
   if time.perf_counter()-start>300 or peak>3000:
    status='watchdog' if peak<=3000 else 'memory_limit';p.kill();p.wait();break
   time.sleep(.25)
  code=p.wait()
 status=status or ('passed' if code==0 else 'failed')
 record=dict(case=case,mode=mode,stage=stage,status=status,peak_mib=peak,wall_seconds=time.perf_counter()-start)
 save(out/f'{stage}-execution.json',record);print(json.dumps(record),flush=True);return record

if __name__=='__main__':
 if len(sys.argv)>1 and sys.argv[1]=='worker':worker(int(sys.argv[2]),sys.argv[3],sys.argv[4])
 else:
  cases=[6,11,59,64,101];modes=['prefix'];rows=[]
  with ThreadPoolExecutor(max_workers=2) as pool:
   futures=[pool.submit(launch,c,m,'search') for c in cases for m in modes]
   for f in as_completed(futures):rows.append(f.result())
  with ThreadPoolExecutor(max_workers=2) as pool:
   futures=[pool.submit(launch,c,m,'decode') for c in cases for m in modes if m=='none' or json.loads((S/'runs'/str(c)/m/'search-execution.json').read_text())['status']=='passed']
   for f in as_completed(futures):rows.append(f.result())
  save(S/'execution.json',rows)
