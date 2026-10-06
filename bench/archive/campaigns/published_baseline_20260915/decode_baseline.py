"""Decode saved baseline searches only. Never launch AAM or competition."""
from pathlib import Path
import sys,os,json,gzip,time,subprocess,signal,concurrent.futures,hashlib
S=Path(__file__).resolve().parent
R=S.parents[1]
W=Path(os.environ['MAPPA_ARCHIVE_ROOT'])
sys.path.insert(0,str(R/'src'))
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def save(p,d):p.write_text(json.dumps(d,separators=(',',':'))+'\n')
def child(case):
 from mappa.artifacts import read_aam_checkpoint
 from mappa.final_branches import FinalBranchCatalogue
 from mappa.event_patterns import SignedEventIndex,extract_path_events
 folder=S/'decoded'/str(case);folder.mkdir(parents=True,exist_ok=True)
 start=time.perf_counter();cpu=time.process_time()
 path=W/f'runs/coordinate/seed1/case{case}/R_to_P/cuts/aam.pkl.gz'
 a=read_aam_checkpoint(path);assert a.config.branch_limit==2000 and a.config.seed_count==1
 cat=FinalBranchCatalogue(a.problem).add_aam(a,'baseline');cat.rebuild_all_symmetries()
 idx=SignedEventIndex(a.problem)
 saved=read(R/'reports/current_validation_20260912/coordinate-decoded.json.gz')[case]
 window=saved['window'];patterns={};done=0
 def checkpoint():save(folder/'result.json',dict(case=case,window=window,paths=cat.path_count,branches=cat.branch_count,families=len(cat.families),complete=done==len(cat.families),completed_families=done,patterns=patterns,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),branch_limit=2000,competition=False))
 checkpoint()
 with (folder/'families.jsonl').open('w',buffering=1) as f:
  for i,fam in enumerate(cat.families):
   fam.validate_representative(a.problem)
   result=extract_path_events(fam.as_path(a.problem),a.problem,idx,max_events=window,seconds=290,max_patterns=None)
   for p in result['patterns']:patterns.setdefault(p['id'],dict(p,family=i))
   done+=result['complete'];f.write(json.dumps(dict(family=i,**result))+'\n')
   checkpoint()
 print(json.dumps(dict(case=case,complete=done==len(cat.families),patterns=len(patterns),cpu=time.process_time()-cpu)),flush=True)
def launch(case):
 import psutil
 folder=S/'decoded'/str(case);folder.mkdir(parents=True,exist_ok=True)
 env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
 start=time.perf_counter();peak=0;status=None
 with (folder/'run.log').open('w') as f:
  p=subprocess.Popen([sys.executable,__file__,'child',str(case)],env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  save(folder/'active.json',dict(pid=p.pid));q=psutil.Process(p.pid)
  try:
   while p.poll() is None:
    try:peak=max(peak,q.memory_info().rss/2**20)
    except psutil.NoSuchProcess:break
    if time.perf_counter()-start>300 or peak>3072:
     status='watchdog';os.killpg(p.pid,signal.SIGKILL);break
    time.sleep(.2)
   code=p.wait()
  finally:
   if p.poll() is None:os.killpg(p.pid,signal.SIGKILL);p.wait()
 result=dict(case=case,status=status or ('finished' if code==0 else 'error'),peak_mib=peak,wall_seconds=time.perf_counter()-start)
 save(folder/'execution.json',result);print(json.dumps(result),flush=True);return result
if __name__=='__main__':
 if sys.argv[1]=='child':child(int(sys.argv[2]))
 else:
  start=time.perf_counter()
  with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:results=list(pool.map(launch,range(140)))
  save(S/'execution.json',dict(results=results,wall_seconds=time.perf_counter()-start,workers=3))
