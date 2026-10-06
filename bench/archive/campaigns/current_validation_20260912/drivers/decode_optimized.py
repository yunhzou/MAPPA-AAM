"""Fresh complete event-window decoding of final unordered family unions."""
import os,sys,time,json,gzip,gc,collections
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parent;sys.path[:0]=[str(ROOT/'optimized/src'),str(ROOT/'optimized/bench')]
from run_event_campaign import MemoryGuard,bounded_process,save,read,sha
REPO=Path('/Users/yunhengz/Desktop/AAM Writing')
def child(case):
 from mappa.artifacts import read_aam_checkpoint
 from mappa.final_branches import FinalBranchCatalogue
 from mappa.event_patterns import SignedEventIndex,extract_path_events
 from validate_event_examples import compare_saved_slap
 start=time.perf_counter();cpu=time.process_time();folder=ROOT/'decoded-optimized'/f'case{case}';folder.mkdir(parents=True,exist_ok=True)
 base=ROOT/'runs/coordinate/seed1'/f'case{case}'/'R_to_P/cuts/aam.pkl.gz'
 assert read(base.parent.parent/'search.json')['archive_sha256']==sha(base)
 assert read(ROOT/'competition'/f'case{case}'/'execution.json')['status']=='passed'
 archives=[base,*sorted((ROOT/'competition'/f'case{case}').glob('takeover*.pkl.gz'))]
 identity=dict(archives={str(p):sha(p) for p in archives},engine_sha256=sha(ROOT/'optimized-engine.json'),driver_sha256=sha(Path(__file__)))
 if (folder/'identity.json').exists():assert read(folder/'identity.json')==identity
 else:save(folder/'identity.json',identity)
 a=read_aam_checkpoint(base);problem=a.problem;idx=SignedEventIndex(problem);cat=FinalBranchCatalogue(problem);cat.add_graph(a.graph,str(base));del a;gc.collect()
 for archive in archives[1:]:
  a=read_aam_checkpoint(archive);cat.add_graph(a.graph,str(archive));del a
 cat.rebuild_all_symmetries()
 metrics=dict(case=case,paths=cat.path_count,branches=cat.branch_count,families=len(cat.families),symmetries=len(cat.symmetries),catalogue_seconds=time.perf_counter()-start)
 # Reference witnesses are loaded only after search and catalogue construction.
 slap={}
 for method,filename,key in [('native_slap','reports/holdout_cap1000_seed1_20260910/slap_scoring.json.gz',str(case)),('slap_sweep','reports/holdout_slap_xyz_sweep_20260910/scored_candidates.json.gz',f'{case}/scored')]:
  rows=read(REPO/filename)[key]['slap'];best=None;patterns={}
  for offset in range(0,len(rows),64):
   vectors=[[dict(row['mapping'])[i] for i in range(idx.n)] for row in rows[offset:offset+64]]
   for vector,count in zip(vectors,idx.counts(vectors).sum(axis=1)):
    total=int(count)
    if best is None or total<best:best=total;patterns={}
    if total==best:
     p=idx.describe(vector);patterns.setdefault(p['id'],p)
  slap[method]=dict(minimum=best,minimum_ids=sorted(patterns),patterns=patterns,saved_witnesses=len(rows))
 window=max(v['minimum'] for v in slap.values())+1
 source_edges={floor:tuple((i,j) for i in range(idx.n) for j in range(i+1,idx.n) if problem.reactant.wbo[i,j]>=floor) for floor in {f.policy[0] for f in cat.families}}
 metrics.update(window=window,event_policy=idx.policy,slap=slap)
 journal=folder/'families.jsonl';done=set();patterns={};reasons=collections.Counter();attempts=[]
 if journal.exists():
  lines=journal.read_bytes().splitlines(keepends=True)
  # Remove only an uncommitted final journal line after a watchdog.
  if lines and not lines[-1].endswith(b'\n'):journal.write_bytes(b''.join(lines[:-1]));lines=lines[:-1]
  for line in lines:
   row=json.loads(line)
   for p in row['patterns']:
    check=idx.describe(p['mapping']);assert check['id']==p['id'];patterns.setdefault(p['id'],p)
   if row['complete']:done.add(row['family'])
   reasons[row['reason']]+=1
 initial=len(done);deadline=start+275;last=time.perf_counter()
 def checkpoint():
  save(folder/'result.json',dict(**metrics,complete_saved_window=len(done)==len(cat.families),completed_families=len(done),patterns=patterns,reasons=dict(reasons),comparisons=compare_saved_slap(slap,patterns,window,len(done)==len(cat.families)),wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,previously_complete=initial))
 checkpoint()
 with journal.open('a',buffering=1) as stream:
  for fid,family in enumerate(cat.families):
   if fid in done:continue
   if time.perf_counter()>deadline:break
   family.validate_representative(problem)
   result=extract_path_events(family.as_path(problem,source_edges[family.policy[0]]),problem,idx,max_events=window,seconds=min(20,max(.01,deadline-time.perf_counter())),max_patterns=None)
   for p in result['patterns']:patterns.setdefault(p['id'],dict(p,family=fid))
   stream.write(json.dumps(dict(family=fid,**result),separators=(',',':'))+'\n')
   if result['complete']:done.add(fid)
   reasons[result['reason']]+=1
   if time.perf_counter()-last>3:checkpoint();last=time.perf_counter()
 checkpoint();print(json.dumps(dict(case=case,complete=len(done)==len(cat.families),done=len(done),families=len(cat.families),patterns=len(patterns),seconds=time.perf_counter()-start)),flush=True)
def main(cases,workers=3):
 guard=MemoryGuard(3072,8192,6144);env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1');start=time.perf_counter()
 def run(c):
  folder=ROOT/'decoded-optimized'/f'case{c}';folder.mkdir(parents=True,exist_ok=True)
  if (folder/'result.json').exists() and read(folder/'result.json')['complete_saved_window']:return dict(case=c,status='previously_complete')
  with (folder/'run.log').open('a') as log:status,peak=bounded_process([sys.executable,__file__,str(c)],env,log,300,guard)
  row=dict(case=c,status=status,peak_mib=peak);save(folder/'execution.json',row);print(json.dumps(row),flush=True);return row
 with ThreadPoolExecutor(max_workers=workers) as pool:rows=list(pool.map(run,cases))
 save(ROOT/'decoded-optimized/execution.json',dict(rows=rows,wall_seconds=time.perf_counter()-start,peak_mib=guard.peak_total_kib/1024))
if __name__=='__main__':
 if sys.argv[1]=='all':main(list(range(140)))
 elif sys.argv[1]=='campaign':main(list(map(int,sys.argv[2].split(','))))
 else:child(int(sys.argv[1]))
