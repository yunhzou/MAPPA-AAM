"""Paired one-seed coordinate cap experiment; separate persisted AAM stages."""
from pathlib import Path
import os,sys,json,gzip,time,hashlib,random,subprocess,resource,signal,uuid,collections
from dataclasses import asdict
S=Path(__file__).resolve().parent
sys.path.insert(0,str(S/'engine/src'))
def read(p):return json.loads(Path(p).read_text())
def save(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(d)+'\n');q.replace(p)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def problem(case):
 from mappa import AAMProblem,MolecularEndpoint
 r=read(S/'inputs'/str(case)/'input.json')
 return AAMProblem(*(MolecularEndpoint(**r[k]) for k in ['reactant','product']),r.get('name',''))
def search(case,cap,out):
 from mappa import search_aam,AAMSearchConfig
 p=problem(case);cfg=AAMSearchConfig(seed_count=1,branch_limit=cap,iso_tolerance=1.0,graph_floor=.2,cut_floor=.2,sweep_cuts=True,random_seed=42,event_threshold=.5,metal_event_threshold=.3)
 c=time.process_time();w=time.perf_counter()
 a=search_aam(p,cfg,execution='reused_native',workers=1,intermediate_dir=out/'cuts',archive_format='checkpoint',resume=False)
 full=sum(len(a.graph.states[t].mapping)==p.source_atom_count for t in a.graph.terminals)
 save(out/'search.json',dict(config=asdict(cfg),cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w,full_terminals=full,terminals=len(a.graph.terminals),capped=a.graph.capped,archive_sha256=sha(out/'cuts/aam.pkl.gz')))
def competition(case,cap,out):
 from mappa.artifacts import read_aam_checkpoint,write_aam_checkpoint
 from mappa.competition import compete_fragments,CompetitionConfig
 c=time.process_time();w=time.perf_counter()
 a=read_aam_checkpoint(out/'cuts/aam.pkl.gz');assert a.config.branch_limit==cap
 cfg=CompetitionConfig(operation_budget=128,seconds=270,parent_limit=8,depth_limit=2,queue_limit=512,dependent_component_limit=8)
 r=compete_fragments(a,cfg)
 (out/'repairs').mkdir(exist_ok=True)
 for i,a in enumerate(r.repairs):write_aam_checkpoint(a,out/f'repairs/repair{i:05d}.pkl.gz')
 save(out/'competition.json',dict(config=asdict(cfg),repairs=len(r.repairs),counts=r.counts,pending=r.pending,offers=r.offers,cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w))
def decode(case,cap,out,passno):
 from mappa.artifacts import read_aam_checkpoint
 from mappa.final_branches import FinalBranchCatalogue
 from mappa.event_patterns import SignedEventIndex,extract_path_events
 c=time.process_time();w=time.perf_counter();deadline=w+275
 a=read_aam_checkpoint(out/'cuts/aam.pkl.gz');p=a.problem;idx=SignedEventIndex(p,threshold=.5,metal_threshold=.3)
 cat=FinalBranchCatalogue(p).add_aam(a,'baseline');del a
 for f in sorted((out/'repairs').glob('*.pkl.gz')):cat.add_aam(read_aam_checkpoint(f),f.name)
 cat.rebuild_all_symmetries()
 K=read(S/'controls.json')['per_case'][case]['window'];repmin=None
 for start in range(0,len(cat.families),64):
  batch=cat.families[start:start+64]
  for f in batch:f.validate_representative(p)
  vs=[[dict(f.mapping)[i] for i in range(idx.n)] for f in batch]
  if vs:
   value=int(min(idx.counts(vs).sum(axis=1)));repmin=value if repmin is None else min(repmin,value)
 catalogue_identity=hashlib.sha256(json.dumps([f.to_record() for f in cat.families],sort_keys=True).encode()).hexdigest()
 identity=dict(catalogue_sha256=catalogue_identity,window=K,input_sha256=sha(S/'inputs'/str(case)/'input.json'),source_manifest_sha256=sha(S/'manifest.json'))
 if (out/'decode-identity.json').exists():assert read(out/'decode-identity.json')==identity
 else:save(out/'decode-identity.json',identity)
 journal=out/'families.jsonl';done=set();patterns={};reasons=collections.Counter()
 if journal.exists():
  lines=journal.read_bytes().splitlines(keepends=True)
  if lines and not lines[-1].endswith(b'\n'):lines=lines[:-1];journal.write_bytes(b''.join(lines))
  for line in lines:
   r=json.loads(line)
   for item in r['patterns']:
    check=idx.describe(item['mapping']);assert check['id']==item['id'];patterns[item['id']]=item
   if r['complete']:done.add(r['family'])
   reasons[r['reason']]+=1
 edges={floor:tuple((i,j) for i in range(idx.n) for j in range(i+1,idx.n) if p.reactant.wbo[i,j]>=floor) for floor in {f.policy[0] for f in cat.families}}
 def snapshot():
  minimum=min((v['total'] for v in patterns.values()),default=None)
  complete=len(done)==len(cat.families)
  save(out/'decoded.json',dict(case=case,cap=cap,window=K,complete_saved_full_window=complete,families=len(cat.families),completed_families=len(done),branches=cat.branch_count,paths=cat.path_count,incomplete_paths=cat.incomplete_paths,patterns=patterns,minimum=minimum,minimum_proven=complete and minimum is not None,minimum_ids=sorted(k for k,v in patterns.items() if v['total']==minimum),representative_upper_bound=repmin,minimum_lower_bound=K+1 if complete and not patterns and cat.families else None,no_full_families=not cat.families,reasons=dict(reasons),last_pass=passno,last_pass_cpu_seconds=time.process_time()-c,last_pass_wall_seconds=time.perf_counter()-w))
 snapshot();last=time.perf_counter()
 with journal.open('a',buffering=1) as stream:
  for i,f in enumerate(cat.families):
   if i in done:continue
   remaining=deadline-time.perf_counter()
   if remaining<=0:break
   r=extract_path_events(f.as_path(p,edges[f.policy[0]]),p,idx,max_events=K,seconds=min([20,60,180,260][min(passno,3)],remaining),max_patterns=None)
   for item in r['patterns']:
    assert idx.describe(item['mapping'])['id']==item['id'];patterns.setdefault(item['id'],item)
   stream.write(json.dumps(dict(family=i,**r))+'\n')
   if r['complete']:done.add(i)
   reasons[r['reason']]+=1
   if time.perf_counter()-last>3:snapshot();last=time.perf_counter()
 snapshot()
def execute(case,cap,out,phase,passno=0):
 import psutil
 label=phase if phase!='decode' else f'decode{passno}'
 c=resource.getrusage(resource.RUSAGE_CHILDREN);w=time.perf_counter();peak=0;status='passed'
 with (out/f'{label}.log').open('w') as log:
  q=subprocess.Popen([sys.executable,__file__,'child',str(case),str(cap),str(out),phase,str(passno)],stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
  while q.poll() is None:
   try:
    proc=psutil.Process(q.pid);peak=max(peak,(proc.memory_info().rss+sum(x.memory_info().rss for x in proc.children(recursive=True)))/2**20)
   except (psutil.NoSuchProcess,psutil.AccessDenied):pass
   if time.perf_counter()-w>=300:status='timeout'
   elif peak>=6144:status='memory_limit'
   if status!='passed':
    try:os.killpg(q.pid,signal.SIGKILL)
    except ProcessLookupError:pass
    break
   time.sleep(.1)
  code=q.wait()
 if code and status=='passed':status='failed'
 used=resource.getrusage(resource.RUSAGE_CHILDREN);cpu=used.ru_utime+used.ru_stime-c.ru_utime-c.ru_stime
 r=dict(phase=phase,passno=passno,status=status,exit_code=code,worker_cpu_seconds=cpu,wall_seconds=time.perf_counter()-w,peak_mib=peak)
 save(out/f'{label}-execution.json',r);return r

def run_case(case):
 if hasattr(os,'sched_getaffinity'):os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
 root=S/'results'/str(case);root.mkdir(parents=True,exist_ok=True)
 if (root/'selected.json').exists():raise RuntimeError('Selected attempt already exists; explicit recovery required')
 attempt=root/('attempt-'+uuid.uuid4().hex[:12]);attempt.mkdir()
 cpuinfo=subprocess.check_output(['lscpu','-J'],text=True) if sys.platform=='linux' else '{}'
 save(attempt/'machine.json',dict(cpu_info=json.loads(cpuinfo),threads={k:os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']},affinity=sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else [],job_id=os.environ.get('SLURM_JOB_ID')))
 caps=[100,2000];random.Random(20260915+case).shuffle(caps);summaries={}
 for cap in caps:
  out=attempt/str(cap);out.mkdir()
  stages=[]
  for phase in ['search','competition']:
   r=execute(case,cap,out,phase);stages.append(r)
   if r['status']!='passed':break
  if all(r['status']=='passed' for r in stages) and len(stages)==2:
   for n in range(4):
    r=execute(case,cap,out,'decode',n);stages.append(r)
    if (out/'decoded.json').exists() and read(out/'decoded.json')['complete_saved_full_window']:break
    if r['status']=='failed':break
  summaries[str(cap)]=dict(stages=stages,decoded=read(out/'decoded.json') if (out/'decoded.json').exists() else None)
  save(attempt/'progress.json',dict(case=case,results=summaries))
  print(json.dumps(dict(case=case,cap=cap,statuses=[r['status'] for r in stages],complete=bool(summaries[str(cap)]['decoded'] and summaries[str(cap)]['decoded']['complete_saved_full_window']))),flush=True)
 save(attempt/'done.json',dict(case=case,cap_order=caps,source_manifest_sha256=sha(S/'manifest.json')))
 save(root/'selected.json',dict(attempt=attempt.name))

def verify():
 m=read(S/'manifest.json')
 for n,h in m['source_sha256'].items():assert sha(S/n)==h,n
 for i,h in m['input_sha256'].items():assert sha(S/'inputs'/i/'input.json')==h,i
 assert sha(S/'controls.json')==m['controls_sha256']
if __name__=='__main__':
 if sys.argv[1]=='child':
  case,cap,out,phase,n=int(sys.argv[2]),int(sys.argv[3]),Path(sys.argv[4]),sys.argv[5],int(sys.argv[6])
  if phase=='decode':decode(case,cap,out,n)
  else:globals()[phase](case,cap,out)
 elif sys.argv[1]=='case':verify();run_case(int(sys.argv[2]))
 elif sys.argv[1]=='verify':verify()
