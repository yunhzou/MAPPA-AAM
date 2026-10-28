"""Portable rerun of the frozen SLAP Golden edge-sweep protocol."""
import sys,os,json,time,collections
from pathlib import Path
from types import SimpleNamespace
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parent;ENGINE=ROOT/'engine';RUN=ROOT/'slap-golden'
sys.path[:0]=[str(ROOT/'slap-upstream/src'),str(ENGINE/'src'),str(ENGINE/'bench')]
from resources import read,save,sha,MemoryGuard,bounded_process

def prepare():
 from rdkit import Chem
 rows=[json.loads(line) for line in (ROOT/'golden-data/audit.jsonl').read_text().splitlines()]
 save(RUN/'inputs.json',[dict(index=r['index'],input_reaction=r['input_reaction']) for r in rows])
 tasks=[]
 for row in rows:
  mols=[Chem.AddHs(Chem.MolFromSmiles(s)) for s in row['input_reaction'].split('>>')]
  for reverse in (False,True):
   for mode in ('binary','weighted'):tasks.append(dict(slot=len(tasks),case=row['index'],direction='P_to_R' if reverse else 'R_to_P',mode=mode,expected_cuts=1+mols[int(reverse)].GetNumBonds()))
 save(RUN/'tasks.json',tasks)
 save(RUN/'manifest.json',dict(upstream_commit='ea248fd9494f52f4865193e87a98cc92c62b5f9e',audit_sha256=sha(ROOT/'golden-data/audit.jsonl'),source_sha256={str(p.relative_to(ROOT/'slap-upstream')):sha(p) for p in (ROOT/'slap-upstream/src').rglob('*.py')},adapter_sha256=sha(Path(__file__)),cases=1851,workers_max=8,watchdog=300,scope='Uncut plus every source edge; both directions and binary/weighted modes; explicit H; heavy symmetry breaking.'))

def child(slot):
 import slap_edge_sweep as m
 # Linux affinity was only machine metadata in the archived adapter.
 if not hasattr(os,'sched_getaffinity'):os.sched_getaffinity=lambda pid:set(range(os.cpu_count()))
 m.case(SimpleNamespace(run=RUN,slot=slot))

def evaluate(case):
 from golden_competitors import signatures
 from slap_edge_sweep import records,folder
 ref=json.loads((ROOT/'golden-data/audit.jsonl').read_text().splitlines()[case]);expected_ep,expected=signatures(ref['mapped_reaction'])
 attempts=[];found=False;complete=True;cpu=0.;errors=0;invalid=0;classes=set();witness=None
 for task in read(RUN/'tasks.json')[4*case:4*case+4]:
  rows=records(folder(RUN,task)/'records.jsonl');hits=[]
  complete &= len(rows)==task['expected_cuts']
  for row in rows:
   cpu+=row['mapping_cpu']+row['graph_cpu']+row['export_cpu'];errors+=row['status']!='mapped'
   for ordinal,c in enumerate(row['candidates']):
    try:
     ep,sig=signatures(c['mapped_rxn']);assert ep==expected_ep
    except Exception:invalid+=1;continue
    classes.add(sig)
    if sig==expected:
     found=True;hits.append(dict(ordinal=row['ordinal'],candidate=ordinal));witness=c['mapped_rxn']
  attempts.append(dict(**task,completed_cuts=len(rows),hits=hits))
 result=dict(case=case,reference_recovery='recovered' if found else 'not_recovered' if complete and not errors and not invalid else 'unknown',complete=bool(complete and not errors and not invalid),workflow_cpu_excluding_io=cpu,mapping_errors=errors,invalid=invalid,unique_classes=len(classes),attempts=attempts,witness=witness)
 save(RUN/'evaluations'/f'{case}.json',result)

def campaign(cases,workers,label):
 tasks=[4*c+d for c in cases for d in range(4)];guard=MemoryGuard(3072,8192,6144)
 env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1');start=time.perf_counter()
 def run(slot):
  d=RUN/'execution'/str(slot);d.mkdir(parents=True,exist_ok=True)
  if (d/'status.json').exists():return read(d/'status.json')
  with (d/'run.log').open('w') as log:state,peak=bounded_process([sys.executable,__file__,'child',str(slot)],env,log,300,guard)
  row=dict(slot=slot,status=state,peak_mib=peak);save(d/'status.json',row);return row
 rows=[]
 with ThreadPoolExecutor(max_workers=workers) as pool:
  for future in as_completed([pool.submit(run,s) for s in tasks]):
   rows.append(future.result());save(RUN/f'{label}-progress.json',dict(completed=len(rows),total=len(tasks),wall_seconds=time.perf_counter()-start,peak_mib=guard.peak_total_kib/1024))
 # Evaluation in separately guarded child processes.
 def score(c):
  d=RUN/'evaluations';d.mkdir(exist_ok=True)
  with (d/f'{c}.log').open('w') as log:state,peak=bounded_process([sys.executable,__file__,'evaluate',str(c)],env,log,300,guard)
  return dict(case=c,status=state,peak_mib=peak)
 with ThreadPoolExecutor(max_workers=workers) as pool:scores=list(pool.map(score,cases))
 save(RUN/f'{label}-execution.json',dict(rows=rows,scores=scores,wall_seconds=time.perf_counter()-start,peak_mib=guard.peak_total_kib/1024))
if __name__=='__main__':
 if sys.argv[1]=='prepare':prepare()
 elif sys.argv[1]=='child':child(int(sys.argv[2]))
 elif sys.argv[1]=='evaluate':evaluate(int(sys.argv[2]))
 else:campaign(read(ROOT/'inputs.json')['pilot'] if sys.argv[2]=='pilot' else list(range(1851)),int(sys.argv[3]),sys.argv[2])
