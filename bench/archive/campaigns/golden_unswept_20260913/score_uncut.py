"""Score only saved uncut graphs from the final one-seed campaign; no search."""
from pathlib import Path
import sys,json,hashlib,time,signal,gc,os
from concurrent.futures import ProcessPoolExecutor,as_completed
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent/'current-validation'
sys.path[:0]=[str(ROOT),str(ROOT/'engine/src'),str(ROOT/'engine/bench')]
OUT=HERE/'uncut';OUT.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def expired(*args):raise TimeoutError('300-second uncut evaluation watchdog')
def score(c):
 from score_checkpoints import inputs,finalized
 from mappa.artifacts import read_raw_cut
 from mappa.aam import checkpoint_manifest
 from mappa.domain import AAMResult,AAMSearchMetrics
 from golden_evaluation import evaluate_planned
 signal.signal(signal.SIGALRM,expired)
 dest=OUT/f'{c}.json'
 if dest.exists():return json.loads(dest.read_text())
 rows=[]
 for direction in ['R_to_P','P_to_R']:
  start=time.perf_counter();cpu=time.process_time();signal.alarm(300)
  row=dict(case=c,direction=direction)
  try:
   plan,config,folder=inputs(1,c,direction)
   assert json.loads((folder/'cuts/manifest.json').read_text())==checkpoint_manifest(plan.problem,config)
   path=folder/'cuts/cut_00000.raw.pkl.gz'
   graph=read_raw_cut(path)
   assert all(not ctx.cuts for ctx in graph.contexts), 'Non-uncut context'
   graph=finalized(graph,plan.problem,config)
   reference=json.loads((ROOT/f'golden-inputs/{c}/reference.json').read_text())
   aam=AAMResult(plan.problem,config,graph,AAMSearchMetrics.from_record({},0))
   result=evaluate_planned(aam,plan,reference['features'],reference['mapping'],seconds=220,query_timeout_ms=5000)
   row.update(result=result,checkpoint_sha256=sha(path),reference_sha256=sha(ROOT/f'golden-inputs/{c}/reference.json'),status='passed')
   del graph,aam,result
  except Exception as e:row.update(status='error',error=repr(e))
  finally:signal.alarm(0)
  row.update(cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start);rows.append(row);gc.collect()
 outcomes=[r.get('result',{}).get('reference_recovery','unknown') for r in rows]
 outcome='recovered' if 'recovered' in outcomes else 'not_recovered' if outcomes==['not_recovered']*2 else 'unknown'
 record=dict(case=c,outcome=outcome,directions=rows)
 dest.write_text(json.dumps(record,indent=2)+'\n');return record
if __name__=='__main__':
 cases=[int(s) for s in sys.argv[1:]] or list(range(1851));start=time.perf_counter()
 with ProcessPoolExecutor(max_workers=4) as pool:
  for i,f in enumerate(as_completed([pool.submit(score,c) for c in cases]),1):
   r=f.result()
   if len(cases)<10 or i%100==0:print(json.dumps(dict(completed=i,total=len(cases),last_case=r['case'],last_outcome=r['outcome'],elapsed=time.perf_counter()-start)),flush=True)
 print(json.dumps(dict(done=len(cases),elapsed=time.perf_counter()-start)),flush=True)
