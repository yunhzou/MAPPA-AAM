"""Cap effects on minimum events, pattern coverage, completeness and stage CPU."""
from pathlib import Path
import json,collections,statistics
import run
S=Path(__file__).resolve().parent

def aggregate():
 controls=run.read(S/'controls.json')['per_case'];cases=[];models=collections.Counter()
 for case in range(140):
  folder=S/'results'/str(case)
  if not (folder/'selected.json').exists():continue
  attempt=folder/run.read(folder/'selected.json')['attempt']
  m=run.read(attempt/'machine.json');model=next((r['data'] for r in m['cpu_info'].get('lscpu',[]) if r['field']=='Model name:'),'unknown');models[model]+=1
  assert all(v=='1' for v in m['threads'].values())
  rec=dict(case=case,caps={})
  for cap in [100,2000]:
   out=attempt/str(cap);d=run.read(out/'decoded.json') if (out/'decoded.json').exists() else None
   times={}
   for phase in ['search','competition','decode']:
    paths=sorted(out.glob('decode*-execution.json')) if phase=='decode' else [out/f'{phase}-execution.json']
    es=[run.read(p) for p in paths if p.exists()]
    times[phase]=dict(worker_cpu_seconds=sum(x['worker_cpu_seconds'] for x in es),wall_seconds=sum(x['wall_seconds'] for x in es),statuses=[x['status'] for x in es])
   rec['caps'][str(cap)]=dict(decoded=d,timing=times)
  a,b=[rec['caps'][str(cap)]['decoded'] for cap in [100,2000]]
  rec['both_resolved']=bool(a and b and a['minimum_proven'] and b['minimum_proven'])
  if rec['both_resolved']:
   rec['minimum_delta_cap100_minus_cap2000']=a['minimum']-b['minimum']
   rec['cap2000_minimum_patterns_missing_at_cap100']=sorted(set(b['minimum_ids'])-set(a['patterns']))
   rec['cap2000_window_patterns_missing_at_cap100']=sorted(set(b['patterns'])-set(a['patterns']))
  if b:
   rec['saved_cap2000_minimum']=controls[case]['mappa_minimum']
   rec['fresh_cap2000_minimum']=b['minimum']
   rec['saved_minimum_patterns_missing_from_fresh2000']=sorted(set(controls[case]['mappa_minimum_ids'])-set(b['patterns']))
  cases.append(rec)
 summary=dict(status='all_cases_processed' if len(cases)==140 else 'in_progress',processed=len(cases),cpu_models=dict(models),paired_minima_resolved=sum(r['both_resolved'] for r in cases),cap100_higher_minimum=[r['case'] for r in cases if r.get('minimum_delta_cap100_minus_cap2000',0)>0],cap100_lower_minimum=[r['case'] for r in cases if r.get('minimum_delta_cap100_minus_cap2000',0)<0],same_minimum=sum(r.get('minimum_delta_cap100_minus_cap2000')==0 for r in cases),cap100_missing_minimum_patterns=[r['case'] for r in cases if r.get('cap2000_minimum_patterns_missing_at_cap100')],unresolved=[r['case'] for r in cases if not r['both_resolved']],per_case=cases)
 run.save(S/'comparison.json',summary)
 compact={k:v for k,v in summary.items() if k!='per_case'};run.save(S/'progress.json',compact);print(json.dumps(compact),flush=True)
 return summary
if __name__=='__main__':aggregate()
