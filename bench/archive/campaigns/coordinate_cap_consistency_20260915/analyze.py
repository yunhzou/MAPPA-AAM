"""Compare every certified event window, including certified empty searches."""
import json,statistics,collections,hashlib
from pathlib import Path
S=Path(__file__).resolve().parent

def stats(v):
 return dict(n=len(v),mean=statistics.mean(v),median=statistics.median(v),total=sum(v)) if v else dict(n=0)

def analyze():
 data=json.loads((S/'comparison.json').read_text());controls=json.loads((S/'controls.json').read_text())['per_case']
 out=dict(processed=data['processed'],complete_pairs=0,unresolved=[],same_minimum=0,higher_minimum=[],lower_minimum=[],no_full_mapping={'100':[],'2000':[]},minimum_patterns_lost=[],window_patterns_lost=[],window_patterns_added=[],fresh2000_disagrees_with_saved=[],per_case=[],timing={})
 for r in data['per_case']:
  case=r['case'];a,b=[r['caps'][str(c)]['decoded'] for c in [100,2000]]
  for c,d in [('100',a),('2000',b)]:
   if d and d['complete_saved_full_window'] and d['no_full_families']:out['no_full_mapping'][c].append(case)
  if not a or not b or not a['complete_saved_full_window'] or not b['complete_saved_full_window']:
   out['unresolved'].append(case);continue
  out['complete_pairs']+=1
  detail=dict(case=case,minimum100=a['minimum'],minimum2000=b['minimum'],families100=a['families'],families2000=b['families'],patterns100=len(a['patterns']),patterns2000=len(b['patterns']),minimum_patterns_lost=sorted(set(b['minimum_ids'])-set(a['patterns'])),window_patterns_lost=sorted(set(b['patterns'])-set(a['patterns'])),window_patterns_added=sorted(set(a['patterns'])-set(b['patterns'])))
  for k in ['minimum_patterns_lost','window_patterns_lost','window_patterns_added']:
   if detail[k]:out[k].append(dict(case=case,count=len(detail[k])))
  if a['minimum_proven'] and b['minimum_proven']:
   if a['minimum']==b['minimum']:out['same_minimum']+=1
   elif a['minimum']>b['minimum']:out['higher_minimum'].append(case)
   else:out['lower_minimum'].append(case)
  if b['minimum']!=controls[case]['mappa_minimum'] or set(controls[case]['mappa_minimum_ids'])-set(b['patterns']):out['fresh2000_disagrees_with_saved'].append(case)
  out['per_case'].append(detail)
 common=[r for r in data['per_case'] if all(r['caps'][str(c)]['decoded'] and r['caps'][str(c)]['decoded']['minimum_proven'] and all(x=='passed' for phase in r['caps'][str(c)]['timing'].values() for x in phase['statuses']) for c in [100,2000])]
 out['timing']['paired_cases']=[r['case'] for r in common]
 out['timing']['scope']='Same completed-mapping cohort, one pinned CPU thread, both caps serial in shuffled order; worker CPU includes imports and persistence. Decoder includes all passes. Search/competition internal stage CPU excludes imports.'
 for cap in ['100','2000']:
  out['timing'][cap]={phase:stats([r['caps'][cap]['timing'][phase]['worker_cpu_seconds'] for r in common]) for phase in ['search','competition','decode']}
  out['timing'][cap]['total_worker']=stats([sum(t['worker_cpu_seconds'] for t in r['caps'][cap]['timing'].values()) for r in common])
  for phase in ['search','competition']:
   out['timing'][cap][phase+'_internal']=stats([r['caps'][cap][phase]['cpu_seconds'] for r in common if r['caps'][cap].get(phase)])
 out['all_processed_worker_cpu']={cap:stats([sum(p['worker_cpu_seconds'] for p in row['caps'][cap]['timing'].values()) for row in data['per_case']]) for cap in ['100','2000']}
 out['scope']='Minimum/event-window coverage within saved full families, not exhaustive search or verified mapping accuracy. Empty searches are search losses, not unknown decoder results.'
 saved=json.loads((S/'saved-windows.json').read_text())['per_case']
 out['fresh2000_window_differences']=[]
 for row in data['per_case']:
  d=row['caps']['2000']['decoded'];old=saved[str(row['case'])]
  if d and d['complete_saved_full_window']:
   assert d['window']==old['window'] and old['complete']
   lost=set(old['pattern_ids'])-set(d['patterns']);added=set(d['patterns'])-set(old['pattern_ids'])
   if lost or added:out['fresh2000_window_differences'].append(dict(case=row['case'],lost=sorted(lost),added=sorted(added)))
 out['comparison_sha256']=hashlib.sha256((S/'comparison.json').read_bytes()).hexdigest()
 (S/'analysis.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k not in ['per_case','timing']}))
 print(json.dumps(out['timing']))
 return out
if __name__=='__main__':analyze()
