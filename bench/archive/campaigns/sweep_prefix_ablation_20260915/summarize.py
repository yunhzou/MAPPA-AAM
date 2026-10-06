from pathlib import Path
import json,collections,os
PACKAGE=Path(__file__).resolve().parent
S=Path(os.environ.get("MAPPA_EXPERIMENT_WORK",PACKAGE/"work")).resolve()
controls=json.loads((PACKAGE/'controls.json').read_text())['per_case'];cases=json.loads((PACKAGE/'manifest.json').read_text())['cases'];rows=[];paired=[]
for c in cases:
 values={}
 for sweep in [0,1]:
  f=S/f'runs/{c}/{sweep}/decode.json'
  if not f.exists():continue
  x=json.loads(f.read_text())
  for mode,v in x['arms'].items():
   targets=controls[c]['mappa_minimum_ids'];slap=controls[c]['comparators']['slap_sweep']['ids'];values[sweep,mode]=v
   search=json.loads((S/f'runs/{c}/{sweep}/search.json').read_text())
   comp=json.loads((S/f'runs/{c}/{sweep}/{mode}.json').read_text()) if mode!='none' else {}
   rows.append(dict(case=c,sweep=sweep,mode=mode,complete=v['complete'],minimum=v['minimum'],no_full=v['no_full_mapping'],above_window=v['minimum_lower_bound'],patterns=len(v['patterns']),mappa_targets=len(targets),mappa_recovered=len(set(targets)&v['patterns'].keys()),slap_targets=len(slap),slap_recovered=len(set(slap)&v['patterns'].keys()),search_comp_cpu=search['cpu_seconds']+comp.get('cpu_seconds',0)))
 if len(values)==6:
  paired.append(c)
  print(c,[(s,m,values[s,m]['minimum'] if values[s,m]['minimum'] is not None else ('empty' if values[s,m]['no_full_mapping'] else 'above')) for s in [0,1] for m in ['none','local','prefix']])
print('paired cases',len(paired),paired)
for sweep in [0,1]:
 for mode in ['none','local','prefix']:
  rs=[r for r in rows if r['case'] in paired and r['sweep']==sweep and r['mode']==mode]
  print(sweep,mode,'full',sum(not r['no_full'] for r in rs),'mappa',sum(r['mappa_recovered'] for r in rs),'/',sum(r['mappa_targets'] for r in rs),'slap',sum(r['slap_recovered'] for r in rs),'/',sum(r['slap_targets'] for r in rs),'cpu',round(sum(r['search_comp_cpu'] for r in rs),3))
(S/'summary.json').write_text(json.dumps(dict(rows=rows,paired_cases=paired),indent=2)+'\n')
