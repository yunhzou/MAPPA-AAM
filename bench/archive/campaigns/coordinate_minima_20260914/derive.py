"""Audit saved coordinate event minima and measured stage CPU; no mapping runs."""
from pathlib import Path
import json,gzip,hashlib,statistics,math,collections
R=Path('/Users/yunhengz/Desktop/AAM Writing');W=Path('/Users/yunhengz/Documents/Codex/2026-09-11/aa/work/current-validation');S=Path(__file__).resolve().parent
load=lambda p:json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sources=['reports/current_validation_20260912/coordinate-decoded.json.gz','reports/current_validation_20260912/coordinate-slap-comparison.json.gz','manuscript/evidence/timing_comparison.json']
mappas=load(R/sources[0]);comps=load(R/sources[1]);timing=load(R/sources[2]);rows=[];stages=collections.defaultdict(list);hashes={}
def record(p):hashes[str(p.relative_to(W))]=sha(p);return load(p)
for g,c in zip(mappas,comps['rows']):
 i=g['case'];assert i==c['case'] and g['complete_saved_window'];m=min(p['total'] for p in g['patterns'].values());ids=sorted(k for k,p in g['patterns'].items() if p['total']==m);assert m<=g['window']
 row=dict(case=i,mappa_minimum=m,mappa_minimum_ids=ids,window=g['window'],comparators={})
 for name,d in c['methods'].items():
  other=set(d['patterns']);assert all(p['total']==d['minimum'] for p in d['patterns'].values());common=sorted(set(ids)&other)
  row['comparators'][name]=dict(minimum=d['minimum'],ids=sorted(other),count_relation='lower' if m<d['minimum'] else 'equal' if m==d['minimum'] else 'higher',shared_minimum_ids=common,covered_in_mappa_window=sorted(set(g['patterns'])&other))
 rows.append(row)
 sr=record(W/f'runs/coordinate/seed1/case{i}/R_to_P/search.json');se=record(W/f'runs/coordinate/seed1/case{i}/R_to_P/search_execution.json');assert se['status']=='passed' and se.get('host','Mac')=='Mac' and sr['config']['seed_count']==1 and sr['config']['branch_limit']==2000
 cp=record(W/f'competition/case{i}/result.json');assert cp['operation_budget']==128
 sl=record(W/f'slap-coordinate/results/{i}.json');assert sl['complete']
 sc=record(W/f'slap-coordinate/slap_sweep/slap_scoring/{i}.refined.json')
 vals=dict(mappa_search=sr['cpu_seconds'],mappa_competition=cp['total_cpu_seconds'],mappa_decoding=timing['coordinate_decoding']['per_case'][i]['cpu_seconds'],slap_mapping=sl['mapping_cpu_seconds'],slap_h_refinement=sc['scoring_cpu_seconds'])
 for k,v in vals.items():assert v>=0;stages[k].append(dict(case=i,cpu_seconds=v))
assert [r['case'] for r in rows]==list(range(140))
def stats(rows):
 xs=sorted(r['cpu_seconds'] for r in rows);q=(len(xs)-1)*.95;l=math.floor(q);u=math.ceil(q)
 return dict(cases=len(xs),mean=statistics.mean(xs),median=statistics.median(xs),p95=xs[l]+(xs[u]-xs[l])*(q-l),total=sum(xs))
summary={}
for name in ['native_slap','slap_sweep']:
 rs=[r['comparators'][name] for r in rows];eq=[r for r in rs if r['count_relation']=='equal'];summary[name]=dict(count_relations=dict(collections.Counter(r['count_relation'] for r in rs)),minimum_patterns=sum(len(r['ids']) for r in rs),shared_minimum_patterns=sum(len(r['shared_minimum_ids']) for r in rs),equal_minimum_patterns=sum(len(r['ids']) for r in eq),equal_cases_complete=sum(set(r['ids'])<=set(r['shared_minimum_ids']) for r in eq),covered_in_mappa_window=sum(len(r['covered_in_mappa_window']) for r in rs),lower_cases=[row['case'] for row in rows if row['comparators'][name]['count_relation']=='lower'])
out=dict(cases=140,scope='Unverified coordinate mappings: comparison of minimum signed bond-event counts and classes, not a mapping-accuracy benchmark. Minima refer to returned candidate sets/saved families, not all atom mappings or proven chemistry.',mapping_runs=0,mappa_minimum_patterns=sum(len(r['mappa_minimum_ids']) for r in rows),comparisons=summary,per_case=rows,timing={k:dict(stats=stats(v),per_case=v) for k,v in stages.items()},timing_scope='Recorded successful stages on Mac, all140cases; MAPPA search includes checkpointI/O, competition includes loading/proposals/admission/persistence, fullwindowdecoder includes continuations. SLAP binarybidirectionalmaptimer includes cut perturbation and mappercall, excludes preparation/export; Hrefinement/initialscoring separately timed. Final shared signed-event rescoring and failed/replaced attempts not included. No end-to-end or minimum-only speed ratio.',source_hashes={f:sha(R/f) for f in sources})
(S/'coordinate_minima.json').write_text(json.dumps(out,separators=(',',':'))+'\n');(S/'timing-source-hashes.json').write_text(json.dumps(hashes,separators=(',',':'))+'\n')
print(json.dumps(dict(mappa_minimum_patterns=out['mappa_minimum_patterns'],comparisons=summary,timing={k:v['stats'] for k,v in out['timing'].items()}),indent=2))
