"""Audit saved per-direction and per-call timings without running mapping."""
from pathlib import Path
import json,gzip,hashlib,statistics,math
R=Path('/Users/yunhengz/Desktop/AAM Writing');W=Path('/Users/yunhengz/Documents/Codex/2026-09-11/aa/work/current-validation');OUT=Path(__file__).resolve().parent
load=lambda p:json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
seed=load(R/'manuscript/evidence/seed_comparison.json');orient=load(R/'manuscript/evidence/direction_recovery.json');common=seed['common_case_indices'];assert len(common)==1821
raw=load(R/'reports/current_validation_20260912/golden-direction-records.json.gz');index={(r['seed'],r['case'],r['direction']):r['files'] for r in raw};slap=load(R/'reports/current_validation_20260912/slap-golden-evaluations.json.gz');slap={r['case']:r for r in slap}
sizes={r['case']:r['counts'] for r in orient['methods']['seeds1']['per_case']}
def direction(c):return ['R_to_P','P_to_R'] if sizes[c][0]<=sizes[c][1] else ['P_to_R','R_to_P']
def stats(xs):
 xs=sorted(xs);q=(len(xs)-1)*.95;l=math.floor(q);u=math.ceil(q)
 return dict(n=len(xs),mean=statistics.mean(xs),median=statistics.median(xs),p95=xs[l]+(xs[u]-xs[l])*(q-l),total=sum(xs))
methods={};records=[];hashes={}
for sn in [1,2]:
 for policy in ['smaller_first','larger_first','bidirectional']:
  rows=[]
  for c in common:
   dirs=direction(c);ds=dirs if policy=='bidirectional' else [dirs[0 if policy=='smaller_first' else 1]]
   fs=[index[(sn,c,d)] for d in ds];assert all(f['search_execution.json']['status']=='passed' for f in fs)
   assert seed['methods'][f'seeds{sn}']['per_case'][c]['search_hosts']==['Mac']  # direction host can refer to later verification
   assert all(f['search_execution.json'].get('host','Mac')=='Mac' for f in fs)
   cpu=sum(f['search.json']['cpu_seconds'] for f in fs)
   wall=sum(f['search.json']['wall_seconds'] for f in fs)
   rows.append(dict(case=c,cpu_seconds=cpu,call_wall_seconds=wall))
   if policy=='bidirectional':assert abs(cpu-seed['methods'][f'seeds{sn}']['per_case'][c]['search_cpu_including_io'])<1e-7
  methods[f'mappa{sn}_{policy}']=dict(method='MAPPA',seeds=sn,policy=policy,recovered=orient['methods'][f'seeds{sn}']['groups']['all'][policy]['recovered'],recovery_denominator=1851,stats=stats([r['cpu_seconds'] for r in rows]),per_case=rows)
  if policy=='bidirectional':assert abs(methods[f'mappa{sn}_{policy}']['stats']['mean']-seed['methods'][f'seeds{sn}']['common_mean_cpu_seconds'])<1e-8
# Read only scalar timers from the saved SLAP rows; validate cut completeness and sums.
slap_times={};slap_uncut={}
for i,c in enumerate(common):
 for d in ['R_to_P','P_to_R']:
  timers=dict(cpu_seconds=0.,call_wall_seconds=0.);uncut=dict(cpu_seconds=0.,call_wall_seconds=0.)
  for mode in ['binary','weighted']:
   folder=W/f'slap-golden/outputs/{c}/{d}/{mode}';path=folder/'records.jsonl';data=path.read_bytes();hashes[str(path.relative_to(W))]=hashlib.sha256(data).hexdigest();rs=[json.loads(l) for l in data.splitlines()];complete=load(folder/'complete.json');assert complete['complete'] and len(rs)==complete['expected_cuts'];assert [r['ordinal'] for r in rs]==list(range(len(rs)))
   for r in rs:
    assert r['status']=='mapped';v=sum(r[k+'_cpu'] for k in ['graph','mapping','export']);w=sum(r[k+'_wall'] for k in ['graph','mapping','export']);timers['cpu_seconds']+=v;timers['call_wall_seconds']+=w
    if r['ordinal']==0:uncut['cpu_seconds']+=v;uncut['call_wall_seconds']+=w
  slap_times[c,d]=timers;slap_uncut[c,d]=uncut
 assert abs(sum(slap_times[c,d]['cpu_seconds'] for d in ['R_to_P','P_to_R'])-slap[c]['workflow_cpu_excluding_io'])<1e-7
 if i%400==0:print('Audited SLAP timings:',i+1,flush=True)
for policy in ['smaller_first','larger_first','bidirectional']:
 rows=[];hits=[]
 for c in range(1851):
  ds=direction(c);ds=ds if policy=='bidirectional' else [ds[0 if policy=='smaller_first' else 1]]
  hit=any(a.get('hits') for a in slap[c].get('attempts',[]) if a['direction'] in ds)
  if policy=='bidirectional':assert hit==(slap[c]['reference_recovery']=='recovered'),c
  if hit:hits.append(c)
  if c in common:rows.append(dict(case=c,cpu_seconds=sum(slap_times[c,d]['cpu_seconds'] for d in ds),call_wall_seconds=sum(slap_times[c,d]['call_wall_seconds'] for d in ds)))
 methods[f'slap_{policy}']=dict(method='SLAP sweep',policy=policy,recovered=len(hits),recovery_denominator=1851,stats=stats([r['cpu_seconds'] for r in rows]),per_case=rows)
assert abs(methods['slap_bidirectional']['stats']['mean']-seed['same_host_timing']['metrics']['slap']['mean_seconds'])<1e-8
uncut=stats([sum(slap_uncut[c,d]['cpu_seconds'] for d in ['R_to_P','P_to_R']) for c in common])
comp=load(R/'manuscript/evidence/competitors.json');defaults=[]
for d in comp['methods']:
 t=d['original_mapping_timing'];n=d['statuses']['mapped'];defaults.append(dict(method=d['name'],key=d['method'],calls=n,failed=1851-n,mean_wall_seconds=t['successful_mapping_wall_sum_seconds']/n,median_wall_seconds=t['successful_mapping_median_seconds'],mean_cpu_seconds=t['successful_mapping_cpu_hours']*3600/n,scope='Original cluster; completed calls including invalid predictions. Startup, queue, evaluation and saving excluded; not the Mac timing cohort.'))
result=dict(status='passed',searches_rerun=0,common_case_indices=common,methods=methods,slap_unswept_bidirectional=uncut,default_comparators=defaults,source_hashes={f:sha(R/f) for f in ['manuscript/evidence/seed_comparison.json','manuscript/evidence/direction_recovery.json','manuscript/evidence/competitors.json','reports/current_validation_20260912/golden-direction-records.json.gz','reports/current_validation_20260912/slap-golden-evaluations.json.gz']},scope='CPU seconds per reaction on identical 1821 completed Mac cases; recovery on all1851. MAPPA includes checkpoint I/O and symmetry finalization; SLAP includes graph/mapping/export and excludes I/O. Reference verification excluded. Opposite directions add CPU cost; summed call wall is not parallel latency. Equal explicit-atom sizes use stored orientation for smaller-first and reverse for larger-first.')
(OUT/'timing_comparison.json').write_text(json.dumps(result,separators=(',',':'))+'\n');(OUT/'slap-timing-input-hashes.json').write_text(json.dumps(hashes,separators=(',',':'))+'\n')
print(json.dumps({k:{'recovered':v['recovered'],**v['stats']} for k,v in methods.items()},indent=2));print('SLAP unswept',uncut)
