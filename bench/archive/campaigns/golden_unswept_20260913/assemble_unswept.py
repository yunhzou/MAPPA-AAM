from pathlib import Path
import json,gzip,hashlib,sys,collections
V=Path(__file__).resolve().parent;ROOT=V.parent/'current-validation';REPO=Path('/Users/yunhengz/Desktop/AAM Writing');OUT=V/'report';OUT.mkdir(exist_ok=True)
sys.path[:0]=[str(ROOT/'engine/src'),str(ROOT/'engine/bench')]
from golden_competitors import signatures
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
save=lambda p,d:p.write_text(json.dumps(d,indent=2)+'\n')
def load(p):return json.loads(p.read_text())
mappa=sorted([load(p) for p in (V/'uncut').glob('*.json')],key=lambda r:r['case']);assert len(mappa)==1851
assert all(r['status']=='passed' for c in mappa for r in c['directions'])
source=REPO/'reports/current_validation_20260912/slap-golden-evaluations.json.gz';original=json.load(gzip.open(source,'rt'));slap=[]
refs=[json.loads(s) for s in (ROOT/'golden-data/audit.jsonl').read_text().splitlines()]
for r in original:
 c=r['case'];hit=any(h['ordinal']==0 for a in r['attempts'] for h in a['hits']);row=dict(case=c,outcome='recovered' if hit else 'not_recovered',attempts=[dict(direction=a['direction'],mode=a['mode'],hits=[h for h in a['hits'] if h['ordinal']==0]) for a in r['attempts']])
 if not hit and not r['complete']:
  ep,expected=signatures(refs[c]['mapped_reaction']);checks=[]
  for a in r['attempts']:
   p=ROOT/f"slap-golden/outputs/{c}/{a['direction']}/{a['mode']}/records.jsonl"
   with p.open('rb') as f:line=f.readline()
   check=dict(direction=a['direction'],mode=a['mode'],record_sha256=hashlib.sha256(line).hexdigest())
   if not line.endswith(b'\n'):check.update(status='missing')
   else:
    raw=json.loads(line);assert raw['ordinal']==0 and raw['cut'] is None
    invalid=[];hits=[]
    for i,candidate in enumerate(raw['candidates']):
     try:
      actual_ep,actual=signatures(candidate['mapped_rxn']);assert actual_ep==ep,'Endpoint chemistry changed'
      if actual==expected:hits.append(i)
     except Exception as e:invalid.append(dict(candidate=i,error=str(e)))
    check.update(status=raw['status'],invalid=invalid,hits=hits)
   checks.append(check)
  assert not any(a.get('hits') for a in checks), 'Uncut outcomes disagree with archived evaluation'
  row.update(rechecked_uncut=checks,outcome='not_recovered' if all(a['status']=='mapped' and not a['invalid'] for a in checks) else 'unknown')
 slap.append(row)
 print('',end='',flush=True)
seed=load(REPO/'manuscript/evidence/seed_comparison.json')['methods']['seeds1']['per_case'];sweep=load(REPO/'manuscript/evidence/slap_sweep.json')['per_case']
methods={}
for name,rows,full in [('mappa',mappa,seed),('slap',slap,sweep)]:
 hits={r['case'] for r in rows if r['outcome']=='recovered'};fullhits={r['case'] for r in full if r['outcome']=='recovered'};assert hits<=fullhits
 methods[name]=dict(counts=dict(collections.Counter(r['outcome'] for r in rows)),recovered_cases=sorted(hits),unknown_cases=[r['case'] for r in rows if r['outcome']=='unknown'],sweep_recovered=len(fullhits),added_by_sweep=sorted(fullhits-hits))
 with gzip.open(OUT/f'{name}-uncut-records.json.gz','wt') as f:json.dump(rows,f,separators=(',',':'))
result=dict(denominator=1851,seed_count=1,branch_cap=100,bidirectional=True,competition=False,algorithm_commit='3a9ef9a8aec3e129ea5f4ee67df47ce2c25ade3c',methods=methods,scope='Paired uncut controls taken from the final campaign; MAPPA saved families freshly finalized and evaluated. SLAP unions uncut binary/weighted outputs in both directions; incomplete sweep records are rechecked only at ordinal zero. No search reruns.',sources={'slap_evaluations_sha256':sha(source),'audit_sha256':sha(ROOT/'golden-data/audit.jsonl'),'mappa_driver_sha256':sha(V/'score_uncut.py')},mappa_rescore_wall_seconds=50.87213141703978)
save(OUT/'unswept.json',result);save(V/'manuscript/evidence/unswept.json',result)
print(json.dumps({k:dict(counts=m['counts'],sweep_added=len(m['added_by_sweep']),unknown=m['unknown_cases']) for k,m in methods.items()}),flush=True)
