"""Derive size-ordered Golden recovery from saved inputs and verdicts; no search."""
from pathlib import Path
import json,hashlib,collections,argparse
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
src=a.repo/'manuscript/evidence/seed_comparison.json';manifest=a.repo/'reports/current_validation_20260912/inputs.json';seed=json.loads(src.read_text());sizes={}
for r in json.loads(manifest.read_text())['records']:
 p=a.inputs/str(r['case'])/'input.json';assert sha(p)==r['input_sha256'];d=json.loads(p.read_text());sizes[r['case']]=dict(case=r['case'],input_sha256=r['input_sha256'],counts=[len(d[k]['elements']) for k in ('reactant','product')])
assert set(sizes)==set(range(1851))
methods={}
for key in ['seeds1','seeds2']:
 rows=[];groups={}
 for r in seed['methods'][key]['per_case']:
  c=r['case'];n,m=sizes[c]['counts'];small='R_to_P' if n<=m else 'P_to_R';large='P_to_R' if n<=m else 'R_to_P';ds={x['direction']:x['outcome'] for x in r['directions']}
  assert set(ds)=={'R_to_P','P_to_R'} and all(v in ('recovered','not_recovered') for v in ds.values())
  union='recovered' if 'recovered' in ds.values() else 'not_recovered';assert union==r['outcome']
  rows.append(dict(**sizes[c],cohort='equal' if n==m else 'unequal',smaller_first=ds[small],larger_first=ds[large],bidirectional=union))
 for group in ['all','unequal','equal']:
  subset=[r for r in rows if group=='all' or r['cohort']==group];groups[group]={'n':len(subset)}
  for policy in ['smaller_first','larger_first','bidirectional']:
   hits=[r['case'] for r in subset if r[policy]=='recovered'];groups[group][policy]=dict(recovered=len(hits),not_recovered=len(subset)-len(hits),unknown=0,percent=100*len(hits)/len(subset),cases=hits)
 methods[key]=dict(groups=groups,per_case=rows)
out=dict(denominator=1851,definition='Endpoint size is total explicit-atom count, including hydrogen. Smaller-first uses the original input orientation on equal-size endpoints; larger-first uses its reverse. Bidirectional is their union. No correspondence labels used to select orientation.',source_hashes={str(src.relative_to(a.repo)):sha(src),str(manifest.relative_to(a.repo)):sha(manifest),'src/mappa/search_orientation.py':sha(a.repo/'src/mappa/search_orientation.py')},methods=methods,no_search_reruns=True)
a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n')
for k,v in methods.items():print(k, {g:{p:(x['recovered'],round(x['percent'],2)) for p,x in d.items() if p!='n'} for g,d in v['groups'].items()})
