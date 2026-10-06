from pathlib import Path
import sys,json,hashlib,time,os
S=Path(os.environ['MAPPA_EXPERIMENT_WORK']).resolve();sys.path.insert(0,str(S/'engine/src'))
BEFORE=Path(os.environ['MAPPA_REFERENCE_WORK']).resolve()
from mappa.artifacts import read_aam_checkpoint

def digest(graph):return hashlib.sha256(json.dumps(graph.to_record(copy=False),sort_keys=True,separators=(',',':')).encode()).hexdigest()
def reference(case,mode):return BEFORE/f'runs/{case}/{mode}'
rows=[]
for c,m in [(c,'local') for c in [6,11,59,64,101]]+[(c,'prefix') for c in [6,11,59,64,101]]+[(11,'extended')]:
 old=reference(c,m);new=S/f'runs/{c}/{m}';a=json.loads((old/'search.json').read_text());b=json.loads((new/'search.json').read_text())
 ca={k:v for k,v in a['counts'].items() if v};cb={k:v for k,v in b['counts'].items() if v and k not in ['symmetry_reuse','validation_reuse']};assert ca==cb,(c,m,'counts')
 assert a['repairs']==b['repairs'] and a['pending']==b['pending'];assert len(a['offers'])==len(b['offers']);offers=[]
 for x,y in zip(a['offers'],b['offers']):
  y=dict(y)
  if 'prefix' not in x:y.pop('prefix',None)
  assert x==y,(c,m,'offer');offers.append(y)
 hashes=[]
 for i in range(a['repairs']):
  aa=read_aam_checkpoint(old/f'takeover{i}.pkl.gz');bb=read_aam_checkpoint(new/f'takeover{i}.pkl.gz')
  assert aa.config==bb.config
  ah,bh=digest(aa.graph),digest(bb.graph);assert ah==bh,(c,m,i,'graph')
  hashes.append(ah)
 row=dict(case=c,mode=m,repairs=a['repairs'],graphs_identical=True,counts_identical=True,offers_identical=True,repair_sha256=hashes,old_search_cpu=a['cpu_seconds'],new_search_cpu=b['cpu_seconds'])
 if (new/'decode.json').exists():
  x=json.loads((old/'decode.json').read_text());y=json.loads((new/'decode.json').read_text())
  assert x['complete'] and y['complete'];assert x['patterns']==y['patterns'];assert (x['branches'],x['families'])==(y['branches'],y['families'])
  row.update(patterns_identical=True,patterns=len(y['patterns']),families=y['families'],old_decode_cpu=x['cpu_seconds'],new_decode_cpu=y['cpu_seconds'])
 rows.append(row);(S/'verification.json').write_text(json.dumps(rows,indent=2)+'\n');print(c,m,'identical',a['repairs'],'search',round(a['cpu_seconds'],3),'->',round(b['cpu_seconds'],3),flush=True)
