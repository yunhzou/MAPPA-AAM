"""Functional gate before the full timing array; no reference-guided search."""
import sys,json,hashlib,subprocess,platform,socket,os
from pathlib import Path
import run
S=Path(__file__).resolve().parent
run.verify_sources()
prior=run.read(S/'prior-cap100-verdicts.json')
# The complete schedule and isolated verification are exercised at both caps.
for c in [9,1001]:
 subprocess.run([sys.executable,str(S/'run.py'),'case',str(c)],check=True)
 for n in [1,2,3,10]:
  for d in ['R_to_P','P_to_R']:
   old=prior.get(f'{c}/{n}/{d}')
   new=run.read(run.folder(c,f'mappa_c100_s{n}_sweep',d)/'evaluation.json')['reference_recovery']
   if old is not None:assert old==new,(c,n,d,old,new)
 for cfg in run.configurations():
  for d in ['R_to_P','P_to_R']:
   p=run.folder(c,cfg['key'],d)
   assert run.read(p/'search-execution.json')['status']=='passed',(c,cfg,d)
   assert run.read(p/'score-execution.json')['status']=='passed',(c,cfg,d)
   assert run.read(p/'search.json')['search_complete']
   assert run.read(p/'evaluation.json')['reference_recovery']!='unknown'
   if not cfg['sweep'] and cfg['method']=='mappa':
    v=run.read(p/'evaluation.json');assert v['checkpoint_verification']['expected_cuts']==1
run.save(S/'pilot-proof.json',dict(status='passed',cases=[9,1001],configurations=12,directions=2,cap100_prior_verdict_parity=True,no_sweep_expected_cut_verified=True,host=socket.gethostname(),python=sys.version,source_manifest_sha256=run.sha(S/'manifest.json')))
print('Controlled timing pilot passed',flush=True)
