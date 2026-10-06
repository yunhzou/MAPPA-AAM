from pathlib import Path
import sys,json,gzip,time,hashlib
from run import S,PACKAGE
sys.path.insert(0,str(S/'engine/src'))
from run import save
from mappa.artifacts import read_aam_checkpoint
from mappa.paired_competition import compete_fragments,CompetitionConfig
from mappa.final_branches import FinalBranchCatalogue
from mappa.event_patterns import SignedEventIndex,extract_path_events
base=PACKAGE.parent/'unrestricted_competition_20260915'
def family_key(f):
 r=f.to_record();r.pop('provenance',None);return hashlib.sha256(json.dumps(r,sort_keys=True).encode()).hexdigest()
def run(case):
 trial=json.loads((PACKAGE/'fixed_trials.json').read_text())[str(case)]
 a=read_aam_checkpoint(base/f'inputs/case{case}.pkl.gz');results={};sets={};patterns={};raw={}
 for mode in ['extended','prefix']:
  t=time.process_time();r=compete_fragments(a,CompetitionConfig(proposal_mode=mode,operation_budget=100000,depth_limit=1,seconds=270),trial=trial);cpu=time.process_time()-t
  assert r.pending==0
  cat=FinalBranchCatalogue(a.problem)
  for i,repair in enumerate(r.repairs):cat.add_aam(repair,f'{mode}:{i}')
  sets[mode]={family_key(f) for f in cat.families};idx=SignedEventIndex(a.problem);pats={};cert=[];dt=time.process_time()
  for i,f in enumerate(cat.families):
   f.validate_representative(a.problem);z=extract_path_events(f.as_path(a.problem),a.problem,idx,max_events=8 if case==64 else 5,seconds=20,max_patterns=None);cert.append(z['complete'])
   for p in z['patterns']:pats[p['id']]=p
  assert all(cert)
  if mode=='prefix':
   assert all(x['available_prefix_sizes']==x['offered_prefix_sizes'] for x in r.counts['prefix_audit'])
  patterns[mode]=set(pats)
  results[mode]=dict(search_cpu=cpu,decode_cpu=time.process_time()-dt,families=len(cat.families),patterns=pats,counts=r.counts,repairs=len(r.repairs),offers=r.offers,complete=True,family_hashes=sorted(sets[mode]))
  raw[mode]=cat.to_record()
 assert sets['extended']<=sets['prefix'];assert patterns['extended']<=patterns['prefix']
 out=S/f'paired-control/{case}';save(out/'result.json',dict(case=case,trial=trial,results=results,maximal_families_retained=True))
 (out/'catalogues.json.gz').write_bytes(gzip.compress(json.dumps(raw).encode(),mtime=0))
 print(case,{m:{k:v[k] for k in ['search_cpu','decode_cpu','families','repairs']} for m,v in results.items()},flush=True)
if __name__=='__main__':
 for case in [6,11,59,64,101]:run(case)
