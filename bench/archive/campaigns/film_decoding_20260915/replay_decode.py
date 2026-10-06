"""Reproduce the film's two canonical event classes from its saved final catalogue."""
import gzip,json,sys,time
from pathlib import Path
S=Path(__file__).resolve().parent
# Run in an installed checkout, or from reports/film_decoding_20260915.
repo=S.parents[1];sys.path.insert(0,str(repo/'src'))
from mappa import AAMProblem,MolecularEndpoint
from mappa.final_branches import FinalBranchCatalogue
from mappa.event_patterns import SignedEventIndex,extract_path_events
source=json.loads((S/'input.json').read_text())
problem=AAMProblem(*(MolecularEndpoint(**source[k]) for k in ['reactant','product']))
cat=FinalBranchCatalogue.from_record(problem,json.loads(gzip.decompress((S/'final-catalogue.json.gz').read_bytes())))
index=SignedEventIndex(problem);ids=set();c=time.process_time()
for i,f in enumerate(cat.families):
 f.validate_representative(problem)
 result=extract_path_events(f.as_path(problem),problem,index,max_events=5,seconds=20,max_patterns=None)
 assert result['complete'],f'Incomplete family {i}'
 for p in result['patterns']:
  assert index.describe(p['mapping'])['id']==p['id'];ids.add(p['id'])
expected=json.loads((S/'final-decoding.json').read_text())
assert ids==set(expected['patterns']) and len(ids)==2
print(json.dumps(dict(complete=True,branches=cat.branch_count,families=len(cat.families),unique_event_classes=len(ids),cpu_seconds=time.process_time()-c)))
