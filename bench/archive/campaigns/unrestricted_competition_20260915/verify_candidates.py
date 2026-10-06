"""Check each new representative against its exact saved family and event ID."""
from pathlib import Path
import gzip,json,os,sys
package=Path(__file__).resolve().parent
work=Path(os.environ.get('MAPPA_EXPERIMENT_WORK',package/'work')).resolve()
sys.path.insert(0,str(work/'engine/src'))
from mappa.artifacts import read_aam_checkpoint
from mappa.final_branches import FinalFamily
from mappa.family_query import query_path
from mappa.event_patterns import SignedEventIndex
proofs=json.loads(gzip.decompress((package/'novel-candidates.json.gz').read_bytes()))
for proof in proofs:
 a=read_aam_checkpoint(package/f"inputs/case{proof['case']}.pkl.gz")
 f=FinalFamily.from_record(proof['family']);p=proof['pattern']
 f.validate_representative(a.problem)
 status,_=query_path(f.as_path(a.problem),a.problem,dict(enumerate(p['mapping'])),source_atoms=range(a.problem.source_atom_count),timeout_ms=10000)
 assert status=='recovered'
 assert SignedEventIndex(a.problem).describe(p['mapping'])['id']==p['id']
 print(proof['case'],proof['mode'],p['total'],status,p['id'])
print(f'All {len(proofs)} new representatives verified.')
