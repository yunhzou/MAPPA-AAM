"""Audit minimum witnesses from saved decoder artifacts without rerunning decoding."""
import argparse,gzip,json,time
from pathlib import Path
import numpy as np
from mappa import AAMProblem,MolecularEndpoint
from mappa.final_branches import FinalBranchCatalogue
from mappa.event_patterns import SignedEventIndex
from mappa.postprocessing import DecodedEvents,EventCandidate,EventDecodeConfig


def main():
    p=argparse.ArgumentParser();p.add_argument('--result',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);a=p.parse_args();start=time.perf_counter()
    ends=[]
    for name in ['reference_arrays.npz','prediction_arrays.npz']:
        with np.load(a.inputs/name,allow_pickle=False) as d:ends.append(MolecularEndpoint(tuple(d['elements'].tolist()),d['coordinates'],d['wbo']))
    problem=AAMProblem(*ends)
    c=FinalBranchCatalogue.from_record(problem,json.load(gzip.open(a.result/'catalogue.json.gz','rt')))
    data=json.load(gzip.open(a.result/'decoded.json.gz','rt'))
    candidates=[EventCandidate(row['id'],{int(k):v for k,v in row['mapping'].items()},row['events'],row['total'],row['family_ids']) for row in data['candidates']]
    decoded=DecodedEvents(candidates,c,SignedEventIndex(problem),data['family_reports'],EventDecodeConfig(**data['config']))
    results=[]
    for candidate in decoded.minimum_candidates:
        q=decoded.query(candidate,candidate.mapping,seconds=120)
        results.append(dict(candidate=candidate.id,status=q.status,reason=q.reason,family_id=q.family_id))
    output=dict(all_allowed=all(row['status']=='allowed' for row in results),witness_count=len(results),results=results,seconds=time.perf_counter()-start)
    (a.result/'membership-audit.json').write_text(json.dumps(output,indent=2));print(json.dumps(output,indent=2))
    assert output['all_allowed'] and results

if __name__=='__main__':main()
