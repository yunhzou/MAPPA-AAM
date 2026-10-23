"""Audit each decoded witness against its recorded producing family."""
import argparse,gzip,json,time
from pathlib import Path
from decode_shard import load
from mappa.final_branches import FinalBranchCatalogue
from mappa.event_patterns import SignedEventIndex
from mappa.postprocessing import DecodedEvents,EventCandidate,EventDecodeConfig


def main():
    p=argparse.ArgumentParser();p.add_argument('--result',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);a=p.parse_args();started=time.perf_counter()
    cat=load(a.result,a.inputs);data=json.load(gzip.open(a.result/'decoded.json.gz','rt'));checked=0;minimum=min(c['total'] for c in data['candidates']);queries=[]
    for row in data['candidates']:
        mapping={int(k):v for k,v in row['mapping'].items()};fid=row['family_ids'][0];family=cat.families[fid]
        # Invariance outputs in these cases are the retained family representatives.
        assert mapping==dict(family.mapping), 'Witness requires an action proof rather than the identity proof'
        family.validate_representative(cat.problem);checked+=1
        if row['total']==minimum:
            # Positive membership is checked in its recorded supporting family,
            # without spending time rejecting thousands of unrelated families.
            single=FinalBranchCatalogue(cat.problem);single.families=[family]
            candidate=EventCandidate(row['id'],mapping,row['events'],row['total'],[0])
            decoded=DecodedEvents([candidate],single,SignedEventIndex(cat.problem),[],EventDecodeConfig(seconds_per_family=120))
            q=decoded.query(candidate,mapping,seconds=120)
            queries.append(dict(candidate=row['id'],producing_family_id=fid,status=q.status,reason=q.reason))
            assert q.status=='allowed'
    result=dict(status='passed',witnesses=checked,all_identity_action_witnesses_valid=True,minimum_witness_queries=queries,minimum_witness_count=len(queries),seconds=time.perf_counter()-started)
    (a.result/'producing-family-audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

if __name__=='__main__':main()
