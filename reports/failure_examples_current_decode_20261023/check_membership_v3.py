"""Audit each decoded witness against its recorded producing family."""
import argparse,gzip,json,time
from pathlib import Path
from decode_shard import load
from mappa.final_branches import FinalBranchCatalogue
from mappa.event_patterns import SignedEventIndex
from mappa.postprocessing import DecodedEvents,EventCandidate,EventDecodeConfig


def main():
    p=argparse.ArgumentParser();p.add_argument('--result',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);a=p.parse_args();started=time.perf_counter()
    cat=load(a.result,a.inputs);data=json.load(gzip.open(a.result/'decoded.json.gz','rt'));checked=0;identity=0;solver=0;minimum=min(c['total'] for c in data['candidates']);queries=[]
    for row in data['candidates']:
        mapping={int(k):v for k,v in row['mapping'].items()};fid=row['family_ids'][0];family=cat.families[fid]
        family.validate_representative(cat.problem)
        is_identity=mapping==dict(family.mapping)
        identity+=int(is_identity)
        checked+=1
        if not is_identity or row['total']==minimum:
            # Positive membership is checked in its recorded supporting family,
            # without spending time rejecting thousands of unrelated families.
            single=FinalBranchCatalogue(cat.problem);single.families=[family]
            candidate=EventCandidate(row['id'],mapping,row['events'],row['total'],[0])
            decoded=DecodedEvents([candidate],single,SignedEventIndex(cat.problem),[],EventDecodeConfig(seconds_per_family=120))
            q=decoded.query(candidate,mapping,same_events=False,seconds=120)
            solver+=1
            queries.append(dict(candidate=row['id'],producing_family_id=fid,status=q.status,reason=q.reason,is_minimum=row['total']==minimum))
            assert q.status=='allowed'
    result=dict(status='passed',witnesses=checked,identity_action_witnesses=identity,exact_family_queries=solver,queries=queries,minimum_witness_count=sum(q['is_minimum'] for q in queries),seconds=time.perf_counter()-started)
    (a.result/'producing-family-audit-v3.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='queries'},indent=2))

if __name__=='__main__':main()
