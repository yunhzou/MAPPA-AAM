"""Require exact family coverage and merge current-decoder reports in family order."""
import argparse,dataclasses,gzip,hashlib,json,os,resource,time
from pathlib import Path
from decode_shard import load
from mappa.artifacts import decoded_events_record
from mappa.event_patterns import SignedEventIndex
from mappa.postprocessing import EventCandidate,DecodedEvents,EventDecodeConfig


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);a=p.parse_args();r=a.run;start=time.perf_counter()
    source=r/'results/search_memory';out=r/'results/search_memory_partitioned';out.mkdir(exist_ok=False)
    for name in ['catalogue.json.gz','search.json','configuration.json']:os.link(source/name,out/name)
    cat=load(source,a.inputs);reports={};metrics=[];hashes={}
    for shard in range(32):
        metrics.append(json.loads((r/f'shards/result-{shard:03d}.json').read_text()))
        path=r/f'shards/shard-{shard:03d}.jsonl';hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest();ids=[]
        with path.open() as f:
            for line in f:
                report=json.loads(line);i=report['family_id'];assert i not in reports;reports[i]=report;ids.append(i)
        assert ids==list(range(shard,len(cat.families),32))
    assert set(reports)==set(range(len(cat.families)))
    candidates={};final_reports=[]
    for i in range(len(cat.families)):
        report=reports.pop(i)
        assert report['max_events'] is None
        for pattern in report.pop('patterns'):
            c=candidates.setdefault(pattern['id'],EventCandidate(pattern['id'],dict(enumerate(pattern['mapping'])),pattern['events'],pattern['total']))
            c.family_ids.append(i)
        final_reports.append(report)
    decoded=DecodedEvents(sorted(candidates.values(),key=lambda c:(c.total,c.id)),cat,SignedEventIndex(cat.problem),final_reports,EventDecodeConfig(seconds_per_family=120))
    with gzip.open(out/'decoded.json.gz','wt') as f:json.dump(decoded_events_record(decoded),f)
    record=dict(complete=decoded.complete,families=len(final_reports),candidates=len(candidates),minimum_events=decoded.observed_minimum,minimum_proven=decoded.minimum_proven,event_window=decoded.event_window,shards=metrics,shard_sha256=hashes,merge_seconds=time.perf_counter()-start,merge_peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
    (out/'merge-validation.json').write_text(json.dumps(record,indent=2));print(json.dumps({k:v for k,v in record.items() if k not in ['shards','shard_sha256']},indent=2))
    assert decoded.complete, 'Family decoding remains incomplete; never a negative verdict'

if __name__=='__main__':main()
