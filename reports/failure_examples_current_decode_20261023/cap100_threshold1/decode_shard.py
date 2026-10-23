"""Partition the current decoder's independent family loop; retain each result."""
import argparse,gzip,hashlib,json,os,resource,sys,time
from pathlib import Path
import numpy as np
from mappa import AAMProblem,MolecularEndpoint
from mappa.final_branches import FinalBranchCatalogue
from mappa.event_patterns import SignedEventIndex,extract_path_events
from mappa.postprocessing import decode_events,EventDecodeConfig


def load(root,inputs):
    cfg=json.loads((root/'configuration.json').read_text());ends=[]
    for name in ['reference_arrays.npz','prediction_arrays.npz']:
        path=inputs/name
        assert hashlib.sha256(path.read_bytes()).hexdigest()==cfg['inputs'][name]
        with np.load(path,allow_pickle=False) as d:ends.append(MolecularEndpoint(tuple(d['elements'].tolist()),d['coordinates'],d['wbo']))
    problem=AAMProblem(*ends,name=cfg['case'])
    return FinalBranchCatalogue.from_record(problem,json.load(gzip.open(root/'catalogue.json.gz','rt')))


def family_reports(cat,ids,seconds):
    index=SignedEventIndex(cat.problem,threshold=.5,metal_threshold=.3);source_edges={}
    for i in ids:
        family=cat.families[i];floor=family.policy[0]
        if floor not in source_edges:
            source_edges[floor]=tuple((int(a),int(b)) for a,b in zip(index.a,index.b) if index.r[a,b]>=floor)
        report=extract_path_events(family.as_path(cat.problem,source_edges[floor]),cat.problem,index,max_events=None,seconds=seconds,max_patterns=None)
        yield dict(family_id=i,**report)


def main():
    p=argparse.ArgumentParser();p.add_argument('--result',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path);p.add_argument('--shard',type=int,default=0);p.add_argument('--shards',type=int,default=32);p.add_argument('--validate',action='store_true');a=p.parse_args()
    if sys.platform=='linux':resource.setrlimit(resource.RLIMIT_AS,(5*1024**3,5*1024**3))
    start=time.perf_counter();cat=load(a.result,a.inputs)
    if a.validate:
        original=decode_events(cat,EventDecodeConfig(seconds_per_family=120))
        reports=list(family_reports(cat,range(len(cat.families)),120))
        patterns={}
        for report in reports:
            for r in report['patterns']:patterns.setdefault(r['id'],(r['total'],r['events'],r['mapping']))
        expected={c.id:(c.total,c.events,[c.mapping[i] for i in range(cat.problem.atom_count)]) for c in original.candidates}
        assert patterns==expected
        for lhs,rhs in zip(reports,original.family_reports):
            assert all(lhs[k]==rhs[k] for k in ['family_id','complete','reason','solver_queries','max_events'])
        print(json.dumps(dict(validation='passed',families=len(reports),candidates=len(patterns))),flush=True);return
    assert a.output is not None
    a.output.mkdir(parents=True,exist_ok=True);path=a.output/f'shard-{a.shard:03d}.jsonl'
    ids=list(range(a.shard,len(cat.families),a.shards));complete=0;incomplete=0
    with path.open('x') as f:
        for report in family_reports(cat,ids,120):
            f.write(json.dumps(report)+'\n');f.flush()
            complete+=int(report['complete']);incomplete+=int(not report['complete'])
            progress=dict(shard=a.shard,total=len(ids),saved=complete+incomplete,complete=complete,incomplete=incomplete,seconds=time.perf_counter()-start,last_family=report['family_id'],last_reason=report['reason'])
            (a.output/f'progress-{a.shard:03d}.json').write_text(json.dumps(progress))
    progress.update(status='complete' if not incomplete else 'incomplete',peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024 if sys.platform=='linux' else 1024**2))
    (a.output/f'result-{a.shard:03d}.json').write_text(json.dumps(progress,indent=2));print(json.dumps(progress),flush=True)

if __name__=='__main__':main()
