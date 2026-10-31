"""One exact cut on macOS; select historical transport only in this diagnostic."""
import argparse
from dataclasses import asdict
import gzip
import importlib.abc
import importlib.util
import json
from pathlib import Path
import pickle
import resource
import sys
import tempfile
import time


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input'); parser.add_argument('output')
    parser.add_argument('--historical-build')
    parser.add_argument('--downstream',action='store_true')
    args=parser.parse_args()
    output=Path(args.output); output.mkdir(parents=True,exist_ok=True)
    if args.historical_build:
        before=Path(args.historical_build)
        paths={'mappa._engine':next((before/'src/mappa').glob('_engine*.so')),
               'mappa.mappa_core.core.search.graph':before/'graph.py',
               'mappa.mappa_core.core.search.branching.native':before/'adapter.py'}
        class HistoricalTransport(importlib.abc.MetaPathFinder):
            def find_spec(self,fullname,path=None,target=None):
                if fullname in paths:
                    return importlib.util.spec_from_file_location(fullname,paths[fullname])
        sys.meta_path.insert(0,HistoricalTransport())
    from mappa.mappa_core.core.scheduler.runtime import Workspace
    from mappa.mappa_core.core.scheduler import search
    from mappa.mappa_core.core.search.session import SearchSession
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    from mappa import decode_events, query_path
    from mappa.mappa_core.utils.io.records import decoded_events_record
    with gzip.open(args.input,'rb') as stream:
        saved=pickle.load(stream)
    source=saved['source']; task=saved['task']
    record=dict(task=task.key,atoms=[len(source.problem.reactant.elements),len(source.problem.product.elements)],
                config=repr(source.config),historical_build=args.historical_build,runs=[])
    search_cut=SearchSession.search_cut
    def measured(self,*values,**options):
        started=time.perf_counter()
        result=search_cut(self,*values,**options)
        record['runs'][-1]['search_seconds']=time.perf_counter()-started
        return result
    SearchSession.search_cut=measured
    with tempfile.TemporaryDirectory(prefix='mappa-point-') as scratch:
        workspace=Workspace(scratch,str(output/'bundles'))
        payload={**task.payload,'input':workspace.writer.put(source)}
        for repetition in range(2):
            run=dict(repetition=repetition); record['runs'].append(run)
            started=time.perf_counter(); cpu=time.process_time()
            result=search.search_complete_cut(payload,workspace)
            run.update(compute_seconds=time.perf_counter()-started,compute_cpu_seconds=time.process_time()-cpu,
                peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                states=len(result.graph.states),transitions=len(result.graph.transitions),terminals=len(result.graph.terminals))
            reference=workspace.writer.put(result)
            run.update(total_seconds=time.perf_counter()-started,reference=asdict(reference),
                       peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
            print(json.dumps(run),flush=True)
            del result
        if args.downstream:
            graph=workspace.load(reference).graph
            started=time.perf_counter()
            catalogue=FinalBranchCatalogue(source.problem)
            catalogue.add_graph(graph,archive='mac-point')
            reduced=catalogue.to_record()
            record['deduplication_seconds']=time.perf_counter()-started
            started=time.perf_counter()
            decoded=decoded_events_record(decode_events(catalogue))
            record['decoding_seconds']=time.perf_counter()-started
            reference_mapping=dict(saved['reference']['mapping'])
            statuses=[]
            started=time.perf_counter()
            for family in catalogue.families:
                status,evidence=query_path(family.as_path(source.problem),source.problem,
                    reference_mapping,source_atoms=tuple(sorted(reference_mapping)),
                    complete_reference=True,timeout_ms=1000)
                statuses.append(status)
                if status=='recovered':
                    break
            record['verification_seconds']=time.perf_counter()-started
            record['verification_status']='recovered' if 'recovered' in statuses else 'unknown' if 'unknown' in statuses else 'not_recovered'
            record['families']=len(catalogue.families)
            with gzip.open(output/'scientific.pkl.gz','wb') as stream:
                pickle.dump(dict(catalogue=reduced,decoded=decoded,verification=record['verification_status']),stream,protocol=5)
            print(json.dumps({k:v for k,v in record.items() if k!='runs'}),flush=True)
    (output/'measurement.json').write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':
    main()
