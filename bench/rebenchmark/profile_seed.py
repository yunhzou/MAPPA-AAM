"""Replay one saved seed for native sampling; verify its exact graph afterward."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile
import tempfile
import time


def run(archive_path,seed,output):
    from mappa import AAMProblem,MolecularEndpoint,AAMSearchConfig,aam
    from mappa.artifacts import read_raw_cut
    from mappa.native_search import find_islands_native
    from mappa.matcher import _nauty_orbits
    from rebenchmark.search import select_contexts
    with tarfile.open(archive_path) as archive:
        manifest=json.load(archive.extractfile('archive-manifest.json'))['files']
        def checked(name):
            data=archive.extractfile(name).read()
            assert hashlib.sha256(data).hexdigest()==manifest[name]['sha256']
            return data
        raw=json.loads(checked('input.json'));task=json.loads(checked('task.json'))
        with tempfile.TemporaryDirectory(prefix='mappa-seed-profile-',dir='/tmp') as td:
            path=Path(td)/'cut.raw.pkl.gz'
            path.write_bytes(checked('cuts/cut_00000.raw.pkl.gz'))
            expected=select_contexts(read_raw_cut(path),[seed])
    endpoints=[MolecularEndpoint(**raw[k]) for k in ('reactant','product')]
    if task['direction']=='P_to_R':endpoints.reverse()
    problem=AAMProblem(*endpoints,raw['name']);config=AAMSearchConfig(**task['config'])
    context=expected.contexts[0]
    aam._initialize_search(problem,config,'reused_native')
    view=aam._SEARCH_REPAIR.for_cut(context.cuts)
    source_orbits=_nauty_orbits(view.source,wbo_tol=config.iso_tolerance)
    profile=[];cpu=time.process_time();wall=time.perf_counter()
    graph=find_islands_native(view.source,view.target,context.seed_order,
        graph_floor=config.graph_floor,iso_tol=config.iso_tolerance,max_branches=config.branch_limit,
        p_orbits=aam._SEARCH_CONTEXT[3],r_orbits=source_orbits,cuts=context.cuts,
        growth_replay=view,random_seed=config.random_seed,profile=profile)
    result=dict(cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        seed=seed,direction=task['direction'],profile=profile,graph_equal=graph==expected,
        states=len(graph.states),stops=[vars(s) for s in graph.stops],
        repair_stats=aam._SEARCH_REPAIR.stats())
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    assert result['graph_equal']
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive');parser.add_argument('output');parser.add_argument('--seed',type=int,required=True)
    args=parser.parse_args();run(args.archive,args.seed,args.output)
