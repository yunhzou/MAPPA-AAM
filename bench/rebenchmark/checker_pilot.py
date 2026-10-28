"""Replay a saved successful membership query without repeating AAM search."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile
import tempfile
import time

from .evidence import read,save


def replay(archive,cut,seed,query,output,seconds):
    from .campaign import load_plan
    from .search import select_contexts
    from mappa.artifacts import read_graph_checkpoint
    from mappa.search_graph import SearchPath
    from mappa.family_query import compile_path
    from golden_evaluation import endpoint_generators,project,colored_graph
    import pynauty
    import z3
    root=f'analysis/{cut:05d}/{seed-1:02d}'
    with tempfile.TemporaryDirectory(prefix='mappa-checker-') as td:
        folder=Path(td)
        with tarfile.open(archive) as tar:
            manifest=json.load(tar.extractfile('archive-manifest.json'))['files']
            names=['input.json','reference.json','task.json',f'finalized/{cut:05d}.pkl.gz',
                   root+'/evaluation.json',root+'/queries/'+query+'.json']
            for name in names:
                data=tar.extractfile(name).read()
                assert hashlib.sha256(data).hexdigest()==manifest[name]['sha256']
                target=folder/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
        plan=load_plan(folder);reference=read(folder/'reference.json')
        evaluation=read(folder/root/'evaluation.json');previous=read(folder/root/'queries'/f'{query}.json')
        assert previous['outcome']=='recovered'
        assert previous['result']['mapping']==evaluation['witness_actions']['mapping']
        graph=select_contexts(read_graph_checkpoint(folder/f'finalized/{cut:05d}.pkl.gz'),[seed-1])
        path=SearchPath(graph,evaluation['witness_terminal'],tuple(evaluation['witness_path']))
        features=list(reversed(reference['features'])) if plan.reversed else reference['features']
        def lifted(feature,degree):
            result=[]
            for permutation in endpoint_generators(feature):
                image=list(range(degree))
                for i,a in enumerate(feature['heavy']):image[a]=feature['heavy'][permutation[i]]
                result.append(tuple(image))
            return tuple(result)
        started,cpu=time.monotonic(),time.process_time()
        compiled=compile_path(path,plan.problem,plan.to_search_mapping(reference['mapping']),
            source_atoms=features[0]['heavy'],
            source_generators=lifted(features[0],plan.problem.source_atom_count),
            target_generators=lifted(features[1],plan.problem.target_atom_count),complete_reference=True)
        encoding_cpu=time.process_time()-cpu
        query_text=compiled.solver.to_smt2()
        output.mkdir(parents=True,exist_ok=True)
        (output/'query.smt2').write_text(query_text)
        compiled.solver.set(timeout=max(1,int(1000*(seconds-(time.monotonic()-started)))))
        solve_cpu=time.process_time();status=compiled.solver.check();solve_cpu=time.process_time()-solve_cpu
        result=dict(status=str(status),encoding_cpu_seconds=encoding_cpu,solve_cpu_seconds=solve_cpu,
                    total_wall_seconds=time.monotonic()-started,previous_solve_cpu_seconds=previous['result']['cpu_seconds'],
                    query_bytes=len(query_text),cut=cut,seed=seed,archive=str(archive),previous_query=query,
                    preserved_input_members={n:manifest[n] for n in names},search_rerun=False)
        if status==z3.sat:
            witness=compiled.realize(compiled.solver.model())
            actual=project(plan.to_input_mapping(witness['mapping']),reference['features'])
            expected=project(reference['mapping'],reference['features'])
            assert pynauty.certificate(colored_graph(reference['features'],actual))==pynauty.certificate(colored_graph(reference['features'],expected))
            result.update(witness=witness,reference_verified=True)
        elif status==z3.unknown:result['reason']=compiled.solver.reason_unknown()
        else:raise AssertionError('Saved SAT query became UNSAT')
        save(output/'result.json',result)
        print(json.dumps({k:v for k,v in result.items() if k not in ('witness','preserved_input_members')}),flush=True)
        return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archive',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--cut',type=int,required=True);p.add_argument('--seed',type=int,required=True)
    p.add_argument('--query',required=True);p.add_argument('--seconds',type=float,default=120)
    a=p.parse_args();replay(a.archive,a.cut,a.seed,a.query,a.output,a.seconds)
