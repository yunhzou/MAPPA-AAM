"""Current public search/checkpoint and event decoder on unchanged supplied inputs."""
import argparse, dataclasses, gc, gzip, hashlib, json, os, resource, sys, time, traceback
from pathlib import Path


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--case', choices=['rp_timeout','search_memory'],required=True)
    p.add_argument('--inputs',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--memory-gib',type=float,default=5)
    args=p.parse_args(); out=args.output;out.mkdir(parents=True,exist_ok=False)
    if sys.platform=='linux':
        limit=int(args.memory_gib*1024**3)
        resource.setrlimit(resource.RLIMIT_AS,(limit,limit))
    import numpy as np
    import mappa
    from mappa import AAMProblem,AAMSearchConfig,AAMResult,MolecularEndpoint,search_aam_checkpoints
    from mappa.growth import native
    from mappa.artifacts import read_raw_cut,decoded_events_record
    from mappa.frag import build_graph,classify_bonds
    from mappa.search_symmetry import finalize_graph_symmetry
    from mappa.final_branches import FinalBranchCatalogue
    from mappa.postprocessing import decode_events,EventDecodeConfig
    assert native.available(), 'Built native engine is required; no Python substitution'
    assert not hasattr(mappa,'select_rp_mappings')
    started=time.perf_counter(); timings={}
    def emit(stage,**kw):
        row=dict(stage=stage,seconds=time.perf_counter()-started,**kw)
        print(json.dumps(row),flush=True)
        (out/'progress.json').write_text(json.dumps(row,indent=2))
    folders={'rp_timeout':'01_RP_timeout_PR22_gpt-6-luna','search_memory':'02_AAM_memory_PMC11894599_10_gpt-6-luna'}
    folder=args.inputs/folders[args.case];endpoints=[]; hashes={}
    for name in ['reference_arrays.npz','prediction_arrays.npz']:
        path=folder/name; hashes[name]=hashlib.sha256(path.read_bytes()).hexdigest()
        with np.load(path,allow_pickle=False) as d:
            endpoints.append(MolecularEndpoint(tuple(d['elements'].tolist()),d['coordinates'],d['wbo']))
    problem=AAMProblem(*endpoints,name=args.case)
    config=AAMSearchConfig(graph_floor=.15,cut_floor=.15,iso_tolerance=1.,event_threshold=1.,metal_event_threshold=1.,seed_count=3,branch_limit=100,sweep_cuts=args.case=='search_memory')
    decoder=EventDecodeConfig(threshold=1.,metal_threshold=1.,seconds_per_family=120)
    (out/'configuration.json').write_text(json.dumps(dict(case=args.case,atoms=problem.atom_count,inputs=hashes,search=dataclasses.asdict(config),decode=dataclasses.asdict(decoder),execution='reference',mappa_path=mappa.__file__,memory_gib=args.memory_gib,python=sys.version),indent=2))
    try:
        emit('search'); begin=time.perf_counter()
        search=search_aam_checkpoints(problem,config,intermediate_dir=out/'cuts',workers=1,execution='reference')
        timings['search']=time.perf_counter()-begin
        (out/'search.json').write_text(json.dumps(dict(metrics=dataclasses.asdict(search.metrics),capped=search.capped,cuts=len(search.cut_paths)),indent=2))
        emit('catalogue',cuts=len(search.cut_paths),capped=search.capped);begin=time.perf_counter()
        cat=FinalBranchCatalogue(problem)
        target=build_graph(problem.product.elements,problem.product.wbo,bond_cut=config.graph_floor)
        for i,path in enumerate(search.cut_paths):
            graph,_=finalize_graph_symmetry(read_raw_cut(path),target,iso_tolerance=config.iso_tolerance)
            cat.add_graph(graph,archive=path.name)
            del graph; gc.collect()
            emit('catalogue',cuts_done=i+1,cuts=len(search.cut_paths),families=len(cat.families))
        timings['catalogue']=time.perf_counter()-begin
        with gzip.open(out/'catalogue.json.gz','wt') as f:json.dump(cat.to_record(),f)
        emit('decode',families=len(cat.families),incomplete_paths=cat.incomplete_paths);begin=time.perf_counter()
        import mappa.postprocessing as postprocessing
        original_extract=postprocessing.extract_path_events
        family_number=0
        report_file=(out/'family-patterns.jsonl').open('x')
        def record_extract(*a,**kw):
            nonlocal family_number
            report=original_extract(*a,**kw)
            report_file.write(json.dumps(dict(family_id=family_number,**report))+'\n');report_file.flush()
            family_number+=1
            return report
        postprocessing.extract_path_events=record_extract
        try:
            decoded=decode_events(cat,decoder)
        finally:
            report_file.close();postprocessing.extract_path_events=original_extract
        timings['decode']=time.perf_counter()-begin
        timings['through_decoding']=time.perf_counter()-started
        record=decoded_events_record(decoded)
        with gzip.open(out/'decoded.json.gz','wt') as f:json.dump(record,f)
        emit('audit',candidates=len(decoded.candidates));begin=time.perf_counter()
        for c in decoded.candidates:
            assert set(c.mapping)==set(range(problem.atom_count))
            assert set(c.mapping.values())==set(range(problem.atom_count))
            assert all(problem.reactant.elements[a]==problem.product.elements[b] for a,b in c.mapping.items())
            broken,formed,*_=classify_bonds(c.mapping,problem.reactant.wbo,problem.product.wbo,dwbo_threshold=1.,elements_R=problem.reactant.elements,elements_P=problem.product.elements,metal_dwbo_threshold=1.)
            norm=lambda bonds:sorted(tuple(sorted(map(int,b))) for b in bonds)
            inverse={b:a for a,b in c.mapping.items()}
            assert norm((b[0],b[1]) for b in broken)==norm(c.events['broken'])
            assert norm((inverse[b[0]],inverse[b[1]]) for b in formed)==norm(c.events['formed'])
            assert c.total==len(broken)+len(formed)
        from audit_actions import verify,program
        wanted={c.id:c for c in decoded.candidates};checked=set();membership=[]
        with (out/'family-patterns.jsonl').open() as f:
            for line in f:
                row=json.loads(line);fid=row['family_id'];family=cat.families[fid];steps=None
                for pattern in row['patterns']:
                    c=wanted.get(pattern['id'])
                    if c is None or c.id in checked or fid!=c.family_ids[0]:continue
                    assert dict(enumerate(pattern['mapping']))==c.mapping
                    if pattern['actions'] and steps is None:steps=program(family,problem)
                    verify(family,problem,pattern,steps)
                    checked.add(c.id)
        assert checked==set(wanted)
        membership=[dict(candidate=c.id,status='allowed',reason='exact constructive certificate replay') for c in decoded.minimum_candidates]
        (out/'independent-audit.json').write_text(json.dumps(dict(status='passed',witnesses=len(checked),raw_wbo_rescore=True,constructive_membership=True,minimum_witnesses=len(membership)),indent=2))
        timings['audit']=time.perf_counter()-begin
        result=dict(status='complete' if decoded.complete else 'incomplete',case=args.case,atoms=problem.atom_count,search_capped=search.capped,cuts=len(search.cut_paths),families=len(cat.families),family_reports=len(decoded.family_reports),incomplete_families=[r for r in decoded.family_reports if not r['complete']],incomplete_paths=cat.incomplete_paths,candidates=len(decoded.candidates),minimum_events=decoded.observed_minimum,minimum_proven=decoded.minimum_proven,event_window=decoded.event_window,witness_rescore_passed=len(decoded.candidates),minimum_membership=membership,timings=timings,wall_seconds=time.perf_counter()-started,peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024 if sys.platform=='linux' else 1024**2),geometry_selection=False)
        (out/'result.json').write_text(json.dumps(result,indent=2));emit('finished',result=result)
    except Exception as exc:
        (out/'failure.json').write_text(json.dumps(dict(error=repr(exc),traceback=traceback.format_exc(),timings=timings),indent=2));raise

if __name__=='__main__':main()
