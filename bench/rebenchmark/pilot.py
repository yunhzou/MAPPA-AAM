"""Small real-case gate before a full cluster campaign."""
from dataclasses import replace
from pathlib import Path
import tempfile
import time

from .evidence import read,save,publish,restore,sha
from .campaign import load_plan,verify_sources,environment
from .search import search,select_contexts
from .analysis import analyze_cut


def run(root,output):
    from mappa import search_aam_checkpoints
    from mappa.artifacts import read_raw_cut
    manifest=verify_sources(root)
    results=[]
    # Fixed small, medium and larger historical timing cases; no selection by outcome.
    for index in (64,604,17):
        for direction in ('R_to_P','P_to_R'):
            task=next(t for t in manifest['tasks'] if t['index']==index and t['direction']==direction)
            with tempfile.TemporaryDirectory(prefix='mappa-pilot-') as td:
                import shutil
                folder=Path(td)/'parallel';folder.mkdir()
                for name in ('input.json','reference.json','source.json'):
                    shutil.copy2(Path(root)/'inputs'/str(index)/name,folder/name)
                save(folder/'task.json',task)
                plan=load_plan(folder);start=time.perf_counter()
                parallel=search(plan.problem,plan.config,folder,8)
                serial=search_aam_checkpoints(plan.problem,plan.config,workers=1,
                    intermediate_dir=Path(td)/'serial',execution='reused_native')
                assert len(parallel.cut_paths)==len(serial.cut_paths)
                for a,b in zip(parallel.cut_paths,serial.cut_paths):assert read_raw_cut(a)==read_raw_cut(b)
                for n in (1,2,3):
                    prefix=search_aam_checkpoints(plan.problem,replace(plan.config,seed_count=n),
                        intermediate_dir=Path(td)/f'prefix{n}',workers=1,execution='reused_native')
                    for a,b in zip(parallel.cut_paths,prefix.cut_paths):
                        assert select_contexts(read_raw_cut(a),range(n))==read_raw_cut(b)
                # One cut through all ten seed evaluations exercises actual scoring.
                analyze_cut((str(folder),0))
                receipt=publish(folder,Path(output)/f'{index}-{direction}.tar')
                restored=Path(td)/'restored';restore(receipt['path'],restored)
                assert read(restored/'search.json')==read(folder/'search.json')
                results.append(dict(index=index,direction=direction,cuts=len(parallel.cut_paths),
                    serial_parallel_equal=True,seed_prefixes_equal=[1,2,3],
                    scoring_seeds=10,archive_verified=True,wall_seconds=time.perf_counter()-start,
                    search=read(folder/'search.json')))
                save(Path(output)/'progress.json',dict(completed=results,environment=environment()))
    save(Path(output)/'validation.json',dict(status='passed',cases=results,environment=environment(),
        manifest_sha256=sha(Path(root)/'manifest.json')))


if __name__=='__main__':
    import sys
    run(Path(sys.argv[1]),Path(sys.argv[2]))
