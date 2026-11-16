"""Build a verified 3D algorithm replay from a reusable archive-selection manifest.

Usage: python tools/build_search_trajectory.py manifest.json output_directory
The capture subprocess can load a frozen engine selected by engine_source.
Rendering uses current shared viewer assets and requires no search or replay.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def saved_selection(selection, base):
    """Load a trusted checkpoint or a checksummed work-pool cut with its input."""
    from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob
    from mappa.mappa_core.core.scheduler.results import CutResult, SearchInput
    from mappa.mappa_core.domain import AAMResult, AAMSearchMetrics
    from mappa.presentation.trajectory import read_archive

    if ('archive' in selection) == ('saved_cut' in selection):
        raise ValueError('Select exactly one of archive or saved_cut')
    selected = dict(selection)
    if 'archive' in selected:
        archive = (base / selected.pop('archive')).resolve()
        result = read_archive(archive)
        provenance = dict(archive=str(archive), sha256=hashlib.sha256(archive.read_bytes()).hexdigest())
    else:
        references = {name: Blob(**{**record, 'path': str((base / record['path']).resolve())})
                      for name, record in selected.pop('saved_cut').items()}
        source = read_blob(references['input'], None)
        cut = read_blob(references['cut'], None)
        if not isinstance(source, SearchInput) or not isinstance(cut, CutResult):
            raise TypeError('saved_cut requires the original SearchInput and CutResult')
        metrics = AAMSearchMetrics.from_record(
            {'cuts': 1, 'raw_result_count': len(cut.graph.terminals)}, 0.0)
        result = AAMResult(source.problem, source.config, cut.graph, metrics)
        provenance = dict(saved_cut={name: vars(ref) for name, ref in references.items()},
                          cut=cut.cut, seeds=cut.seeds, metrics='Unmeasured archive adapter')
    for terminal, expected in selected.pop('expected_transitions', {}).items():
        if list(next(result.graph.paths(int(terminal))).transitions) != expected:
            raise ValueError(f'Selected history differs from saved provenance for terminal {terminal}')
    return {**selected, 'aam': result}, provenance


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifest',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--open',action='store_true',help='Open the finished HTML in the browser')
    p.add_argument('--geometry-cache',type=Path,help='Cache missing 3D display conformers')
    p.add_argument('--capture-only',action='store_true',help=argparse.SUPPRESS)
    a=p.parse_args();a.manifest=a.manifest.resolve();a.output=a.output.resolve()
    manifest=json.loads(a.manifest.read_text())
    resolve=lambda value:str((a.manifest.parent/value).resolve())
    a.output.mkdir(parents=True,exist_ok=True)
    if a.capture_only:
        engine=resolve(manifest['engine_source']) if manifest.get('engine_source') else str(ROOT/'src')
        sys.path.insert(0,engine)
        import mappa
        # Load this pipeline module against the selected engine's matcher APIs.
        import importlib.util
        spec=importlib.util.spec_from_file_location('mappa.presentation.trajectory',ROOT/'src/mappa/presentation/trajectory.py')
        module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
        loaded=[saved_selection(s,a.manifest.parent) for s in manifest['selections']]
        selections=[selection for selection,_ in loaded]
        started=time.perf_counter()
        document=module.build_trajectory(selections,title=manifest.get('title'),key_atoms=manifest.get('key_atoms'),
            watch_targets=manifest.get('watch_targets'),event_tolerance=manifest.get('event_tolerance',.5))
        for run,(_,provenance) in zip(document['runs'],loaded):
            run['provenance']=provenance
        document['engine_source']=str(Path(mappa.__file__).parent)
        document['requested_engine_commit']=manifest.get('engine_commit')
        (a.output/'trace.json').write_text(json.dumps(document,separators=(',',':'),default=lambda x:int(x))+'\n')
        (a.output/'validation.json').write_text(json.dumps(dict(status='passed',checks=document['checks'],
            capture_seconds=time.perf_counter()-started),indent=2)+'\n')
        print(f'Captured {len(document["checks"])} verified saved paths.',flush=True)
        return
    env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'}
    subprocess.run([sys.executable,str(Path(__file__).resolve()),str(a.manifest),str(a.output),'--capture-only'],
                   env=env,check=True)
    sys.path.insert(0,str(ROOT/'src'))
    from mappa.presentation.viewers import growth_trace_html
    document=json.loads((a.output/'trace.json').read_text())
    cache=a.geometry_cache or (Path(resolve(manifest['geometry_cache'])) if manifest.get('geometry_cache') else None)
    if cache is not None:
        from mappa.presentation.trajectory import with_display_geometry
        document=with_display_geometry(document,cache=cache)
        (a.output/'display-trace.json').write_text(json.dumps(document,separators=(',',':'))+'\n')
    started=time.perf_counter()
    output=a.output/'algorithm_trajectory.html'
    output.write_text(growth_trace_html(document))
    validation=json.loads((a.output/'validation.json').read_text())
    validation['render_seconds']=time.perf_counter()-started
    (a.output/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(output)
    if a.open:
        if sys.platform=='darwin':
            subprocess.run(['open',str(output)],check=True)
        else:
            import webbrowser
            if not webbrowser.open(output.as_uri()):
                raise RuntimeError('Browser did not open the trajectory')


if __name__=='__main__':main()
