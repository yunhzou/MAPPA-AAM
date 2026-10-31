"""Replay exact archived failed tasks in fresh processes; leave the campaign untouched."""
import argparse, gc, importlib, json, multiprocessing as mp, os, sqlite3, tempfile, time, traceback
from pathlib import Path
import psutil
from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob
from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint


def child(task, pipe, directory):
    from mappa.mappa_core.core.scheduler.runtime import Workspace
    import mappa.mappa_core.core.scheduler.search as search
    from mappa.mappa_core.core.search.session import SearchSession
    def event(phase, **fields):
        pipe.send(dict(phase=phase, rss=psutil.Process().memory_info().rss, time=time.monotonic(), **fields))
    from mappa.mappa_core.core.matching.growth import native
    from mappa.mappa_core.core.search.graph import AAMSearchGraph
    native_search=native._engine.search_fragments
    def measured_native(*args,**kwargs):
        event('native_search_start')
        record=native_search(*args,**kwargs)
        event('native_record_exported',states=len(record['states']),transitions=len(record['transitions']),
              mapping_pairs=sum(len(s['mapping']) for s in record['states']),
              island_pairs=sum(len(s['islands']) for s in record['states']))
        return record
    native._engine.search_fragments=measured_native
    from_record=AAMSearchGraph.from_record
    def measured_record(cls,*args,**kwargs):
        event('python_graph_start')
        graph=from_record(*args,**kwargs)
        event('python_graph_end')
        return graph
    AAMSearchGraph.from_record=classmethod(measured_record)
    original_search=SearchSession.search_cut
    def measured_search(self,*args,**kwargs):
        event('growth_start')
        result=original_search(self,*args,**kwargs)
        event('growth_end',counts=result[1],terminals=len(result[0].terminals))
        return result
    SearchSession.search_cut=measured_search
    original_finalize=search._finalize
    def measured_finalize(*args,**kwargs):
        event('finalization_start')
        result=original_finalize(*args,**kwargs)
        event('finalization_end')
        return result
    search._finalize=measured_finalize
    workspace=Workspace(directory,directory+'/bundles')
    try:
        record=workspace.load(task.payload['input'])
        event('input_loaded',atoms=[len(record.problem.reactant.elements),len(record.problem.product.elements)],config=repr(record.config),cut=task.payload['cut'])
        del record
        module,name=task.operation.split(':')
        operation=getattr(importlib.import_module(module),name)
        for repetition in range(2):
            event('task_start',repetition=repetition)
            result=operation(task.payload,workspace)
            event('publication_start')
            blob=workspace.writer.put(result)
            event('publication_end',bytes=blob.size)
            del result
            event('output_released')
            gc.collect()
            event('after_gc')
        event('complete')
    except Exception:
        event('error',error=traceback.format_exc())
    finally:pipe.close()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('campaign');parser.add_argument('output')
    args=parser.parse_args();root=Path(args.campaign);out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    summaries=[]
    for bin_number,case,cut in [(0,'05735/forward',111),(1,'06711/forward',190)]:
        with sqlite3.connect(f'file:{root}/dfw160/{bin_number}/result/ledger.sqlite?mode=ro',uri=True) as db:
            row=db.execute('SELECT audit FROM results WHERE case_name=? AND stage=?',(case,'matching')).fetchone()
        audit=read_blob(Blob(**json.loads(row[0])),None)
        tasks=[r for r in audit['tasks'] if f'/cut/{cut:05d}/' in r[0]]
        taskrow,=tasks;task=_loads_checkpoint(taskrow[1]);original=_loads_checkpoint(taskrow[8])
        result=dict(task=task.key,original=original,events=[],samples=[])
        with tempfile.TemporaryDirectory(prefix='mappa-memory-') as temporary:
            parent,connection=mp.get_context('spawn').Pipe(False)
            process=mp.get_context('spawn').Process(target=child,args=(task,connection,temporary));process.start();connection.close()
            started=time.monotonic();result['started']=started;peak=0;phase='startup'
            while process.is_alive() or parent.poll():
                while parent.poll():
                    try:event=parent.recv()
                    except EOFError:break
                    result['events'].append(event);phase=event['phase'];print(task.key,event,flush=True)
                if process.is_alive():
                    try:rss=psutil.Process(process.pid).memory_info().rss
                    except psutil.NoSuchProcess:break
                    peak=max(peak,rss);result['samples'].append([time.monotonic()-started,phase,rss])
                    if rss>12*2**30 or time.monotonic()-started>180:
                        result['stopped']='12GiB guard or 180s';process.terminate();break
                else:break
                time.sleep(.05)
            process.join();result.update(peak_bytes=peak,exitcode=process.exitcode)
        (out/f'{bin_number}-{cut}.json').write_text(json.dumps(result,indent=2,default=str))
        summaries.append(dict(task=task.key,peak_bytes=peak,exitcode=process.exitcode,stopped=result.get('stopped')))
    (out/'summary.json').write_text(json.dumps(summaries,indent=2));print(summaries,flush=True)

if __name__=='__main__':main()
