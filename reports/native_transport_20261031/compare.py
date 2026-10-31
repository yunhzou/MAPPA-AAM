"""Fresh-process replay of archived failures, with full saved-graph comparison."""
import argparse
from dataclasses import asdict
import importlib
import json
import multiprocessing as mp
from pathlib import Path
import sqlite3
import tempfile
import time
import traceback

import psutil
from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob
from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint


def replay(task, directory, pipe):
    from mappa.mappa_core.core.scheduler.runtime import Workspace
    from mappa.mappa_core.core.search.session import SearchSession
    search_cut = SearchSession.search_cut
    def measured(self, *args, **kwargs):
        started = time.perf_counter()
        result = search_cut(self, *args, **kwargs)
        pipe.send(dict(phase='searched', seconds=time.perf_counter()-started))
        return result
    SearchSession.search_cut = measured
    try:
        with tempfile.TemporaryDirectory(prefix='mappa-compare-') as scratch:
            workspace = Workspace(scratch, str(Path(directory)/'bundles'))
            source = workspace.load(task.payload['input'])
            module, name = task.operation.split(':')
            operation = getattr(importlib.import_module(module), name)
            for repetition in range(2):
                started = time.perf_counter()
                pipe.send(dict(phase='started', repetition=repetition))
                result = operation(task.payload, workspace)
                computed = time.perf_counter()-started
                graph = result.graph
                reference = workspace.writer.put(result)
                pipe.send(dict(phase='published', repetition=repetition, compute_seconds=computed,
                    total_seconds=time.perf_counter()-started, reference=asdict(reference),
                    states=len(graph.states), transitions=len(graph.transitions),
                    terminals=len(graph.terminals), mapping_pairs=sum(len(s.mapping) for s in graph.states)))
                del graph, result
            pipe.send(dict(phase='complete'))
    except Exception:
        pipe.send(dict(phase='error', error=traceback.format_exc()))
    finally:
        pipe.close()


def archived_task(campaign, batch, case, cut):
    with sqlite3.connect(f'file:{campaign}/dfw160/{batch}/result/ledger.sqlite?mode=ro', uri=True) as db:
        audit, = db.execute('SELECT audit FROM results WHERE case_name=? AND stage=?', (case, 'matching')).fetchone()
    rows = read_blob(Blob(**json.loads(audit)), None)['tasks']
    row, = [row for row in rows if f'/cut/{cut:05d}/' in row[0]]
    return _loads_checkpoint(row[1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign')
    parser.add_argument('output')
    parser.add_argument('--baseline')
    args = parser.parse_args()
    output = Path(args.output); output.mkdir(parents=True, exist_ok=True)
    summaries = []
    for batch, case, cut in [(0, '05735/forward', 111), (1, '06711/forward', 190)]:
        task = archived_task(args.campaign, batch, case, cut)
        name = f'{batch}-{cut}'
        record = dict(task=task.key, events=[], peaks=[0, 0])
        parent, connection = mp.get_context('spawn').Pipe(False)
        process = mp.get_context('spawn').Process(target=replay, args=(task, str(output/name), connection))
        process.start(); connection.close()
        started = time.monotonic(); repetition = 0
        while process.is_alive() or parent.poll():
            while parent.poll():
                try:
                    event = parent.recv()
                except EOFError:
                    break
                record['events'].append(event)
                repetition = event.get('repetition', repetition)
                print(name, event, flush=True)
            if not process.is_alive():
                break
            try:
                rss = psutil.Process(process.pid).memory_info().rss
            except psutil.NoSuchProcess:
                break
            record['peaks'][repetition] = max(record['peaks'][repetition], rss)
            if rss > 12*2**30 or time.monotonic()-started > 600:
                record['stopped'] = '12 GiB or 600 seconds diagnostic guard'
                process.terminate()
                break
            time.sleep(.02)
        process.join(); record['exitcode'] = process.exitcode
        assert process.exitcode == 0 and record['events'][-1]['phase'] == 'complete', record
        if args.baseline:
            baseline = json.loads((Path(args.baseline)/f'{name}.json').read_text())
            old = [e for e in baseline['events'] if e['phase'] == 'published']
            new = [e for e in record['events'] if e['phase'] == 'published']
            for expected, actual in zip(old, new):
                before = read_blob(Blob(**expected['reference']), None)
                after = read_blob(Blob(**actual['reference']), None)
                assert before.graph == after.graph, name
                assert (before.cut, before.seeds) == (after.cut, after.seeds)
                del before, after
            record['full_graph_equal'] = True
        (output/f'{name}.json').write_text(json.dumps(record, indent=2)+'\n')
        summaries.append(record)
    (output/'summary.json').write_text(json.dumps(summaries, indent=2)+'\n')


if __name__ == '__main__':
    main()
