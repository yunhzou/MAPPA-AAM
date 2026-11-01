"""Measure exact cut phases and retained memory without changing search policy."""
import argparse
from dataclasses import asdict
import gc
import gzip
import json
import multiprocessing as mp
from pathlib import Path
import pickle
import tempfile
import time
import traceback
import sys

import psutil


def replay(input_path, cuts, directory, connection, trace=False, instrumented_engine=None):
    try:
        if instrumented_engine:
            import importlib.util
            spec = importlib.util.spec_from_file_location('mappa._engine', instrumented_engine)
            module = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = module
            spec.loader.exec_module(module)
        from mappa import _engine
        from mappa.mappa_core.core.scheduler.runtime import Workspace
        from mappa.mappa_core.core.scheduler import search
        with gzip.open(input_path, 'rb') as stream:
            saved = pickle.load(stream)
        source = saved['source']
        current = None
        process = psutil.Process()
        if trace:
            import tracemalloc
            tracemalloc.start()
        def phase(name, **details):
            record = dict(phase=name, cut=current, rss=process.memory_info().rss,
                          pid=process.pid, time=time.perf_counter(),
                          cpu_seconds=time.process_time(), **details)
            if trace:
                record['python_bytes'], record['python_peak'] = tracemalloc.get_traced_memory()
            connection.send(record)
        with tempfile.TemporaryDirectory(prefix='mappa-memory-') as scratch:
            workspace = Workspace(scratch, str(Path(directory)/'bundles'))
            tasks, _ = search.plan_search('diagnostic', source.problem, source.config,
                                         workspace.writer, execution=source.execution)
            by_number = {task.payload['number']: task for task in tasks}
            original_native = _engine.search_fragments
            def native(*args, **kwargs):
                from mappa.mappa_core.core.search import graph as records
                originals = {name: getattr(records, name) for name in ('SearchState', 'FragmentTransition', 'SearchStop')}
                def boundary(name, stage):
                    original = originals[name]
                    first = True
                    def construct(*values, **options):
                        nonlocal first
                        if first:
                            first = False
                            phase(stage, symmetry_work=args[1].symmetry_stats())
                        return original(*values, **options)
                    return construct
                for name, stage in [('SearchState', 'growth_finished'), ('FragmentTransition', 'states_exported'),
                                    ('SearchStop', 'transitions_exported')]:
                    setattr(records, name, boundary(name, stage))
                phase('native_start')
                try:
                    result = original_native(*args, **kwargs)
                finally:
                    for name, value in originals.items():
                        setattr(records, name, value)
                phase('native_finished', states=len(result[0].states), transitions=len(result[0].transitions),
                      symmetry_work=args[1].symmetry_stats())
                if hasattr(args[1], 'memory_stats'):
                    phase('target_memory', **args[1].memory_stats())
                matches = [edge.match for edge in result[0].transitions if edge.match is not None]
                deferred = [match['deferred_edges'] for match in matches]
                phase('graph_fields', matches=len(matches),
                      fragments=sum(len(match['fragment']) for match in matches),
                      witness_entries=sum(len(match['symmetry']['witness']) for match in matches),
                      deferred_pairs=sum(map(len, deferred)),
                      deferred_list_bytes=sum(sys.getsizeof(rows)+sum(sys.getsizeof(pair) for pair in rows) for rows in deferred),
                      unique_deferred_values=len({tuple(map(tuple, rows)) for rows in deferred}))
                return result
            _engine.search_fragments = native
            original_finalize = search.finalize_graph_symmetry
            def finalize(*args, **kwargs):
                phase('finalization_start')
                result = original_finalize(*args, **kwargs)
                phase('finalization_finished')
                return result
            search.finalize_graph_symmetry = finalize
            phase('prepared')
            for current in cuts:
                phase('cut_start')
                result = search.search_complete_cut(by_number[current].payload, workspace)
                symmetry = workspace.resources['finalization'][1]
                phase('publication_start', repair=workspace.prepared_value.repair.stats(),
                      symmetry=symmetry.stats(), cache_entries=len(symmetry.cache),
                      generators=len(symmetry.generators), groups=len(symmetry.groups))
                reference = workspace.writer.put(result)
                phase('publication_finished', bundle_bytes=reference.size, reference=asdict(reference))
                del result, symmetry
                phase('output_released')
            phase('before_collect')
            collected = gc.collect()
            phase('after_collect', collected=collected)
            if trace:
                snapshot = tracemalloc.take_snapshot()
                phase('python_allocations', top=[str(s) for s in snapshot.statistics('lineno')[:20]])
            workspace.resources.clear()
            gc.collect()
            phase('symmetry_released')
            workspace.prepared_value = None
            workspace.cache.clear()
            gc.collect()
            phase('session_released')
            phase('complete')
    except Exception:
        connection.send(dict(phase='error', error=traceback.format_exc()))
    finally:
        connection.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input'); parser.add_argument('output')
    parser.add_argument('--cuts', type=int, nargs='+', required=True)
    parser.add_argument('--trace', action='store_true')
    parser.add_argument('--instrumented-engine')
    parser.add_argument('--guard-gib', type=float, default=8.)
    parser.add_argument('--seconds', type=float, default=900.)
    args = parser.parse_args()
    output = Path(args.output); output.mkdir(parents=True, exist_ok=True)
    context = mp.get_context('spawn')
    parent, connection = context.Pipe(False)
    process = context.Process(target=replay, args=(args.input, args.cuts, str(output), connection,
                                                  args.trace, args.instrumented_engine))
    record = dict(cuts=args.cuts, events=[], peaks={}, guard_gib=args.guard_gib,
                  seconds=args.seconds)
    process.start(); connection.close()
    started = time.monotonic()
    current = 'startup'
    while process.is_alive() or parent.poll():
        while parent.poll():
            try:
                event = parent.recv()
            except EOFError:
                break
            record['events'].append(event)
            current = f"{event.get('cut')}:{event['phase']}"
            record['peaks'][current] = event.get('rss', 0)
            print(json.dumps(event), flush=True)
        if not process.is_alive():
            break
        try:
            rss = psutil.Process(process.pid).memory_info().rss
        except psutil.NoSuchProcess:
            break
        record['peaks'][current] = max(record['peaks'].get(current, 0), rss)
        elapsed = time.monotonic()-started
        if rss > args.guard_gib*2**30 or elapsed > args.seconds:
            record['stopped'] = 'memory_guard' if rss > args.guard_gib*2**30 else 'time_guard'
            record['guard_sample'] = dict(rss=rss, elapsed_seconds=elapsed, phase=current)
            process.terminate()
            break
        time.sleep(.02)
    process.join()
    record['elapsed_seconds'] = time.monotonic()-started
    record['exitcode'] = process.exitcode
    (output/'measurements.json').write_text(json.dumps(record, indent=2)+'\n')
    if process.exitcode or record['events'][-1]['phase'] != 'complete':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
