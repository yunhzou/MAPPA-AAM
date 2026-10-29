"""Read-only logging around the public checkpoint search and seed prefixes."""
from dataclasses import replace, asdict
import json
import os
from pathlib import Path
import resource
import time

from .evidence import save, sha


def select_contexts(graph, selected):
    """Reindex independent contexts without changing any represented path."""
    selected = tuple(selected)
    if len(set(selected)) != len(selected) or any(i not in range(len(graph.contexts)) for i in selected):
        raise ValueError('Invalid context selection')
    contexts = {old: new for new, old in enumerate(selected)}
    states = [s for s in graph.states if s.context in contexts]
    ids = {s.id: n for n, s in enumerate(states)}
    edges = [e for e in graph.transitions if e.source in ids]
    assert all(e.target in ids for e in edges), 'Contexts must be independent'
    return type(graph)(tuple(graph.contexts[i] for i in selected),
        tuple(ids[r] for r in graph.roots if r in ids),
        tuple(replace(s, id=ids[s.id], context=contexts[s.context]) for s in states),
        tuple(replace(e, id=n, source=ids[e.source], target=ids[e.target]) for n, e in enumerate(edges)),
        tuple(replace(s, state=ids[s.state]) for s in graph.stops if s.state in ids))


class SeedLog:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)

    def __enter__(self):
        from mappa.mappa_core.core.search import native as native_search
        self.original = native_search.find_islands_native
        original, directory = self.original, self.directory
        def observed(*args, **kwargs):
            cpu, wall = time.process_time(), time.perf_counter()
            graph = original(*args, **kwargs)
            cpu, wall = time.process_time()-cpu, time.perf_counter()-wall
            row = dict(context=asdict(graph.contexts[0]), cpu_seconds=cpu, wall_seconds=wall,
                pid=os.getpid(), states=len(graph.states), transitions=len(graph.transitions),
                terminals=len(graph.terminals), cap_stops=sum(s.reason=='capped' for s in graph.stops),
                process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
            with (directory / f'{os.getpid()}.jsonl').open('a') as stream:
                stream.write(json.dumps(row, separators=(',', ':'))+'\n')
            return graph
        native_search.find_islands_native = observed
        return self

    def __exit__(self, *exception):
        from mappa.mappa_core.core.search import native as native_search
        native_search.find_islands_native = self.original


def search(problem, config, folder, workers):
    from mappa import search_aam_checkpoints
    from publication_timing import SearchProfiler
    folder = Path(folder)
    resume = (folder/'cuts/manifest.json').exists()
    from mappa.artifacts import raw_cut_paths
    resumed_cuts = len(raw_cut_paths(folder/'cuts')) if resume else 0
    timing_directory, seed_directory = folder/'search-timing', folder/'seed-timing'
    if resume:
        attempt = len(list(timing_directory.glob('continuation-*'))) + 1
        timing_directory /= f'continuation-{attempt:04d}'
        seed_directory /= f'continuation-{attempt:04d}'
    with SearchProfiler(timing_directory) as profiler, SeedLog(seed_directory):
        result = search_aam_checkpoints(problem, config, workers=workers,
            execution='reused_native', intermediate_dir=folder/'cuts', resume=resume)
    save(folder/'search.json', dict(complete=True, metrics=asdict(result.metrics),
        capped=result.capped, timing=profiler.summary(),
        resumed_cuts=resumed_cuts,
        timing_directory=str(timing_directory.relative_to(folder)),
        metric_scope='Counts cover all cuts; growth peak and worker/checkpoint durations cover newly searched cuts.',
        cut_hashes={p.name:sha(p) for p in result.cut_paths},
        seed_timing_scope='Seed-call CPU with shared runtime caches; not standalone seed-count latency.'))
    return result
