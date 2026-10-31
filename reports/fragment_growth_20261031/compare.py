"""Compare fragment policies in fresh processes; matching includes final compression."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
CASES = ('ring3', 'ring6', 'ring12', 'ring24', 'ring6_hydrogens')


def worker(case, policy):
    import gc
    import numpy as np
    import psutil
    from mappa import FragmentMatchConfig, match_fragment
    from mappa.mappa_core.core.matching.graphs import build_graph
    from mappa.mappa_core.core.matching.representation.records import FragmentMatchContext
    from mappa.mappa_core.core.matching.symmetry.orbits import _nauty_orbits
    from mappa.mappa_core.core.matching.growth import native

    size = int(case.split('_')[0][4:])
    hydrogens = case.endswith('_hydrogens')
    graphs = []
    for extra in ((1, 2) if hydrogens else (0, 0)):
        elements = ['C'] * size + ['H'] * (size * extra)
        weights = np.zeros((len(elements), len(elements)))
        edges = [(i, (i + 1) % size) for i in range(size)]
        edges += [(i, size + i * extra + h) for i in range(size) for h in range(extra)]
        for a, b in edges:
            weights[a, b] = weights[b, a] = 1.
        graphs.append(build_graph(elements, weights, bond_cut=.2))
    source, target = graphs
    context = FragmentMatchContext(source_orbits=_nauty_orbits(source, wbo_tol=1.),
                                   target_orbits=_nauty_orbits(target, wbo_tol=1.))
    config = FragmentMatchConfig(iso_tolerance=1., growth_policy=policy)
    source_view = native.source_graph(source)
    target_view = native.target_graph(target, context.target_orbits)
    process = psutil.Process()
    gc.collect()
    baseline_rss = process.memory_info().rss

    def match():
        return match_fragment(source, target, seed=0, context=context, config=config)

    started = time.perf_counter()
    result = match()
    cold_seconds = time.perf_counter() - started
    assert not result.capped and len(result.matches) == 1
    scientific_record = [(dict(m), sorted(m.fragment), m.preserved_bonds,
                          sorted(m.deferred_edges), m.symmetry) for m in result.matches]
    digest = hashlib.sha256(json.dumps(scientific_record, sort_keys=True).encode()).hexdigest()

    def measure(call):
        started = time.process_time()
        for _ in range(3):
            call()
        per_call = (time.process_time() - started) / 3
        repetitions = max(3, min(500, int(.06 / max(per_call, 1e-9))))
        cpu, wall = [], []
        for _ in range(7):
            wall_start, cpu_start = time.perf_counter(), time.process_time()
            for _ in range(repetitions):
                call()
            cpu.append((time.process_time() - cpu_start) / repetitions)
            wall.append((time.perf_counter() - wall_start) / repetitions)
        return dict(cpu_seconds=statistics.median(cpu), wall_seconds=statistics.median(wall),
                    cpu_batches=cpu, repetitions_per_batch=repetitions)

    public = measure(match)
    arguments = (source_view.graph, target_view.graph, 0, [-1] * len(source),
                 .2, 1., 1, 100, None, [], False)
    kernel = measure(lambda: native._engine.grow_island(*arguments, growth_policy=policy))
    profile = []
    match_fragment(source, target, seed=0, context=context, config=config, profile=profile)
    key_bytes, value_bytes = target_view.graph.memo_storage_bytes()
    hits, misses, entries = target_view.graph.memo_stats()
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != 'darwin':
        peak *= 1024
    return dict(case=case, policy=policy, atoms=[len(source), len(target)],
                public=public, native_kernel=kernel, cold_seconds=cold_seconds,
                baseline_rss_bytes=baseline_rss, peak_process_rss_bytes=peak,
                ending_rss_bytes=process.memory_info().rss,
                certificate_cache=dict(entries=entries, hits=hits, misses=misses,
                                       key_bytes=key_bytes, value_bytes=value_bytes),
                result_sha256=digest, profile=profile[0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', choices=CASES)
    parser.add_argument('--policy', choices=('incremental', 'existence'))
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('results.json'))
    parser.add_argument('--process-repeats', type=int, default=3)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args.worker, args.policy)), flush=True)
        return
    environment = {**os.environ, 'MAPPA_NATIVE': '1', 'OMP_NUM_THREADS': '1',
                   'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
                   'PYTHONPATH': str(ROOT / 'src')}
    rows = []
    for repeat in range(args.process_repeats):
        policies = ('incremental', 'existence') if repeat % 2 == 0 else ('existence', 'incremental')
        for case in CASES:
            for policy in policies:
                completed = subprocess.run([sys.executable, str(Path(__file__).resolve()),
                    '--worker', case, '--policy', policy], cwd=ROOT, env=environment,
                    capture_output=True, text=True, check=True)
                row = json.loads(completed.stdout)
                row['process_repeat'] = repeat
                rows.append(row)
                print(case, policy, 'CPU ms', round(row['public']['cpu_seconds'] * 1000, 4),
                      'peak MiB', round(row['peak_process_rss_bytes'] / 2**20, 2), flush=True)
    summaries = []
    for case in CASES:
        selected = [r for r in rows if r['case'] == case]
        assert len({r['result_sha256'] for r in selected}) == 1, (case, 'different compressed records')
        summary = dict(case=case, identical_compressed_records=True)
        for policy in ('incremental', 'existence'):
            subset = [r for r in selected if r['policy'] == policy]
            summary[policy] = dict(
                public_cpu_ms=1000*statistics.median(r['public']['cpu_seconds'] for r in subset),
                native_kernel_cpu_ms=1000*statistics.median(r['native_kernel']['cpu_seconds'] for r in subset),
                cold_ms=1000*statistics.median(r['cold_seconds'] for r in subset),
                peak_process_mib=statistics.median(r['peak_process_rss_bytes']/2**20 for r in subset),
                cached_payload_bytes=statistics.median(
                    r['certificate_cache']['key_bytes']+r['certificate_cache']['value_bytes'] for r in subset))
        summary['public_time_ratio'] = summary['existence']['public_cpu_ms']/summary['incremental']['public_cpu_ms']
        summaries.append(summary)
    from mappa import _engine
    report = dict(platform=platform.platform(), processor=platform.processor(),
        native_sha256=hashlib.sha256(Path(_engine.__file__).read_bytes()).hexdigest(),
        process_repeats=args.process_repeats, rows=rows, summaries=summaries,
        scope='Prepared graphs, matching through final compressed results; no decoding. Separate single-thread processes. '
              'Warm measurements retain the normal target cache. Cold measurement starts before the first match, '
              'after graph preparation. Peak RSS includes Python and imported libraries; cache payload excludes allocator overhead.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(summaries, indent=2), flush=True)


if __name__ == '__main__':
    main()
