"""Serial, isolated revision comparisons over identical prepared single-cut inputs."""
from __future__ import annotations

import argparse
import cProfile
import gc
import hashlib
import json
import os
from pathlib import Path
import pickle
import resource
import subprocess
import sys
import tempfile
import time


def source_identity(root):
    package = root / 'src' / 'mappa'
    paths = sorted(set(package.rglob('*.py')) | set(package.glob('_*.so')) |
                   set(package.glob('_*.pyd')))
    if not paths:
        raise ValueError(f'No MAPPA source in {root}')
    digest = hashlib.sha256()
    for path in paths:
        digest.update(str(path.relative_to(package)).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def emit(record):
    print(json.dumps(record), flush=True)


def peak_rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (
        1 if sys.platform == 'darwin' else 1024)


def worker(args, case):
    root = args.revision.resolve()
    if os.environ.get('PYTHONPATH') != str(root / 'src'):
        raise ValueError('Worker PYTHONPATH must contain only the selected revision/src')
    import mappa
    from mappa import AAMSearchConfig, decode_events, _engine, _group_ops
    from mappa.execution.identity import code_identity
    from mappa.mappa_core.core.matching.graphs import build_graph
    from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob
    from mappa.mappa_core.core.scheduler.results import SearchInput
    from mappa.mappa_core.core.scheduler.runtime import Workspace
    from mappa.mappa_core.core.scheduler.search import search_complete_cut
    from mappa.mappa_core.core.scheduler.storage import MemoryStore
    from mappa.mappa_core.core.search.planning import cut_seed, _generate_seed_orders
    from mappa.mappa_core.core.search.session import SearchSession
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    from mappa.mappa_core.utils.io.records import problem_from_record, decoded_events_record

    modules = {module.__name__: str(Path(module.__file__).resolve())
               for module in (mappa, _engine, _group_ops)}
    if any(not Path(path).is_relative_to(root / 'src') for path in modules.values()):
        raise ValueError(f'Imported a different MAPPA revision: {modules}')
    identity = code_identity()
    if identity != args.identity:
        raise ValueError('Selected source or native engine changed before execution')
    emit(dict(kind='identity', source_identity=identity, modules=modules))
    profiler = cProfile.Profile() if args.profile else None
    with tempfile.TemporaryDirectory(prefix='mappa-performance-') as scratch:
        problem = problem_from_record(read_blob(Blob(**case['input']), scratch))
        config = AAMSearchConfig(**case['config'])
        cut = tuple(tuple(edge) for edge in case['cut'])
        phases = {}

        def measure(name, operation):
            emit(dict(kind='stage', repeat=repeat, stage=name, state='started'))
            cpu, wall = time.process_time(), time.perf_counter()
            result = operation()
            phases[name] = dict(cpu_seconds=time.process_time()-cpu,
                                wall_seconds=time.perf_counter()-wall,
                                process_lifetime_peak_rss_bytes=peak_rss())
            emit(dict(kind='stage', repeat=repeat, stage=name, state='finished',
                      **phases[name]))
            return result

        def plan():
            store = MemoryStore()
            workspace = Workspace(scratch, store=store)
            reference = store.put(SearchInput(problem, config, 'reused_native'))
            source = build_graph(problem.reactant.elements, problem.reactant.wbo,
                                 bond_cut=config.graph_floor)
            source.remove_edges_from(cut)
            orders = tuple(map(tuple, _generate_seed_orders(source, config.seed_count,
                rng_seed=cut_seed(cut, config.random_seed), seed_selection=config.seed_selection)))
            payload = dict(input=reference, cut=cut, number=0, seeds=len(orders),
                           indices=tuple(range(len(orders))), orders=orders)
            return workspace, reference, payload

        repeat = 0
        if profiler:
            profiler.enable()
        workspace, reference, payload = measure('planning', plan)
        try:
            for repeat in range(args.warm_repeats + 1):
                gc.collect()
                session = measure('preparation', lambda: workspace.prepared(reference,
                    lambda: SearchSession(problem, config, 'reused_native')))
                cache_before = session.repair.stats()
                completed = measure('matching', lambda: search_complete_cut(payload, workspace))
                catalogue = measure('reduction', lambda:
                    FinalBranchCatalogue(problem).add_graph(completed.graph))
                decoded = measure('decoding', lambda: decode_events(catalogue,
                    threshold=config.event_threshold, metal_threshold=config.metal_event_threshold))
                record = measure('export', lambda: decoded_events_record(decoded))
                binary = measure('serialization', lambda: pickle.dumps(record, protocol=5))
                restored = measure('restoration', lambda: pickle.loads(binary))
                if restored != record:
                    raise ValueError('Decoded record changed during binary round trip')
                graph = completed.graph
                # Explicit scientific audit is outside the timing measurements.
                audit = (hashlib.sha256(json.dumps(record, sort_keys=True,
                    separators=(',', ':')).encode()).hexdigest() if args.audit_output else None)
                emit(dict(kind='result', repeat=repeat, cache_state='cold' if repeat == 0 else 'warm',
                    phases=phases, cpu_seconds=sum(p['cpu_seconds'] for p in phases.values()),
                    wall_seconds=sum(p['wall_seconds'] for p in phases.values()),
                    process_lifetime_peak_rss_bytes=peak_rss(),
                    states=len(graph.states), transitions=len(graph.transitions),
                    terminals=len(graph.terminals), stops=[vars(stop) for stop in graph.stops],
                    families=len(catalogue.families), fragments=len(decoded.inventory.fragments),
                    tables=len(record['tables']), responses=len(record['responses']),
                    output_schema=record['schema'], serialized_bytes=len(binary),
                    output_sha256=audit,
                    cache_before=cache_before, cache_after=session.repair.stats(),
                    search_counts=completed.search_counts))
                del completed, graph, catalogue, decoded, record, binary, restored
                phases = {}
                if profiler:
                    profiler.dump_stats(str(args.profile))
        finally:
            if profiler:
                profiler.disable()
                profiler.dump_stats(str(args.profile))
    loaded = [str(Path(module.__file__).resolve()) for name, module in sys.modules.items()
              if (name == 'mappa' or name.startswith('mappa.')) and module.__file__]
    if any(not Path(path).is_relative_to(root / 'src') for path in loaded):
        raise ValueError('Pipeline imported MAPPA modules outside the selected revision')
    if source_identity(root) != identity:
        raise ValueError('MAPPA source or engine changed during execution')


def save(path, report):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(report, indent=2) + '\n')
    temporary.replace(path)


def compare(args, cases):
    revisions = {name: dict(root=str(getattr(args, name).resolve()),
                           revision=getattr(args, name + '_revision'),
                           source_identity=source_identity(getattr(args, name).resolve()))
                 for name in ('baseline', 'candidate')}
    report = dict(schema='mappa.performance/1', revisions=revisions,
                  runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  manifest_sha256=hashlib.sha256(args.cases.read_bytes()).hexdigest(),
                  cases=cases, repetitions=args.repetitions, warm_repeats=args.warm_repeats,
                  timeout_seconds=args.timeout, profiled=args.profile_directory is not None, runs=[])
    report['audited_output'] = args.audit_output
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        raise FileExistsError('Use a new output file for each comparison')
    if args.profile_directory:
        args.profile_directory.mkdir(parents=True, exist_ok=True)
    save(args.output, report)
    for repetition in range(args.repetitions):
        order = ('baseline', 'candidate') if repetition % 2 == 0 else ('candidate', 'baseline')
        for index, case in enumerate(cases):
            for name in order:
                selected = revisions[name]
                command = [sys.executable, str(Path(__file__).resolve()), '--worker',
                    '--revision', selected['root'], '--identity', selected['source_identity'],
                    '--cases', str(args.cases.resolve()), '--case', case['name'],
                    '--warm-repeats', str(args.warm_repeats)]
                if args.audit_output:
                    command.append('--audit-output')
                if args.profile_directory:
                    command += ['--profile', str((args.profile_directory /
                        f'{repetition}-{index}-{name}.pstats').resolve())]
                environment = dict(os.environ, PYTHONPATH=str(Path(selected['root']) / 'src'),
                    MAPPA_NATIVE='1', PYTHONHASHSEED='0', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                    MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
                run = dict(revision=name, case=case['name'], repetition=repetition)
                started = time.perf_counter()
                try:
                    process = subprocess.run(command, env=environment, cwd=selected['root'],
                        capture_output=True, text=True, timeout=args.timeout)
                    run.update(status='ok' if process.returncode == 0 else 'error',
                        returncode=process.returncode, stdout=process.stdout, stderr=process.stderr)
                except subprocess.TimeoutExpired as error:
                    run.update(status='timeout', stdout=(error.stdout or b'').decode(),
                               stderr=(error.stderr or b'').decode())
                run['process_wall_seconds'] = time.perf_counter()-started
                run['records'] = []
                for line in run['stdout'].splitlines():
                    try:
                        run['records'].append(json.loads(line))
                    except json.JSONDecodeError:
                        run.setdefault('unparsed_stdout', []).append(line)
                        if run['status'] == 'ok':
                            run['status'] = 'invalid_output'
                report['runs'].append(run)
                save(args.output, report)
                emit({key: run[key] for key in ('revision', 'case', 'repetition', 'status')})
    return int(any(run['status'] != 'ok' for run in report['runs']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--candidate', type=Path)
    parser.add_argument('--baseline-revision', help='Explicit revision provenance for an exported source tree')
    parser.add_argument('--candidate-revision')
    parser.add_argument('--cases', type=Path, required=True)
    parser.add_argument('--case', action='append', help='Select named cases; repeat to select several')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--repetitions', type=int, default=3)
    parser.add_argument('--warm-repeats', type=int, default=1)
    parser.add_argument('--timeout', type=float, default=120.)
    parser.add_argument('--profile-directory', type=Path, help='Profile each worker; changes timing overhead')
    parser.add_argument('--audit-output', action='store_true',
                        help='Hash complete scientific records outside stage timings; adds memory use')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--revision', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--identity', help=argparse.SUPPRESS)
    parser.add_argument('--profile', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    cases = json.loads(args.cases.read_text())['cases']
    names = [case['name'] for case in cases]
    if len(set(names)) != len(names) or (args.case and set(args.case)-set(names)):
        parser.error('Case names must be unique and selected cases must exist')
    cases = [case for case in cases if args.case is None or case['name'] in args.case]
    if not cases or args.repetitions < 1 or args.warm_repeats < 0 or args.timeout <= 0:
        parser.error('Require cases, positive repetitions/timeout and nonnegative warm repeats')
    if any('self_tolerance' in case['config'] for case in cases):
        parser.error('Comparison settings must use the common configuration fields')
    if args.worker:
        if len(cases) != 1 or args.revision is None or args.identity is None:
            parser.error('Worker requires one case, --revision and --identity')
        worker(args, cases[0])
        return 0
    if any(value is None for value in (args.baseline, args.candidate, args.output)):
        parser.error('Comparison requires --baseline, --candidate and --output')
    return compare(args, cases)


if __name__ == '__main__':
    raise SystemExit(main())
