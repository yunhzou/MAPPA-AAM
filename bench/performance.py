"""Isolated revision timings and bounded audits of saved decoding outputs."""
from __future__ import annotations

import argparse
import cProfile
from collections import Counter
from contextlib import closing
from dataclasses import asdict
import fcntl
import gc
import hashlib
import json
import os
from pathlib import Path
import pickle
import resource
import signal
import sqlite3
import subprocess
import sys
import tempfile
import time
import uuid


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


def saved_matching_worker(args, case):
    """Compare postprocessing revisions on identical immutable cut results."""
    from mappa import decode_events
    from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    from mappa.mappa_core.utils.io.records import decoded_events_record
    from golden import verify_catalogue, verify_decoded

    phases = {}
    def measure(name, operation):
        cpu, wall = time.process_time(), time.perf_counter()
        result = operation()
        value = phases.setdefault(name, dict(cpu_seconds=0., wall_seconds=0.))
        value['cpu_seconds'] += time.process_time() - cpu
        value['wall_seconds'] += time.perf_counter() - wall
        value['process_lifetime_peak_rss_bytes_after'] = peak_rss()
        return result

    with tempfile.TemporaryDirectory(prefix='mappa-reduction-comparison-') as scratch:
        matching = read_blob(Blob(**case['matching']), scratch)
        source = read_blob(matching.input, scratch)
        if json.loads(json.dumps(asdict(source.config))) != case['config']:
            raise ValueError('Saved matching configuration differs')
        reference = read_blob(Blob(**case['reference']), scratch)
        catalogue = FinalBranchCatalogue(source.problem)
        for number, cut in enumerate(matching.cuts):
            graph = measure('loading', lambda: read_blob(cut, scratch)).graph
            measure('reduction', lambda: catalogue.add_graph(graph, archive=f'cut/{number}'))
            del graph
        emit(dict(kind='stage', stage='reduction', state='finished', **phases['reduction']))
        verified = measure('verification', lambda: verify_catalogue(catalogue, reference,
                                                                   reverse=case['reversed']))
        def prepare_decoding():
            decoded = decode_events(catalogue, **case['decode_config'])
            _ = decoded.window
            return decoded
        decoded = measure('decoding', prepare_decoding)
        selected = measure('decoded_reference', lambda: verify_decoded(decoded, reference,
                                                                      reverse=case['reversed']))
        record = measure('export', lambda: decoded_events_record(decoded))
        binary = measure('serialization', lambda: pickle.dumps(record, protocol=5))
        restored = measure('restoration', lambda: pickle.loads(binary))
        if restored != record:
            raise ValueError('Decoded output changed during transport')
        catalogue_digest = None
        if args.audit_output:
            scientific = dict(record['catalogue'])
            scientific.pop('schema')  # Storage format may change without changing scientific contents.
            scientific.pop('histories', None)
            scientific['families'] = [dict(saved, provenance=list(family.provenance))
                for saved, family in zip(scientific['families'], catalogue.families)]
            catalogue_digest = canonical_digest(scientific)[0]
        emit(dict(kind='result', phases=phases, matching_repeated=False,
            cuts=len(matching.scope), seed_searches=sum(len(c.seeds) for c in matching.scope),
            branches=catalogue.branch_count, families=len(catalogue.families),
            fragments=len(decoded.inventory.fragments), paths=catalogue.path_count,
            branch_composition_sha256=canonical_digest(sorted(catalogue.branches))[0],
            catalogue_reference=verified['status'], decoded_reference=selected,
            catalogue_content_sha256=catalogue_digest,
            minimum_events=decoded.minimum_events, minimum_proven=decoded.minimum_proven,
            serialized_bytes=len(binary), process_lifetime_peak_rss_bytes=peak_rss()))


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
    if 'matching' in case:
        profiler = cProfile.Profile() if args.profile else None
        if profiler:
            profiler.enable()
        try:
            saved_matching_worker(args, case)
        finally:
            if profiler:
                profiler.disable()
                profiler.dump_stats(str(args.profile))
        if source_identity(root) != identity:
            raise ValueError('Source changed during saved-matching comparison')
        return
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
                if run['status'] != 'ok':
                    return 1
            if 'matching' in case:
                paired = {run['revision']: next(r for r in run['records'] if r['kind'] == 'result')
                          for run in report['runs'][-2:]}
                before, after = paired['baseline'], paired['candidate']
                regression = (before['branch_composition_sha256'] != after['branch_composition_sha256']
                    or before['catalogue_content_sha256'] != after['catalogue_content_sha256']
                    or (before['catalogue_reference'] == 'recovered' and after['catalogue_reference'] != 'recovered')
                    or (before['decoded_reference']['status'] == 'recovered'
                        and after['decoded_reference']['status'] != 'recovered')
                    or any(after['phases'][phase]['cpu_seconds'] > max(.1, 1.25 * before['phases'][phase]['cpu_seconds'])
                           for phase in ('reduction', 'decoding')))
                if regression:
                    report['stopped_for_investigation'] = case['name']
                    save(args.output, report)
                    return 1
    return int(any(run['status'] != 'ok' for run in report['runs']))


def canonical_digest(value):
    """Hash every JSON field without constructing a second serialized record."""
    digest = hashlib.sha256()
    size = 0
    encoder = json.JSONEncoder(sort_keys=True, separators=(',', ':'), allow_nan=False)
    for chunk in encoder.iterencode(value):
        raw = chunk.encode()
        digest.update(raw)
        size += len(raw)
    return digest.hexdigest(), size


def audit_decoding_direction(payload, workspace):
    """Release one complete decoded artifact before loading the other revision."""
    from mappa.mappa_core.core.scheduler.bundles import Blob
    from mappa.mappa_core.core.scheduler.results import DecodingResult, MatchingResult

    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != payload['runner_sha256']:
        raise ValueError('Decoding audit implementation changed')
    matching = workspace.load(Blob(**payload['matching']))
    if not isinstance(matching, MatchingResult) or matching.case != payload['case']:
        raise ValueError('Matching scope belongs to another case')
    scope = [asdict(item) for item in matching.scope]
    del matching
    result = dict(case=payload['case'], baseline_state=payload['baseline_state'],
                  candidate_state=payload['candidate_state'])
    for arm in ('baseline', 'candidate'):
        reference = payload[arm]
        if reference is None:
            result[arm] = None
            continue
        wall, cpu = time.perf_counter(), time.process_time()
        saved = workspace.load(Blob(**reference))
        loaded = dict(wall_seconds=time.perf_counter()-wall, cpu_seconds=time.process_time()-cpu)
        if not isinstance(saved, DecodingResult):
            raise TypeError('Decoding audit requires a DecodingResult')
        metadata = dict(case=saved.case, config=asdict(saved.config),
            scope=[asdict(item) for item in saved.scope], stops=saved.stops,
            search_config=asdict(saved.search_config), execution=saved.execution)
        if (saved.case != payload['case'] or metadata['scope'] != scope
                or canonical_digest(metadata['search_config'])[0] != canonical_digest(payload['search_config'])[0]
                or metadata['config'] != payload['decode_config']):
            raise ValueError('Decoded output differs from its declared input or search scope')
        wall, cpu = time.perf_counter(), time.process_time()
        record_hash, size = canonical_digest(saved.record)
        metadata_hash, _ = canonical_digest(metadata)
        result[arm] = dict(record_sha256=record_hash, metadata_sha256=metadata_hash,
            canonical_record_bytes=size, schema=saved.record.get('schema'),
            cuts=len(scope), seed_searches=sum(len(item['seeds']) for item in scope),
            loading=loaded, hashing=dict(wall_seconds=time.perf_counter()-wall,
                                         cpu_seconds=time.process_time()-cpu))
        del saved, metadata
        gc.collect()
    if result['candidate'] is None:
        result['status'] = 'candidate_error'
    elif result['baseline'] is None:
        result['status'] = 'previously_missing'
    else:
        fields = ('record_sha256', 'metadata_sha256')
        result['status'] = ('equal' if all(result['baseline'][key] == result['candidate'][key]
                                          for key in fields) else 'different')
    return result


def decoding_campaign(path, *, candidate):
    """Read immutable campaign provenance and terminal ledger rows, never payloads."""
    from enzymemap import decoding_policy, require_drained

    path = Path(path).resolve()
    files = {}
    def read(relative):
        source = path/relative
        raw = source.read_bytes()
        files[str(source)] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)
    manifest = read('manifest.json')
    launch = read('launch.json')
    identity = launch['source_identity']
    cases = {}
    for batch in manifest['batches']:
        directory = str(batch['batch'])
        index = read(f'{directory}/index.json')
        plan = read(f'{directory}/result/plan.json')
        digest = files[str(path/directory/'index.json')]
        if (digest != batch['index_sha256'] or digest != plan['input_index_sha256']
                or index['config'] != manifest['config'] or plan['config'] != manifest['config']):
            raise ValueError('Campaign input hashes or scientific settings differ')
        if not candidate and (decoding_policy(index) != decoding_policy(manifest)
                              or decoding_policy(plan) != decoding_policy(manifest)):
            raise ValueError('Baseline decoding policies differ')
        if candidate:
            progress = read(f'{directory}/result/progress.json')
            if (plan.get('mode') != 'decode_only'
                    or index.get('decoding_source') != manifest.get('decoding_source')
                    or plan.get('decoding_source') != manifest.get('decoding_source')):
                raise ValueError('Candidate must be a fresh decode-only campaign with shared provenance')
            if plan.get('decoder_source') != identity or progress['source'] != identity:
                raise ValueError('Candidate decoder identities differ')
            if (index.get('decode_config') != manifest.get('decode_config')
                    or plan.get('decode_config') != manifest.get('decode_config')):
                raise ValueError('Candidate decoding policies differ')
        output = path/directory/'result'
        require_drained(output)
        ledger = output/'ledger.sqlite'
        files[str(ledger)] = hashlib.sha256(ledger.read_bytes()).hexdigest()
        with closing(sqlite3.connect(f'{ledger.as_uri()}?mode=ro', uri=True)) as db:
            records = {(row[0], row[1]): row for row in db.execute('SELECT * FROM results')}
        names = {f"{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}"
                 for row in index['directions']}
        if (len(names) != len(index['directions']) or set(cases).intersection(names)
                or {name for name, stage in records if stage == 'decoding'} != names):
            raise ValueError('Campaign has duplicate or missing decoding directions')
        for row in index['directions']:
            name = f"{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}"
            decoded = records[name, 'decoding']
            if (decoded[2] not in ('complete', 'error')
                    or (decoded[2] == 'complete') != bool(decoded[3])):
                raise ValueError('Audit requires terminal decoding and references for completed outputs')
            if not candidate and any(not records.get((name, stage))
                    or records[name, stage][2] != 'complete' or not records[name, stage][3]
                    for stage in ('matching', 'deduplication')):
                raise ValueError('Baseline matching and reduction must be complete')
            cases[name] = dict(input=row, state=decoded[2],
                reference=json.loads(decoded[3]) if decoded[3] else None,
                matching=None if candidate else json.loads(records[name, 'matching'][3]),
                catalogue=row['catalogue'] if candidate else json.loads(records[name, 'deduplication'][3]))
    identities = {(item['input']['case'], item['input']['reversed']) for item in cases.values()}
    if (len(cases) != manifest['directions'] or len(identities) != len(cases)
            or len({case for case, reverse in identities}) != manifest['cases']
            or any((case, not reverse) not in identities for case, reverse in identities)):
        raise ValueError('Campaign manifest does not cover both directions of every reaction')
    return manifest, cases, files, dict(source_identity=identity, source_commit=launch['source_commit'])


def decoding_audit_plan(baseline, candidate, *, memory_bytes, task_seconds):
    from enzymemap import decoding_policy
    from mappa.execution.identity import code_identity

    old, before, old_files, old_source = decoding_campaign(baseline, candidate=False)
    new, after, new_files, new_source = decoding_campaign(candidate, candidate=True)
    if old['config'] != new['config'] or set(before) != set(after):
        raise ValueError('Decoding campaigns have different scientific scope')
    decode_config = asdict(decoding_policy(old))
    if new.get('decode_config') != decode_config:
        raise ValueError('Decoding configuration changed')
    origin = new['decoding_source']
    baseline = Path(baseline).resolve()
    if (Path(origin['campaign']).resolve() != baseline
            or origin['manifest_sha256'] != old_files[str(baseline/'manifest.json')]):
        raise ValueError('Candidate provenance identifies another original campaign')
    source_files = {**old_files, **new_files}
    evidence = [(baseline/'launch.json', origin['launch_sha256'])]
    evidence.extend((Path(item['path']), item['sha256']) for item in origin['source_indexes'])
    for batch in origin['original_batches']:
        directory = baseline/str(batch['batch'])
        evidence.extend((directory/path, batch[key]) for path, key in (
            ('index.json', 'input_index_sha256'), ('result/plan.json', 'plan_sha256'),
            ('result/ledger.sqlite', 'ledger_sha256'), ('result/bundles/queue.sqlite', 'queue_sha256'),
            ('result/progress.json', 'progress_sha256')))
    if {item['batch'] for item in origin['original_batches']} != {item['batch'] for item in old['batches']}:
        raise ValueError('Candidate provenance omits an original batch')
    for path, digest in evidence:
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Original source changed after fresh decoding was prepared')
        source_files[str(path)] = digest
    runner = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    rows = []
    for name in sorted(before):
        left, right = before[name], after[name]
        if (left['catalogue'] != right['catalogue']
                or any(right['input'].get(key) != value for key, value in left['input'].items())):
            raise ValueError('Candidate catalogue or original input lineage differs')
        rows.append(dict(case=name, baseline=left['reference'], candidate=right['reference'],
            baseline_state=left['state'], candidate_state=right['state'], matching=left['matching'],
            search_config=old['config'], decode_config=decode_config, runner_sha256=runner))
    return dict(schema='mappa.decoding_audit/v1', baseline=str(Path(baseline).resolve()),
        candidate=str(Path(candidate).resolve()), source_files=source_files,
        sources=dict(baseline=old_source, candidate=new_source),
        code_identity=code_identity(), runner_sha256=runner, memory_bytes=memory_bytes,
        task_seconds=task_seconds, directions=len(rows),
        baseline_states=dict(Counter(row['baseline_state'] for row in rows)),
        candidate_states=dict(Counter(row['candidate_state'] for row in rows)), rows=rows)


def audit_decoding(baseline, candidate, output, *, workers=1, memory_bytes=8*2**30,
                   task_seconds=1800., run_seconds=9600., resume=False, slurm=False):
    """Audit one bounded queue; resume only unsubmitted work after a clean drain."""
    from enzymemap import Ledger, require_drained
    from mappa.execution import BundleStore, WorkPool
    from mappa.mappa_core.core.scheduler import Task
    from mappa.mappa_core.core.scheduler.bundles import Blob

    if min(workers, memory_bytes, task_seconds, run_seconds) <= 0:
        raise ValueError('Workers, memory and time limits must be positive')
    output = Path(output).resolve()
    if any(output.is_relative_to(Path(path).resolve()) or Path(path).resolve().is_relative_to(output)
           for path in (baseline, candidate)):
        raise ValueError('Audit output must be separate from both source campaigns')
    plan = decoding_audit_plan(baseline, candidate, memory_bytes=memory_bytes, task_seconds=task_seconds)
    output.mkdir(parents=True, exist_ok=True)
    with (output/'.owner.lock').open('a+b') as owner:
        fcntl.flock(owner, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if resume:
            if json.loads((output/'plan.json').read_text()) != plan:
                raise ValueError('Audit sources, implementation or budgets changed')
            if not json.loads((output/'progress.json').read_text()).get('drained'):
                raise ValueError('Interrupted audit ownership requires explicit recovery')
            require_drained(output)
        elif (output/'plan.json').exists():
            raise ValueError('Existing audit requires --resume')
        else:
            save(output/'plan.json', plan)
        started, stop = time.monotonic(), False
        def drain(signum, frame):
            nonlocal stop
            stop = True
        handlers = {sig: signal.signal(sig, drain) for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGUSR1)}
        try:
            with tempfile.TemporaryDirectory(prefix='mappa-decode-audit-', dir=os.environ.get('SLURM_TMPDIR')) as scratch:
                ledger = Ledger(output, scratch, resume=resume)
                store = BundleStore(output/'bundles', scratch=scratch)
                try:
                    factory = WorkPool.slurm if slurm else WorkPool.local
                    with factory(workers=workers, store=store, memory_bytes=memory_bytes,
                                 case_seconds=task_seconds, shared_cache=False) as pool:
                        status = dict(drained=False, expected_directions=plan['directions'],
                            ready_workers=len(pool.queue.ready_workers()), startup_seconds=pool.startup_elapsed_seconds)
                        def checkpoint(active):
                            ledger.checkpoint(); pool.checkpoint()
                            status.update(active=len(active), completed=ledger.counts(),
                                elapsed_seconds=time.monotonic()-started)
                            save(output/'progress.json', status)
                        pending, exhausted, active = iter(plan['rows']), False, {}
                        checkpoint(active); next_checkpoint = time.monotonic()+30
                        while active or not exhausted:
                            stop |= time.monotonic()-started >= run_seconds
                            while not stop and not exhausted and len(active) < 2*status['ready_workers']:
                                row = next(pending, None)
                                if row is None:
                                    exhausted = True
                                    break
                                if ledger.get(row['case'], 'audit') is not None:
                                    continue
                                key = f"{uuid.uuid4().hex}/{row['case']}/audit"
                                pool.queue.add([Task(key, 'performance:audit_decoding_direction', row,
                                    case=key, case_seconds=task_seconds, memory_bytes=memory_bytes)])
                                active[row['case']] = key
                            finished = pool.queue.results(list(active.values()))
                            for name, key in list(active.items()):
                                result = finished[key]
                                if result['state'] in ('pending', 'running'):
                                    continue
                                measured = result['result'] or {}
                                complete = result['state'] == 'complete'
                                audit = pool.queue.archive_case(key, store)
                                ledger.put(name, 'audit', 'complete' if complete else 'error',
                                    measured['blob'] if complete else None, audit,
                                    None if complete else dict(state=result['state'], error=measured.get('error')))
                                del active[name]
                            if time.monotonic() >= next_checkpoint:
                                checkpoint(active); next_checkpoint = time.monotonic()+30
                            if stop and not active:
                                break
                            if active:
                                pool.check_workers(); time.sleep(.1)
                        status['drained'] = True
                        checkpoint(active)
                    records = []
                    for name, state, reference, error in ledger.db.execute(
                            'SELECT case_name,state,reference,error FROM results ORDER BY case_name'):
                        records.append(store.load(Blob(**json.loads(reference))) if state == 'complete'
                                       else dict(case=name, status='audit_error', error=json.loads(error)))
                    counts = dict(Counter(row['status'] for row in records))
                    changed = [path for path, digest in plan['source_files'].items()
                               if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest]
                    report = dict(schema=plan['schema'], expected_directions=plan['directions'],
                        sources=plan['sources'], auditor_identity=plan['code_identity'],
                        baseline_states=plan['baseline_states'], candidate_states=plan['candidate_states'],
                        terminal_complete=len(records) == plan['directions'], status_counts=counts,
                        source_files_unchanged=not changed, changed_source_files=changed, records=records)
                    report['complete'] = (report['terminal_complete'] and not changed
                        and not any(counts.get(key) for key in ('audit_error', 'candidate_error', 'different')))
                    save(output/'summary.json', report)
                    emit({key: report[key] for key in ('terminal_complete', 'complete', 'status_counts')})
                    return report
                finally:
                    ledger.db.close()
        finally:
            for sig, handler in handlers.items():
                signal.signal(sig, handler)


def decoding_audit_main(argv):
    parser = argparse.ArgumentParser(description='Compare complete saved decoded outputs outside decoding timings')
    parser.add_argument('--baseline-campaign', required=True, type=Path)
    parser.add_argument('--candidate-campaign', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--memory-gib', type=float, default=8.)
    parser.add_argument('--task-seconds', type=float, default=1800.)
    parser.add_argument('--run-seconds', type=float, default=9600.)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--slurm', action='store_true')
    args = parser.parse_args(argv)
    report = audit_decoding(args.baseline_campaign, args.candidate_campaign, args.output,
        workers=args.workers, memory_bytes=int(args.memory_gib*2**30), task_seconds=args.task_seconds,
        run_seconds=args.run_seconds, resume=args.resume, slurm=args.slurm)
    return int(not report['complete'])


def decoding_timings(campaign):
    """Summarize a drained decode-only pass, keeping failures in its denominator."""
    import numpy as np
    from mappa.mappa_core.core.scheduler.bundles import Blob, read_blobs
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint

    campaign = Path(campaign).resolve()
    manifest, cases, files, source = decoding_campaign(campaign, candidate=True)
    audits = {}
    for batch in manifest['batches']:
        ledger = campaign/str(batch['batch'])/'result/ledger.sqlite'
        with closing(sqlite3.connect(f'{ledger.as_uri()}?mode=ro', uri=True)) as db:
            for name, raw in db.execute("SELECT case_name,audit FROM results WHERE stage='decoding'"):
                audits[Blob(**json.loads(raw))] = name
            if db.execute('SELECT count(*) FROM retries').fetchone()[0]:
                raise ValueError('Fresh decoding timing summary requires no retries')
    if len(audits) != len(cases):
        raise ValueError('Every direction requires its own archived decoding audit')
    records = []
    for reference, audit in read_blobs(sorted(audits, key=lambda ref: (ref.path, ref.offset)), None):
        name = audits[reference]
        task, = audit['tasks']
        specification = _loads_checkpoint(task[1])
        if (specification.operation != 'mappa.mappa_core.core.scheduler.decoding:decode_catalogue'
                or asdict(specification.payload['catalogue']) != cases[name]['catalogue']):
            raise ValueError('Decoding audit used a different operation or catalogue')
        metrics = _loads_checkpoint(task[8])
        if cases[name]['state'] == 'complete':
            if task[4] != 'complete' or asdict(metrics['blob']) != cases[name]['reference']:
                raise ValueError('Completed ledger output differs from its archived task')
        record = dict(case=name, status=task[4], attempts=len(audit['attempts']),
            consumed_seconds=audit['case'][3], memory_guard_bytes=task[3])
        for field in ('compute_seconds', 'cpu_seconds', 'publish_seconds', 'publish_cpu_seconds',
                      'sampled_cpu_seconds', 'elapsed_seconds', 'peak_tree_rss_bytes', 'phase_at_exit'):
            record[field] = metrics.get(field)
        records.append(record)
    records.sort(key=lambda row: row['case'])

    def distribution(rows, field):
        values = [row[field] for row in rows]
        if not values:
            return None
        slowest = max(rows, key=lambda row: row[field])
        return dict(count=len(rows), total=sum(values),
            **dict(zip(('median', 'p90', 'p95', 'p99'), map(float, np.quantile(values, [.5, .9, .95, .99])))),
            maximum=max(values), maximum_case=slowest['case'])

    complete = [row for row in records if row['status'] == 'complete']
    grouped = {}
    for row in records:
        grouped.setdefault(row['case'].split('/')[0], []).append(row)
    reactions = [dict(case=name, **{field: sum(row[field] for row in pair)
                 for field in ('compute_seconds', 'cpu_seconds')})
                 for name, pair in grouped.items() if all(row['status'] == 'complete' for row in pair)]
    return dict(schema='mappa.decoding_timings/v1', source=source,
        expected_reactions=manifest['cases'], expected_directions=manifest['directions'],
        terminal_directions=len(records), status_counts=dict(Counter(row['status'] for row in records)),
        successful_reactions=len(reactions),
        timing_scope='Catalogue loading, decoder and record construction; publication reported separately',
        conversion=manifest['decoding_source'].get('conversion'), decode_config=manifest['decode_config'],
        direction={field: distribution(complete, field) for field in (
            'compute_seconds', 'cpu_seconds', 'publish_seconds', 'publish_cpu_seconds')},
        reaction={field: distribution(reactions, field) for field in ('compute_seconds', 'cpu_seconds')},
        source_files={str(Path(path).relative_to(campaign)): digest for path, digest in files.items()},
        records=records)


def main():
    if sys.argv[1:2] == ['audit-decoding']:
        return decoding_audit_main(sys.argv[2:])
    if sys.argv[1:2] == ['decoding-timings']:
        parser = argparse.ArgumentParser(description='Report full decoding CPU/wall distributions and failures')
        parser.add_argument('campaign', type=Path)
        parser.add_argument('output', type=Path)
        args = parser.parse_args(sys.argv[2:])
        report = decoding_timings(args.campaign)
        save(args.output, report)
        emit({key: report[key] for key in ('terminal_directions', 'status_counts', 'direction', 'reaction')})
        return 0
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
