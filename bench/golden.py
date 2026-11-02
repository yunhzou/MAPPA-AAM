"""Golden preparation, reference verification and bounded local diagnostics."""
import argparse
from collections import Counter, defaultdict
from contextlib import contextmanager
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import platform
import signal
import sys
import time

from mappa import AAMSearchConfig, query_path
from mappa.mappa_core.core.scheduler.bundles import Blob, BundleWriter, read_blob
from mappa.mappa_core.utils.io.records import problem_from_record


REFERENCE_POLICY = 'exact heavy-atom relation including unmatched atoms, modulo endpoint chemical symmetry'


def endpoint_generators(feature, degree):
    """Lift chemically colored heavy-atom symmetries into original atom indexing."""
    import pynauty
    heavy = feature['heavy']
    colors, adjacency = defaultdict(set), defaultdict(set)
    for atom, color in enumerate(feature['colors']):
        colors[('atom', tuple(color))].add(atom)
    count = len(heavy)
    for left, right, color in feature['bonds']:
        colors[('bond', tuple(color))].add(count)
        adjacency[count].update((left, right))
        adjacency[left].add(count); adjacency[right].add(count)
        count += 1
    graph = pynauty.Graph(count, adjacency_dict={i: sorted(adjacency[i]) for i in range(count)},
        vertex_coloring=[colors[key] for key in sorted(colors, key=repr)])
    result = []
    for generator in pynauty.autgrp(graph)[0]:
        images = list(range(degree))
        for atom, original in enumerate(heavy):
            images[original] = heavy[generator[atom]]
        result.append(tuple(images))
    return tuple(result)


def verify(payload, workspace):
    saved = workspace.load(payload['catalogue'])
    reference = read_blob(payload['reference'], workspace.scratch)
    return verify_catalogue(saved.catalogue, reference, reverse=payload['reversed'])


def verify_catalogue(catalogue, reference, *, reverse=False):
    """The same full heavy-atom membership query for campaigns and local studies."""
    mapping = dict(reference['mapping'])
    features, generators = reference['features'], reference['generators']
    if reverse:
        mapping = {product: reactant for reactant, product in mapping.items()}
        features, generators = features[::-1], generators[::-1]
    unknown = False
    for number, family in enumerate(catalogue.families):
        status, evidence = query_path(family.as_path(catalogue.problem), catalogue.problem,
            mapping, source_atoms=tuple(features[0]['heavy']), source_generators=generators[0],
            target_generators=generators[1], complete_reference=True, timeout_ms=10000)
        if status == 'recovered':
            return dict(status=status, family=number, evidence=evidence,
                        heavy_atoms=len(mapping), policy=REFERENCE_POLICY)
        unknown |= status == 'unknown'
    return dict(status='unknown' if unknown else 'not_recovered',
                heavy_atoms=len(mapping), policy=REFERENCE_POLICY)


def prepare(inputs, manifest, output, *, batches=5):
    """Copy verified original endpoints into aggregated bundles; never rerun old code."""
    inputs, manifest, output = Path(inputs), Path(manifest), Path(output)
    records = json.loads(manifest.read_text())['records']
    if sorted(row['case'] for row in records) != list(range(1851)):
        raise ValueError('Golden must contain exactly the original 1,851 records')
    output.mkdir(parents=True, exist_ok=False)
    writer = BundleWriter(output/'inputs', scratch=None)
    rows, evidence = [], []
    for row in sorted(records, key=lambda row: row['case']):
        folder = inputs/str(row['case'])
        raw = (folder/'input.json').read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if digest != row['input_sha256']:
            raise ValueError(f'Golden input checksum differs for {row["case"]}')
        source = json.loads(raw)
        problem = problem_from_record(source)
        raw_reference = (folder/'reference.json').read_bytes()
        reference = json.loads(raw_reference)
        reference['generators'] = [endpoint_generators(feature, len(endpoint.elements))
            for feature, endpoint in zip(reference['features'], (problem.reactant, problem.product))]
        reference_blob = asdict(writer.put(reference))
        for reverse in (False, True):
            record = dict(source, reactant=source['product'], product=source['reactant']) if reverse else source
            rows.append(dict(case=row['case'], reversed=reverse,
                             input=asdict(writer.put(record)), reference=reference_blob))
        evidence.append(dict(case=row['case'], input_sha256=digest,
            reference_sha256=hashlib.sha256(raw_reference).hexdigest()))
    config = asdict(AAMSearchConfig(seed_count=1, branch_limit=100, iso_tolerance=1.,
                                   random_seed=42))
    for batch in range(batches):
        folder = output/str(batch); folder.mkdir()
        subset = [row for row in rows if row['case'] % batches == batch]
        (folder/'index.json').write_text(json.dumps(dict(config=config, directions=subset,
            verification='golden:verify', reference_policy=REFERENCE_POLICY), indent=2)+'\n')
    (output/'manifest.json').write_text(json.dumps(dict(cases=1851, directions=len(rows), batches=batches,
        input_manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest(),
        config=config, reference_policy=REFERENCE_POLICY, records=evidence), indent=2)+'\n')
    print(json.dumps(dict(cases=1851, directions=len(rows), batches=batches, output=str(output))))


def _publish(path, value):
    path = Path(path)
    pending = path.with_suffix(path.suffix+'.pending')
    pending.write_text(json.dumps(value, indent=2)+'\n')
    pending.replace(path)


@contextmanager
def _measure(metrics, phase):
    cpu, wall = time.process_time(), time.perf_counter()
    yield
    values = metrics.setdefault(phase, dict(cpu_seconds=0., wall_seconds=0.))
    values['cpu_seconds'] += time.process_time()-cpu
    values['wall_seconds'] += time.perf_counter()-wall


def _digest(value):
    encoded = json.dumps(value, sort_keys=True, separators=(',', ':')).encode()
    return hashlib.sha256(encoded).hexdigest()


def _catalogue_signature(catalogue):
    """Compare saved scientific records without expanding shared search histories."""
    families = [_digest((f.policy, f.mapping, f.required_edges, f.actions)) for f in catalogue.families]
    branches = [_digest((key, sorted(families[i] for i in ids)))
                for key, ids in catalogue.branches.items()]
    coupled = [_digest((families[i], sorted(partitions)))
               for i, partitions in catalogue.coupled_branches.items()]
    return _digest(('catalogue-signature/v2', sorted(families),
                    sorted(branches), sorted(coupled)))


def diagnose_direction(payload, workspace):
    """One cold process, one complete direction; retain its cache across cuts."""
    from mappa.mappa_core.core.matching.symmetry.workspace import ConditionedSymmetryWorkspace
    from mappa.mappa_core.core.search.branching.reuse import FragmentRepair
    from mappa.mappa_core.core.search.branching.symmetry import finalize_graph_symmetry
    from mappa.mappa_core.core.search.planning import cut_sweep_items
    from mappa.mappa_core.core.search.session import SearchSession
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue

    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != payload['runner_sha256']:
        raise ValueError('Diagnostic runner changed after the run was planned')
    row, metrics = payload['row'], {}
    with _measure(metrics, 'preparation'):
        problem = problem_from_record(read_blob(Blob(**row['input']), workspace.scratch))
        reference = read_blob(Blob(**row['reference']), workspace.scratch)
        config = AAMSearchConfig(**payload['config'])
        session = SearchSession(problem, config)
        if not payload['cache']:
            repair = session.repair
            session.repair = FragmentRepair(repair.source, repair.target, session.context[3], cache_bytes=0)
            del repair
        symmetry = ConditionedSymmetryWorkspace(session.context[2], config.iso_tolerance)
        cuts = cut_sweep_items(problem.reactant.wbo, config.cut_floor) if config.sweep_cuts else [()]
        catalogue = FinalBranchCatalogue(problem)
    result = dict(key=payload['key'], atoms=[problem.source_atom_count, problem.target_atom_count],
        expected_cuts=len(cuts), completed_cuts=0, expected_seed_searches=len(cuts)*config.seed_count,
        completed_seed_searches=0, metrics=metrics, cuts=[], stops={}, cap_stages={},
        graph_counts=dict(states=0, transitions=0, terminals=0), max_growth_candidates=0)
    scope, stops, cap_stages = hashlib.sha256(), Counter(), Counter()

    def progress(phase):
        result.update(phase=phase, cache=session.repair.stats(), stops=dict(stops),
                      cap_stages=dict(cap_stages))
        _publish(payload['progress'], result)

    for number, cut in enumerate(cuts):
        progress('matching')
        with _measure(metrics, 'matching'):
            graph, counts = session.search_cut(cut)
        progress('finalization')
        with _measure(metrics, 'finalization'):
            graph, _ = finalize_graph_symmetry(graph, session.context[2],
                iso_tolerance=config.iso_tolerance, workspace=symmetry)
        progress('reduction')
        with _measure(metrics, 'reduction'):
            catalogue.add_graph(graph, archive=f'cut/{number}')
        with _measure(metrics, 'audit'):
            scope.update(json.dumps([cut, [context.seed_order for context in graph.contexts]],
                                    separators=(',', ':')).encode())
            cut_stops = Counter()
            for stop in graph.stops:
                stops[stop.reason] += stop.occurrences
                cut_stops[stop.reason] += stop.occurrences
                if stop.reason == 'capped':
                    cap_stages[stop.stage] += stop.occurrences
            result['cuts'].append(dict(number=number, bonds=cut, states=len(graph.states),
                transitions=len(graph.transitions), terminals=len(graph.terminals), stops=dict(cut_stops)))
            result['graph_counts']['states'] += len(graph.states)
            result['graph_counts']['transitions'] += len(graph.transitions)
            result['graph_counts']['terminals'] += len(graph.terminals)
            result['completed_seed_searches'] += len(graph.contexts)
            result['max_growth_candidates'] = max(result['max_growth_candidates'], counts['max_growth_candidates'])
        result['completed_cuts'] = number+1
        del graph
    result['catalogue_audited'] = payload['audit_catalogue']
    if payload['audit_catalogue']:
        progress('audit')
        with _measure(metrics, 'audit'):
            result['catalogue_sha256'] = _catalogue_signature(catalogue)
    result.update(scope_sha256=scope.hexdigest(), families=len(catalogue.families),
                  atom_set_branches=len(catalogue.branches), retained_paths=catalogue.path_count)
    progress('verification')
    with _measure(metrics, 'verification'):
        result['verification'] = verify_catalogue(catalogue, reference, reverse=row['reversed'])
    progress('complete')
    return result


def diagnose(indices, output, *, cases, seeds, cache, branch_limit, direction_seconds,
             run_seconds, memory_bytes, cpu_model=None, directions=('forward', 'reverse'),
             audit_catalogue=False):
    """Run only explicitly selected cases, with supervised budgets and cold controls."""
    from mappa.execution import WorkPool
    from mappa.execution.identity import code_identity
    from mappa.mappa_core.core.scheduler import Task

    if not cases or len(set(cases)) != len(cases) or min(cases) < 0:
        raise ValueError('Select distinct nonnegative Golden case IDs')
    if not seeds or len(set(seeds)) != len(seeds) or min(seeds) < 1:
        raise ValueError('Seed counts must be distinct positive integers')
    if len(set(cache)) != len(cache) or not set(cache) <= {'on', 'off'} or not cache:
        raise ValueError('Cache settings must be on and/or off, without repetitions')
    if not directions or len(set(directions)) != len(directions) or not set(directions) <= {'forward', 'reverse'}:
        raise ValueError('Directions must be forward and/or reverse, without repetitions')
    if min(branch_limit, direction_seconds, run_seconds, memory_bytes) <= 0:
        raise ValueError('The cap, time budgets and memory guard must be positive')
    audit_catalogue = bool(audit_catalogue or len(cache) > 1)
    selected_directions = {direction == 'reverse' for direction in directions}
    indexes = [json.loads(Path(path).read_text()) for path in indices]
    if any(index['config'] != indexes[0]['config'] or index.get('verification') != 'golden:verify'
           for index in indexes):
        raise ValueError('Use Golden indexes with the same scientific configuration and verifier')
    rows = [row for index in indexes for row in index['directions']
            if row['case'] in cases and row['reversed'] in selected_directions]
    wanted = {(case, reverse) for case in cases for reverse in selected_directions}
    if len(rows) != len(wanted) or {(r['case'], r['reversed']) for r in rows} != wanted:
        raise ValueError('Indexes must supply each selected case/direction exactly once')
    rows.sort(key=lambda row: (cases.index(row['case']), row['reversed']))
    output = Path(output).resolve(); output.mkdir(parents=True, exist_ok=False)
    config = AAMSearchConfig(**indexes[0]['config'])
    report = dict(schema='mappa.golden_local_diagnostic/v1', code_identity=code_identity(),
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        platform=dict(system=platform.system(), release=platform.release(), machine=platform.machine(),
                      cpu_model=cpu_model),
        cases=cases, directions=list(directions), seeds=seeds, cache=cache,
        config=asdict(replace(config, branch_limit=branch_limit)),
        direction_seconds=direction_seconds, run_seconds=run_seconds, memory_guard_bytes=memory_bytes,
        reference_policy=REFERENCE_POLICY, complete=False, results=[], comparisons=[],
        catalogue_signature='catalogue-signature/v2' if audit_catalogue else None,
        input_indexes_sha256=[hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in indices],
        input_sha256=[dict(case=r['case'], reversed=r['reversed'], input=r['input']['sha256'],
                           reference=r['reference']['sha256']) for r in rows],
        timing_scope='Separate preparation, matching, symmetry finalization, family reduction, '
                     'reference verification and audit CPU/wall; no decoding. Worker totals include publication; '
                     'process startup is separate. Cache off disables completed-fragment retention only.')
    started, stop = time.monotonic(), False
    def cancel(signum, frame):
        nonlocal stop
        stop = True
    previous = {sig: signal.signal(sig, cancel) for sig in (signal.SIGINT, signal.SIGTERM)}
    plans = [(row, count, setting) for case in cases for count in seeds for setting in cache
             for row in rows if row['case'] == case]
    report['expected_directions'] = len(plans)
    try:
        for row, count, setting in plans:
            if stop or time.monotonic()-started >= run_seconds:
                break
            direction = 'reverse' if row['reversed'] else 'forward'
            key = f"{row['case']:05d}/seeds-{count}/cache-{setting}/{direction}"
            report.update(active=key, elapsed_seconds=time.monotonic()-started)
            _publish(output/'results.json', report)
            # A new pool gives this direction a new process and uncontaminated RSS peak.
            with WorkPool.local(workers=1, memory_bytes=memory_bytes, shared_cache=False) as pool:
                remaining = run_seconds-(time.monotonic()-started)
                if stop or remaining <= 0:
                    break
                selected = replace(config, seed_count=count, branch_limit=branch_limit)
                task = Task(key, 'golden:diagnose_direction', dict(key=key, row=row,
                    audit_catalogue=audit_catalogue,
                    config=asdict(selected), cache=setting == 'on', progress=str(output/'progress.json'),
                    runner_sha256=report['runner_sha256']),
                    case=key, case_seconds=min(direction_seconds, remaining), memory_bytes=memory_bytes)
                pool.queue.add([task])
                next_report = 0.
                while True:
                    status = pool.queue.results([key])[key]
                    if status['state'] not in ('pending', 'running'):
                        break
                    pool.check_workers()
                    if time.monotonic() >= next_report:
                        report.update(active_state=status['state'], elapsed_seconds=time.monotonic()-started,
                                      cancellation_requested=stop)
                        _publish(output/'results.json', report)
                        print(json.dumps(dict(key=key, state=status['state'], elapsed_seconds=report['elapsed_seconds'])), flush=True)
                        next_report = time.monotonic()+10
                    time.sleep(.1)
                measured = status['result'] or {}
                record = dict(case=row['case'], reversed=row['reversed'], seeds=count, cache=setting,
                    config=asdict(selected), state=status['state'], startup_seconds=pool.startup_elapsed_seconds,
                    worker_metrics={name: value for name, value in measured.items()
                        if name.endswith(('_seconds', '_bytes', '_samples'))})
                if status['state'] == 'complete':
                    record['result'] = pool.store.load(measured['blob'])
                else:
                    record['error'] = measured.get('error')
                    if (output/'progress.json').exists():
                        partial = json.loads((output/'progress.json').read_text())
                        if partial['key'] == key:
                            record['partial'] = partial
                report['results'].append(record)
            peers = [r for r in report['results'] if r['case'] == row['case']
                     and r['seeds'] == count and r['reversed'] == row['reversed']]
            if len(peers) == 2:
                comparable = all(r['state'] == 'complete' for r in peers)
                equal = (all(peers[0]['result'][field] == peers[1]['result'][field]
                             for field in ('scope_sha256', 'catalogue_sha256', 'cuts', 'cap_stages'))
                         if comparable else None)
                report['comparisons'].append(dict(case=row['case'], reversed=row['reversed'],
                    seeds=count, complete=comparable, catalogue_and_scope_equal=equal))
                if equal is False:
                    raise AssertionError(f'Cache control changed the catalogue or search scope: {key}')
            report.update(active=None, elapsed_seconds=time.monotonic()-started)
            _publish(output/'results.json', report)
            print(json.dumps(dict(key=key, state=record['state'], completed=len(report['results']),
                                  expected=len(plans))), flush=True)
        report.update(complete=len(report['results']) == len(plans)
                      and all(r['state'] == 'complete' for r in report['results']),
                      stopped='cancelled' if stop else ('run_budget' if len(report['results']) < len(plans) else None))
    finally:
        report.update(elapsed_seconds=time.monotonic()-started, active=None)
        _publish(output/'results.json', report)
        for sig, handler in previous.items():
            signal.signal(sig, handler)
    return report


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ['diagnose']:
        parser = argparse.ArgumentParser(description='Bounded, sequential Golden diagnostics on selected cases')
        parser.add_argument('--index', action='append', required=True)
        parser.add_argument('--output', required=True)
        parser.add_argument('--cases', type=int, nargs='+', required=True)
        parser.add_argument('--directions', choices=['forward', 'reverse'], nargs='+',
                            default=['forward', 'reverse'])
        parser.add_argument('--seeds', type=int, nargs='+', default=[1, 2, 5, 10])
        parser.add_argument('--cache', choices=['on', 'off'], nargs='+', default=['on', 'off'])
        parser.add_argument('--audit-catalogue', action='store_true',
                            help='Keep catalogue fingerprints for comparison across separate runs')
        parser.add_argument('--branch-limit', type=int, default=2000)
        parser.add_argument('--direction-seconds', type=float, required=True)
        parser.add_argument('--run-seconds', type=float, required=True)
        parser.add_argument('--memory-gib', type=float, required=True)
        parser.add_argument('--cpu-model', help='Verified CPU model for the timing record; omit if unknown')
        args = parser.parse_args(argv[1:])
        report = diagnose(args.index, args.output, cases=args.cases, seeds=args.seeds, cache=args.cache,
            directions=args.directions, audit_catalogue=args.audit_catalogue,
            branch_limit=args.branch_limit, direction_seconds=args.direction_seconds,
            run_seconds=args.run_seconds, memory_bytes=int(args.memory_gib*2**30), cpu_model=args.cpu_model)
        return 0 if report['complete'] else 75
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs'); parser.add_argument('manifest'); parser.add_argument('output')
    parser.add_argument('--batches', type=int, default=5)
    args = parser.parse_args(argv)
    prepare(args.inputs, args.manifest, args.output, batches=args.batches)
    return 0


if __name__ == '__main__':
    sys.exit(main())
