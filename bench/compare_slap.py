"""Paired Golden timings: current MAPPA and pinned upstream SLAP on one machine."""
import argparse
from collections import Counter, defaultdict
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import signal
import tempfile
import time

from golden import _match_direction, _measure, _publish, verify_catalogue
from mappa import AAMSearchConfig
from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob
from mappa.mappa_core.utils.io.records import problem_from_record


SLAP_REVISION = 'ea248fd9494f52f4865193e87a98cc92c62b5f9e'
SLAP_CORE_SHA256 = '62964b12866e691698d594bfdc65f84ae714fb049df4267a5f5964b5c7e5fa3c'
CONFIG = AAMSearchConfig(seed_count=1, branch_limit=100, iso_tolerance=2., random_seed=42)
DECODE = dict(event_slack=5, threshold=.5, metal_threshold=.3, fragment_seconds=5.)


def hardware():
    if platform.system() == 'Linux':
        return next(line.split(':', 1)[1].strip() for line in
                    Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name'))
    import subprocess
    return subprocess.check_output(['sysctl', '-n', 'machdep.cpu.brand_string'], text=True).strip()


def slap_source():
    import slapmapper
    root = Path(slapmapper.__file__).parent
    hashes = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(root.rglob('*.py'))}
    if hashes['core.py'] != SLAP_CORE_SHA256:
        raise ValueError('SLAP core differs from the pinned published baseline')
    return hashes


def slap_graphs(problem):
    """Preserve prepared atom indices, using upstream weights and isolated balancing atoms."""
    from rdkit import Chem
    from slapmapper.core import LabeledGraph
    table = Chem.GetPeriodicTable()
    ends = (problem.reactant, problem.product)
    labels = [[table.GetAtomicNumber(str(e)) for e in end.elements] for end in ends]
    difference = Counter(labels[1]); difference.subtract(labels[0])
    for element, count in difference.items():
        labels[0 if count > 0 else 1].extend([element] * abs(count))
    graphs = []
    for end, colored in zip(ends, labels):
        n = len(end.elements)
        adjacency = {i: {} for i in range(len(colored))}
        for i in range(n):
            for j in range(i + 1, n):
                if end.wbo[i, j] > .2:
                    weight = int(2 * end.wbo[i, j] + .1)
                    adjacency[i][j] = adjacency[j][i] = weight
        graphs.append(LabeledGraph(adjacency, colored))
    return graphs


def cut_graphs(base, cut):
    from slapmapper.core import LabeledGraph
    adjacency = {i: dict(neighbors) for i, neighbors in base[0].graph.items()}
    for a, b in cut:
        del adjacency[a][b]; del adjacency[b][a]
    # Rebuild initial refinement labels after the cut, as in the audited adapter.
    return [LabeledGraph(adjacency, list(base[0].labels)), base[1].copy()]


def heavy_mapping(labels, problem):
    """Read native heavy labels; repeated labels are invalid, never free assignments."""
    ends = (problem.reactant, problem.product)
    indices = []
    for values, end in zip(labels, ends):
        numbered = {}
        for i, element in enumerate(end.elements):
            if element == 'H':
                continue
            label = int(values[i])
            if label in numbered:
                raise ValueError('SLAP returned a repeated heavy-atom label')
            numbered[label] = i
        indices.append(numbered)
    mapping = tuple(sorted((i, indices[1][label]) for label, i in indices[0].items()
                           if label in indices[1]))
    if any(ends[0].elements[a] != ends[1].elements[b] for a, b in mapping):
        raise ValueError('SLAP returned a non-element-preserving mapping')
    return mapping


def match(payload, workspace):
    import numpy as np
    from mappa.mappa_core.core.search.planning import cut_sweep_items
    from mappa.mappa_core.core.scheduler.results import CatalogueResult, CutScope
    from slapmapper.core import SlapMapper
    metrics = {}
    with _measure(metrics, 'loading'):
        problem = problem_from_record(read_blob(Blob(**payload['input']), workspace.scratch))
    result = dict(method=payload['method'], cpu_model=hardware(), metrics=metrics)
    if payload['method'] == 'mappa':
        catalogue, search = _match_direction(problem, CONFIG, key=payload['name'], metrics=metrics)
        scope = tuple(CutScope(c['number'], tuple(map(tuple, c['bonds'])), (0,))
                      for c in search.pop('cuts'))
        saved = CatalogueResult(payload['name'], catalogue, scope,
            tuple(sorted(search['stops'].items())), CONFIG, 'reused_native')
        result['search'] = search
    else:
        with _measure(metrics, 'preparation'):
            base = slap_graphs(problem)
            cuts = cut_sweep_items(problem.reactant.wbo, .2)
            targets = [i for i, label in enumerate(base[0].labels) if label > 1]
        records, mappings = [], set()
        for mode in ('binary', 'weighted'):
            random.seed(20260910); np.random.seed(20260910)
            mapper = SlapMapper(binary=mode == 'binary')
            for cut in cuts:
                with _measure(metrics, 'preparation'):
                    graphs = cut_graphs(base, cut)
                with _measure(metrics, 'matching'):
                    mapper.get_maps(graphs, break_sym_targets=targets)
                with _measure(metrics, 'reduction'):
                    for found in mapper.results:
                        labels = tuple(tuple(map(int, g.labels)) for g in found['lgp'])
                        mapping = heavy_mapping(labels, problem)
                        mappings.add(mapping)
                        records.append(dict(mode=mode, cut=cut, labels=labels,
                                            mapping=mapping, cost=float(found['val'])))
        saved = dict(records=records, mappings=tuple(sorted(mappings)))
        result['search'] = dict(expected_cuts=len(cuts), completed_cuts=len(cuts),
            calls=2 * len(cuts), candidates=len(records), unique_heavy_mappings=len(mappings))
    with _measure(metrics, 'publication'):
        result['artifact'] = asdict(workspace.writer.put(saved))
    return result


def decode(payload, workspace):
    from mappa import decode_events
    from mappa.mappa_core.utils.io.records import decoded_events_record
    metrics = {}
    with _measure(metrics, 'loading'):
        matched = workspace.load(payload['matching'])
        saved = workspace.load(Blob(**matched['artifact']))
    with _measure(metrics, 'decoding'):
        interpreted = decode_events(saved.catalogue, **DECODE)
        window = (interpreted.minimum_events, interpreted.event_window)
        record = decoded_events_record(interpreted)
        del interpreted
    with _measure(metrics, 'publication'):
        artifact = workspace.writer.put(record)
    return dict(artifact=asdict(artifact), metrics=metrics, window=window, cpu_model=hardware())


def reference_certificate(features, mapping):
    """Canonical paired heavy graph, preserving chemistry and all unmatched atoms."""
    import pynauty
    colors, adjacency = defaultdict(set), defaultdict(set)
    locations, count = [], 0
    for side, feature in enumerate(features):
        locations.append({atom: count + i for i, atom in enumerate(feature['heavy'])})
        for color in feature['colors']:
            colors[('atom', side, tuple(color))].add(count); count += 1
    edges = []
    for side, feature in enumerate(features):
        for a, b, color in feature['bonds']:
            edges.append((locations[side][feature['heavy'][a]],
                          locations[side][feature['heavy'][b]], ('bond', side, tuple(color))))
    edges.extend((locations[0][a], locations[1][b], ('mapping',)) for a, b in mapping)
    for a, b, color in edges:
        colors[color].add(count)
        adjacency[count].update((a, b)); adjacency[a].add(count); adjacency[b].add(count)
        count += 1
    return pynauty.certificate(pynauty.Graph(count,
        adjacency_dict={i: sorted(adjacency[i]) for i in range(count)},
        vertex_coloring=[colors[k] for k in sorted(colors, key=repr)]))


def verify(payload, workspace):
    from decoded_reference import verify_decoded
    from mappa.mappa_core.utils.io.records import decoded_events_from_record
    metrics = {}
    with _measure(metrics, 'loading'):
        matched = workspace.load(payload['matching'])
        saved = workspace.load(Blob(**matched['artifact']))
        reference = read_blob(Blob(**payload['reference']), workspace.scratch)
    with _measure(metrics, 'verification'):
        if payload['method'] == 'mappa':
            answer = verify_catalogue(saved.catalogue, reference, reverse=payload['reversed'])
        else:
            features = reference['features'][::-1] if payload['reversed'] else reference['features']
            mapping = [(b, a) for a, b in reference['mapping']] if payload['reversed'] else reference['mapping']
            expected = reference_certificate(features, mapping)
            hit = next((m for m in saved['mappings'] if reference_certificate(features, m) == expected), None)
            answer = dict(status='recovered' if hit is not None else 'not_recovered', witness=hit)
    result = dict(verification=answer, metrics=metrics)
    if payload.get('decoding'):
        with _measure(metrics, 'decoded_loading'):
            decoded = workspace.load(payload['decoding'])
            interpreted = decoded_events_from_record(workspace.load(Blob(**decoded['artifact'])))
        with _measure(metrics, 'decoded_verification'):
            result['decoded_verification'] = verify_decoded(interpreted, reference,
                reverse=payload['reversed'], preferred_family=answer.get('family'), seconds=10.)
    return result


def summary(plan, ledger, store):
    import numpy as np
    rows = []
    for name, stage, state, reference, audit, error in ledger.db.execute('SELECT * FROM results'):
        value = store.load(Blob(**json.loads(reference))) if state == 'complete' else {}
        if stage == 'matching':
            value.get('search', {}).pop('metrics', None)
        if stage == 'verification':
            # Full replayable witnesses remain in the aggregate bundle.
            value = {key: (dict(status=v['status']) if key.endswith('verification') else v)
                     for key, v in value.items()}
        rows.append(dict(name=name, stage=stage, state=state, result=value,
                         audit=json.loads(audit) if audit else None, error=json.loads(error) if error else None))
    indexed = {(r['name'], r['stage']): r for r in rows}
    def completed(method, case, stage):
        return [indexed.get((f'{method}/{case:05d}/{direction}', stage))
                for direction in ('forward', 'reverse')]
    def success(records):
        return all(r is not None and r['state'] == 'complete' for r in records)
    common = [c for c in plan['cases'] if all(success(completed(m, c, 'matching'))
                                            for m in ('mappa', 'slap'))]
    costs, recovery = {}, {}
    for method in ('mappa', 'slap'):
        times = [sum(sum(r['result']['metrics'].get(phase, {}).get('cpu_seconds', 0.)
                         for phase in ('preparation', 'matching', 'reduction'))
                     for r in completed(method, c, 'matching')) for c in common]
        costs[method] = dict(reactions=len(times), total_cpu_seconds=sum(times),
            mean_cpu_seconds=float(np.mean(times)) if times else None,
            median_cpu_seconds=float(np.median(times)) if times else None,
            p95_cpu_seconds=float(np.percentile(times, 95)) if times else None)
        outcomes = []
        for case in plan['cases']:
            records = completed(method, case, 'verification')
            statuses = [r['result']['verification']['status'] if r and r['state'] == 'complete'
                        else 'unknown' for r in records]
            outcomes.append('recovered' if 'recovered' in statuses else
                            'unknown' if 'unknown' in statuses else 'not_recovered')
        recovery[method] = dict(Counter(outcomes))
    return dict(plan_sha256=hashlib.sha256(json.dumps(plan, sort_keys=True).encode()).hexdigest(),
                expected_reactions=len(plan['cases']), common_completed_cases=common,
                paired_search=costs, catalogue_recovery=recovery, rows=rows)


def run(args):
    import fcntl
    from enzymemap import Ledger, require_drained
    from mappa.execution import BundleStore, WorkPool
    from mappa.execution.identity import code_identity
    from mappa.mappa_core.core.scheduler import Task
    from importlib.metadata import version

    index = json.loads(args.index.read_bytes())
    cases = list(range(1851)) if args.cases is None else sorted(set(args.cases))
    rows = sorted((r for r in index['directions'] if r['case'] in cases),
                  key=lambda r: (r['case'], r['reversed']))
    if len(rows) != 2 * len(cases) or {(r['case'], r['reversed']) for r in rows} != {
            (c, rev) for c in cases for rev in (False, True)}:
        raise ValueError('Require exactly both directions of every selected Golden case')
    for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        if os.environ.get(variable) != '1':
            raise ValueError(f'Set {variable}=1 before starting Python')
    plan = dict(schema='mappa.paired_slap/v1', cases=cases, inputs=rows,
        input_sha256=hashlib.sha256(args.index.read_bytes()).hexdigest(),
        mappa_identity=code_identity(), mappa_config=asdict(CONFIG), decode_config=DECODE,
        slap_revision=SLAP_REVISION, slap_hashes=slap_source(),
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        helper_sha256={name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                       for name in ('golden.py', 'decoded_reference.py', 'enzymemap.py')},
        cpu_model=hardware(), python=platform.python_version(),
        packages={name: version(name) for name in ('numpy', 'scipy', 'rdkit', 'pynauty')},
        memory_gib=args.memory_gib, stage_seconds=args.seconds,
        timing='Prepared identical graphs; preparation + search + reduction; input loading, publication, decoding and verification separate',
        reference_policy='Full heavy mapping including unmatched atoms modulo chemical endpoint symmetry; decoded audit separately uses exact indices',
        slap_policy='Binary and weighted, explicit H, heavy symmetry breaking, isolated upstream element balancing, uncut plus all source edges',
        method_order='Alternate first admitted method by case parity; one persistent pool on one CPU model')
    plan = json.loads(json.dumps(plan))
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output/'.owner.lock').open('a+b') as owner:
        fcntl.flock(owner, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if args.resume:
            if json.loads((args.output/'plan.json').read_text()) != plan:
                raise ValueError('Resume must preserve inputs, source, settings and budgets')
            if not json.loads((args.output/'progress.json').read_text()).get('drained'):
                raise ValueError('Previous owner did not drain')
            require_drained(args.output)
        elif (args.output/'plan.json').exists():
            raise ValueError('Existing run requires explicit --resume')
        else:
            _publish(args.output/'plan.json', plan)
        stop, started = False, time.monotonic()
        def drain(signum, frame):
            nonlocal stop
            stop = True
        handlers = {s: signal.signal(s, drain) for s in (signal.SIGTERM, signal.SIGINT)}
        try:
            with tempfile.TemporaryDirectory(prefix='mappa-paired-') as scratch:
                ledger = Ledger(args.output, scratch, resume=args.resume)
                store = BundleStore(args.output/'bundles', scratch=scratch)
                factory = WorkPool.slurm if args.slurm else WorkPool.local
                with factory(workers=args.workers, store=store, memory_bytes=int(args.memory_gib*2**30),
                             case_seconds=args.seconds, shared_cache=False) as pool:
                    status = dict(ready_workers=len(pool.queue.ready_workers()), cpu_model=hardware(),
                                  startup_seconds=pool.startup_elapsed_seconds, drained=False)
                    def checkpoint(stage, active):
                        ledger.checkpoint(); pool.checkpoint()
                        status.update(stage=stage, counts=ledger.counts(), active=len(active),
                                      elapsed_seconds=time.monotonic()-started)
                        _publish(args.output/'progress.json', status)
                        print(json.dumps(status), flush=True)
                    for stage in ('matching', 'decoding', 'verification'):
                        pending = iter((method, row) for row in rows for method in
                            (('mappa', 'slap') if row['case'] % 2 == 0 else ('slap', 'mappa'))
                            if stage != 'decoding' or method == 'mappa')
                        active, exhausted, next_save = {}, False, 0.
                        while active or not exhausted:
                            stop |= time.monotonic()-started >= args.run_seconds
                            while not stop and not exhausted and len(active) < args.workers*2:
                                try:
                                    method, row = next(pending)
                                except StopIteration:
                                    exhausted = True; break
                                name = f"{method}/{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}"
                                if ledger.get(name, stage):
                                    continue
                                payload = dict(name=name, method=method, input=row['input'])
                                if stage != 'matching':
                                    prior = ledger.get(name, 'matching')
                                    if prior is None:
                                        raise RuntimeError('Matching barrier is incomplete')
                                    if prior[0] != 'complete':
                                        ledger.put(name, stage, 'blocked', error='Matching failed'); continue
                                    payload['matching'] = Blob(**json.loads(prior[1]))
                                if stage == 'verification':
                                    payload.update(reference=row['reference'], reversed=row['reversed'])
                                    prior = ledger.get(name, 'decoding')
                                    if prior and prior[0] == 'complete':
                                        payload['decoding'] = Blob(**json.loads(prior[1]))
                                key = f'{name}/{stage}'
                                operation = dict(matching='match', decoding='decode', verification='verify')[stage]
                                pool.queue.add([Task(key, f'compare_slap:{operation}', payload,
                                    case=key, case_seconds=args.seconds)])
                                active[name] = key
                            for name, key in list(active.items()):
                                state = pool.queue.results([key])[key]
                                if state['state'] in ('pending', 'running'):
                                    continue
                                measured = state['result'] or {}
                                success = state['state'] == 'complete'
                                ledger.put(name, stage, 'complete' if success else 'error',
                                    measured['blob'] if success else None, pool.queue.archive_case(key, store),
                                    None if success else dict(state=state['state'], error=measured.get('error')))
                                del active[name]
                                if not success:
                                    stop = True  # Preserve evidence and investigate before continuing.
                            if time.monotonic() >= next_save or not active:
                                checkpoint(stage, active); next_save = time.monotonic()+30
                            if stop and not active:
                                break
                            if active:
                                pool.check_workers(); time.sleep(.1)
                        if stop:
                            break
                    status['drained'] = True
                    checkpoint(stage, {})
                report = summary(plan, ledger, store)
                report['complete'] = not stop
                report['elapsed_seconds'] = time.monotonic()-started
                _publish(args.output/'results.json', report)
                ledger.db.close()
        finally:
            for s, handler in handlers.items():
                signal.signal(s, handler)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('index', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--cases', type=int, nargs='+')
    parser.add_argument('--workers', type=int, default=16)
    parser.add_argument('--memory-gib', type=float, default=4.)
    parser.add_argument('--seconds', type=float, default=600.)
    parser.add_argument('--run-seconds', type=float, default=9600.)
    parser.add_argument('--slurm', action='store_true')
    parser.add_argument('--resume', action='store_true')
    arguments = parser.parse_args()
    result = run(arguments)
    raise SystemExit(0 if result['complete'] else 75)
