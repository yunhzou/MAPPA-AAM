"""Explicit event-count experiments over immutable completed catalogues."""
import argparse
from collections import Counter
from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from mappa.execution import BundleStore, WorkPool
from mappa.execution.identity import code_identity
from mappa.mappa_core.core.scheduler import Task
from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob
from mappa.mappa_core.core.scheduler.results import DecodingResult
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex
from mappa.mappa_core.post_processing.decoder.api import _validate_atom_identity
from mappa.mappa_core.utils.io.records import problem_from_record
from mappa.mappa_workflow.mapping_analytics.event_query import _event_model, _event_count_constraint
from mappa.mappa_core.core.search.queries.family import compile_path


def heavy_event_index(problem, *, threshold=.5, metal_threshold=.3):
    """Project event counting only; keep full atom indices for original constraints."""
    endpoints = []
    for endpoint in (problem.reactant, problem.product):
        heavy = np.asarray([element != 'H' for element in endpoint.elements])
        weights = np.where(heavy[:, None] & heavy[None, :], endpoint.wbo, 0.)
        endpoints.append(replace(endpoint, wbo=weights))
    counted = replace(problem, reactant=endpoints[0], product=endpoints[1])
    index = SignedEventIndex(counted, threshold=threshold, metal_threshold=metal_threshold)
    # Omitted pairs have identically zero responses. Retain the original atom
    # numbering so the existing symbolic family expressions need no translation.
    for endpoint, names in ((counted.reactant, ('a', 'b', 'tau')),
                            (counted.product, ('pa', 'pb', 'ptau'))):
        a, b, thresholds = (getattr(index, name) for name in names)
        heavy = np.asarray([element != 'H' for element in endpoint.elements])
        keep = heavy[a] & heavy[b]
        for name, values in zip(names, (a, b, thresholds)):
            setattr(index, name, values[keep])
    return index


def saved_response_bounds(record, problem, deadline):
    """Admissible heavy-event bounds, scoped to this complete decoded record.

    Each partition contains each source bond once. Alternative partitions are
    separate decompositions, so their bounds combine by maximum, never addition.
    Product-boundary occupancy is not proved here and contributes a safe zero.
    Correlation and required edges can restrict choices further, never invalidate
    an event whose saved local responses are all nonzero.
    """
    from mappa.mappa_core.utils.io.records import decoded_events_from_record
    decoded = decoded_events_from_record(record)
    heavy = [element != 'H' for element in problem.reactant.elements]
    table_bounds = {}
    bounds = [0] * len(record['families'])
    checked = 0
    for fid, number, table in decoded._export_tables():
        if time.perf_counter() >= deadline:
            break
        if table is not None:
            table_bounds[number] = sum(
                heavy[term.atoms[0]] and heavy[term.atoms[1]]
                and bool(len(term.response.categories)) and bool(all(term.response.categories))
                for group in (*table.internal.values(), *table.boundary.values()) for term in group)
            checked += 1
        bounds[fid] = max(bounds[fid], table_bounds[number])
    return bounds, checked


def minimum_events(decoded, reference, *, seconds=30.):
    """Certify a retained-catalogue minimum; unknown never becomes a point estimate.

    Reference and candidate counts use heavy-heavy original-WBO events. Matching
    constraints still use the original full-atom problem. The ordinary decoder's
    mapped/unmapped boundary rule is unchanged; witness coverage is reported.
    This is a soft encoding/solving budget. The runner supplies a process guard.
    """
    import z3

    if not math.isfinite(seconds) or seconds < 0:
        raise ValueError('seconds must be finite and nonnegative')
    if decoded['schema'] != 'mappa.decoded_events/v14':
        raise ValueError('Require original-space decoded responses')
    started, cpu = time.perf_counter(), time.process_time()
    deadline = started + seconds
    problem = problem_from_record(decoded['problem'])
    catalogue = decoded['catalogue']
    families = FinalBranchCatalogue.from_record(problem, catalogue).families
    if (FinalBranchCatalogue(problem).problem_sha256 != catalogue['problem_sha256'] or
            len(decoded['families']) != len(families)):
        raise ValueError('Decoded responses and catalogue must belong to the same problem')
    event_problem = problem
    if 'event_problem' in decoded:
        event_problem = problem_from_record(decoded['event_problem'])
        _validate_atom_identity(problem, event_problem)
        if (decoded['matching_problem_sha256'] != catalogue['problem_sha256'] or
                decoded['event_problem_sha256'] != FinalBranchCatalogue(event_problem).problem_sha256):
            raise ValueError('Saved matching and event identities differ from their matrices')
    threshold = decoded['event_policy']['threshold']
    metal_threshold = decoded['event_policy']['metal_threshold']
    index = heavy_event_index(event_problem, threshold=threshold, metal_threshold=metal_threshold)
    reference = {int(a): int(b) for a, b in reference.items()}
    def count(mapping):
        return int(index.counts([index.vector(mapping)])[0].sum())

    reference_count = count(reference)
    upper, witness, witness_family = None, None, None
    scores = []
    lower, checked_tables = saved_response_bounds(decoded, problem, deadline)
    saved_bounds = tuple(lower)
    proofs, unknown = Counter(), []
    queries = 0
    initial_global_lower = min(lower, default=0)

    def finish():
        exact = upper is not None and (upper == 0 or min(lower, default=0) >= upper)
        selected = {} if witness is None else witness
        return dict(status='exact' if exact else 'empty_catalogue' if not families else 'unknown',
            reference_events=reference_count, minimum_events=upper if exact else None,
            minimum_lower_bound=min(lower, default=0), minimum_upper_bound=upper,
            difference=reference_count-upper if exact else None,
            witness_family=witness_family, witness=sorted(selected.items()),
            reference_mapped_heavy_atoms=sum(problem.reactant.elements[a] != 'H' for a in reference),
            witness_mapped_heavy_atoms=sum(problem.reactant.elements[a] != 'H' for a in selected),
            same_heavy_coverage=(set(a for a in selected if problem.reactant.elements[a] != 'H') ==
                                 set(a for a in reference if problem.reactant.elements[a] != 'H') and
                                 set(b for b in selected.values() if problem.product.elements[b] != 'H') ==
                                 set(b for b in reference.values() if problem.product.elements[b] != 'H')),
            families=len(families), witnesses_scored=len(scores), proof_counts=dict(proofs),
            saved_response_tables_checked=checked_tables,
            saved_response_certified_families=sum(upper is not None and value >= upper
                                                 for value in saved_bounds),
            unresolved_family_count=sum(upper is None or value < upper for value in lower),
            solver_queries=queries, unresolved_families=[row for row in unknown
                if upper is None or lower[row['family']] < upper],
            incomplete_paths=catalogue['incomplete_paths'],
            problem_sha256=catalogue['problem_sha256'], catalogue_schema=catalogue['schema'],
            event_problem_sha256=(decoded['event_problem_sha256']
                if 'event_problem' in decoded else catalogue['problem_sha256']),
            decoded_schema=decoded['schema'],
            cpu_seconds=time.process_time()-cpu, wall_seconds=time.perf_counter()-started,
            event_policy=dict(threshold=threshold, metal_threshold=metal_threshold,
                              atom_scope='heavy-heavy', comparison='absolute difference > threshold'),
            scope='Minimum over this saved direction catalogue; no claim of exhaustive search coverage')

    for fid, family in enumerate(families):
        if time.perf_counter() >= deadline:
            return finish()
        value = count(dict(family.mapping))
        scores.append((value, fid))
        if lower[fid] > value:
            raise AssertionError('Saved response bound exceeds a valid witness count')
        if upper is None or value < upper:
            family.validate_representative(problem)
            upper, witness, witness_family = value, dict(family.mapping), fid
        if value == 0 or initial_global_lower >= upper:
            return finish()

    for score, fid in sorted(scores):
        if upper == 0:
            break
        if time.perf_counter() >= deadline:
            unknown.append(dict(family=fid, reason='time_budget'))
            break
        if lower[fid] >= upper:
            continue
        family = families[fid]
        if not family.actions:
            lower[fid] = score
            proofs['fixed_mapping'] += 1
            continue
        path = family.as_path(problem)
        if index.invariant_path(path):
            lower[fid] = score
            proofs['event_invariance'] += 1
            continue
        compiled = compile_path(path, problem, {}, source_atoms=(), complete_reference=False)
        try:
            terms, objective, metrics = _event_model(compiled, index, deadline)
        except TimeoutError:
            unknown.append(dict(family=fid, reason='event_encoding_budget'))
            continue
        lower[fid] = max(lower[fid], metrics['event_lower_bound'])
        while lower[fid] < upper:
            remaining = deadline-time.perf_counter()
            if remaining <= 0:
                unknown.append(dict(family=fid, reason='time_budget'))
                break
            compiled.solver.add(_event_count_constraint(terms, upper-1, at_most=True))
            compiled.solver.set(timeout=max(1, int(remaining*1000)))
            status = compiled.solver.check()
            queries += 1
            if status == z3.unsat:
                lower[fid] = upper
                proofs['no_improving_assignment'] += 1
                break
            if status == z3.unknown:
                unknown.append(dict(family=fid, reason=compiled.solver.reason_unknown()))
                break
            model = compiled.solver.model()
            mapping = dict(compiled.realize(model)['mapping'])
            value = count(mapping)
            if value != model.eval(objective).as_long() or value >= upper:
                raise AssertionError('Exact symbolic and direct event counts differ')
            upper, witness, witness_family = value, mapping, fid
        if lower[fid] >= upper and lower[fid] == metrics['event_lower_bound']:
            proofs['event_lower_bound'] += 1
    return finish()


def event_minimum(payload, workspace):
    """Query the exact catalogue already contained in the immutable decoder output."""
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != payload['runner_sha256']:
        raise ValueError('Analysis runner changed after planning')
    prepared = time.process_time()
    saved = workspace.load(Blob(**payload['row']['decoded']))
    if not isinstance(saved, DecodingResult):
        raise TypeError('Event minima require a saved original-space DecodingResult')
    if (saved.record['event_policy']['threshold'] != saved.config.threshold or
            saved.record['event_policy']['metal_threshold'] != saved.config.metal_threshold):
        raise ValueError('Saved event configuration and response policy differ')
    reference = read_blob(Blob(**payload['row']['reference']), workspace.scratch)
    mapping = dict(reference['mapping'])
    if payload['row']['reversed']:
        mapping = {b: a for a, b in mapping.items()}
    preparation_cpu = time.process_time()-prepared
    result = minimum_events(saved.record, mapping, seconds=payload['seconds'])
    result['preparation_cpu_seconds'] = preparation_cpu
    return result


def _publish(path, report):
    pending = path.with_suffix('.pending')
    pending.write_text(json.dumps(report, indent=2)+'\n')
    pending.replace(path)


def detected_position(decoded, reference, *, seconds=10.):
    """Count heavy events in one reference-independent saved-minimum witness.

    This retrieves a witness at the already saved model minimum. It does not
    optimize a new heavy-only objective or use the reference to choose a mapping.
    """
    started = time.perf_counter()
    minimum = decoded.minimum_events
    result = dict(status='unknown', saved_model_minimum=minimum,
                  model_minimum_proven=decoded.minimum_proven)
    if minimum is None:
        return dict(result, status='empty')
    for key in sorted(decoded.window.contexts, key=lambda key: (decoded.window.bounds[key], key)):
        if decoded.window.bounds[key] > minimum:
            continue
        remaining = seconds - (time.perf_counter() - started)
        if remaining <= 0:
            break
        table = decoded._prepare_table(*key[:2], decoded._table_key(*key))
        answer = table.program.select(minimum, minimum).find(seconds=remaining, max_nodes=100000)
        if answer.status != 'found':
            del table
            continue
        witness = dict(answer.witness)
        if (decoded.realize(key[0], answer.choices) != witness or
                len(set(witness.values())) != len(witness) or
                int(table.evaluate(decoded.index.vector(witness))['total']) != minimum):
            raise AssertionError('Saved-minimum witness failed replay or counting')
        index = heavy_event_index(decoded.index.problem, threshold=decoded.index.threshold,
                                  metal_threshold=decoded.index.metal_threshold)
        counts = index.counts([index.vector(reference), index.vector(witness)]).sum(axis=1)
        heavy = lambda endpoint, atoms: {a for a in atoms if endpoint.elements[a] != 'H'}
        return dict(result, status='complete', reference_events=int(counts[0]),
            detected_events=int(counts[1]), difference=int(counts[0]-counts[1]),
            same_heavy_coverage=(heavy(index.problem.reactant, reference) == heavy(index.problem.reactant, witness)
                and heavy(index.problem.product, reference.values()) == heavy(index.problem.product, witness.values())),
            witness=sorted(witness.items()), context=key,
            query_seconds=time.perf_counter()-started)
    return dict(result, query_seconds=time.perf_counter()-started)


def event_position(payload, workspace):
    """One reaction; choose the lower saved model minimum, ties forward first."""
    from mappa.mappa_core.utils.io.records import decoded_events_from_record
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != payload['runner_sha256']:
        raise ValueError('Analysis runner changed after planning')
    candidates = []
    for row in payload['directions']:
        minimum = row.get('saved_minimum')
        if minimum is None:
            saved = workspace.load(Blob(**row['decoded']))
            minimum = saved.record['minimum_events']
            del saved
        if minimum is not None:
            candidates.append((minimum, row['reversed'], row))
    if not candidates:
        return dict(status='empty')
    minimum, reverse, row = min(candidates, key=lambda item: item[:2])
    saved = workspace.load(Blob(**row['decoded']))
    if saved.record['minimum_events'] != minimum:
        raise ValueError('Saved minimum metadata differs from decoded output')
    decoded = decoded_events_from_record(saved.record)
    reference = dict(workspace.load(Blob(**row['reference']))['mapping'])
    if reverse:
        reference = {b: a for a, b in reference.items()}
    result = detected_position(decoded, reference, seconds=payload['seconds'])
    return dict(result, reversed=reverse, decoded_sha256=row['decoded']['sha256'],
                reference_sha256=row['reference']['sha256'])


def run_positions(manifest, output, *, workers=1, seconds=10., task_seconds=120.,
                  memory_bytes=2**33, slurm=False):
    """Read-only reaction-level positions; keep incomplete queries outside bins."""
    from collections import defaultdict
    manifest, output = Path(manifest), Path(output)
    source = json.loads(manifest.read_bytes())
    sources = dict(source['source_files'])
    sources[str(manifest.resolve())] = hashlib.sha256(manifest.read_bytes()).hexdigest()
    def unchanged():
        return all(hashlib.sha256(Path(p).read_bytes()).hexdigest() == digest for p, digest in sources.items())
    if not unchanged():
        raise ValueError('Saved analysis inputs changed')
    grouped = defaultdict(list)
    for row in source['directions']:
        grouped[(row['dataset'], row['case'])].append(row)
    if any(sorted(row['reversed'] for row in pair) != [False, True] for pair in grouped.values()):
        raise ValueError('Require both directions exactly once per reaction')
    output.mkdir(parents=True, exist_ok=False)
    digest = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    report = dict(schema='mappa.detected_event_positions/v1', code_identity=code_identity(),
        runner_sha256=digest, original_hashes=sources, expected=len(grouped), results=[],
        policy='Heavy-heavy physical events: reference minus one saved-model-minimum witness; no heavy-only optimization',
        direction_policy='Lower saved model minimum across both directions; forward wins ties',
        resources=dict(workers=workers, query_seconds=seconds, task_seconds=task_seconds, memory_bytes=memory_bytes))
    store = BundleStore(output/'bundles')
    factory = WorkPool.slurm if slurm else WorkPool.local
    start = time.perf_counter()
    with factory(workers=workers, store=store, memory_bytes=memory_bytes,
                 case_seconds=task_seconds, shared_cache=False) as pool:
        report['ready_workers'] = len(pool.queue.ready_workers())
        from itertools import islice
        waiting = iter((f'{dataset}/{case:05d}', (dataset, case, rows))
                       for (dataset, case), rows in grouped.items())
        pending, exhausted = {}, False
        published = 0.
        while pending or not exhausted:
            capacity = 512-len(pending)
            admitted = dict(islice(waiting, capacity))
            if capacity:
                exhausted = not admitted
            pending.update(admitted)
            if admitted:
                pool.queue.add([Task(key, 'analytics:event_position',
                    dict(directions=rows, seconds=seconds, runner_sha256=digest), case=key,
                    case_seconds=task_seconds, memory_bytes=memory_bytes)
                    for key, (_, _, rows) in admitted.items()])
            pool.check_workers()
            for key, state in pool.queue.results(tuple(pending)).items():
                if state['state'] in ('pending', 'running'):
                    continue
                dataset, case, _ = pending.pop(key)
                measured = state['result'] or {}
                row = dict(dataset=dataset, case=case, state=state['state'],
                    metrics={k: v for k, v in measured.items() if k.endswith(('_seconds', '_bytes'))})
                if state['state'] == 'complete':
                    value = store.load(measured['blob'])
                    row.update({k: v for k, v in value.items() if k != 'witness'})
                    from dataclasses import asdict
                    row['result_blob'] = asdict(measured['blob'])
                else:
                    row.update(status='unknown', error=measured.get('error'))
                    if state['state'] == 'error':
                        pool.queue.finish_input(cancel_pending=True)
                report['results'].append(row)
                if state['state'] == 'error':
                    _publish(output/'results.json', report)
                    raise RuntimeError('Saved-minimum analysis failed: '+str(row))
            if time.perf_counter() >= published or not pending:
                report['elapsed_seconds'] = time.perf_counter()-start
                _publish(output/'results.json', report)
                print(json.dumps(dict(completed=len(report['results']), expected=len(grouped))), flush=True)
                published = time.perf_counter()+30
            if pending:
                time.sleep(.1)
    if len(report['results']) != len(grouped):
        raise AssertionError('Analysis did not cover the complete reaction scope')
    if not unchanged():
        raise ValueError('Saved sources changed during analysis')
    report.update(source_unchanged=True, elapsed_seconds=time.perf_counter()-start,
                  complete=all(row['status'] == 'complete' for row in report['results']))
    _publish(output/'results.json', report)
    return report


def plot_positions(results, output):
    """Publication figure from checked per-reaction counts; zero bins stay empty."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogLocator, NullFormatter, MaxNLocator
    data = json.loads(Path(results).read_bytes())
    rows = data['results']
    datasets = sorted({r['dataset'] for r in rows}, key=lambda name: (name != 'Golden', name))
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    fig, axes = plt.subplots(1, len(datasets), figsize=(11.5, 4.5), squeeze=False, layout='constrained')
    for ax, dataset, color in zip(axes[0], datasets, ('#3568a8', '#258a78')):
        cohort = [r for r in rows if r['dataset'] == dataset]
        measured = [r for r in cohort if r['status'] == 'complete']
        counts = Counter(r['difference'] for r in measured)
        xs = sorted(counts)
        ax.bar(xs, [counts[x] for x in xs], width=.82, color=color, zorder=3)
        ax.set_yscale('log')
        ax.set_ylim(.7, max(counts.values(), default=1)*2.5)
        ax.set_xlim(min(min(xs, default=0), -1)-.8, max(max(xs, default=0), 5)+.8)
        ax.axvline(0, color='#273543', lw=1, alpha=.65)
        ax.set_axisbelow(True)
        ax.grid(axis='y', which='major', color='#dce1e6', lw=.7)
        ax.yaxis.set_major_locator(LogLocator(base=10))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=12))
        ax.spines[['top', 'right']].set_visible(False)
        below = sum(v for x, v in counts.items() if x < 0)
        above = sum(v for x, v in counts.items() if x > 0)
        ax.set_title(f'{dataset}  |  {len(measured):,} / {len(cohort):,} reactions', loc='left', fontweight='bold', pad=15)
        ax.text(.98, .96, f'Below: {below:,}   Equal: {counts[0]:,}   Above: {above:,}',
                transform=ax.transAxes, ha='right', va='top', fontsize=9, color='#43505d')
        ax.set_xlabel(r'$n$ = ground-truth events − detected events', labelpad=10)
        ax.set_ylabel('Number of reactions (log scale)')
    fig.suptitle('Where ground truth falls relative to MAPPA’s detected mapping', fontsize=14, fontweight='bold')
    output = Path(output)
    for suffix in ('png', 'pdf', 'svg'):
        fig.savefig(output.with_suffix('.'+suffix), dpi=220, bbox_inches='tight')
    plt.close(fig)


def run(manifest, output, *, workers=1, seconds=30., task_seconds=120., memory_bytes=2**31):
    """A bounded pilot manifest, one independent task per saved direction."""
    if workers < 1 or not 0 < seconds < task_seconds or memory_bytes < 1:
        raise ValueError('Require positive workers/memory and 0 < seconds < task_seconds')
    manifest, output = Path(manifest), Path(output)
    raw = manifest.read_bytes()
    rows = json.loads(raw)['directions']
    keys = [f"{r['dataset']}/{r['case']}/{'reverse' if r['reversed'] else 'forward'}" for r in rows]
    if len(keys) != len(set(keys)):
        raise ValueError('Each dataset/case/direction must occur exactly once')
    output.mkdir(parents=True, exist_ok=False)
    report = dict(schema='mappa.event_minimum_experiment/v1', code_identity=code_identity(),
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        manifest_sha256=hashlib.sha256(raw).hexdigest(), expected=len(rows), terminal_complete=False,
        successful_analysis=False, minimum_scope_complete=False,
        counting='Heavy-heavy original-WBO changes; ordinary mapped/unmapped boundary policy',
        statistic='reference event count minus certified saved-direction minimum',
        resources=dict(workers=workers, query_seconds=seconds, task_seconds=task_seconds,
                       memory_guard_bytes=memory_bytes), results=[])
    _publish(output/'results.json', report)
    store = BundleStore(output/'bundles')
    with WorkPool.local(workers=workers, store=store, memory_bytes=memory_bytes,
                        case_seconds=task_seconds, shared_cache=False) as pool:
        pending = dict(zip(keys, rows))
        pool.queue.add([Task(key, 'analytics:event_minimum',
            dict(row=row, seconds=seconds, runner_sha256=report['runner_sha256']),
            case=key, case_seconds=task_seconds, memory_bytes=memory_bytes)
            for key, row in pending.items()])
        while pending:
            pool.check_workers()
            for key, state in pool.queue.results(tuple(pending)).items():
                if state['state'] in ('pending', 'running'):
                    continue
                row = pending.pop(key)
                measured = state['result'] or {}
                result = dict(dataset=row['dataset'], case=row['case'], reversed=row['reversed'],
                    state=state['state'], decoded_sha256=row['decoded']['sha256'],
                    reference_sha256=row['reference']['sha256'],
                    source_provenance=row.get('source_provenance'),
                    worker_metrics={k:v for k,v in measured.items()
                                    if k.endswith(('_seconds', '_bytes', '_samples'))})
                if state['state'] == 'complete':
                    result['result'] = store.load(measured['blob'])
                else:
                    result['error'] = measured.get('error')
                report['results'].append(result)
                _publish(output/'results.json', report)
                print(json.dumps(dict(key=key, state=state['state'], completed=len(report['results']),
                                      expected=len(rows))), flush=True)
            time.sleep(.1)
    report['terminal_complete'] = True
    report['successful_analysis'] = all(row['state'] == 'complete' for row in report['results'])
    report['histograms'] = {dataset: dict(sorted(Counter(
        row['result']['difference'] for row in report['results']
        if row['dataset'] == dataset and row.get('result', {}).get('status') == 'exact').items()))
        for dataset in sorted({row['dataset'] for row in rows})}
    report['exact_directions'] = sum(row.get('result', {}).get('status') == 'exact'
                                     for row in report['results'])
    report['minimum_scope_complete'] = report['exact_directions'] == len(rows)
    report['status_counts'] = dict(Counter(row.get('result', {}).get('status', 'task_error')
                                         for row in report['results']))
    _publish(output/'results.json', report)
    print(json.dumps({key: report[key] for key in (
        'terminal_complete', 'successful_analysis', 'minimum_scope_complete', 'status_counts')}), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--seconds', type=float, default=30.)
    parser.add_argument('--task-seconds', type=float, default=120.)
    parser.add_argument('--memory-gib', type=float, default=2.)
    parser.add_argument('--positions', action='store_true', help='Heavy-event positions at the already saved model minimum')
    parser.add_argument('--slurm', action='store_true', help='Use the existing allocation for --positions')
    args = parser.parse_args()
    if args.positions:
        run_positions(args.manifest, args.output, workers=args.workers, seconds=args.seconds,
            task_seconds=args.task_seconds, memory_bytes=int(args.memory_gib*2**30), slurm=args.slurm)
    else:
        run(args.manifest, args.output, workers=args.workers, seconds=args.seconds,
            task_seconds=args.task_seconds, memory_bytes=int(args.memory_gib*2**30))
