"""Complete supplied reference pairs with one uncut seed in each direction.

This reference-conditioned experiment does not measure independent recovery and
does not replace minimum queries over the original benchmark's saved branches.
"""
import argparse
from collections import Counter
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time

from mappa import AAMSearchConfig, search_aam, deduplicate, decode_events
from mappa.execution import BundleStore, WorkPool
from mappa.execution.identity import code_identity
from mappa.mappa_core.core.scheduler import Task
from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob
from mappa.mappa_core.utils.io.records import problem_from_record
from decoded_reference import _publish, audit_decoded_witness


SETTINGS = dict(seed_count=1, branch_limit=100, iso_tolerance=2,
                random_seed=42, sweep_cuts=False)


def reference_in_input_frame(reference):
    """Translate separately indexed annotations before using them as anchors.

    Preserve element, charge, isotope, hydrogen count and exact bond order.
    Stereochemical labels do not enter MAPPA's endpoint graph representation.
    """
    if reference.get('input_reference_index_agrees', True):
        return reference, None
    from slap_audit import paired_graph, equivalent

    def structural(features):
        return [dict(heavy=f['heavy'], colors=[tuple(c[:4]) for c in f['colors']],
                     bonds=[(a, b, (color[0],)) for a, b, color in f['bonds']])
                for f in features]

    transport = equivalent(paired_graph(structural(reference['features']), []),
                           paired_graph(structural(reference['input_features']), []))
    if transport is None:
        raise ValueError('Reference endpoint graph differs from the prepared input')
    mapping = [(transport[0, a][1], transport[1, b][1]) for a, b in reference['mapping']]
    normalized = dict(reference, mapping=mapping, features=reference['input_features'])
    return normalized, dict(original_mapping=reference['mapping'], input_mapping=mapping,
                            endpoint_transports=sorted(transport.items()))


def complete_reference(problem, reference, *, reverse=False, query_seconds=10.):
    """Preserve annotated pairs and audit an actual minimum-window witness."""
    anchors = tuple((b, a) if reverse else (a, b) for a, b in reference['mapping'])
    timing = {}

    def measure(stage, operation):
        cpu, wall = time.process_time(), time.perf_counter()
        result = operation()
        timing[stage] = dict(cpu_seconds=time.process_time()-cpu,
                             wall_seconds=time.perf_counter()-wall)
        return result

    def matching():
        raw = search_aam(problem, AAMSearchConfig(anchors=anchors, **SETTINGS), workers=1)
        assert raw.metrics.cut_count == 1
        return deduplicate(raw), asdict(raw.metrics)

    catalogue, metrics = measure('matching', matching)
    result = dict(timing=timing, search_metrics=metrics, families=len(catalogue.families),
                  anchor_count=len(anchors), endpoint_atoms=[len(p.elements)
                      for p in (problem.reactant, problem.product)])
    if not catalogue.families:
        return dict(result, status='no_retained_family')

    def decoding():
        decoded = decode_events(catalogue, event_slack=5)
        _ = decoded.window
        return decoded

    decoded = measure('decoding', decoding)
    result.update(model_minimum=decoded.minimum_events,
                  model_lower_bound=decoded.minimum_lower_bound,
                  model_minimum_proven=decoded.minimum_proven)

    def witness():
        deadline = time.perf_counter()+query_seconds
        for context in decoded.window.contexts:
            if decoded.window.bounds[context] > decoded.minimum_events:
                continue
            remaining = deadline-time.perf_counter()
            if remaining <= 0:
                break
            family, partition, cut = context
            answer = decoded.select(family, decoded.minimum_events,
                                    partition=partition, cut=cut).find(seconds=remaining)
            if answer.status == 'found':
                # Exact anchors, full joint replay and scalar model/physical
                # event counts are checked by the existing independent audit.
                return audit_decoded_witness(decoded, reference, context, answer,
                                             reverse=reverse, unmatched_heavy=False)
        return None

    evidence = measure('witness_audit', witness)
    if evidence is None:
        return dict(result, status='minimum_witness_not_found_within_allowance')
    mapping = dict(evidence['mapping'])
    assert all(mapping.get(a) == b for a, b in anchors)
    assert len(mapping) == len(set(mapping.values()))
    assert all(problem.reactant.elements[a] == problem.product.elements[b]
               for a, b in mapping.items())
    extra = sorted((a, b) for a, b in mapping.items()
                   if (a, b) not in anchors and problem.reactant.elements[a] != 'H')
    if reverse:
        mapping = {b: a for a, b in mapping.items()}
        extra = sorted((b, a) for a, b in extra)
    result.update(status='completed', all_anchors_preserved=True,
                  mapped_atoms=len(mapping), complete_endpoints=(len(mapping) ==
                      len(problem.reactant.elements) == len(problem.product.elements)),
                  mapped_hydrogens=sum(problem.reactant.elements[a] == 'H'
                      for a, _ in evidence['mapping']),
                  additional_heavy_pairs=extra, mapping=sorted(mapping.items()),
                  physical_event_count=evidence['physical_event_count'],
                  model_event_count=evidence['event_count'], evidence=evidence)
    return result


def complete(payload, workspace):
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != payload['runner_sha256']:
        raise ValueError('Anchored completion runner changed after launch')
    row = payload['row']
    problem = problem_from_record(read_blob(Blob(**row['input']), workspace.scratch))
    reference = read_blob(Blob(**row['reference']), workspace.scratch)
    reference, transport = reference_in_input_frame(reference)
    result = complete_reference(problem, reference, reverse=row['reversed'],
                                query_seconds=payload['query_seconds'])
    result['reference_index_translation'] = transport
    return result


def run(manifest, output, *, workers=8, slurm=False, task_seconds=180., memory_gib=4.):
    source = json.loads(manifest.read_bytes())
    rows = source['directions']
    keys = [f"{r['dataset']}/{r['case']:05d}/{'reverse' if r['reversed'] else 'forward'}" for r in rows]
    assert len(keys) == len(set(keys))
    counts = Counter((r['dataset'], r['case']) for r in rows)
    assert all(n == 2 for n in counts.values())
    assert all({r['reversed'] for r in rows if (r['dataset'], r['case']) == key} == {False, True}
               for key in counts)
    output.mkdir(parents=True, exist_ok=False)
    runner = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    report = dict(schema='mappa.anchored_completion/v1', code_identity=code_identity(),
        runner_sha256=runner, manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest(),
        settings=dict(SETTINGS, event_slack=5), reactions=len(counts), directions=len(rows),
        scope='Reference-conditioned completion; not independent benchmark recovery or an original-branch offset',
        resources=dict(workers=workers, task_seconds=task_seconds, memory_gib=memory_gib),
        terminal_complete=False, results=[])
    _publish(output/'results.json', report)
    store = BundleStore(output/'bundles')
    factory = WorkPool.slurm if slurm else WorkPool.local
    started = time.perf_counter()
    with factory(workers=workers, store=store, memory_bytes=int(memory_gib*2**30),
                 case_seconds=task_seconds, shared_cache=False) as pool:
        report['ready_workers'] = len(pool.queue.ready_workers())
        pending = dict(zip(keys, rows))
        pool.queue.add([Task(key, 'anchored_completion:complete',
            dict(row=row, runner_sha256=runner, query_seconds=10.), case=key,
            case_seconds=task_seconds, memory_bytes=int(memory_gib*2**30))
            for key, row in pending.items()])
        published = 0.
        while pending:
            pool.check_workers()
            for key, state in pool.queue.results(tuple(pending)).items():
                if state['state'] in ('pending', 'running'):
                    continue
                original = pending.pop(key)
                measured = state['result'] or {}
                row = dict(dataset=original['dataset'], case=original['case'],
                    reversed=original['reversed'], prior_status=original['prior_status'],
                    input_sha256=original['input']['sha256'], reference_sha256=original['reference']['sha256'],
                    state=state['state'], worker_metrics={k: v for k, v in measured.items()
                        if k.endswith(('_seconds', '_bytes'))})
                if state['state'] == 'complete':
                    row.update(store.load(measured['blob']))
                else:
                    row.update(status='task_failed', error=measured.get('error'))
                report['results'].append(row)
                if state['state'] == 'error':
                    pool.queue.finish_input(cancel_pending=True)
                    _publish(output/'results.json', report)
                    raise RuntimeError('Anchored completion failed: '+key)
            if time.perf_counter() >= published or not pending:
                report['elapsed_seconds'] = time.perf_counter()-started
                _publish(output/'results.json', report)
                print(json.dumps(dict(completed=len(report['results']), expected=len(rows))), flush=True)
                published = time.perf_counter()+30
            if pending:
                time.sleep(.1)
    report.update(terminal_complete=True, elapsed_seconds=time.perf_counter()-started,
                  status_counts=dict(Counter(r['status'] for r in report['results'])))
    _publish(output/'results.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--slurm', action='store_true')
    parser.add_argument('--task-seconds', type=float, default=180.)
    parser.add_argument('--memory-gib', type=float, default=4.)
    args = parser.parse_args()
    run(args.manifest, args.output, workers=args.workers, slurm=args.slurm,
        task_seconds=args.task_seconds, memory_gib=args.memory_gib)
