"""Independent reference membership and certificate replay for saved decoded outputs."""
import argparse
from collections import Counter, defaultdict
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time

from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob
from mappa.mappa_core.utils.io.records import problem_from_record


def _publish(path, report):
    pending = path.with_suffix('.pending')
    pending.write_text(json.dumps(report, indent=2)+'\n')
    pending.replace(path)


def verifier_identity():
    digest = hashlib.sha256()
    for name in ('decoded_reference.py', 'decoded_constraints.py'):
        digest.update(name.encode())
        digest.update(Path(__file__).with_name(name).read_bytes())
    return digest.hexdigest()


def verify_decoded(decoded, reference, *, reverse=False, seconds=10., preferred_family=None,
                   unmatched_heavy=True, max_conflicts=10000):
    """Test exact annotated indices in the retained event window, separately from catalogue membership.

    This diagnostic does not relax reference indices by endpoint symmetry.
    The symbolic query enforces saved ordered symmetry, occupancy and the event window.
    """
    from mappa.mappa_core.post_processing.decoder.events import AssignedEventBound, prepare_boundary_events
    from decoded_constraints import ProjectedChoices

    mapping = dict(reference['mapping'])
    if reverse:
        mapping = {p: r for r, p in mapping.items()}
    # A prior catalogue certificate is only an ordering hint. Its mapping is
    # never accepted without the independent saved-decoder query and replay.
    contexts = sorted(decoded.window.contexts, key=lambda c: c[0] != preferred_family)
    started = time.perf_counter()
    bound = AssignedEventBound(decoded.index, mapping)
    filtering = dict(contexts=0, rejected=0, symmetry_queries=0, seconds=time.perf_counter()-started,
                     solver_limited_queries=0)
    unknown = False
    for family, partition, cut in contexts:
        mapped = dict(decoded.catalogue.families[family].mapping)
        if unmatched_heavy and set(mapped).intersection(reference['features'][int(reverse)]['heavy']) != set(mapping):
            continue
        remaining = seconds - (time.perf_counter() - started)
        if remaining <= 0:
            unknown = True
            break
        filtering['contexts'] += 1
        filter_started = time.perf_counter()
        boundary = prepare_boundary_events(decoded.index, mapped,
            decoded.inventory.uses[family][partition], decoded.cut_contexts[family][partition][cut])
        outside = bound.count(boundary) > decoded.event_window
        filtering['seconds'] += time.perf_counter() - filter_started
        if outside:
            filtering['rejected'] += 1
            continue
        # A diagnostic scan owns one working table at a time. Public select()
        # intentionally retains tables for interactive reuse, not this scan.
        key = decoded._table_key(family, partition, cut)
        table = decoded._prepare_table(family, partition, key, boundary=boundary)
        remaining = seconds - (time.perf_counter() - started)
        if remaining <= 0:
            unknown = True
            break
        filtering['symmetry_queries'] += 1
        try:
            query = ProjectedChoices(table, (decoded.minimum_lower_bound, decoded.event_window),
                                     started + seconds, assignments=mapping)
            answer = query.find(max_conflicts=max_conflicts)
        except TimeoutError:
            filtering['solver_limited_queries'] += 1
            unknown = True
            continue
        if answer is not None:
            evidence = audit_decoded_witness(decoded, reference, (family, partition, cut),
                                             answer, reverse=reverse, unmatched_heavy=unmatched_heavy)
            return dict(status='recovered', event_count=answer.event_count, evidence=evidence,
                        event_filter=filtering,
                        policy='exact annotated indices in selected decoded window')
        del table
    return dict(status='unknown' if unknown else 'not_recovered',
                event_filter=filtering,
                policy='exact annotated indices in selected decoded window')


def audit_decoded_witness(decoded, reference, context, answer, *, reverse=False, unmatched_heavy=True):
    """Replay a positive query proof and count its events from the original matrices.

    Catalogue membership alone is insufficient: every local discovery must be
    the saved decoder start. Replay checks the correlated binary symmetries,
    anchors and ordered placement groups without asking another solver.
    """
    from itertools import combinations
    from mappa.mappa_core.utils.bonds import bond_event_threshold

    def require(condition, reason):
        if not condition:
            raise ValueError('Invalid decoded verification certificate: ' + reason)

    require(answer.status == 'found', 'query did not return a witness')
    require(context in decoded.window.contexts, 'context is outside the retained window')
    family_id, partition, cut = context
    family = decoded.catalogue.families[family_id]
    require(len(answer.choices) == len(family.actions), 'action count differs')
    for action, choice in zip(family.actions, answer.choices):
        if action[0] == 'relation':
            selected = decoded.local.prepare(action[1]).optimization.discovery_witness
            require(dict(choice['discovery']) == dict(selected), 'different starting orbit')
    mapping = decoded.realize(family_id, answer.choices)
    require(tuple(sorted(mapping.items())) == tuple(sorted(answer.witness)), 'replay differs')
    require(len(set(mapping.values())) == len(mapping), 'duplicate target occupancy')
    expected = dict(reference['mapping'])
    if reverse:
        expected = {p: r for r, p in expected.items()}
    if unmatched_heavy:
        features = reference['features'][::-1] if reverse else reference['features']
        require(set(mapping).intersection(features[0]['heavy']) == set(expected), 'unmatched source atoms differ')
        require(set(mapping.values()).intersection(features[1]['heavy']) == set(expected.values()),
                'unmatched target atoms differ')
    require(all(mapping.get(a) == b for a, b in expected.items()), 'annotated correspondence differs')

    # Scalar counting is independent of compressed response tables and query
    # cost projections. Keep model counts distinct from whole-reaction counts.
    index = decoded.index
    source, target = index.problem.reactant, index.problem.product
    def threshold(endpoint, a, b):
        return bond_event_threshold(endpoint.elements, a, b,
            default_threshold=index.threshold, metal_threshold=index.metal_threshold)
    def changed(a, b):
        product_weight = target.wbo[mapping[a], mapping[b]] if a in mapping and b in mapping else 0.
        return bool(abs(float(product_weight) - float(source.wbo[a, b])) > threshold(source, a, b))
    pieces = decoded.fragment_partitions[family_id][partition]
    owners = {a: number for number, atoms in enumerate(pieces) for a in atoms}
    require(len(owners) == sum(map(len, pieces)) and set(owners) == set(mapping), 'invalid fragment partition')
    target_owners = {p: owners[a] for a, p in family.mapping}
    source_boundary = {(a, b) for a, b in combinations(range(index.n), 2)
        if owners.get(a, -1) != owners.get(b, -1) and source.wbo[a, b] > threshold(source, a, b)}
    target_boundary = {(a, b) for a, b in combinations(range(index.m), 2)
        if target_owners.get(a, -1) != target_owners.get(b, -1)
        and target.wbo[a, b] > threshold(target, a, b)}
    model_count = len(source_boundary) + len(target_boundary)
    model_count += sum(changed(a, b) for atoms in pieces for a, b in combinations(atoms, 2))
    cuts = decoded.cut_contexts[family_id][partition][cut]
    for a, b in cuts:
        if a in owners and b in owners and owners[a] != owners[b]:
            model_count += int(changed(a, b)) - int(tuple(sorted((a, b))) in source_boundary)
            model_count -= int(tuple(sorted((mapping[a], mapping[b]))) in target_boundary)
    images = set(mapping.values())
    physical_count = sum(changed(a, b) for a, b in combinations(range(index.n), 2)
                         if a in mapping or b in mapping)
    physical_count += sum(target.wbo[a, b] > threshold(target, a, b)
        for a, b in combinations(range(index.m), 2) if (a in images) != (b in images))
    require(model_count == answer.event_count, 'query and scalar model counts differ')
    require(decoded.minimum_lower_bound <= model_count <= decoded.event_window, 'outside event window')
    return dict(family=family_id, partition=partition, cut=cut, mapping=tuple(sorted(mapping.items())),
        choices=answer.choices, event_count=int(model_count), physical_event_count=int(physical_count),
        event_window=(decoded.minimum_lower_bound, decoded.event_window),
        certificate='replayed_saved_decoder_choices_with_scalar_event_count')


def replay_catalogue_witness(decoded, reference, certificate, *, reverse=False, unmatched_heavy=True):
    """Check an existing candidate against the saved starts, symmetries and window.

    Rebase the correspondence on either side only when the resulting self
    operation preserves that side's binary graph and anchors.
    Failure to convert this one witness says nothing about other allowed choices.
    """
    from types import SimpleNamespace
    from mappa.mappa_core.core.symmetry.actions import decision_blocks

    expected = dict(reference['mapping'])
    if reverse:
        expected = {p: r for r, p in expected.items()}
    candidate = dict(certificate['evidence']['mapping'])
    # Catalogue recovery allows endpoint-equivalent labels; this audit requires
    # the exact annotated pairs. An unsuitable hint cannot stop that query.
    if any(candidate.get(a) != b for a, b in expected.items()):
        return None
    if unmatched_heavy:
        features = reference['features'][::-1] if reverse else reference['features']
        if (set(candidate).intersection(features[0]['heavy']) != set(expected)
                or set(candidate.values()).intersection(features[1]['heavy']) != set(expected.values())):
            return None
    family_id = certificate['family']
    contexts = [key for key in decoded.window.contexts if key[0] == family_id]
    if not contexts:
        return None
    family = decoded.catalogue.families[family_id]
    evidence = certificate['evidence']
    actions = {(a['transition'], a['kind']): a for a in evidence['actions']}
    blocks = decision_blocks(family.actions)
    transitions = [len(blocks)-1-i for i, block in enumerate(blocks) for _ in block]
    if len(actions) != len(evidence['actions']) or set(actions) - set(zip(transitions, (a[0] for a in family.actions))):
        raise ValueError('Catalogue certificate action identifiers differ')
    choices = []
    for (kind, fragment), transition in zip(family.actions, transitions):
        action = actions.get((transition, kind))
        if kind == 'relation':
            selected = dict(decoded.local.prepare(fragment).optimization.discovery_witness)
            inverse = {p: r for r, p in selected.items()}
            positions = {p: i for i, p in enumerate(fragment.target_atoms)}
            # Certificates omit identity operations. For a relation, identity
            # means its original discovery and identity self permutations.
            discovery = dict(fragment.witness if action is None else action['discovery'])
            beta = tuple(range(len(fragment.target_atoms))) if action is None else action['target_permutation']
            alpha = tuple(range(len(fragment.source_atoms))) if action is None else action['source_permutation']
            adjusted = tuple(beta[positions[discovery[inverse[p]]]] for p in fragment.target_atoms)
            if fragment.relation.target.contains(adjusted):
                beta = adjusted  # beta' = beta o discovery o selected^-1
            else:
                inverse = {p: r for r, p in discovery.items()}
                positions = {r: i for i, r in enumerate(fragment.source_atoms)}
                adjusted = tuple(alpha[positions[inverse[selected[r]]]] for r in fragment.source_atoms)
                if not fragment.relation.source.contains(adjusted):
                    return None
                alpha = adjusted  # alpha' = alpha o discovery^-1 o selected
            choices.append(dict(discovery=selected, source_permutation=alpha, target_permutation=beta))
        elif kind == 'group':
            choices.append(dict(permutation=tuple(range(decoded.index.m)) if action is None
                                else tuple(p for _, p in action['permutation'])))
        elif kind == 'pool':
            choices.append(tuple(fragment) if action is None else tuple(p for _, p in action['permutation']))
        else:
            raise ValueError('Unsupported catalogue certificate action')
    mapping = decoded.realize(family_id, choices)
    if mapping != dict(evidence['mapping']):
        raise ValueError('Rebased catalogue witness changed its mapping')
    vector = decoded.index.vector(mapping)
    started = time.perf_counter()
    for examined, (family_id, partition, cut) in enumerate(contexts, 1):
        table = decoded._prepare_table(family_id, partition, decoded._table_key(family_id, partition, cut))
        count = int(table.evaluate(vector)['total'])
        del table
        if decoded.minimum_lower_bound <= count <= decoded.event_window:
            answer = SimpleNamespace(status='found', witness=tuple(sorted(mapping.items())),
                                     choices=tuple(choices), event_count=count)
            proof = audit_decoded_witness(decoded, reference, (family_id, partition, cut), answer,
                                          reverse=reverse, unmatched_heavy=unmatched_heavy)
            return dict(status='recovered', event_count=count, evidence=proof,
                event_filter=dict(contexts=examined, rejected=0, symmetry_queries=0, proof_replays=1,
                    seconds=time.perf_counter()-started, solver_limited_queries=0),
                policy='exact annotated indices in selected decoded window')
    return None


def audit_saved_decoding(payload, workspace):
    """External verification of immutable output; never rerun mapping or decoding."""
    from mappa.mappa_core.utils.io.records import decoded_events_from_record
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    if verifier_identity() != payload['verifier_sha256']:
        raise ValueError('Verifier changed beneath the audit pool')
    row = payload['row']
    started = time.perf_counter()
    saved = read_blob(Blob(**row['decoded']), workspace.scratch)
    reference = read_blob(Blob(**row['reference']), workspace.scratch)
    decoded = decoded_events_from_record(saved.record)
    original = problem_from_record(read_blob(Blob(**row['input']), workspace.scratch))
    if FinalBranchCatalogue(original).problem_sha256 != decoded.matching_problem_sha256:
        raise ValueError('Decoded record belongs to different prepared endpoints')
    certificate = None
    if row['verification'] is not None:
        certificate = read_blob(Blob(**row['verification']), workspace.scratch)
    loaded = time.perf_counter()
    result = None
    if certificate is not None and certificate.get('status') == 'recovered' and 'evidence' in certificate:
        result = replay_catalogue_witness(decoded, reference, certificate, reverse=row['reversed'],
                                         unmatched_heavy=payload['unmatched_heavy'])
    if result is None:
        result = verify_decoded(decoded, reference, reverse=row['reversed'], seconds=payload['seconds'],
            preferred_family=certificate.get('family') if certificate is not None else None,
            unmatched_heavy=payload['unmatched_heavy'], max_conflicts=payload['max_conflicts'])
    result.update(decoded_sha256=row['decoded']['sha256'], reference_sha256=row['reference']['sha256'],
                  input_sha256=row['input']['sha256'], loading_seconds=loaded-started,
                  query_seconds=time.perf_counter()-loaded)
    return result


def audit_inputs(campaign):
    """Resolve a completed batch or an explicit manifest of recovered outputs."""
    import sqlite3
    if campaign.is_file():
        value = json.loads(campaign.read_bytes())
        if value['schema'] != 'mappa.decoded_reference_inputs/v1':
            raise ValueError('Unsupported decoded audit manifest')
        files = dict(value['source_files'])
        files[str(campaign)] = hashlib.sha256(campaign.read_bytes()).hexdigest()
        for path, digest in files.items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
                raise ValueError('Audit source changed: ' + path)
        return value['directions'], files, value['unmatched_heavy']
    original = {str(campaign/name): hashlib.sha256((campaign/name).read_bytes()).hexdigest()
                for name in ('index.json', 'result/plan.json', 'result/ledger.sqlite')}
    index = json.loads((campaign/'index.json').read_text())
    plan = json.loads((campaign/'result/plan.json').read_text())
    if original[str(campaign/'index.json')] != plan['input_index_sha256']:
        raise ValueError('Campaign input index changed')
    with sqlite3.connect(f'file:{campaign}/result/ledger.sqlite?mode=ro', uri=True) as db:
        records = {(name, stage): (state, reference) for name, stage, state, reference in
                   db.execute('SELECT case_name,stage,state,reference FROM results')}
    rows = []
    for row in index['directions']:
        name = f"{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}"
        state, reference = records[name, 'decoding']
        if state != 'complete':
            raise ValueError('Audit requires a saved successful decoding for ' + name)
        verification_state, verification = records.get((name, 'verification'), ('missing', None))
        decoded = json.loads(reference)
        rows.append(dict(case=row['case'], reversed=row['reversed'], reference=row['reference'], input=row['input'],
                         decoded=decoded, verification=json.loads(verification) if verification_state == 'complete'
                         else None, name=name))
    return rows, original, True


def audit_decoding(campaign, output, *, cases=None, workers=1, seconds=10., task_seconds=60.,
                   memory_bytes=2**31, unknown_from=None, slurm=False, max_conflicts=10000):
    """Audit immutable saved outputs in one guarded pool and aggregated evidence bundle."""
    from mappa.execution import BundleStore, WorkPool
    from mappa.execution.identity import code_identity
    from mappa.mappa_core.core.scheduler import Task

    if workers < 1 or memory_bytes < 1 or max_conflicts < 1 or not 0 <= seconds < task_seconds:
        raise ValueError('Require positive workers/memory and 0 <= query seconds < task seconds')
    campaign, output = Path(campaign).resolve(), Path(output).resolve()
    if output == campaign or campaign in output.parents:
        raise ValueError('Use an audit output outside the immutable campaign')
    rows, original, unmatched_heavy = audit_inputs(campaign)
    if any(Path(name).is_relative_to(output) for name in original):
        raise ValueError('Audit output overlaps immutable source files')
    previous = None
    if unknown_from is not None:
        previous = json.loads(Path(unknown_from).read_bytes())
        if previous['original_hashes'] != original or previous['certificate_failures']:
            raise ValueError('Continuation requires the same unchanged source and no certificate failures')
        unresolved = {r['name']: r for r in previous['results'] if r['status'] == 'unknown'}
        rows = [row for row in rows if row['name'] in unresolved]
        if any(unresolved[row['name']][key+'_sha256'] != row[key]['sha256']
               for row in rows for key in ('decoded', 'reference', 'input')):
            raise ValueError('Continuation input identity differs')
    if cases is not None:
        rows = [row for row in rows if row['case'] in cases]
    if not rows or len({r['name'] for r in rows}) != len(rows):
        raise ValueError('Require unique, nonempty direction scope')
    if cases is not None and set(cases) != {r['case'] for r in rows}:
        raise ValueError('Requested cases are missing')
    output.mkdir(parents=True, exist_ok=False)
    report = dict(schema='mappa.decoded_reference_audit/v1', code_identity=code_identity(),
        verifier_sha256=verifier_identity(),
        original_hashes=original, expected=len(rows), complete=False,
        policy='exact annotated indices in the saved decoded window; no endpoint-symmetry relaxation',
        unmatched_heavy=unmatched_heavy,
        hydrogen_policy='unannotated; must admit one complete injective witness',
        event_count_policy='saved independent-fragment model; physical count reported separately',
        resources=dict(workers=workers, query_seconds=seconds, task_seconds=task_seconds,
                       memory_bytes=memory_bytes, query_max_conflicts=max_conflicts), results=[])
    if previous is not None:
        report['previous_audit_sha256'] = hashlib.sha256(Path(unknown_from).read_bytes()).hexdigest()
    store = BundleStore(output/'bundles')
    started = time.perf_counter()
    factory = WorkPool.slurm if slurm else WorkPool.local
    with factory(workers=workers, store=store, memory_bytes=memory_bytes,
                        case_seconds=task_seconds, shared_cache=False) as pool:
        report.update(ready_workers=len(pool.queue.ready_workers()),
                      startup_seconds=pool.startup_elapsed_seconds)
        pending = {row['name']: row for row in rows}
        pool.queue.add([Task(name, 'decoded_reference:audit_saved_decoding',
            dict(row=row, seconds=seconds, max_conflicts=max_conflicts, unmatched_heavy=unmatched_heavy,
                 verifier_sha256=report['verifier_sha256']),
            case=name, case_seconds=task_seconds, memory_bytes=memory_bytes) for name, row in pending.items()])
        next_publish = 0.
        while pending:
            pool.check_workers()
            for name, state in pool.queue.results(tuple(pending)).items():
                if state['state'] in ('pending', 'running'):
                    continue
                row = pending.pop(name)
                measured = state['result'] or {}
                result = dict(name=name, case=row['case'], reversed=row['reversed'], state=state['state'],
                    attempt_audit=asdict(pool.queue.archive_case(name, store)),
                    input_sha256=row['input']['sha256'],
                    decoded_sha256=row['decoded']['sha256'], reference_sha256=row['reference']['sha256'],
                    metrics={k: v for k, v in measured.items() if k.endswith(('_seconds', '_bytes', '_samples'))})
                if state['state'] == 'complete':
                    value = store.load(measured['blob'])
                    result.update(status=value['status'], certificate_blob=asdict(measured['blob']),
                                  event_filter=value['event_filter'], loading_seconds=value['loading_seconds'],
                                  query_seconds=value['query_seconds'])
                    if value['status'] == 'recovered':
                        proof = value['evidence']
                        result.update({key: proof[key] for key in (
                            'family', 'partition', 'cut', 'event_count', 'physical_event_count', 'event_window')})
                else:
                    result.update(status='unknown', error=measured.get('error'))
                    if state['state'] == 'error':
                        pool.queue.finish_input(cancel_pending=True)
                report['results'].append(result)
            if time.perf_counter() >= next_publish or not pending:
                report.update(elapsed_seconds=time.perf_counter()-started,
                    status_counts=dict(Counter(r['status'] for r in report['results'])))
                _publish(output/'results.json', report)
                pool.checkpoint()
                print(json.dumps(dict(completed=len(report['results']), expected=len(rows),
                                      status_counts=report['status_counts'])), flush=True)
                next_publish = time.perf_counter()+30
            if pending:
                time.sleep(.1)
    if any(hashlib.sha256(Path(name).read_bytes()).hexdigest() != digest for name, digest in original.items()):
        raise ValueError('Source campaign changed during the audit')
    by_case = defaultdict(list)
    for row in report['results']:
        by_case[row['case']].append(row['status'])
    report.update(complete=all(r['state'] == 'complete' for r in report['results']), source_unchanged=True,
        reaction_counts=dict(Counter('recovered' if 'recovered' in statuses else
            'unknown' if 'unknown' in statuses else 'not_recovered' for statuses in by_case.values())),
        certificate_failures=sum(r['state'] == 'error' for r in report['results']))
    _publish(output/'results.json', report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', type=Path, help='Completed Golden batch or explicit decoded-output manifest')
    parser.add_argument('output', type=Path)
    parser.add_argument('--cases', nargs='+', type=int)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--seconds', type=float, default=10.,
                        help='Search allowance after candidate checking; zero checks existing candidates only')
    parser.add_argument('--task-seconds', type=float, default=60.)
    parser.add_argument('--max-conflicts', type=int, default=10000,
                        help='Explicit solver-conflict allowance per decoded query; reaching it remains unknown')
    parser.add_argument('--memory-gib', type=float, default=2.)
    parser.add_argument('--unknown-from', type=Path)
    parser.add_argument('--slurm', action='store_true')
    args = parser.parse_args()
    report = audit_decoding(args.inputs, args.output, cases=args.cases, workers=args.workers,
        seconds=args.seconds, task_seconds=args.task_seconds, memory_bytes=int(args.memory_gib*2**30),
        unknown_from=args.unknown_from, slurm=args.slurm, max_conflicts=args.max_conflicts)
    return int(report['certificate_failures'] > 0)


if __name__ == '__main__':
    raise SystemExit(main())

