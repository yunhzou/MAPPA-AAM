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


def reference_view_cases(index, cases, *, geometry_cache):
    """Load only selected original inputs/references; never load scientific results."""
    import numpy as np
    from mappa.presentation.geometry import display_endpoint
    from mappa.presentation.trajectory import bond_events

    rows = {row['case']: row for row in json.loads(Path(index).read_text())['directions']
            if not row['reversed']}
    documents = []
    for case in dict.fromkeys(cases):
        row = rows[case]
        raw = read_blob(Blob(**row['input']), None)
        reference = read_blob(Blob(**row['reference']), None)
        endpoints, indices, geometry = [], [], []
        for side, feature in zip(('reactant', 'product'), reference['features']):
            endpoint = display_endpoint(raw[side], cache=Path(geometry_cache))
            heavy = feature['heavy']
            endpoints.append(dict(elements=[endpoint['elements'][i] for i in heavy],
                coordinates=[endpoint['coordinates'][i] for i in heavy],
                wbo=np.asarray(endpoint['wbo'])[np.ix_(heavy, heavy)].tolist(),
                atom_indices=heavy))
            indices.append({original: dense for dense, original in enumerate(heavy)})
            geometry.append(endpoint['geometry'])
        mapping = {indices[0][a]: indices[1][b] for a, b in reference['mapping']}
        events = bond_events(dict(reactant=endpoints[0], product=endpoints[1]), mapping)
        documents.append(dict(name=f'Golden {case}', index=case, endpoints=endpoints,
            note='Ground-truth heavy-atom mapping, in original forward orientation. '
                 'Unmatched atoms are retained; hydrogen correspondence is unannotated. '
                 'Generated geometry, when needed, uses ETKDGv3 + UFF for display only; '
                 'unspecified stereochemistry is illustrative. Separate molecules are spaced for display.',
            records=[dict(name='Ground-truth mapping', mapping=list(mapping.items()), events=events,
                provenance=dict(kind='ground_truth_reference', geometry=geometry,
                    input_sha256=row['input']['sha256'], reference_sha256=row['reference']['sha256']))]))
    return documents


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
        evidence = _atom_set_reference(family, catalogue.problem, mapping, features[0]['heavy'])
        if evidence is not None:
            return dict(status='recovered', family=number, evidence=evidence,
                        heavy_atoms=len(mapping), policy=REFERENCE_POLICY)
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


def _atom_set_reference(family, problem, reference, heavy):
    """Construct and validate a direct proof on disjoint completed atom sets.

    Unspecified atoms are filled injectively by element. A failed construction
    says nothing about membership through global placement or endpoint symmetry;
    those choices still require the ordinary joint query.
    """
    if set(dict(family.mapping)).intersection(heavy) != set(reference):
        return None
    discoveries, occupied = {}, set()
    for fragment in family.fragments:
        if fragment.accepted_edges or occupied.intersection(fragment.source_atoms):
            return None
        occupied.update(fragment.source_atoms)
        mapping = dict(fragment.fixed_assignments)
        for atom in fragment.source_atoms:
            if atom in reference:
                if atom in mapping and mapping[atom] != reference[atom]:
                    return None
                mapping[atom] = reference[atom]
        if len(set(mapping.values())) != len(mapping) or not set(mapping.values()) <= set(fragment.target_atoms):
            return None
        remaining = {}
        for atom in reversed(fragment.target_atoms):
            if atom not in mapping.values():
                remaining.setdefault(problem.product.elements[atom], []).append(atom)
        for atom in fragment.source_atoms:
            if atom not in mapping:
                options = remaining.get(problem.reactant.elements[atom])
                if not options:
                    return None
                mapping[atom] = options.pop()
        if not fragment.relation.allows_discovery(mapping):
            return None
        discoveries[fragment] = dict(discovery=mapping)
    path = family.as_path(problem)
    choices = {edge: discoveries[path.graph.transitions[edge].match] for edge in path.transitions
               if path.graph.transitions[edge].match is not None}
    realized = path.realize(choices)
    actual = dict(realized.mapping)
    if (len(set(actual.values())) != len(actual) or any(actual.get(r) != p for r, p in reference.items())):
        return None
    return dict(mapping=realized.mapping, fragment_choices=realized.fragment_choices,
                placement_choices=realized.placement_choices, scope='full_explicit',
                certificate='validated_atom_set_correspondence')


def verify_decoded(decoded, reference, *, reverse=False, seconds=10., preferred_family=None):
    """Test exact annotated indices in the retained event window, separately from catalogue membership.

    This diagnostic does not relax reference indices by endpoint symmetry.
    The shared query enforces local symmetry, occupancy and the event window.
    """
    from mappa.mappa_core.post_processing.decoder.events import AssignedEventBound, prepare_boundary_events

    mapping = dict(reference['mapping'])
    feature = reference['features'][int(reverse)]
    if reverse:
        mapping = {p: r for r, p in mapping.items()}
    # A prior catalogue certificate is only an ordering hint. Its mapping is
    # never accepted without the independent saved-decoder query and replay.
    contexts = sorted(decoded.window.contexts, key=lambda c: c[0] != preferred_family)
    started = time.perf_counter()
    bound = AssignedEventBound(decoded.index, mapping)
    filtering = dict(contexts=0, rejected=0, symmetry_queries=0, seconds=time.perf_counter()-started)
    unknown = False
    for family, partition, cut in contexts:
        mapped = dict(decoded.catalogue.families[family].mapping)
        if set(mapped).intersection(feature['heavy']) != set(mapping):
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
        table = decoded._prepare_table(family, partition, key)
        remaining = seconds - (time.perf_counter() - started)
        if remaining <= 0:
            unknown = True
            break
        filtering['symmetry_queries'] += 1
        answer = table.program.select(decoded.minimum_lower_bound, decoded.event_window).find(
            assignments=mapping, seconds=remaining)
        if answer.status == 'found':
            evidence = audit_decoded_witness(decoded, reference, (family, partition, cut),
                                             answer, reverse=reverse)
            return dict(status='recovered', event_count=answer.event_count, evidence=evidence,
                        event_filter=filtering,
                        policy='exact annotated indices in selected decoded window')
        del table
        unknown |= answer.status == 'unknown'
    return dict(status='unknown' if unknown else 'not_recovered',
                event_filter=filtering,
                policy='exact annotated indices in selected decoded window')


def audit_decoded_witness(decoded, reference, context, answer, *, reverse=False):
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
    features = reference['features']
    if reverse:
        expected = {p: r for r, p in expected.items()}
        features = features[::-1]
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


def audit_saved_decoding(payload, workspace):
    """External verification of immutable output; never rerun mapping or decoding."""
    from mappa.mappa_core.utils.io.records import decoded_events_from_record
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != payload['verifier_sha256']:
        raise ValueError('Verifier changed beneath the audit pool')
    row = payload['row']
    saved = read_blob(Blob(**row['decoded']), workspace.scratch)
    reference = read_blob(Blob(**row['reference']), workspace.scratch)
    decoded = decoded_events_from_record(saved.record)
    original = problem_from_record(read_blob(Blob(**row['input']), workspace.scratch))
    if FinalBranchCatalogue(original).problem_sha256 != decoded.matching_problem_sha256:
        raise ValueError('Decoded record belongs to different prepared endpoints')
    catalogue_check = read_blob(Blob(**row['verification']), workspace.scratch)
    result = verify_decoded(decoded, reference, reverse=row['reversed'], seconds=payload['seconds'],
                            preferred_family=catalogue_check.get('family'))
    result.update(decoded_sha256=row['decoded']['sha256'], reference_sha256=row['reference']['sha256'],
                  input_sha256=row['input']['sha256'])
    return result


def audit_decoding(campaign, output, *, cases=None, workers=1, seconds=10., task_seconds=60.,
                   memory_bytes=2**31, unknown_from=None):
    """Audit a completed single batch in a separate pool and evidence bundle."""
    import sqlite3
    from mappa.execution import BundleStore, WorkPool
    from mappa.execution.identity import code_identity
    from mappa.mappa_core.core.scheduler import Task

    if workers < 1 or not 0 < seconds < task_seconds:
        raise ValueError('Require positive workers and 0 < query seconds < task seconds')
    campaign, output = Path(campaign).resolve(), Path(output).resolve()
    if output == campaign or campaign in output.parents:
        raise ValueError('Use an audit output outside the immutable campaign')
    original = {name: hashlib.sha256((campaign/name).read_bytes()).hexdigest()
                for name in ('index.json', 'result/plan.json', 'result/ledger.sqlite')}
    index = json.loads((campaign/'index.json').read_text())
    plan = json.loads((campaign/'result/plan.json').read_text())
    if original['index.json'] != plan['input_index_sha256']:
        raise ValueError('Campaign input index changed')
    with sqlite3.connect(f'file:{campaign}/result/ledger.sqlite?mode=ro', uri=True) as db:
        records = {(name, stage): (state, reference) for name, stage, state, reference in
                   db.execute('SELECT case_name,stage,state,reference FROM results')}
    previous = None
    if unknown_from is not None:
        previous = json.loads(Path(unknown_from).read_bytes())
        if previous['original_hashes'] != original or previous['certificate_failures']:
            raise ValueError('Continuation requires the same unchanged source and no certificate failures')
        unresolved = {r['name']: r for r in previous['results'] if r['status'] == 'unknown'}
    rows = []
    for row in index['directions']:
        if cases is not None and row['case'] not in cases:
            continue
        name = f"{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}"
        if previous is not None and name not in unresolved:
            continue
        state, reference = records[name, 'decoding']
        if state != 'complete':
            raise ValueError('Audit requires a saved successful decoding for ' + name)
        verification_state, verification = records[name, 'verification']
        if verification_state != 'complete':
            raise ValueError('Audit requires a saved catalogue check for ' + name)
        decoded = json.loads(reference)
        if previous is not None and any(unresolved[name][key+'_sha256'] != value['sha256']
            for key, value in (('decoded', decoded), ('reference', row['reference']), ('input', row['input']))):
            raise ValueError('Continuation input identity differs')
        rows.append(dict(case=row['case'], reversed=row['reversed'], reference=row['reference'], input=row['input'],
                         decoded=decoded, verification=json.loads(verification), name=name))
    if not rows or len({r['name'] for r in rows}) != len(rows):
        raise ValueError('Require unique, nonempty direction scope')
    if cases is not None and set(cases) != {r['case'] for r in rows}:
        raise ValueError('Requested cases are missing')
    output.mkdir(parents=True, exist_ok=False)
    report = dict(schema='mappa.decoded_reference_audit/v1', code_identity=code_identity(),
        verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        original_hashes=original, expected=len(rows), complete=False,
        policy='exact annotated heavy-atom indices, including unmatched atoms, in the saved decoded window',
        hydrogen_policy='unannotated; must admit one complete injective witness',
        event_count_policy='saved independent-fragment model; physical count reported separately',
        resources=dict(workers=workers, query_seconds=seconds, task_seconds=task_seconds,
                       memory_bytes=memory_bytes), results=[])
    if previous is not None:
        report['previous_audit_sha256'] = hashlib.sha256(Path(unknown_from).read_bytes()).hexdigest()
    store = BundleStore(output/'bundles')
    started = time.perf_counter()
    with WorkPool.local(workers=workers, store=store, memory_bytes=memory_bytes,
                        case_seconds=task_seconds, shared_cache=False) as pool:
        pending = {row['name']: row for row in rows}
        pool.queue.add([Task(name, 'golden:audit_saved_decoding',
            dict(row=row, seconds=seconds, verifier_sha256=report['verifier_sha256']),
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
                    input_sha256=row['input']['sha256'],
                    decoded_sha256=row['decoded']['sha256'], reference_sha256=row['reference']['sha256'],
                    metrics={k: v for k, v in measured.items() if k.endswith(('_seconds', '_bytes', '_samples'))})
                if state['state'] == 'complete':
                    value = store.load(measured['blob'])
                    result.update(status=value['status'], certificate_blob=asdict(measured['blob']),
                                  event_filter=value['event_filter'])
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
                print(json.dumps(dict(completed=len(report['results']), expected=len(rows),
                                      status_counts=report['status_counts'])), flush=True)
                next_publish = time.perf_counter()+30
            if pending:
                time.sleep(.1)
    if any(hashlib.sha256((campaign/name).read_bytes()).hexdigest() != digest
           for name, digest in original.items()):
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
    families = [_digest((f.policy, f.mapping,
                         tuple((kind, data.to_record() if kind == 'relation' else data)
                               for kind, data in f.actions)))
                for f in catalogue.families]
    branches = [_digest((key, sorted(families[i] for i in ids)))
                for key, ids in catalogue.branches.items()]
    coupled = [_digest((families[i], sorted(partitions)))
               for i, partitions in catalogue.coupled_branches.items()]
    return _digest(('catalogue-signature/v4', sorted(families),
                    sorted(branches), sorted(coupled)))


def _match_direction(problem, config, *, key, cache=True, progress=None,
                     audit_catalogue=False, metrics=None):
    """One direction owns its cache and reduces completed cuts as they finish."""
    from mappa.mappa_core.core.search.branching.reuse import FragmentRepair
    from mappa.mappa_core.core.search.planning import cut_sweep_items
    from mappa.mappa_core.core.search.session import SearchSession
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue

    metrics = {} if metrics is None else metrics
    with _measure(metrics, 'preparation'):
        session = SearchSession(problem, config)
        if not cache:
            repair = session.repair
            session.repair = FragmentRepair(repair.source, repair.target, cache_bytes=0)
            del repair
        cuts = cut_sweep_items(problem.reactant.wbo, config.cut_floor) if config.sweep_cuts else [()]
        catalogue = FinalBranchCatalogue(problem)
    result = dict(key=key, atoms=[problem.source_atom_count, problem.target_atom_count],
        expected_cuts=len(cuts), completed_cuts=0, expected_seed_searches=len(cuts)*config.seed_count,
        completed_seed_searches=0, metrics=metrics, cuts=[], stops={}, cap_stages={},
        graph_counts=dict(states=0, transitions=0, terminals=0), max_growth_candidates=0)
    scope, stops, cap_stages = hashlib.sha256(), Counter(), Counter()

    def checkpoint(phase):
        result.update(phase=phase, cache=session.repair.stats(), stops=dict(stops),
                      cap_stages=dict(cap_stages))
        if progress is not None:
            progress(result)

    for number, cut in enumerate(cuts):
        checkpoint('matching')
        with _measure(metrics, 'matching'):
            graph, counts = session.search_cut(cut)
        checkpoint('reduction')
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
    result['catalogue_audited'] = audit_catalogue
    if audit_catalogue:
        checkpoint('audit')
        with _measure(metrics, 'audit'):
            result['catalogue_sha256'] = _catalogue_signature(catalogue)
    result.update(scope_sha256=scope.hexdigest(), families=len(catalogue.families),
                  atom_set_branches=len(catalogue.branches), retained_paths=catalogue.path_count)
    checkpoint('matching_complete')
    return catalogue, result


def diagnose_direction(payload, workspace):
    """One cold process, one complete direction; retain its cache across cuts."""
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != payload['runner_sha256']:
        raise ValueError('Diagnostic runner changed after the run was planned')
    row, metrics = payload['row'], {}
    with _measure(metrics, 'preparation'):
        problem = problem_from_record(read_blob(Blob(**row['input']), workspace.scratch))
        reference = read_blob(Blob(**row['reference']), workspace.scratch)
        config = AAMSearchConfig(**payload['config'])
    catalogue, result = _match_direction(problem, config, key=payload['key'],
        cache=payload['cache'], audit_catalogue=payload['audit_catalogue'], metrics=metrics,
        progress=lambda value: _publish(payload['progress'], value))
    result['phase'] = 'verification'
    _publish(payload['progress'], result)
    with _measure(metrics, 'verification'):
        result['verification'] = verify_catalogue(catalogue, reference, reverse=row['reversed'])
    result['phase'] = 'complete'
    _publish(payload['progress'], result)
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
        catalogue_signature='catalogue-signature/v4' if audit_catalogue else None,
        input_indexes_sha256=[hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in indices],
        input_sha256=[dict(case=r['case'], reversed=r['reversed'], input=r['input']['sha256'],
                           reference=r['reference']['sha256']) for r in rows],
        timing_scope='Separate preparation, matching, family reduction, '
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


COMPARISON_ARMS = ('binary_iso1', 'weighted_iso10')


def comparison_config(arm):
    if arm not in COMPARISON_ARMS:
        raise ValueError('Unknown Golden comparison arm')
    return AAMSearchConfig(seed_count=1, branch_limit=2000, random_seed=42,
                          iso_tolerance=1. if arm == 'binary_iso1' else 10.)


def comparison_problem(original, arm):
    """Change matching weights only; indexed original endpoints remain separate."""
    comparison_config(arm)
    if arm == 'weighted_iso10':
        return original
    def binary(endpoint):
        import numpy as np
        weights = (endpoint.wbo >= .2).astype(float)
        np.fill_diagonal(weights, 0.)
        return replace(endpoint, wbo=weights)
    return replace(original, reactant=binary(original.reactant), product=binary(original.product))


def _comparison_runner(payload):
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != payload['runner_sha256']:
        raise ValueError('Golden comparison runner changed after planning')


def compare_match(payload, workspace):
    """Publish one reduced catalogue, without transporting individual cut graphs."""
    from mappa.mappa_core.core.scheduler.results import CatalogueResult, CutScope
    _comparison_runner(payload)
    metrics = {}
    with _measure(metrics, 'preparation'):
        original = problem_from_record(read_blob(Blob(**payload['row']['input']), workspace.scratch))
        problem = comparison_problem(original, payload['arm'])
        config = comparison_config(payload['arm'])
    catalogue, measured = _match_direction(problem, config, key=payload['name'], metrics=metrics)
    if (measured['completed_cuts'] != measured['expected_cuts'] or
            measured['completed_seed_searches'] != measured['expected_seed_searches']):
        raise ValueError('Cannot publish an incomplete comparison catalogue')
    scope = tuple(CutScope(cut['number'], tuple(map(tuple, cut['bonds'])), (0,))
                  for cut in measured.pop('cuts'))
    saved = CatalogueResult(payload['name'], catalogue, scope,
                            tuple(sorted(measured['stops'].items())), config, 'reused_native')
    with _measure(metrics, 'catalogue_publication'):
        reference = workspace.writer.put(saved)
    measured.pop('metrics')
    return dict(catalogue=asdict(reference), original_input=payload['row']['input'],
                matching_problem_sha256=catalogue.problem_sha256, metrics=metrics,
                search=measured)


def compare_decode(payload, workspace):
    from mappa import decode_events
    from mappa.mappa_core.core.scheduler.results import DecodingResult
    from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig
    from mappa.mappa_core.utils.io.records import decoded_events_record
    _comparison_runner(payload)
    metrics = {}
    with _measure(metrics, 'loading'):
        matched = workspace.load(payload['matching'])
        saved = workspace.load(Blob(**matched['catalogue']))
        original = problem_from_record(read_blob(Blob(**matched['original_input']), workspace.scratch))
    config = EventDecodeConfig()
    with _measure(metrics, 'decoding'):
        interpreted = decode_events(saved.catalogue, config, event_problem=original)
        record = decoded_events_record(interpreted)
        identities = dict(matching_problem_sha256=interpreted.matching_problem_sha256,
                          event_problem_sha256=interpreted.event_problem_sha256)
        del interpreted
    result = DecodingResult(saved.case, config, saved.scope, saved.stops,
                            record, saved.search_config, saved.execution)
    del saved, original, matched
    with _measure(metrics, 'decoded_publication'):
        reference = workspace.writer.put(result)
    return dict(decoded=asdict(reference), metrics=metrics, schema=record['schema'],
                families=len(record['families']), tables=len(record['tables']),
                fragment_responses=len(record['fragment_responses']), **identities)


def compare_verify(payload, workspace):
    _comparison_runner(payload)
    metrics = {}
    with _measure(metrics, 'loading'):
        matched = workspace.load(payload['matching'])
        saved = workspace.load(Blob(**matched['catalogue']))
        reference = read_blob(Blob(**payload['row']['reference']), workspace.scratch)
    with _measure(metrics, 'verification'):
        verdict = verify_catalogue(saved.catalogue, reference, reverse=payload['row']['reversed'])
    return dict(verification=verdict, metrics=metrics)


def _case_verdict(directions):
    if 'recovered' in directions:
        return 'recovered'
    if all(value == 'not_recovered' for value in directions):
        return 'not_recovered'
    return 'unknown' if all(value in ('not_recovered', 'unknown') for value in directions) else 'incomplete'


def _comparison_report(plan, ledger, store, baseline):
    """Read small stage descriptors and attempt audits, never molecular catalogues."""
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint
    records, totals, intervals = {}, {}, []
    for name, stage, state, raw_reference, raw_audit, error in ledger.db.execute('SELECT * FROM results'):
        row = dict(state=state)
        if raw_reference:
            value = store.load(Blob(**json.loads(raw_reference)))
            row.update(value)
            if stage == 'verification' and 'evidence' in row['verification']:
                evidence = row['verification'].pop('evidence')
                row['verification']['evidence_sha256'] = _digest(evidence)
        if error:
            row['error'] = json.loads(error)
        arm = name.split('/')[0]
        aggregate = totals.setdefault(arm, {}).setdefault(stage, dict(states=Counter(),
            measured_phases={}, total_worker_cpu_seconds=0., failed_sampled_cpu_seconds=0.,
            maximum_sampled_rss_bytes=0, attempts=0))
        aggregate['states'][state] += 1
        for phase, values in row.get('metrics', {}).items():
            combined = aggregate['measured_phases'].setdefault(phase, dict(cpu_seconds=0., wall_seconds=0.))
            for key in combined:
                combined[key] += values[key]
        if raw_audit:
            audit = store.load(Blob(**json.loads(raw_audit)))
            row['spent_seconds'] = audit['case'][3]
            row['attempts'] = len(audit['attempts'])
            row['worker_metrics'] = []
            for token, task, worker, start, end, outcome, raw_metrics in audit['attempts']:
                measured = _loads_checkpoint(raw_metrics) if raw_metrics is not None else {}
                values = {key: value for key, value in measured.items()
                          if key.endswith(('_seconds', '_bytes', '_samples'))}
                row['worker_metrics'].append(dict(state=outcome, **values))
                aggregate['attempts'] += 1
                aggregate['total_worker_cpu_seconds'] += measured.get('total_cpu_seconds', measured.get('cpu_seconds', 0.))
                if 'total_cpu_seconds' not in measured and 'cpu_seconds' not in measured:
                    aggregate['failed_sampled_cpu_seconds'] += measured.get('sampled_cpu_seconds', 0.)
                aggregate['maximum_sampled_rss_bytes'] = max(aggregate['maximum_sampled_rss_bytes'], measured.get('peak_tree_rss_bytes', 0))
                if start is not None and end is not None:
                    intervals.extend(((start, 1), (end, -1)))
        records.setdefault(name, {})[stage] = row
    cases = []
    for case in plan['cases']:
        value = dict(case=case, baseline=baseline[case]['status'], arms={})
        for arm in COMPARISON_ARMS:
            directions = {}
            for direction in ('forward', 'reverse'):
                stages = records.get(f'{arm}/{case:05d}/{direction}', {})
                verdict = stages.get('verification', {}).get('verification', {}).get('status', 'incomplete')
                directions[direction] = dict(verdict=verdict, stages=stages)
            value['arms'][arm] = dict(status=_case_verdict([d['verdict'] for d in directions.values()]),
                                      directions=directions)
        cases.append(value)
    active = overlap = 0
    for _, change in sorted(intervals):
        active += change
        overlap = max(overlap, active)
    comparisons = dict(shared_recovered=[], binary_recovered_weighted_not_recovered=[],
                       weighted_recovered_binary_not_recovered=[], both_not_recovered=[], inconclusive=[])
    for row in cases:
        binary, weighted = (row['arms'][arm]['status'] for arm in COMPARISON_ARMS)
        category = ({('recovered', 'recovered'): 'shared_recovered',
                     ('recovered', 'not_recovered'): 'binary_recovered_weighted_not_recovered',
                     ('not_recovered', 'recovered'): 'weighted_recovered_binary_not_recovered',
                     ('not_recovered', 'not_recovered'): 'both_not_recovered'}.get((binary, weighted), 'inconclusive'))
        comparisons[category].append(row['case'])
    baseline_comparison = {arm: dict(
        gained=[r['case'] for r in cases if r['baseline'] == 'not_recovered' and r['arms'][arm]['status'] == 'recovered'],
        lost=[r['case'] for r in cases if r['baseline'] == 'recovered' and r['arms'][arm]['status'] == 'not_recovered'],
        unresolved=[r['case'] for r in cases if r['arms'][arm]['status'] in ('unknown', 'incomplete')])
        for arm in COMPARISON_ARMS}
    expected = len(plan['cases'])*2*len(COMPARISON_ARMS)*3
    terminal = ledger.db.execute('SELECT count(*) FROM results').fetchone()[0]
    return dict(schema='mappa.golden_comparison_results/v1', plan=plan,
        terminal_stage_rows=terminal, terminal_complete=terminal == expected,
        expected_stage_rows=expected,
        complete=len(records) == len(plan['cases'])*4 and all(
            len(stages) == 3 and all(row['state'] == 'complete' for row in stages.values())
            for stages in records.values()),
        case_outcomes={arm: dict(Counter(row['arms'][arm]['status'] for row in cases)) for arm in COMPARISON_ARMS},
        stage_metrics=totals, peak_overlapping_attempts=overlap,
        attempt_span_seconds=(max(point[0] for point in intervals)-min(point[0] for point in intervals)
                              if intervals else 0.),
        paired_comparison=comparisons, baseline_comparison=baseline_comparison, cases=cases,
        measurement_scope='Stage measurements include explicit large-artifact publication; supervised worker CPU also includes this publication. Failed sampled CPU is separate from measured worker CPU. Summed wall times are not parallel elapsed time. Process-lifetime RSS can include earlier tasks.',
        baseline_scope='Archived 1830/1851 Golden result: one seed, cap100, original WBO tolerance1. The two new arms use cap2000; differences are not attributed to tolerance alone.')


def compare(index_path, output, *, workers=4, window=None, memory_bytes=8*2**30,
            stage_seconds=1800., run_seconds=86400., cases=None, resume=False):
    """Two declared matching arms, one persistent pool and a global pass barrier."""
    import fcntl
    import gzip
    import tempfile
    import uuid
    from enzymemap import Ledger, require_drained
    from mappa.execution import BundleStore, WorkPool
    from mappa.execution.identity import code_identity
    from mappa.mappa_core.core.scheduler import Task

    window = workers*2 if window is None else window
    if min(workers, window, memory_bytes, stage_seconds, run_seconds) <= 0:
        raise ValueError('Workers, admission window, memory and time budgets must be positive')
    index_path, output = Path(index_path), Path(output)
    raw = index_path.read_bytes(); index = json.loads(raw)
    if index.get('verification') != 'golden:verify' or index.get('reference_policy') != REFERENCE_POLICY:
        raise ValueError('Use the immutable Golden input index and its exact reference policy')
    selected = list(range(1851)) if cases is None else list(cases)
    if not selected or len(set(selected)) != len(selected) or any(case not in range(1851) for case in selected):
        raise ValueError('Select distinct Golden cases in 0..1850')
    rows = [row for row in index['directions'] if row['case'] in selected]
    wanted = {(case, reverse) for case in selected for reverse in (False, True)}
    if len(rows) != len(wanted) or {(r['case'], r['reversed']) for r in rows} != wanted:
        raise ValueError('Each selected Golden case requires exactly one input for each direction')
    if cases is None and len(index['directions']) != 3702:
        raise ValueError('The full Golden index must contain exactly 3702 directions')
    rows.sort(key=lambda row: (row['case'], row['reversed']))
    baseline_path = Path(__file__).resolve().parents[1]/'reports/golden_retest_20261101/reference-results.json.gz'
    baseline = {row['case']: row for row in json.loads(gzip.decompress(baseline_path.read_bytes()))}
    plan = dict(schema='mappa.golden_comparison/v1', cases=sorted(selected),
        full_dataset=cases is None, directions_per_arm=len(rows),
        arms={arm: asdict(comparison_config(arm)) for arm in COMPARISON_ARMS},
        binary_rule='One iff original WBO >= 0.2, otherwise zero; diagonal zero; elements, coordinates and atom order unchanged',
        event_weights='Original endpoint WBOs in both arms', reference_policy=REFERENCE_POLICY,
        input_index_sha256=hashlib.sha256(raw).hexdigest(), inputs=rows,
        code_identity=code_identity(), runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        baseline_sha256=hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
        memory_guard_bytes=memory_bytes, stage_seconds=stage_seconds,
        cache_budget_bytes=64*2**20, automatic_retries=False)
    plan = json.loads(json.dumps(plan))
    output.mkdir(parents=True, exist_ok=True)
    with (output/'.owner.lock').open('a+b') as owner:
        try:
            fcntl.flock(owner, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Golden comparison output already has an active owner') from None
        plan_path = output/'plan.json'
        if resume:
            if json.loads(plan_path.read_text()) != plan:
                raise ValueError('Resume input, scientific source, runner or stage budget differs')
            if not json.loads((output/'progress.json').read_text()).get('drained'):
                raise ValueError('Resume requires a completed drain; interrupted ownership needs an explicit audit')
            require_drained(output)
        elif plan_path.exists():
            raise ValueError('Existing comparison requires --resume')
        else:
            _publish(plan_path, plan)
        _publish(output/'progress.json', dict(drained=False, stage='starting',
                                             code_identity=plan['code_identity']))
        started, stop = time.monotonic(), False
        def drain(signum, frame):
            nonlocal stop
            stop = True
        handlers = {sig: signal.signal(sig, drain) for sig in (signal.SIGINT, signal.SIGTERM)}
        try:
            with tempfile.TemporaryDirectory(prefix='mappa-golden-comparison-') as scratch:
                ledger = Ledger(output, scratch, resume=resume)
                store = BundleStore(output/'bundles', scratch=scratch)
                try:
                    with WorkPool.local(workers=workers, store=store, memory_bytes=memory_bytes,
                                        case_seconds=stage_seconds, shared_cache=False) as pool:
                        status = dict(schema='mappa.golden_comparison_progress/v1',
                            code_identity=plan['code_identity'], cases=len(selected),
                            directions_per_arm=len(rows), expected_directions=len(rows)*2,
                            ready_workers=len(pool.queue.ready_workers()),
                            startup_seconds=pool.startup_elapsed_seconds, window=window)
                        def checkpoint(stage, active):
                            ledger.checkpoint(); pool.checkpoint()
                            status.update(stage=stage, active_directions=len(active), completed=ledger.counts(),
                                queue=pool.queue.snapshot(), elapsed_seconds=time.monotonic()-started)
                            _publish(output/'progress.json', status)
                            print(json.dumps(status), flush=True)
                        for stage in ('matching', 'decoding', 'verification'):
                            pending = iter((arm, row) for row in rows for arm in COMPARISON_ARMS)
                            exhausted, active = False, {}
                            checkpoint(stage, active)
                            next_checkpoint = time.monotonic()+30
                            while active or not exhausted:
                                stop |= time.monotonic()-started >= run_seconds
                                while not stop and not exhausted and len(active) < window:
                                    try:
                                        arm, row = next(pending)
                                    except StopIteration:
                                        exhausted = True
                                        break
                                    name = f"{arm}/{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}"
                                    if ledger.get(name, stage) is not None:
                                        continue
                                    payload = dict(name=name, arm=arm, row=row, runner_sha256=plan['runner_sha256'])
                                    if stage != 'matching':
                                        prior = ledger.get(name, 'matching')
                                        if prior is None:
                                            raise RuntimeError('Matching pass is incomplete')
                                        if prior[0] != 'complete':
                                            ledger.put(name, stage, 'blocked', error='Matching did not complete')
                                            continue
                                        payload['matching'] = Blob(**json.loads(prior[1]))
                                    key = f'{uuid.uuid4().hex}/{name}/{stage}'
                                    operation = dict(matching='compare_match', decoding='compare_decode', verification='compare_verify')[stage]
                                    pool.queue.add([Task(key, f'golden:{operation}', payload, case=key,
                                                         case_seconds=stage_seconds, memory_bytes=memory_bytes)])
                                    active[name] = key
                                finished = pool.queue.results(list(active.values()))
                                for name, key in list(active.items()):
                                    result = finished[key]
                                    if result['state'] in ('pending', 'running'):
                                        continue
                                    succeeded = result['state'] == 'complete'
                                    measured = result['result'] or {}
                                    audit = pool.queue.archive_case(key, store)
                                    ledger.put(name, stage, 'complete' if succeeded else 'error',
                                        measured['blob'] if succeeded else None, audit,
                                        None if succeeded else dict(state=result['state'], error=measured.get('error')))
                                    del active[name]
                                if time.monotonic() >= next_checkpoint:
                                    checkpoint(stage, active); next_checkpoint = time.monotonic()+30
                                if stop and not active:
                                    status['drained'] = True
                                    checkpoint(stage, active)
                                    break
                                if active:
                                    pool.check_workers(); time.sleep(.1)
                            checkpoint(stage, active)
                            if stop:
                                break
                        status['drained'] = True
                        checkpoint(stage, {})
                    report = _comparison_report(plan, ledger, store, baseline)
                    report.update(invocation_elapsed_seconds=time.monotonic()-started, stopped=stop,
                                  workers=workers, window=window)
                    _publish(output/'summary.json', report)
                    status.update(terminal_complete=report['terminal_complete'], complete=report['complete'])
                    _publish(output/'progress.json', status)
                    return report
                finally:
                    ledger.db.close()
        finally:
            for sig, handler in handlers.items():
                signal.signal(sig, handler)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ['audit-decoding']:
        parser = argparse.ArgumentParser(description='Audit saved decoded reference membership with replayable witnesses')
        parser.add_argument('campaign'); parser.add_argument('output')
        parser.add_argument('--cases', type=int, nargs='+')
        parser.add_argument('--workers', type=int, default=1)
        parser.add_argument('--seconds', type=float, default=10.)
        parser.add_argument('--task-seconds', type=float, default=60.)
        parser.add_argument('--memory-gib', type=float, default=2.)
        parser.add_argument('--unknown-from', help='Explicit additional audit of only unresolved directions in this prior report')
        args = parser.parse_args(argv[1:])
        report = audit_decoding(args.campaign, args.output, cases=args.cases, workers=args.workers,
            seconds=args.seconds, task_seconds=args.task_seconds, memory_bytes=int(args.memory_gib*2**30),
            unknown_from=args.unknown_from)
        return 1 if report['certificate_failures'] else 0
    if argv[:1] == ['compare']:
        parser = argparse.ArgumentParser(description='Paired full Golden study with original-WBO decoding')
        parser.add_argument('--index', required=True)
        parser.add_argument('--output', required=True)
        parser.add_argument('--cases', type=int, nargs='+', help='Bounded pilot; omit for all 1851 cases')
        parser.add_argument('--workers', type=int, default=4)
        parser.add_argument('--window', type=int)
        parser.add_argument('--memory-gib', type=float, default=8)
        parser.add_argument('--stage-seconds', type=float, default=1800.)
        parser.add_argument('--run-seconds', type=float, required=True)
        parser.add_argument('--resume', action='store_true')
        args = parser.parse_args(argv[1:])
        report = compare(args.index, args.output, workers=args.workers, window=args.window,
            memory_bytes=int(args.memory_gib*2**30), stage_seconds=args.stage_seconds,
            run_seconds=args.run_seconds, cases=args.cases, resume=args.resume)
        return 0 if report['complete'] else 75
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
