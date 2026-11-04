"""Golden keeps chemical endpoint symmetry and explicit unmatched-atom semantics."""
from dataclasses import asdict
import hashlib
from itertools import permutations
import json
from pathlib import Path
import sqlite3
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'bench'))
from golden import REFERENCE_POLICY, endpoint_generators
from enzymemap import run
from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob
from mappa.execution import BundleStore
from test_enzymemap_campaign import prepared_index


def test_reference_view_preserves_selected_native_atoms_and_unmatched_boundaries(tmp_path, monkeypatch):
    import mappa
    import golden
    from mappa.presentation.viewers import comparison_document

    def endpoint(elements, coordinates, bonds):
        weights = np.zeros((len(elements), len(elements)))
        for a, b, value in bonds:
            weights[a, b] = weights[b, a] = value
        return dict(elements=elements, coordinates=coordinates, wbo=weights.tolist())

    reactant = endpoint(['H', 'C', 'O', 'H', 'N', 'F'],
        [[-.2, 0., .3], [0., 0., 0.], [1.1, .2, -.1],
         [2.2, .4, 0.], [.1, 1.2, .9], [1.5, 1.1, -.4]],
        [(0, 1, 1.), (1, 2, 1.), (2, 3, 1.), (1, 4, 1.), (2, 5, 1.)])
    product = endpoint(['O', 'H', 'C', 'H', 'Cl', 'N'],
        [[5., 1.3, .5], [2., -.4, .1], [3., .2, .4],
         [6., 1., -.2], [4., 1.2, 1.5], [2.3, 1.1, -.6]],
        [(1, 2, 1.), (0, 3, 1.), (2, 4, 1.), (2, 5, 2.), (0, 4, 1.)])
    source = tmp_path/'prepared'
    store = BundleStore(source)
    input_blob = store.put(dict(reactant=reactant, product=product))
    reference_blob = store.put(dict(mapping=[(1, 2), (4, 5)],
        features=[dict(heavy=[1, 2, 4, 5]), dict(heavy=[0, 2, 4, 5])]))
    unreadable = dict(path=str(source/'not-selected.bundle'), offset=0, size=1, sha256='a'*64)
    index = source/'index.json'
    index.write_text(json.dumps(dict(directions=[
        dict(case=99, reversed=False, input=unreadable, reference=unreadable),
        dict(case=17, reversed=True, input=unreadable, reference=unreadable),
        dict(case=17, reversed=False, input=asdict(input_blob), reference=asdict(reference_blob))])))
    original = {p: p.read_bytes() for p in source.rglob('*') if p.is_file()}

    def forbidden(*args, **kwargs):
        raise AssertionError('Viewing a saved reference must not run mapping, decoding or verification')

    monkeypatch.setattr(mappa, 'search_aam', forbidden)
    monkeypatch.setattr(mappa, 'decode_events', forbidden)
    monkeypatch.setattr(golden, 'query_path', forbidden)
    cases = golden.reference_view_cases(index, [17], geometry_cache=tmp_path/'geometry')
    assert [case['index'] for case in cases] == [17]
    case = cases[0]
    for shown, native, heavy in zip(case['endpoints'], (reactant, product),
                                    ([1, 2, 4, 5], [0, 2, 4, 5])):
        assert shown['atom_indices'] == heavy
        assert shown['elements'] == [native['elements'][i] for i in heavy]
        assert shown['coordinates'] == [native['coordinates'][i] for i in heavy]
        assert shown['wbo'] == np.asarray(native['wbo'])[np.ix_(heavy, heavy)].tolist()
    record, = case['records']
    assert dict(record['mapping']) == {0: 1, 2: 3}
    assert {e['kind']: (e['r'], e['p'], e['wbo']) for e in record['events']} == {
        'broken': ([0, 1], [1], [1., None]),
        'formed': ([0], [1, 2], [None, 1.]),
        'strengthened': ([0, 2], [1, 3], [1., 2.]),
    }
    assert len(record['events']) == 3  # Wholly unmatched bonds create no correspondence.
    assert record['provenance']['input_sha256'] == input_blob.sha256
    assert record['provenance']['reference_sha256'] == reference_blob.sha256
    assert record['provenance']['geometry'] == [dict(source='native', display_only=False)]*2
    shown = comparison_document(case)
    assert shown['reactant']['atom_indices'] == [1, 2, 4, 5]
    assert shown['product']['atom_indices'] == [0, 2, 4, 5]
    assert shown['mechanisms'][0]['product_xyz_in_R_aligned'] is None
    assert shown['mechanisms'][0]['broken_bonds_R'] == [[0, 1]]
    assert shown['mechanisms'][0]['formed_bonds_P'] == [[1, 2]]
    from tools.render_mapping_comparison import render

    case['note'] = 'Saved reference note; no inferred search settings.'
    case['recovery'] = [dict(seeds=10, status='represented')]
    comparison = tmp_path/'comparison.json'
    comparison.write_text(json.dumps(cases))
    html = tmp_path/'view.html'
    assert render(comparison, html) == 1
    assert case['note'] in html.read_text()
    assert {p: p.read_bytes() for p in source.rglob('*') if p.is_file()} == original


def decoding_campaign(tmp_path):
    """Archived stage metadata; catalogue files deliberately do not exist."""
    from prepare_batches import prepare_fresh
    from mappa import AAMSearchConfig

    config = asdict(AAMSearchConfig(seed_count=1, branch_limit=2000))
    rows = [dict(case=case, reversed=reverse,
                 input=dict(path=str(tmp_path/'unopened-input.bundle'), offset=case*20,
                            size=10, sha256='a'*64),
                 reference=dict(path=str(tmp_path/'unopened-reference.bundle'), offset=0,
                                size=10, sha256='b'*64))
            for case in range(3) for reverse in (False, True)]
    index = tmp_path/'original-index.json'
    index.write_text(json.dumps(dict(config=config, directions=rows)))
    source = tmp_path/'source'
    prepare_fresh([index], source, workers=[1, 1])
    (source/'launch.json').write_text(json.dumps(dict(source_commit='a'*40, source_identity='c'*64)))
    for batch in range(2):
        folder = source/str(batch)
        data = json.loads((folder/'index.json').read_text())
        output = folder/'result'; output.mkdir()
        store = BundleStore(output/'bundles')
        (output/'bundles/cache').mkdir()
        (output/'bundles/cache/.owner.lock').touch()
        with sqlite3.connect(output/'bundles/queue.sqlite') as queue:
            queue.execute('CREATE TABLE tasks (key TEXT)')
            queue.execute('CREATE TABLE attempts (key TEXT)')
        plan = dict(input_index=str(folder/'index.json'),
            input_index_sha256=hashlib.sha256((folder/'index.json').read_bytes()).hexdigest(),
            config=config, cases=len({r['case'] for r in data['directions']}),
            directions=len(data['directions']))
        (output/'plan.json').write_text(json.dumps(plan))
        (output/'progress.json').write_text(json.dumps(dict(source='c'*64)))
        with sqlite3.connect(output/'ledger.sqlite') as ledger:
            ledger.execute('CREATE TABLE results (case_name TEXT,stage TEXT,state TEXT,reference TEXT,audit TEXT,error TEXT,PRIMARY KEY(case_name,stage))')
            ledger.execute('CREATE TABLE retries (utc REAL,reason TEXT,case_name TEXT,stage TEXT,state TEXT,reference TEXT,audit TEXT,error TEXT)')
            for row in data['directions']:
                name = f"{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}"
                catalogue = dict(path=str(tmp_path/'unopened-catalogue.bundle'),
                    offset=100*row['case']+int(row['reversed']), size=100+row['case'], sha256='d'*64)
                for stage in ('matching', 'deduplication'):
                    ledger.execute('INSERT INTO results VALUES (?,?,?,?,?,?)',
                                   (name, stage, 'complete', json.dumps(catalogue), None, None))
                state = 'memory_limit' if name == '00000/forward' else 'complete'
                attempts = [(None, None, None, None, None, state)]
                if name == '00002/reverse':
                    attempts.insert(0, (None, None, None, None, None, 'memory_limit'))
                audit = store.put(dict(tasks=[(name, None, None, None, state)], attempts=attempts))
                ledger.execute('INSERT INTO results VALUES (?,?,?,?,?,?)',
                    (name, 'decoding', 'error' if state == 'memory_limit' else 'complete',
                     None if state == 'memory_limit' else json.dumps(catalogue),
                     json.dumps(asdict(audit)), json.dumps({name: state}) if state == 'memory_limit' else None))
                if name == '00001/reverse':
                    previous = store.put(dict(tasks=[(name, None, None, None, 'memory_limit')], attempts=[]))
                    ledger.execute('INSERT INTO retries VALUES (?,?,?,?,?,?,?,?)',
                        (1., 'resource retry', name, 'decoding', 'error', None,
                         json.dumps(asdict(previous)), json.dumps({name: 'memory_limit'})))
    return source


def test_fresh_decoding_batches_preserve_full_scope_and_prior_memory_failures(tmp_path):
    from prepare_batches import prepare_decoding
    from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig

    source = decoding_campaign(tmp_path)
    before = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in source.rglob('*') if p.is_file()}
    output = tmp_path/'decoded'
    manifest = prepare_decoding(source, output, [2, 1], memory_failures_workers=[1],
                                decode_config=EventDecodeConfig(.7, .4))
    assert manifest['cases'] == 3 and manifest['directions'] == 6
    assert manifest['imported_stage_results'] == 0 and manifest['reused_catalogues'] == 6
    assert manifest['prior_memory_failure_directions'] == ['00000/forward', '00001/reverse', '00002/reverse']
    assert manifest['decode_config'] == dict(threshold=.7, metal_threshold=.4)
    assert manifest['decoding_source']['source_commit'] == 'a'*40
    assert manifest['decoding_source']['source_identity'] == 'c'*64
    batches = [json.loads((output/str(i)/'index.json').read_text()) for i in range(3)]
    rows = [row for batch in batches for row in batch['directions']]
    assert {(r['case'], r['reversed']) for r in rows} == {(i, r) for i in range(3) for r in (False, True)}
    assert len(rows) == 6
    assert {(r['case'], r['reversed']) for r in batches[2]['directions']} == {(0, False), (1, True), (2, True)}
    assert [r['resource_group'] for r in manifest['batches']] == ['ordinary', 'ordinary', 'memory_failures']
    original = { (r['case'],r['reversed']):r for r in json.loads((tmp_path/'original-index.json').read_text())['directions'] }
    for row in rows:
        assert {key: value for key, value in row.items() if key != 'catalogue'} == original[row['case'],row['reversed']]
        assert row['catalogue']['sha256'] == 'd'*64
    for batch in manifest['batches']:
        raw = (output/str(batch['batch'])/'index.json').read_bytes()
        assert hashlib.sha256(raw).hexdigest() == batch['index_sha256']
    assert not list(output.rglob('*.sqlite')) and not list(output.rglob('*.bundle'))
    assert {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in source.rglob('*') if p.is_file()} == before
    with pytest.raises(FileExistsError):
        prepare_decoding(source, output, [1])
    with pytest.raises(ValueError, match='separate output'):
        prepare_decoding(source, source/'nested', [1])


@pytest.mark.parametrize('corruption', ['matching', 'deduplication', 'source_index', 'batch_index',
                                        'config', 'identity', 'missing_batch', 'retry_scope'])
def test_fresh_decoding_rejects_changed_or_incomplete_source_before_output(tmp_path, corruption):
    from prepare_batches import prepare_decoding

    source = decoding_campaign(tmp_path)
    result = source/'0/result'
    if corruption in ('matching', 'deduplication'):
        with sqlite3.connect(result/'ledger.sqlite') as db:
            db.execute("UPDATE results SET state='error' WHERE stage=?", (corruption,))
    elif corruption == 'source_index':
        with (tmp_path/'original-index.json').open('a') as stream:
            stream.write(' ')
    elif corruption == 'batch_index':
        with (source/'0/index.json').open('a') as stream:
            stream.write(' ')
    elif corruption == 'config':
        path = result/'plan.json'; value = json.loads(path.read_text())
        value['config']['seed_count'] = 10; path.write_text(json.dumps(value))
    elif corruption == 'identity':
        (result/'progress.json').write_text(json.dumps(dict(source='wrong')))
    elif corruption == 'missing_batch':
        path = source/'manifest.json'; value = json.loads(path.read_text())
        value['batches'].pop(); path.write_text(json.dumps(value))
    elif corruption == 'retry_scope':
        with sqlite3.connect(result/'ledger.sqlite') as db:
            db.execute('INSERT INTO retries VALUES (?,?,?,?,?,?,?,?)',
                       (1., 'invalid', '99999/forward', 'decoding', 'error', None, None, None))
    output = tmp_path/'decoded'
    with pytest.raises(ValueError):
        prepare_decoding(source, output, [1])
    assert not output.exists()


@pytest.mark.parametrize('ownership', ['lock', 'tasks', 'attempts'])
def test_fresh_decoding_rejects_active_or_unarchived_source(tmp_path, ownership):
    import fcntl
    from prepare_batches import prepare_decoding

    source = decoding_campaign(tmp_path)
    result = source/'0/result/bundles'
    with (result/'cache/.owner.lock').open('rb') as owner:
        if ownership == 'lock':
            fcntl.flock(owner, fcntl.LOCK_EX | fcntl.LOCK_NB)
        else:
            with sqlite3.connect(result/'queue.sqlite') as queue:
                queue.execute(f'INSERT INTO {ownership} VALUES (?)', ('unarchived',))
        with pytest.raises((RuntimeError, ValueError), match='active coordinator|drained saved queue'):
            prepare_decoding(source, tmp_path/'decoded', [1])
    assert not (tmp_path/'decoded').exists()


def comparison_index(tmp_path):
    from mappa.mappa_core.utils.io.records import endpoint_record
    from test_fragment_conversion import endpoint

    store = BundleStore(tmp_path/'prepared')
    left = endpoint('CO', [(0, 1, 3.)])
    right = endpoint('CO', [(0, 1, 1.)])
    reference = store.put(dict(mapping=[(0, 0), (1, 1)],
        features=[dict(heavy=[0, 1])]*2, generators=[(), ()]))
    rows = []
    for reverse in (False, True):
        source, target = (right, left) if reverse else (left, right)
        original = store.put(dict(name='triple-to-single',
            reactant=endpoint_record(source), product=endpoint_record(target)))
        row = dict(case=0, reversed=reverse, input=asdict(original), reference=asdict(reference))
        rows.append(row)
    index = tmp_path/'index.json'
    index.write_text(json.dumps(dict(directions=rows, verification='golden:verify',
                                    reference_policy=REFERENCE_POLICY)))
    return index


def test_binary_comparison_uses_the_matching_floor_without_mutating_original_weights():
    from golden import comparison_problem
    from mappa import AAMProblem
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    from test_fragment_conversion import endpoint

    below, above = np.nextafter(.2, 0.), np.nextafter(.2, 1.)
    source = endpoint('CONH', [(0, 1, below), (0, 2, .2), (1, 2, above), (2, 3, 3.)])
    target = endpoint('CONH', [(0, 1, 3.), (0, 2, above), (1, 2, .2), (2, 3, below)])
    original = AAMProblem(source, target)
    original_hash = FinalBranchCatalogue(original).problem_sha256
    before = [ep.wbo.copy() for ep in (source, target)]
    binary = comparison_problem(original, 'binary_iso1')
    weighted = comparison_problem(original, 'weighted_iso10')
    expected_source = np.asarray([[0, 0, 1, 0], [0, 0, 1, 0], [1, 1, 0, 1], [0, 0, 1, 0]])
    expected_target = np.asarray([[0, 1, 1, 0], [1, 0, 1, 0], [1, 1, 0, 0], [0, 0, 0, 0]])
    assert np.array_equal(binary.reactant.wbo, expected_source)
    assert np.array_equal(binary.product.wbo, expected_target)
    for raw, transformed, unchanged, weights in zip(
            (source, target), (binary.reactant, binary.product),
            (weighted.reactant, weighted.product), before):
        assert transformed.elements == raw.elements
        assert np.array_equal(transformed.coordinates, raw.coordinates)
        assert np.array_equal(raw.wbo, weights)
        assert np.array_equal(unchanged.wbo, weights)
    assert FinalBranchCatalogue(original).problem_sha256 == original_hash
    assert FinalBranchCatalogue(weighted).problem_sha256 == original_hash
    assert FinalBranchCatalogue(binary).problem_sha256 != original_hash


def test_two_arm_pipeline_preserves_event_weights_barrier_and_terminal_failure_budgets(tmp_path):
    from golden import COMPARISON_ARMS, compare
    from mappa.mappa_core.post_processing.decoder.events import EVENT_CATEGORIES
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint
    from mappa.mappa_core.utils.io.records import problem_from_record

    index = comparison_index(tmp_path)
    inputs = json.loads(index.read_text())['directions']
    input_files = {Path(row[field]['path']) for row in inputs for field in ('input', 'reference')}
    original_hashes = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in input_files}
    output = tmp_path/'comparison'
    options = dict(workers=2, window=2, memory_bytes=2*2**30, stage_seconds=60., cases=[0])
    report = compare(index, output, **options)
    assert report['complete']
    assert report['terminal_stage_rows'] == report['expected_stage_rows'] == 12
    assert report['case_outcomes'] == {arm: {'recovered': 1} for arm in COMPARISON_ARMS}
    store = BundleStore(output/'bundles')
    with sqlite3.connect(output/'ledger.sqlite') as db:
        rows = {(name, stage): (state, reference, audit) for name, stage, state, reference, audit
                in db.execute('SELECT case_name,stage,state,reference,audit FROM results')}
    attempts = {stage: [] for stage in ('matching', 'decoding', 'verification')}
    for (name, stage), (state, _, raw_audit) in rows.items():
        assert state == 'complete'
        audit = store.load(Blob(**json.loads(raw_audit)))
        assert len(audit['attempts']) == 1
        attempt, = audit['attempts']
        assert attempt[3] <= attempt[4] and attempt[5] == 'complete'
        attempts[stage].append(attempt)
        task, = audit['tasks']
        specification = _loads_checkpoint(task[1])
        assert specification.memory_bytes == options['memory_bytes']
        assert specification.case_seconds == options['stage_seconds']
    assert max(row[4] for row in attempts['matching']) <= min(row[3] for row in attempts['decoding'])
    assert max(row[4] for row in attempts['decoding']) <= min(row[3] for row in attempts['verification'])

    for row in inputs:
        reverse = row['reversed']
        original = problem_from_record(read_blob(Blob(**row['input']), None))
        event_hash = FinalBranchCatalogue(original).problem_sha256
        direction = 'reverse' if reverse else 'forward'
        category = 'strengthening' if reverse else 'weakening'
        for arm in COMPARISON_ARMS:
            name = f'{arm}/00000/{direction}'
            def stage_value(stage):
                return store.load(Blob(**json.loads(rows[name, stage][1])))
            matched, decoded, verified = (stage_value(stage) for stage in attempts)
            saved = store.load(Blob(**matched['catalogue']))
            assert matched['original_input'] == row['input']
            assert saved.search_config.seed_count == 1
            assert saved.search_config.branch_limit == 2000
            assert saved.search_config.iso_tolerance == (1. if arm == 'binary_iso1' else 10.)
            assert saved.scope and all(cut.seeds == (0,) for cut in saved.scope)
            assert verified['verification']['status'] == 'recovered'
            assert verified['verification']['policy'] == REFERENCE_POLICY
            record = store.load(Blob(**decoded['decoded'])).record
            assert decoded['event_problem_sha256'] == event_hash
            assert decoded['matching_problem_sha256'] == saved.catalogue.problem_sha256
            assert record['schema'] == ('mappa.decoded_events/v8' if arm == 'binary_iso1'
                                         else 'mappa.decoded_events/v7')
            event_record = record['event_problem'] if arm == 'binary_iso1' else record['problem']
            assert np.array_equal(event_record['reactant']['wbo'], original.reactant.wbo)
            assert np.array_equal(event_record['product']['wbo'], original.product.wbo)
            assert record['catalogue'] == saved.catalogue.to_record()
            family = verified['verification']['family']
            # With two distinct elements the recovered mapping is fixed. Read
            # its saved local response directly, without invoking another decoder.
            for number in record['families'][family]['tables']:
                table = record['tables'][number]
                terms = [term for region in (*table['internal'], *table['boundary'])
                         for term in region['terms'] if tuple(term['atoms']) == (0, 1)]
                assert len(terms) == 1
                response = record['responses'][terms[0]['response']]
                target_weight = float(original.product.wbo[0, 1])
                code = response['categories'][response['target_weights'].index(target_weight)]
                assert EVENT_CATEGORIES[code - 1] == category and code != 0
                assert not table['product_boundary_pairs']
    assert {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in input_files} == original_hashes

    # An empty checkpoint may predate admitted work. Only a deliberate drain
    # permits resume; refusal must leave the saved queue and ledger untouched.
    progress_path = output/'progress.json'
    progress = json.loads(progress_path.read_text())
    assert progress['drained']
    progress_path.write_text(json.dumps(dict(progress, drained=False)))
    before_refusal = {path: path.read_bytes() for path in
                      (output/'ledger.sqlite', output/'bundles/queue.sqlite', progress_path)}
    with pytest.raises(ValueError, match='drain'):
        compare(index, output, resume=True, **options)
    assert {path: path.read_bytes() for path in before_refusal} == before_refusal
    progress_path.write_text(json.dumps(progress))

    # Seed terminal failures in archived evidence instead of exhausting real
    # memory or time. Resume must retain their attempts and spent budgets.
    failures = [('binary_iso1/00000/forward', 'decoding', 'memory_limit', 17.),
                ('weighted_iso10/00000/reverse', 'verification', 'timeout', 60.)]
    with sqlite3.connect(output/'ledger.sqlite') as db:
        for name, stage, failure, spent in failures:
            audit = store.load(Blob(**json.loads(rows[name, stage][2])))
            audit['tasks'] = [tuple(failure if i == 4 else v for i, v in enumerate(task))
                              for task in audit['tasks']]
            audit['attempts'] = [tuple(failure if i == 5 else v for i, v in enumerate(attempt))
                                 for attempt in audit['attempts']]
            budget = list(audit['case']); budget[3] = spent
            audit['case'] = tuple(budget)
            db.execute("UPDATE results SET state='error',reference=NULL,audit=?,error=? WHERE case_name=? AND stage=?",
                       (json.dumps(asdict(store.put(audit))), json.dumps({'state': failure}), name, stage))
        db.execute("DELETE FROM results WHERE case_name='binary_iso1/00000/forward' AND stage='verification'")
        before = {(row[0], row[1]): row for row in db.execute('SELECT * FROM results')}
    resumed = compare(index, output, resume=True, **options)
    assert not resumed['complete'] and resumed['terminal_stage_rows'] == 12
    with sqlite3.connect(output/'ledger.sqlite') as db:
        after = {(row[0], row[1]): row for row in db.execute('SELECT * FROM results')}
        assert db.execute('SELECT count(*) FROM retries').fetchone()[0] == 0
    assert {key: after[key] for key in before} == before
    assert set(after)-set(before) == {('binary_iso1/00000/forward', 'verification')}
    verified = store.load(Blob(**json.loads(after['binary_iso1/00000/forward', 'verification'][3])))
    assert verified['verification']['status'] == 'recovered'
    for name, stage, failure, spent in failures:
        audit = store.load(Blob(**json.loads(after[name, stage][4])))
        assert audit['case'][3] == spent
        assert len(audit['attempts']) == 1 and audit['attempts'][0][5] == failure
    with sqlite3.connect(output/'bundles/queue.sqlite') as db:
        assert db.execute('SELECT count(*) FROM tasks').fetchone()[0] == 0
        assert db.execute('SELECT count(*) FROM attempts').fetchone()[0] == 0


def test_endpoint_group_matches_independent_colored_graph_permutations():
    heavy = [0, 2, 4, 5]
    feature = dict(heavy=heavy, colors=[['C'], ['O'], ['C'], ['O']],
        bonds=[[0, 1, [1]], [1, 2, [1]], [2, 3, [1]], [3, 0, [1]]])
    generators = endpoint_generators(feature, 6)
    identity = tuple(range(6))
    group, pending = {identity}, [identity]
    while pending:
        current = pending.pop()
        for generator in generators:
            moved = tuple(generator[a] for a in current)
            if moved not in group:
                group.add(moved); pending.append(moved)
    edges = {frozenset((a, b)) for a, b, _ in feature['bonds']}
    expected = set()
    for permutation in permutations(range(4)):
        if any(feature['colors'][a] != feature['colors'][permutation[a]] for a in range(4)):
            continue
        if {frozenset((permutation[a], permutation[b])) for a, b, _ in feature['bonds']} != edges:
            continue
        lifted = list(identity)
        for a, b in enumerate(permutation):
            lifted[heavy[a]] = heavy[b]
        expected.add(tuple(lifted))
    assert group == expected
    assert all(g[1] == 1 and g[3] == 3 for g in group)


@pytest.mark.parametrize('omit_atom', [False, True])
def test_direct_inputs_and_golden_unmatched_policy_across_complete_pipeline(tmp_path, omit_atom):
    index = prepared_index(tmp_path)
    data = json.loads(index.read_text())
    store = BundleStore(tmp_path/'prepared')
    mapping = [(i, i) for i in range(1 if omit_atom else 0, 4)]
    reference = store.put(dict(mapping=mapping, features=[dict(heavy=list(range(4)))]*2,
                              generators=[(), ()]))
    for row in data['directions']:
        preparation = read_blob(Blob(**row.pop('blob')), None)
        row['input'] = asdict(preparation[0].payload['input'])
        row['reference'] = asdict(reference)
    data.update(verification='golden:verify', reference_policy=REFERENCE_POLICY)
    index.write_text(json.dumps(data))
    output = tmp_path/'result'
    assert run(index, output, local=True, workers=2, memory_bytes=2*2**30)['finished']
    with sqlite3.connect(output/'ledger.sqlite') as db:
        rows = db.execute("SELECT reference FROM results WHERE stage='verification'").fetchall()
    for raw, in rows:
        result = store.load(Blob(**json.loads(raw)))
        assert result['status'] == ('not_recovered' if omit_atom else 'recovered')
        assert result['policy'] == REFERENCE_POLICY


@pytest.mark.parametrize('change,expected', [
    (None, 'equal'), ('record', 'different'), ('metadata', 'different'),
    ('missing_baseline', 'previously_missing'), ('missing_candidate', 'candidate_error')])
def test_decoding_audit_hashes_complete_records_and_releases_old_before_loading_new(change, expected):
    import weakref
    import performance
    from mappa import AAMSearchConfig
    from mappa.mappa_core.core.scheduler.results import CutScope, DecodingResult, MatchingResult
    from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig

    references = {name: Blob('unused.bundle', i, 1, str(i)*64)
                  for i, name in enumerate(('matching', 'baseline', 'candidate'))}
    scope = (CutScope(0, (), (0, 1)),)
    config = AAMSearchConfig(seed_count=2)
    record = dict(schema='audit-test', tables=[dict(boundary=[[0, 1, .5], [2, 3, .7]])])
    payload = dict(case='00000/forward', matching=asdict(references['matching']),
        search_config=json.loads(json.dumps(asdict(config))), decode_config=asdict(EventDecodeConfig()),
        runner_sha256=hashlib.sha256(Path(performance.__file__).read_bytes()).hexdigest())
    for arm in ('baseline', 'candidate'):
        missing = change == 'missing_'+arm
        payload[arm] = None if missing else asdict(references[arm])
        payload[arm+'_state'] = 'error' if missing else 'complete'

    class Workspace:
        old = None
        def load(self, reference):
            if reference == references['matching']:
                return MatchingResult(payload['case'], reference, scope, ())
            candidate = reference == references['candidate']
            if candidate and self.old is not None:
                assert self.old() is None, 'Old decoded artifact is still resident'
            value = json.loads(json.dumps(record))
            if candidate and change == 'record':
                value['tables'][0]['boundary'][1][2] = .8
            stops = (('capped', 1),) if candidate and change == 'metadata' else ()
            saved = DecodingResult(payload['case'], EventDecodeConfig(), scope, stops,
                                   value, config, 'reused_native')
            if not candidate:
                self.old = weakref.ref(saved)
            return saved

    result = performance.audit_decoding_direction(payload, Workspace())
    assert result['status'] == expected
    if result['baseline']:
        canonical = json.dumps(record, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
        assert result['baseline']['record_sha256'] == hashlib.sha256(canonical).hexdigest()
        assert result['baseline']['canonical_record_bytes'] == len(canonical)
    if change == 'metadata':
        assert result['baseline']['record_sha256'] == result['candidate']['record_sha256']
        assert result['baseline']['metadata_sha256'] != result['candidate']['metadata_sha256']


def test_decoding_audit_real_stage_outputs_resume_and_changed_source_rejection(tmp_path):
    from performance import audit_decoding
    from prepare_batches import prepare_decoding, prepare_fresh
    from mappa.execution.identity import code_identity

    index = prepared_index(tmp_path)
    baseline = tmp_path/'baseline'
    prepare_fresh([index], baseline, [1])
    (baseline/'launch.json').write_text(json.dumps(dict(source_commit='a'*40, source_identity=code_identity())))
    assert run(baseline/'0/index.json', baseline/'0/result', local=True, workers=1,
               memory_bytes=2**30, stop_after='decoding')['pass_finished'] == 'decoding'
    candidate = tmp_path/'candidate'
    prepare_decoding(baseline, candidate, [1])
    (candidate/'launch.json').write_text(json.dumps(dict(source_commit='b'*40, source_identity=code_identity())))
    assert run(candidate/'0/index.json', candidate/'0/result', local=True, workers=1,
               memory_bytes=2**30, decode_only=True)['pass_finished'] == 'decoding'
    frozen = {path: hashlib.sha256(path.read_bytes()).hexdigest()
              for source in (baseline, candidate) for path in source.rglob('*') if path.is_file()}
    output = tmp_path/'audit'
    report = audit_decoding(baseline, candidate, output, workers=1,
                            memory_bytes=2**30, task_seconds=60.)
    assert report['complete'] and report['terminal_complete']
    assert report['status_counts'] == {'equal': 2}
    assert report['sources']['candidate'] == dict(source_commit='b'*40, source_identity=code_identity())
    assert report['baseline_states'] == report['candidate_states'] == {'complete': 2}
    assert {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in frozen} == frozen
    with sqlite3.connect(output/'ledger.sqlite') as db:
        rows = db.execute('SELECT * FROM results').fetchall()
    bundles = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in (output/'bundles').glob('*.bundle')}
    assert audit_decoding(baseline, candidate, output, workers=1, memory_bytes=2**30,
                          task_seconds=60., resume=True)['complete']
    with sqlite3.connect(output/'ledger.sqlite') as db:
        assert db.execute('SELECT * FROM results').fetchall() == rows
    assert {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in bundles} == bundles
    with (candidate/'0/result/ledger.sqlite').open('ab') as stream:
        stream.write(b'changed-after-audit')
    with pytest.raises(ValueError, match='sources, implementation or budgets changed'):
        audit_decoding(baseline, candidate, output, memory_bytes=2**30, task_seconds=60., resume=True)


def test_decoding_audit_plans_missing_baseline_outputs_and_rejects_incomplete_candidate(tmp_path):
    from performance import decoding_audit_plan
    from prepare_batches import prepare_decoding

    source = decoding_campaign(tmp_path)
    target = tmp_path/'target'
    manifest = prepare_decoding(source, target, [1])
    (target/'launch.json').write_text(json.dumps(dict(source_commit='b'*40, source_identity='e'*64)))
    folder = target/'0'; index = json.loads((folder/'index.json').read_text())
    output = folder/'result'; output.mkdir(); (output/'bundles').mkdir()
    with sqlite3.connect(output/'bundles/queue.sqlite') as queue:
        queue.execute('CREATE TABLE tasks (key TEXT)'); queue.execute('CREATE TABLE attempts (key TEXT)')
    (output/'plan.json').write_text(json.dumps(dict(mode='decode_only', config=index['config'],
        decoder_source='e'*64, decode_config=manifest['decode_config'],
        decoding_source=manifest['decoding_source'],
        input_index_sha256=manifest['batches'][0]['index_sha256'])))
    (output/'progress.json').write_text(json.dumps(dict(source='e'*64)))
    with sqlite3.connect(output/'ledger.sqlite') as db:
        db.execute('CREATE TABLE results (case_name TEXT,stage TEXT,state TEXT,reference TEXT,audit TEXT,error TEXT)')
        for row in index['directions']:
            name = f"{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}"
            db.execute('INSERT INTO results VALUES (?,?,?,?,?,?)',
                       (name, 'decoding', 'complete', json.dumps(row['catalogue']), None, None))
    plan = decoding_audit_plan(source, target, memory_bytes=2**30, task_seconds=60.)
    assert plan['directions'] == 6
    assert plan['baseline_states'] == {'complete': 5, 'error': 1}
    assert plan['candidate_states'] == {'complete': 6}
    assert sum(row['baseline'] is None for row in plan['rows']) == 1
    original_index = (folder/'index.json').read_text()
    original_manifest = (target/'manifest.json').read_text()
    original_plan = (output/'plan.json').read_text()
    for field, value, message in (
            ('decoder_source', 'f'*64, 'decoder identities differ'),
            ('decode_config', dict(threshold=.7, metal_threshold=.4), 'decoding policies differ')):
        changed_plan = json.loads(original_plan); changed_plan[field] = value
        (output/'plan.json').write_text(json.dumps(changed_plan))
        with pytest.raises(ValueError, match=message):
            decoding_audit_plan(source, target, memory_bytes=2**30, task_seconds=60.)
    (output/'plan.json').write_text(original_plan)
    index['directions'][0]['catalogue']['offset'] += 1
    (folder/'index.json').write_text(json.dumps(index))
    digest = hashlib.sha256((folder/'index.json').read_bytes()).hexdigest()
    manifest['batches'][0]['index_sha256'] = digest
    (target/'manifest.json').write_text(json.dumps(manifest))
    changed_plan = json.loads(original_plan); changed_plan['input_index_sha256'] = digest
    (output/'plan.json').write_text(json.dumps(changed_plan))
    with pytest.raises(ValueError, match='catalogue or original input lineage differs'):
        decoding_audit_plan(source, target, memory_bytes=2**30, task_seconds=60.)
    (folder/'index.json').write_text(original_index)
    (target/'manifest.json').write_text(original_manifest)
    (output/'plan.json').write_text(original_plan)
    with sqlite3.connect(output/'ledger.sqlite') as db:
        db.execute("DELETE FROM results WHERE case_name='00000/forward'")
    with pytest.raises(ValueError, match='duplicate or missing decoding directions'):
        decoding_audit_plan(source, target, memory_bytes=2**30, task_seconds=60.)
