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
            assert record['schema'] == ('mappa.decoded_events/v6' if arm == 'binary_iso1'
                                         else 'mappa.decoded_events/v5')
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
