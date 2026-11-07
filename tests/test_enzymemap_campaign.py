"""Fresh benchmark orchestration uses public stages and bounded archived history."""
from dataclasses import asdict
import json
from pathlib import Path
import sqlite3
import sys

import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'bench'))
from enzymemap import run
from mappa import AAMSearchConfig
from mappa.mappa_core.core.scheduler import Task
from mappa.mappa_core.core.scheduler.storage import BundleStore
from mappa.mappa_core.utils.io.records import endpoint_record
from test_reused_aam_api import problem


def prepared_index(tmp_path, *, corrupt=False):
    store=BundleStore(tmp_path/'prepared')
    case=problem()
    reference=store.put(dict(mapping=[(i,i) for i in range(4)]))
    rows=[]
    for reverse in (False,True):
        left,right=(case.product,case.reactant) if reverse else (case.reactant,case.product)
        source=store.put(dict(name="test",reactant=endpoint_record(left),product=endpoint_record(right)))
        saved=store.put([Task('retired','retired.must_never_execute:run',dict(input=source))])
        row=dict(case=0,reversed=reverse,blob=asdict(saved),reference=asdict(reference))
        if corrupt: row['blob']['sha256']='0'*64
        rows.append(row)
    path=tmp_path/'index.json'
    path.write_text(json.dumps(dict(config=asdict(AAMSearchConfig(seed_count=1,branch_limit=10)),directions=rows)))
    return path


def test_full_passes_and_resume_preserve_results_without_old_operations(tmp_path):
    index=prepared_index(tmp_path)
    data = json.loads(index.read_text())
    data['decode_config'] = dict(threshold=.75, metal_threshold=None, event_slack=5, fragment_seconds=5.)
    index.write_text(json.dumps(data))
    output=tmp_path/'result'
    status=run(index,output,local=True,workers=2,window=2)
    assert status['finished']
    with sqlite3.connect(output/'ledger.sqlite') as db:
        assert db.execute("SELECT count(*) FROM results WHERE state='complete'").fetchone()[0] == 8
        record=db.execute("SELECT reference FROM results WHERE stage='verification'").fetchone()[0]
    from mappa.mappa_core.core.scheduler.bundles import Blob
    assert BundleStore(output/'bundles').load(Blob(**json.loads(record)))['status']=='recovered'
    assert json.loads((output/'plan.json').read_text())['decode_config'] == data['decode_config']
    with sqlite3.connect(output/'ledger.sqlite') as db:
        for raw, in db.execute("SELECT reference FROM results WHERE stage='decoding'"):
            result = BundleStore(output/'bundles').load(Blob(**json.loads(raw)))
            assert asdict(result.config) == data['decode_config']
            assert result.record['event_slack'] == 5
            assert result.record['event_window'] == result.record['minimum_events'] + 5
    with sqlite3.connect(output/'bundles/queue.sqlite') as db:
        assert db.execute('SELECT count(*) FROM tasks').fetchone()[0] == 0
    before=len(list((output/'bundles').glob('*.bundle')))
    assert run(index,output,local=True,workers=1,resume=True)['finished']
    assert len(list((output/'bundles').glob('*.bundle'))) == before
    data['decode_config']['event_slack'] = 6
    index.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='Resume input/configuration differs'):
        run(index, output, local=True, workers=1, resume=True)


def test_reduction_only_preserves_full_results_without_search_or_decoding(tmp_path, monkeypatch):
    from mappa.execution import WorkPool
    from mappa.mappa_core.core.scheduler.bundles import Blob
    index = prepared_index(tmp_path)
    source = tmp_path/'original'
    assert run(index, source, local=True, workers=1, stop_after='deduplication')['pass_finished'] == 'deduplication'
    data = json.loads(index.read_text())
    original_bytes = (source/'ledger.sqlite').read_bytes()
    with sqlite3.connect(source/'ledger.sqlite') as db:
        refs = {(name, stage): json.loads(reference) for name, stage, reference in
                db.execute('SELECT case_name,stage,reference FROM results')}
    data['directions'] = [dict(case=0, reversed=reverse,
        matching=refs[f"00000/{'reverse' if reverse else 'forward'}", 'matching'],
        catalogue=refs[f"00000/{'reverse' if reverse else 'forward'}", 'deduplication'])
        for reverse in (False, True)]
    data['reduction_source'] = dict(source_identity='original', ledger_sha256='recorded')
    index = tmp_path/'reduce-index.json'; index.write_text(json.dumps(data))
    def forbidden(*args, **kwargs):
        pytest.fail('Reduction-only execution must not run matching or decoding')
    monkeypatch.setattr(WorkPool, 'search', forbidden)
    monkeypatch.setattr(WorkPool, 'decode', forbidden)
    output = tmp_path/'reduced'; options = dict(local=True, workers=2, window=2,
        reduce_only=True, memory_bytes=2**30, case_seconds=123.)
    assert run(index, output, **options)['pass_finished'] == 'reduction_validation'
    with sqlite3.connect(output/'ledger.sqlite') as db:
        records = db.execute('SELECT * FROM results ORDER BY case_name,stage').fetchall()
        assert db.execute('SELECT stage,state,count(*) FROM results GROUP BY stage,state').fetchall() == [
            ('deduplication', 'complete', 2), ('reduction_validation', 'complete', 2), ('snapshot', 'complete', 2)]
    store = BundleStore(output/'bundles')
    for name, stage, state, reference, audit, error in records:
        if stage == 'reduction_validation':
            assert store.load(Blob(**json.loads(reference)))['equal']
        history = store.load(Blob(**json.loads(audit)))
        assert len(history['tasks']) == len(history['attempts']) == 1
        assert history['case'][1] == 123.
    assert (source/'ledger.sqlite').read_bytes() == original_bytes
    assert run(index, output, resume=True, **options)['pass_finished'] == 'reduction_validation'
    with sqlite3.connect(output/'ledger.sqlite') as db:
        assert db.execute('SELECT * FROM results ORDER BY case_name,stage').fetchall() == records
    # An atom-set-equal catalogue with a changed permission must fail the audit.
    from dataclasses import replace
    from enzymemap import compare_reduction
    from mappa.mappa_core.core.scheduler.runtime import Workspace
    original = Blob(**data['directions'][0]['catalogue'])
    saved = store.load(original)
    family = saved.catalogue.families[0]
    saved.catalogue.families[0] = replace(family, actions=())
    with pytest.raises(ValueError, match='changed the saved scientific catalogue'):
        compare_reduction(dict(original=original, catalogue=store.put(saved)), Workspace(tmp_path, store=store))


def test_input_failures_remain_in_denominator_and_stop_next_pass(tmp_path):
    index=prepared_index(tmp_path,corrupt=True)
    output=tmp_path/'result'
    with pytest.raises(RuntimeError,match='matching contains failed'):
        run(index,output,local=True,workers=1)
    with sqlite3.connect(output/'ledger.sqlite') as db:
        assert db.execute('SELECT stage,state,count(*) FROM results GROUP BY stage,state').fetchall()==[('matching','error',2)]


@pytest.mark.parametrize('corrupt', [False, True])
def test_saved_matching_postprocessing_preserves_source_and_stage_barrier(tmp_path, monkeypatch, corrupt):
    from mappa.execution import WorkPool
    from mappa.mappa_core.core.scheduler.bundles import Blob
    index = prepared_index(tmp_path)
    source = tmp_path/'original'
    run(index, source, local=True, workers=1, stop_after='matching')
    original = (source/'ledger.sqlite').read_bytes()
    data = json.loads(index.read_text())
    with sqlite3.connect(source/'ledger.sqlite') as db:
        matching = {name: json.loads(reference) for name, reference in
                    db.execute("SELECT case_name,reference FROM results WHERE stage='matching'")}
    for row in data['directions']:
        name = f"00000/{'reverse' if row['reversed'] else 'forward'}"
        row['matching'] = matching[name]
    if corrupt:
        data['directions'][0]['matching']['sha256'] = '0'*64
    data['reduction_source'] = dict(source_identity='recorded', ledger_sha256='recorded')
    data['decode_config'] = dict(event_slack=5)
    index = tmp_path/'postprocess.json'
    index.write_text(json.dumps(data))
    output = tmp_path/'postprocessed'
    def forbidden(*args, **kwargs):
        pytest.fail('Saved matching must never be searched again')
    decode = WorkPool.decode
    def after_reduction(pool, *args, **kwargs):
        with sqlite3.connect(output/'ledger.sqlite') as db:
            assert db.execute("SELECT count(*) FROM results WHERE stage='deduplication' AND state='complete'").fetchone()[0] == 2
        return decode(pool, *args, **kwargs)
    monkeypatch.setattr(WorkPool, 'search', forbidden)
    monkeypatch.setattr(WorkPool, 'decode', after_reduction)
    options = dict(local=True, workers=1, window=1, postprocess_only=True,
                   memory_bytes=2**30, case_seconds=123.)
    result = run(index, output, **options)
    assert (source/'ledger.sqlite').read_bytes() == original
    with sqlite3.connect(output/'ledger.sqlite') as db:
        records = db.execute('SELECT * FROM results ORDER BY case_name,stage').fetchall()
    if corrupt:
        assert result['stopped_on_failure'] and result['allocation_drained']
        assert len(records) == 1 and records[0][1:3] == ('snapshot', 'error')
        return
    assert result['finished']
    assert {row[1] for row in records} == {'snapshot', 'deduplication', 'decoding', 'verification'}
    assert len(records) == 8 and all(row[2] == 'complete' for row in records)
    store = BundleStore(output/'bundles')
    for name, stage, state, reference, audit, error in records:
        saved = store.load(Blob(**json.loads(reference)))
        if stage == 'decoding':
            assert saved.config.event_slack == 5
        if stage == 'verification':
            assert saved['status'] == 'recovered'
        assert store.load(Blob(**json.loads(audit)))['case'][1] == 123.
    assert run(index, output, resume=True, **options)['finished']
    with sqlite3.connect(output/'ledger.sqlite') as db:
        assert db.execute('SELECT * FROM results ORDER BY case_name,stage').fetchall() == records


def test_seed_batch_is_saved_and_changes_task_count_without_changing_scope(tmp_path):
    from mappa.mappa_core.core.scheduler.bundles import Blob
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint
    index = prepared_index(tmp_path)
    data = json.loads(index.read_text())
    data['config']['seed_count'] = 10
    index.write_text(json.dumps(data))
    output = tmp_path/'result'
    run(index, output, local=True, workers=2, seed_batch=10, stop_after='matching')
    assert json.loads((output/'plan.json').read_text())['seed_batch'] == 10
    store = BundleStore(output/'bundles')
    with sqlite3.connect(output/'ledger.sqlite') as db:
        for reference, audit in db.execute("SELECT reference,audit FROM results WHERE stage='matching'"):
            matched = store.load(Blob(**json.loads(reference)))
            history = store.load(Blob(**json.loads(audit)))
            assert len(history['tasks']) == len(matched.scope)
            for task in history['tasks']:
                specification = _loads_checkpoint(task[1])
                assert specification.operation.endswith(':search_complete_cut')
                assert specification.payload['indices'] == tuple(range(10))
            for attempt in history['attempts']:
                measured = _loads_checkpoint(attempt[-1])['measurements']['search']
                assert [r['index'] for r in measured['seeds']] == list(range(10))
            for ref in matched.cuts:
                cut = store.load(ref)
                assert cut.seeds == 10
                assert [r['index'] for r in cut.search_counts[0]['seeds']] == list(range(10))
    with pytest.raises(ValueError, match='Resume input/configuration differs'):
        run(index, output, local=True, workers=1, resume=True, seed_batch=1)


@pytest.mark.parametrize('verify_only, stop_after, verification_failed', [
    (True, None, False),
    (True, None, True),
    (False, None, False),
    (False, 'decoding', False),
    (False, 'verification', False),
    (False, None, True),
])
def test_verification_preserves_decode_failures_and_existing_results(
        tmp_path, monkeypatch, verify_only, stop_after, verification_failed):
    from mappa.execution import WorkPool
    from mappa.mappa_core.core.scheduler.bundles import Blob
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint

    index = prepared_index(tmp_path)
    output = tmp_path/'result'
    memory = 2*2**30
    run(index, output, local=True, workers=2, memory_bytes=memory)
    store = BundleStore(output/'bundles')
    name = '00000/forward'
    with sqlite3.connect(output/'ledger.sqlite') as db:
        # Seed terminal resource failures without allocating GiBs or waiting an hour.
        failures = [('decoding', 'memory_limit')]
        if verification_failed:
            failures.append(('verification', 'timeout'))
        for stage, state in failures:
            raw, = db.execute('SELECT audit FROM results WHERE case_name=? AND stage=?',
                              (name, stage)).fetchone()
            audit = store.load(Blob(**json.loads(raw)))
            audit['tasks'] = [tuple(state if i == 4 else value for i, value in enumerate(task))
                              for task in audit['tasks']]
            audit['attempts'] = [tuple(state if i == 5 else value for i, value in enumerate(attempt))
                                 for attempt in audit['attempts']]
            budget = list(audit['case'])
            budget[3] = budget[1] if stage == 'verification' else 37.
            audit['case'] = tuple(budget)
            db.execute("UPDATE results SET state='error',reference=NULL,audit=?,error=? WHERE case_name=? AND stage=?",
                       (json.dumps(asdict(store.put(audit))), json.dumps(state), name, stage))
        if not verification_failed:
            db.execute("DELETE FROM results WHERE case_name=? AND stage='verification'", (name,))
        before = {(row[0], row[1]): row for row in db.execute('SELECT * FROM results')}
        assert db.execute('SELECT count(*) FROM retries').fetchone()[0] == 0

    entered_pools = []
    original_enter = WorkPool.__enter__
    def enter_pool(pool):
        entered_pools.append(pool.identity)
        return original_enter(pool)
    monkeypatch.setattr(WorkPool, '__enter__', enter_pool)
    options = dict(local=True, workers=1, memory_bytes=memory, resume=True,
                   verify_only=verify_only, stop_after=stop_after)
    if verification_failed or not verify_only:
        failed_stage = 'verification' if verification_failed else 'decoding'
        with pytest.raises(RuntimeError, match=f'{failed_stage} contains failed'):
            run(index, output, **options)
    else:
        status = run(index, output, **options)
        assert status['pass_finished'] == 'verification'
        assert not status.get('finished', False)
    assert len(entered_pools) == 1

    with sqlite3.connect(output/'ledger.sqlite') as db:
        after = {(row[0], row[1]): row for row in db.execute('SELECT * FROM results')}
        assert db.execute('SELECT count(*) FROM retries').fetchone()[0] == 0
    assert {key: after[key] for key in before} == before
    # Failed decode attempts and their consumed budget retain their original artifact.
    decode_audit = store.load(Blob(**json.loads(after[name, 'decoding'][4])))
    assert decode_audit['case'][3] == 37.
    assert decode_audit['tasks'][0][4] == 'memory_limit'
    if verification_failed:
        assert after == before
        failed_audit = store.load(Blob(**json.loads(after[name, 'verification'][4])))
        assert failed_audit['case'][3] == failed_audit['case'][1] == 3600.
    elif stop_after == 'decoding':
        assert after == before
    else:
        assert set(after)-set(before) == {(name, 'verification')}
        result = store.load(Blob(**json.loads(after[name, 'verification'][3])))
        assert result['status'] == 'recovered'
        assert result['policy'] == 'exact annotated heavy-atom indices; no endpoint-symmetry relaxation'
        audit = store.load(Blob(**json.loads(after[name, 'verification'][4])))
        task, = audit['tasks']
        specification = _loads_checkpoint(task[1])
        assert specification.operation == 'enzymemap:verify'
        assert specification.memory_bytes == memory
        assert specification.case_seconds == 3600.
    progress = json.loads((output/'progress.json').read_text())
    assert progress.get('verification_only', False) == verify_only
    assert not progress.get('finished', False)
    if not verify_only:
        assert progress['failed_stages'] == (['decoding', 'verification']
            if verification_failed else ['decoding'])


def test_verify_only_preserves_failed_reduction_and_verifies_available_catalogue(tmp_path):
    from mappa.mappa_core.core.scheduler.bundles import Blob
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint

    index = prepared_index(tmp_path)
    output = tmp_path/'result'
    memory = 2*2**30
    run(index, output, local=True, workers=1, memory_bytes=memory, stop_after='deduplication')
    store = BundleStore(output/'bundles')
    failed, available = '00000/forward', '00000/reverse'
    with sqlite3.connect(output/'ledger.sqlite') as db:
        raw, = db.execute("SELECT audit FROM results WHERE case_name=? AND stage='deduplication'",
                          (failed,)).fetchone()
        audit = store.load(Blob(**json.loads(raw)))
        audit['tasks'] = [tuple('memory_limit' if i == 4 else value for i, value in enumerate(task))
                          for task in audit['tasks']]
        audit['attempts'] = [tuple('memory_limit' if i == 5 else value for i, value in enumerate(attempt))
                             for attempt in audit['attempts']]
        audit['case'] = tuple(37. if i == 3 else value for i, value in enumerate(audit['case']))
        saved_audit = store.put(audit)
        db.execute("UPDATE results SET state='error',reference=NULL,audit=?,error=? "
                   "WHERE case_name=? AND stage='deduplication'",
                   (json.dumps(asdict(saved_audit)), json.dumps('memory_limit'), failed))
        before = {(row[0], row[1]): row for row in db.execute('SELECT * FROM results')}

    options = dict(local=True, workers=1, memory_bytes=memory, resume=True, verify_only=True)
    with pytest.raises(RuntimeError, match='verification contains failed directions'):
        run(index, output, **options)
    with sqlite3.connect(output/'ledger.sqlite') as db:
        after = {(row[0], row[1]): row for row in db.execute('SELECT * FROM results')}
        assert db.execute('SELECT count(*) FROM retries').fetchone()[0] == 0
    assert {key: after[key] for key in before} == before
    assert store.load(saved_audit) == audit
    assert set(after)-set(before) == {(failed, 'verification'), (available, 'verification')}
    missing = after[failed, 'verification']
    assert missing[2:5] == ('error', None, None)
    assert 'Previous pass did not complete' in json.loads(missing[5])
    completed = after[available, 'verification']
    assert completed[2] == 'complete'
    assert store.load(Blob(**json.loads(completed[3])))['status'] == 'recovered'
    verification_audit = store.load(Blob(**json.loads(completed[4])))
    task, = verification_audit['tasks']
    specification = _loads_checkpoint(task[1])
    assert specification.operation == 'enzymemap:verify'
    assert specification.memory_bytes == memory
    assert specification.case_seconds == 3600.
    progress = json.loads((output/'progress.json').read_text())
    assert progress['verification_only'] and progress['failed_stages'] == ['verification']
    assert not progress.get('finished', False)

    with pytest.raises(RuntimeError, match='verification contains failed directions'):
        run(index, output, **options)
    with sqlite3.connect(output/'ledger.sqlite') as db:
        assert {(row[0], row[1]): row for row in db.execute('SELECT * FROM results')} == after
        assert db.execute('SELECT count(*) FROM retries').fetchone()[0] == 0


def test_verify_only_requires_original_campaign_and_terminal_reduction(tmp_path):
    with pytest.raises(ValueError, match='requires --resume'):
        run(tmp_path/'missing.json', tmp_path/'unused', verify_only=True)
    assert not (tmp_path/'unused').exists()

    index = prepared_index(tmp_path)
    output = tmp_path/'result'
    run(index, output, local=True, workers=1, stop_after='matching')
    progress_path = output/'progress.json'
    progress = json.loads(progress_path.read_text())
    progress_path.write_text(json.dumps(dict(progress, source='another-engine')))
    with pytest.raises(ValueError, match='source or compiled engine differs'):
        run(index, output, local=True, resume=True, verify_only=True)
    progress_path.write_text(json.dumps(progress))
    for state in (None, 'pending', 'running', 'complete'):
        if state is not None:
            with sqlite3.connect(output/'ledger.sqlite') as db:
                db.executemany('INSERT OR REPLACE INTO results VALUES (?,?,?,?,?,?)',
                               [(f'00000/{direction}', 'deduplication', state, None, None, None)
                                for direction in ('forward', 'reverse')])
        ledger = (output/'ledger.sqlite').read_bytes()
        with pytest.raises(ValueError, match='terminal reduction for every direction'):
            run(index, output, local=True, resume=True, verify_only=True)
        assert (output/'ledger.sqlite').read_bytes() == ledger
        assert json.loads(progress_path.read_text()) == progress


@pytest.mark.parametrize('unarchived', ['tasks', 'attempts'])
def test_verify_only_rejects_unarchived_queue_evidence_without_writing(tmp_path, unarchived):
    from mappa.mappa_core.core.scheduler.bundles import Blob

    index = prepared_index(tmp_path)
    output = tmp_path/'result'
    run(index, output, local=True, workers=1, stop_after='deduplication')
    with sqlite3.connect(output/'ledger.sqlite') as db:
        raw, = db.execute("SELECT audit FROM results WHERE case_name='00000/forward' AND stage='deduplication'").fetchone()
    audit = BundleStore(output/'bundles').load(Blob(**json.loads(raw)))
    with sqlite3.connect(output/'bundles/queue.sqlite') as db:
        # Restore persisted work evidence as if its coordinator stopped before retirement.
        db.execute('INSERT INTO cases VALUES (?,?,?,?,?)', audit['case'])
        record = audit[unarchived][0]
        marks = ','.join('?' for _ in record)
        db.execute(f'INSERT INTO {unarchived} VALUES ({marks})', record)
    before = {path.relative_to(output): path.read_bytes()
              for path in output.rglob('*') if path.is_file()}

    with pytest.raises(ValueError, match='unarchived tasks or attempts remain'):
        run(index, output, local=True, workers=1, resume=True, verify_only=True)

    after = {path.relative_to(output): path.read_bytes()
             for path in output.rglob('*') if path.is_file()}
    assert after == before


def test_retry_ledger_retains_original_failure_and_audit(tmp_path):
    from enzymemap import Ledger
    output=tmp_path/'result'; output.mkdir()
    scratch=tmp_path/'scratch'; scratch.mkdir()
    store=BundleStore(output/'bundles')
    audit=store.put({'attempts': ['original memory failure']})
    ledger=Ledger(output,scratch,resume=False)
    ledger.put('case','matching','error',audit=audit,error='memory_limit')
    assert ledger.retry_audit('case','matching')==audit
    ledger.retire_failure('case','matching','diagnosed native export defect')
    ledger.put('case','matching','complete',reference=store.put('replacement'))
    ledger.checkpoint(); ledger.db.close()
    with sqlite3.connect(output/'ledger.sqlite') as db:
        assert db.execute('SELECT state FROM results').fetchone()==('complete',)
        reason,raw,error=db.execute('SELECT reason,audit,error FROM retries').fetchone()
        assert reason=='diagnosed native export defect'
        assert json.loads(raw)==asdict(audit)
        assert json.loads(error)=='memory_limit'


def test_resource_retry_recovers_old_timeout_but_keeps_same_memory_failure(tmp_path):
    from mappa.mappa_core.core.scheduler.bundles import Blob
    index = prepared_index(tmp_path)
    output = tmp_path/'result'
    run(index, output, local=True, workers=2, stop_after='matching')
    store = BundleStore(output/'bundles')
    originals = {}
    with sqlite3.connect(output/'ledger.sqlite') as db:
        rows = db.execute('SELECT case_name,audit FROM results ORDER BY case_name').fetchall()
        for (case, raw), memory in zip(rows, (2*2**30, 4*2**30)):
            audit = store.load(Blob(**json.loads(raw)))
            originals[case] = audit
            task = list(audit['tasks'][-1])
            task[3], task[4] = memory, 'timeout' if memory == 2*2**30 else 'memory_limit'
            audit['tasks'][-1] = tuple(task)
            budget = list(audit['case']); budget[3] = budget[1]
            audit['case'] = tuple(budget)
            db.execute("UPDATE results SET state='error',reference=NULL,audit=? WHERE case_name=?",
                (json.dumps(asdict(store.put(audit))), case))
    with pytest.raises(RuntimeError, match='matching contains failed'):
        run(index, output, local=True, workers=2, resume=True, stop_after='matching',
            retry_memory_limits=True, case_seconds=60.)
    with sqlite3.connect(output/'ledger.sqlite') as db:
        assert db.execute('SELECT case_name,state FROM results ORDER BY case_name').fetchall() == [
            ('00000/forward', 'complete'), ('00000/reverse', 'error')]
        assert db.execute('SELECT case_name FROM retries').fetchall() == [('00000/forward',)]
        raw, = db.execute("SELECT audit FROM results WHERE case_name='00000/forward'").fetchone()
    recovered = store.load(Blob(**json.loads(raw)))
    original = originals['00000/forward']
    assert len(recovered['attempts']) == len(original['attempts'])+1
    assert recovered['tasks'][:-1] == original['tasks'][:-1]
    assert recovered['retries'][0][4:] == (60., 4*2**30)


def test_batches_are_disjoint_reuse_results_and_stop_before_decoding(tmp_path):
    from prepare_batches import prepare
    index=prepared_index(tmp_path)
    data=json.loads(index.read_text())
    data['decode_config'] = dict(threshold=.5, metal_threshold=.3, event_slack=5, fragment_seconds=5.)
    data['directions'] += [dict(row,case=1) for row in data['directions']]
    index.write_text(json.dumps(data))
    baseline=tmp_path/'baseline'
    assert run(index,baseline,local=True,workers=2,stop_after='matching')['pass_finished']=='matching'
    prepare(index,baseline,tmp_path/'batches',count=2,
            additional_ledgers=[baseline/'ledger.sqlite'])
    for batch in range(2):
        root=tmp_path/'batches'/str(batch)
        assert {r['case'] for r in json.loads((root/'index.json').read_text())['directions']}=={batch}
        with sqlite3.connect(root/'result/ledger.sqlite') as db:
            assert db.execute('SELECT stage,count(*) FROM results GROUP BY stage').fetchall()==[('matching',2)]
        alias=tmp_path/f'alias-{batch}'
        alias.symlink_to(root, target_is_directory=True)
        assert run(alias/'index.json',alias/'result',local=True,workers=1,resume=True)['finished']


def test_fresh_batches_preserve_inputs_and_import_no_old_results(tmp_path):
    from prepare_batches import prepare_fresh
    index = prepared_index(tmp_path)
    data = json.loads(index.read_text())
    data['decode_config'] = dict(threshold=.5, metal_threshold=.3, event_slack=5, fragment_seconds=5.)
    data['directions'] = [dict(row, case=case) for case in range(8)
                          for row in data['directions']]
    index.write_text(json.dumps(data))
    output = tmp_path/'fresh'
    manifest = prepare_fresh([index], output, workers=[3, 1])
    assert manifest['directions'] == 16
    assert 'growth_policy' not in manifest['config']
    assert manifest['reused_scientific_results'] == 0
    assert not list(output.rglob('*.sqlite'))
    assert not list(output.rglob('*.bundle'))
    batches = [json.loads((output/str(i)/'index.json').read_text()) for i in range(2)]
    assert manifest['decode_config'] == data['decode_config']
    assert all(batch['decode_config'] == data['decode_config'] for batch in batches)
    assert [len(batch['directions']) for batch in batches] == [12, 4]
    actual = [row for batch in batches for row in batch['directions']]
    key = lambda row: (row['case'], row['reversed'])
    assert sorted(actual, key=key) == sorted(data['directions'], key=key)
    assert not ({r['case'] for r in batches[0]['directions']} &
                {r['case'] for r in batches[1]['directions']})
    with pytest.raises(FileExistsError):
        prepare_fresh([index], output, workers=[3, 1])


def test_fresh_batches_reject_duplicate_and_incomplete_inputs(tmp_path):
    from prepare_batches import prepare_fresh
    index = prepared_index(tmp_path)
    with pytest.raises(ValueError, match='Duplicate reaction/direction'):
        prepare_fresh([index, index], tmp_path/'duplicate', workers=[1])
    data = json.loads(index.read_text())
    data['directions'].pop()
    index.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='both directions'):
        prepare_fresh([index], tmp_path/'partial', workers=[1])


def failed_decode_campaign(tmp_path, *, bad_catalogue=False, decode_config=None, stage='decoding'):
    """Keep a terminal campaign, with two diagnosed failures and two successes."""
    from prepare_batches import prepare_fresh
    from mappa.mappa_core.core.scheduler.bundles import Blob
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint
    import pickle

    index = prepared_index(tmp_path)
    data = json.loads(index.read_text())
    if decode_config is not None:
        data['decode_config'] = decode_config
    data['directions'] += [dict(row, case=1) for row in data['directions']]
    index.write_text(json.dumps(data))
    source = tmp_path/'original'
    prepare_fresh([index], source, workers=[1, 1])
    selected = {}
    for batch in range(2):
        output = source/str(batch)/'result'
        run(output.parent/'index.json', output, local=True, workers=1, memory_bytes=2*2**30)
        store = BundleStore(output/'bundles')
        with sqlite3.connect(output/'ledger.sqlite') as ledger:
            name, raw = ledger.execute("SELECT case_name,audit FROM results WHERE stage=? ORDER BY case_name LIMIT 1", (stage,)).fetchone()
            record = store.load(Blob(**json.loads(raw)))
            task_row = list(record['tasks'][0])
            if bad_catalogue and batch == 0:
                task = _loads_checkpoint(task_row[1])
                task.payload['catalogue'] = store.put('invalid catalogue type')
                task_row[1] = pickle.dumps(task, protocol=5)
                ledger.execute("UPDATE results SET reference=? WHERE case_name=? AND stage='deduplication'",
                    (json.dumps(asdict(task.payload['catalogue'])), name))
            task_row[4] = 'memory_limit'
            record['tasks'] = [tuple(task_row)]
            record['attempts'] = [tuple('memory_limit' if i == 5 else item for i, item in enumerate(attempt))
                                  for attempt in record['attempts']]
            record['case'] = tuple(37. if i == 3 else item for i, item in enumerate(record['case']))
            audit = store.put(record)
            ledger.execute("UPDATE results SET state='error',reference=NULL,audit=?,error=? WHERE case_name=? AND stage=?",
                (json.dumps(asdict(audit)), json.dumps({task_row[0]: 'memory_limit'}), name, stage))
        progress = output/'progress.json'
        progress.write_text(json.dumps(dict(json.loads(progress.read_text()), source='original-code')))
        selected[name] = dict(audit=audit, record=record, output=output)
    (source/'launch.json').write_text(json.dumps(dict(code_identity='original-code')))
    return source, selected


def fresh_decoding_index(source, destination):
    """Reference all saved catalogues, regardless of the previous decode outcome."""
    import hashlib

    manifest = source/'manifest.json'
    data = json.loads(manifest.read_text())
    rows = []
    for batch in data['batches']:
        with sqlite3.connect(source/str(batch['batch'])/'result/ledger.sqlite') as ledger:
            for name, reference in ledger.execute("SELECT case_name,reference FROM results WHERE stage='deduplication'"):
                case, direction = name.split('/')
                rows.append(dict(case=int(case), reversed=direction == 'reverse', catalogue=json.loads(reference)))
    destination.write_text(json.dumps(dict(config=data['config'], directions=rows,
        decoding_source=dict(source_commit='original-revision', source_identity='original-code',
                             manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest()))))
    return destination


def test_fresh_decoding_uses_all_catalogues_with_new_budgets_and_preserves_source(tmp_path, monkeypatch):
    import enzymemap
    from mappa import EventDecodeConfig
    from mappa.execution import WorkPool
    from mappa.mappa_core.core.scheduler.bundles import Blob
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint

    source, _ = failed_decode_campaign(tmp_path)
    index = fresh_decoding_index(source, tmp_path/'decode-index.json')
    data = json.loads(index.read_text())
    data['decode_config'] = dict(threshold=.75, metal_threshold=None, event_slack=6, fragment_seconds=5.)
    index.write_text(json.dumps(data))
    source_files = {path.relative_to(source): path.read_bytes() for path in source.rglob('*') if path.is_file()}
    def forbidden(*args, **kwargs):
        raise AssertionError('A fresh decoding pass cannot run other stages or restore an old attempt')
    for method in ('search', 'deduplicate', 'retry_stage'):
        monkeypatch.setattr(WorkPool, method, forbidden)
    monkeypatch.setattr(enzymemap, 'prepared_problem', forbidden)
    monkeypatch.setattr(enzymemap, 'read_blob', forbidden)
    output = tmp_path/'new-decoding'
    options = dict(local=True, workers=2, window=2, decode_only=True,
                   memory_bytes=2**30, case_seconds=123.)
    status = run(index, output, **options)
    assert status['pass_finished'] == 'decoding' and status['decoding_only']
    assert not status.get('finished', False)
    plan = json.loads((output/'plan.json').read_text())
    assert plan['mode'] == 'decode_only'
    assert plan['decoding_source'] == data['decoding_source']
    assert plan['decoder_source'] == enzymemap.code_identity()
    assert plan['directions_to_decode'] == data['directions']
    assert plan['decode_config'] == data['decode_config']
    assert plan['case_seconds'] == 123. and plan['memory_bytes'] == 2**30
    with sqlite3.connect(output/'ledger.sqlite') as ledger:
        records = ledger.execute('SELECT * FROM results ORDER BY case_name').fetchall()
        assert ledger.execute('SELECT count(*) FROM retries').fetchone()[0] == 0
    assert len(records) == 4
    store = BundleStore(output/'bundles')
    for name, stage, state, reference, audit, error in records:
        assert (stage, state, error) == ('decoding', 'complete', None)
        result = store.load(Blob(**json.loads(reference)))
        assert result.config == EventDecodeConfig(**data['decode_config'])
        history = store.load(Blob(**json.loads(audit)))
        assert len(history['attempts']) == len(history['tasks']) == 1
        assert history['case'][1] == 123. and history['case'][3] < 123.
        task = _loads_checkpoint(history['tasks'][0][1])
        assert task.operation == 'mappa.mappa_core.core.scheduler.decoding:decode_catalogue'
        assert task.case_seconds == 123. and task.memory_bytes == 2**30
    bundles = sorted((output/'bundles').glob('*.bundle'))
    assert run(index, output, resume=True, **options)['pass_finished'] == 'decoding'
    assert sorted((output/'bundles').glob('*.bundle')) == bundles
    with sqlite3.connect(output/'ledger.sqlite') as ledger:
        assert ledger.execute('SELECT * FROM results ORDER BY case_name').fetchall() == records
    assert {path.relative_to(source): path.read_bytes() for path in source.rglob('*') if path.is_file()} == source_files
    for changed in (dict(case_seconds=124.), dict(memory_bytes=2*2**30)):
        with pytest.raises(ValueError, match='Resume input/configuration differs'):
            run(index, output, resume=True, **dict(options, **changed))
    monkeypatch.setattr(enzymemap, 'code_identity', lambda: 'changed-decoder')
    with pytest.raises(ValueError, match='Resume input/configuration differs'):
        run(index, output, resume=True, **options)


def test_fresh_decoding_preserves_new_failures_without_retry_or_other_stages(tmp_path):
    source, _ = failed_decode_campaign(tmp_path, bad_catalogue=True)
    index = fresh_decoding_index(source, tmp_path/'decode-index.json')
    output = tmp_path/'new-decoding'
    options = dict(local=True, workers=1, window=1, decode_only=True, case_seconds=123.)
    with pytest.raises(RuntimeError, match='decoding contains failed directions'):
        run(index, output, **options)
    with sqlite3.connect(output/'ledger.sqlite') as ledger:
        records = ledger.execute('SELECT * FROM results ORDER BY case_name').fetchall()
    assert len(records) == 4
    assert sorted(row[2] for row in records) == ['complete', 'complete', 'complete', 'error']
    assert {row[1] for row in records} == {'decoding'}
    plan = json.loads((output/'plan.json').read_text())
    assert plan['decode_config'] == dict(threshold=.5, metal_threshold=.3, event_slack=6, fragment_seconds=5.)
    with pytest.raises(RuntimeError, match='decoding contains failed directions'):
        run(index, output, resume=True, **options)
    with sqlite3.connect(output/'ledger.sqlite') as ledger:
        assert ledger.execute('SELECT * FROM results ORDER BY case_name').fetchall() == records


@pytest.mark.parametrize('options', [dict(verify_only=True), dict(recover_stage='decoding'),
    dict(retry_memory_limits=True), dict(stop_after='matching'), dict(stop_after='verification')])
def test_fresh_decoding_rejects_competing_stage_modes(tmp_path, options):
    with pytest.raises(ValueError, match='Fresh decoding cannot select another pass'):
        run(tmp_path/'unused', tmp_path/'output', decode_only=True, **options)
    assert not (tmp_path/'output').exists()


@pytest.mark.parametrize('stage', ['decoding', 'verification'])
def test_stage_recovery_explicit_resources_preserve_exhausted_attempts(tmp_path, stage):
    from mappa.mappa_core.core.scheduler.bundles import Blob
    source, selected = failed_decode_campaign(tmp_path, stage=stage)
    for name, entry in selected.items():
        record = entry['record']
        record['case'] = tuple(record['case'][1] if i == 3 else value
                               for i, value in enumerate(record['case']))
        record['tasks'] = [tuple('timeout' if i == 4 else value
                                for i, value in enumerate(record['tasks'][0]))]
        record['attempts'] = [tuple('timeout' if i == 5 else value for i, value in enumerate(attempt))
                              for attempt in record['attempts']]
        store = BundleStore(entry['output']/'bundles')
        with sqlite3.connect(entry['output']/'ledger.sqlite') as ledger:
            ledger.execute('UPDATE results SET audit=? WHERE case_name=? AND stage=?',
                           (json.dumps(asdict(store.put(record))), name, stage))
    output = tmp_path/'recovery'
    options = dict(local=True, workers=2, window=2, memory_bytes=4*2**30,
                   recover_stage=stage, retry_seconds=120.)
    assert run(source, output, **options)['pass_finished'] == stage
    with sqlite3.connect(output/'ledger.sqlite') as ledger:
        rows = ledger.execute('SELECT * FROM results').fetchall()
    assert len(rows) == 2 and all(row[2] == 'complete' for row in rows)
    store = BundleStore(output/'bundles')
    for name, _, _, _, audit, _ in rows:
        before = selected[name]['record']
        after = store.load(Blob(**json.loads(audit)))
        assert after['tasks'][0][:3] == before['tasks'][0][:3]
        assert after['tasks'][0][3] == 4*2**30
        assert after['case'][1] == before['case'][3] + 120.
        assert after['case'][3] >= before['case'][3]
        assert after['attempts'][:-1] == before['attempts']
        assert after['retries'][-1][1:] == (before['case'][0], before['case'][1],
                                           before['case'][3], 120., 4*2**30)
    assert run(source, output, resume=True, **options)['pass_finished'] == stage


def test_decoding_recovery_is_separate_budgeted_and_resumable(tmp_path, monkeypatch):
    import enzymemap
    from mappa.execution import WorkPool
    from mappa.mappa_core.core.scheduler.bundles import Blob
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint

    policy = dict(threshold=.5, metal_threshold=.3, event_slack=5, fragment_seconds=5.)
    source, selected = failed_decode_campaign(tmp_path, decode_config=policy)
    original_files = {path.relative_to(source): path.read_bytes() for path in source.rglob('*') if path.is_file()}
    output = tmp_path/'recovery'
    options = dict(local=True, workers=2, window=2, memory_bytes=2*2**30, recover_stage='decoding')
    original_enter = WorkPool.__enter__
    def enter_recovery(pool):
        with pytest.raises(RuntimeError, match='active coordinator or stage recovery'):
            run(source, tmp_path/'duplicate', **options)
        assert not (tmp_path/'duplicate').exists()
        for entry in selected.values():
            source_pool = WorkPool.local(store=BundleStore(entry['output']/'bundles'))
            with pytest.raises(RuntimeError, match='already has an active coordinator'):
                original_enter(source_pool)
        return original_enter(pool)
    monkeypatch.setattr(WorkPool, '__enter__', enter_recovery)
    status = run(source, output, **options)
    assert status['pass_finished'] == 'decoding'
    assert not status.get('finished', False)
    assert status['expected_directions'] == 2
    plan = json.loads((output/'plan.json').read_text())
    assert plan['original_source'] == 'original-code'
    assert plan['decoder_source'] == enzymemap.code_identity()
    assert plan['memory_bytes'] == 2*2**30
    registration = source/'decoding-recovery.json'
    assert json.loads(registration.read_text())['output'] == str(output.resolve())
    assert {row['task'] for row in plan['directions_to_recover']} == {
        entry['record']['tasks'][0][0] for entry in selected.values()}
    with sqlite3.connect(output/'ledger.sqlite') as db:
        rows = db.execute('SELECT * FROM results ORDER BY case_name').fetchall()
        assert db.execute('SELECT count(*) FROM retries').fetchone()[0] == 0
    assert len(rows) == 2
    store = BundleStore(output/'bundles')
    for name, stage, state, reference, audit, error in rows:
        assert (stage, state, error) == ('decoding', 'complete', None)
        previous = selected[name]['record']
        recovered = store.load(Blob(**json.loads(audit)))
        assert recovered['case'][0:2] == previous['case'][0:2]
        assert previous['case'][3] <= recovered['case'][3] < previous['case'][1]
        assert recovered['attempts'][0] == previous['attempts'][0]
        assert len(recovered['attempts']) == len(previous['attempts'])+1
        assert recovered['tasks'][0][0:4] == previous['tasks'][0][0:4]
        assert not recovered['retries']
        task = _loads_checkpoint(recovered['tasks'][0][1])
        assert store.load(Blob(**json.loads(reference))).config == task.payload['config']
        assert asdict(task.payload['config']) == policy
    assert run(source, output, resume=True, **options)['pass_finished'] == 'decoding'
    with sqlite3.connect(output/'ledger.sqlite') as db:
        assert db.execute('SELECT * FROM results ORDER BY case_name').fetchall() == rows
    after = {path.relative_to(source): path.read_bytes() for path in source.rglob('*') if path.is_file()}
    assert after.pop(Path('decoding-recovery.json')) == registration.read_bytes()
    assert after == original_files
    with pytest.raises(ValueError, match='registered stage recovery'):
        run(source, tmp_path/'later-duplicate', **options)
    assert not (tmp_path/'later-duplicate').exists()
    with pytest.raises(ValueError, match='registered stage recovery'):
        run(source, output, **options)

    with sqlite3.connect(output/'bundles/queue.sqlite') as queue:
        queue.execute('INSERT INTO attempts VALUES (?,?,?,?,?,?,?)', next(iter(selected.values()))['record']['attempts'][0])
    before = {path.relative_to(output): path.read_bytes() for path in output.rglob('*') if path.is_file()}
    with pytest.raises(ValueError, match='unarchived tasks or attempts remain'):
        run(source, output, resume=True, **options)
    assert {path.relative_to(output): path.read_bytes() for path in output.rglob('*') if path.is_file()} == before
    with sqlite3.connect(output/'bundles/queue.sqlite') as queue:
        queue.execute('DELETE FROM attempts')
    before = {path.relative_to(output): path.read_bytes() for path in output.rglob('*') if path.is_file()}
    monkeypatch.setattr(enzymemap, 'code_identity', lambda: 'changed-decoder')
    with pytest.raises(ValueError, match='registered stage recovery'):
        run(source, output, resume=True, **options)
    assert {path.relative_to(output): path.read_bytes() for path in output.rglob('*') if path.is_file()} == before


def test_decoding_recovery_attempts_all_selected_and_does_not_retry_new_failures(tmp_path):
    source, selected = failed_decode_campaign(tmp_path, bad_catalogue=True)
    output = tmp_path/'recovery'
    options = dict(local=True, workers=1, window=1, memory_bytes=2*2**30, recover_stage='decoding')
    with pytest.raises(RuntimeError, match='decoding contains failed directions'):
        run(source, output, **options)
    with sqlite3.connect(output/'ledger.sqlite') as db:
        before = db.execute('SELECT * FROM results ORDER BY case_name').fetchall()
    assert len(before) == len(selected) == 2
    assert sorted(row[2] for row in before) == ['complete', 'error']
    with pytest.raises(RuntimeError, match='decoding contains failed directions'):
        run(source, output, resume=True, **options)
    with sqlite3.connect(output/'ledger.sqlite') as db:
        assert db.execute('SELECT * FROM results ORDER BY case_name').fetchall() == before
    assert json.loads((output/'progress.json').read_text())['failed_stages'] == ['decoding']


def failed_recovery(tmp_path, *, retry_seconds=None, failure='memory_limit', stage='decoding'):
    """One second-attempt memory failure beside an already recovered direction."""
    from mappa.mappa_core.core.scheduler.bundles import Blob

    source, selected = failed_decode_campaign(tmp_path, stage=stage)
    output = tmp_path/'first-recovery'
    options = dict(local=True, workers=1, window=1, memory_bytes=2*2**30, recover_stage=stage)
    if retry_seconds is not None:
        options.update(memory_bytes=4*2**30, retry_seconds=retry_seconds)
    run(source, output, **options)
    store = BundleStore(output/'bundles')
    with sqlite3.connect(output/'ledger.sqlite') as ledger:
        rows = ledger.execute('SELECT * FROM results ORDER BY case_name').fetchall()
        failed, completed = rows
        record = store.load(Blob(**json.loads(failed[4])))
        task = list(record['tasks'][0])
        task[4] = failure
        record['tasks'] = [tuple(task)]
        record['attempts'] = [tuple(failure if i == 5 else value for i, value in enumerate(attempt))
                              if attempt[0] == task[6] else attempt for attempt in record['attempts']]
        record['case'] = tuple(73. if i == 3 else value for i, value in enumerate(record['case']))
        # Archive order is not a scientific contract; preserve attempts by token.
        record['attempts'].reverse()
        audit = store.put(record)
        ledger.execute("UPDATE results SET state='error',reference=NULL,audit=?,error=? WHERE case_name=?",
            (json.dumps(asdict(audit)), json.dumps({task[0]: failure}), failed[0]))
    return source, output, failed[0], completed, record, options


@pytest.mark.parametrize('stage', ['decoding', 'verification'])
def test_execution_error_retry_is_explicit_and_preserves_remaining_budget(tmp_path, stage):
    from mappa.mappa_core.core.scheduler.bundles import Blob

    source, previous, name, completed, before, options = failed_recovery(
        tmp_path, failure='error', stage=stage)
    output = tmp_path/'repaired-launch'
    with pytest.raises(ValueError, match='only archived resource failures'):
        run(previous, output, **options)
    assert not output.exists()
    reason = 'Corrected the input relocation configuration after diagnosis'
    original = (previous/'ledger.sqlite').read_bytes()
    assert run(previous, output, retry_error_reason=reason, **options)['pass_finished'] == stage
    plan = json.loads((output/'plan.json').read_text())
    assert plan['retry_error_reason'] == reason and plan['retry_seconds'] is None
    with sqlite3.connect(output/'ledger.sqlite') as ledger:
        rows = ledger.execute('SELECT * FROM results').fetchall()
    assert len(rows) == 1 and rows[0][:3] == (name, stage, 'complete')
    after = BundleStore(output/'bundles').load(Blob(**json.loads(rows[0][4])))
    assert after['case'][:2] == before['case'][:2]
    assert after['case'][3] >= before['case'][3]
    assert after['retries'] == before['retries']
    assert len(after['attempts']) == len(before['attempts']) + 1
    attempts = {row[0]: row for row in after['attempts']}
    assert all(attempts[row[0]] == row for row in before['attempts'])
    assert (previous/'ledger.sqlite').read_bytes() == original
    assert run(previous, output, resume=True, retry_error_reason=reason, **options)['pass_finished'] == stage
    with pytest.raises(ValueError, match='registered stage recovery'):
        run(previous, output, resume=True, retry_error_reason='Different reason', **options)


def test_recovery_chain_preserves_explicit_resource_grants(tmp_path):
    from mappa.mappa_core.core.scheduler.bundles import Blob
    source, previous, name, completed, old_audit, options = failed_recovery(tmp_path, retry_seconds=120.)
    output = tmp_path/'next-recovery'
    options.update(memory_bytes=8*2**30, retry_seconds=240.)
    assert run(previous, output, **options)['pass_finished'] == 'decoding'
    with sqlite3.connect(output/'ledger.sqlite') as ledger:
        rows = ledger.execute('SELECT * FROM results').fetchall()
    assert len(rows) == 1 and rows[0][:3] == (name, 'decoding', 'complete')
    audit = BundleStore(output/'bundles').load(Blob(**json.loads(rows[0][4])))
    assert audit['retries'][:-1] == old_audit['retries']
    assert audit['attempts'][:-1] == old_audit['attempts']
    assert audit['tasks'][0][3] == 8*2**30
    assert audit['case'][1] == old_audit['case'][3] + 240.


def test_decoding_recovery_chain_preserves_latest_budget_history_and_successes(tmp_path, monkeypatch):
    import enzymemap
    import shutil
    from mappa.execution import WorkPool
    from mappa.mappa_core.core.scheduler.bundles import Blob

    source, previous, name, completed, old_audit, options = failed_recovery(tmp_path)
    copy = tmp_path/'copied-predecessor'
    shutil.copytree(previous, copy)
    with pytest.raises(ValueError, match='registered stage recovery'):
        run(copy, tmp_path/'copied-successor', **options)
    assert not (tmp_path/'copied-successor').exists()
    with pytest.raises(ValueError, match='registered stage recovery'):
        run(source, tmp_path/'stale-original', **options)

    output = tmp_path/'second-recovery'
    before = (previous/'ledger.sqlite').read_bytes()
    original_enter = WorkPool.__enter__
    def enter_recovery(pool):
        predecessor_pool = WorkPool.local(store=BundleStore(previous/'bundles'))
        with pytest.raises(RuntimeError, match='already has an active coordinator'):
            original_enter(predecessor_pool)
        return original_enter(pool)
    monkeypatch.setattr(WorkPool, '__enter__', enter_recovery)
    monkeypatch.setattr(enzymemap, 'code_identity', lambda: 'repaired-decoder')
    status = run(previous, output, **options)
    assert status['pass_finished'] == 'decoding'
    assert status['expected_directions'] == 1 and not status.get('finished', False)
    plan = json.loads((output/'plan.json').read_text())
    assert plan['source_kind'] == 'stage_recovery'
    assert plan['recovery_source'] == str(previous.resolve())
    assert plan['original_campaign'] == str(source.resolve())
    assert plan['original_source'] == 'original-code'
    assert plan['decoder_source'] == 'repaired-decoder'
    assert plan['directions_to_recover'][0]['spent_seconds'] == 73.
    assert len(plan['predecessors']) == 1
    with sqlite3.connect(output/'ledger.sqlite') as db:
        rows = db.execute('SELECT * FROM results').fetchall()
    assert len(rows) == 1 and rows[0][:3] == (name, 'decoding', 'complete')
    audit = BundleStore(output/'bundles').load(Blob(**json.loads(rows[0][4])))
    attempts = {attempt[0]: attempt for attempt in audit['attempts']}
    assert len(attempts) == 3 and len(old_audit['attempts']) == 2
    assert all(attempts[attempt[0]] == attempt for attempt in old_audit['attempts'])
    assert audit['tasks'][0][:4] == old_audit['tasks'][0][:4]
    assert audit['case'][:2] == old_audit['case'][:2]
    assert 73. <= audit['case'][3] < audit['case'][1]
    assert audit['retries'] == old_audit['retries'] == []
    assert (previous/'ledger.sqlite').read_bytes() == before
    with sqlite3.connect(previous/'ledger.sqlite') as db:
        assert db.execute('SELECT * FROM results WHERE case_name=?', (completed[0],)).fetchone() == completed
    assert run(previous, output, resume=True, **options)['pass_finished'] == 'decoding'
    with sqlite3.connect(output/'ledger.sqlite') as db:
        assert db.execute('SELECT * FROM results').fetchall() == rows
    with pytest.raises(ValueError, match='registered successor'):
        run(source, previous, resume=True, **options)
    with pytest.raises(ValueError, match='registered stage recovery'):
        run(previous, tmp_path/'stale-predecessor', **options)


@pytest.mark.parametrize('change, message', [
    ('scope', 'incomplete terminal scope'),
    ('budget', 'consumed budget'),
    ('attempts', 'attempts or consumed budget'),
    ('result', 'completed result'),
])
def test_decoding_recovery_chain_rejects_changed_predecessor(tmp_path, change, message):
    source, previous, name, completed, record, options = failed_recovery(tmp_path)
    store = BundleStore(previous/'bundles')
    with sqlite3.connect(previous/'ledger.sqlite') as ledger:
        if change == 'scope':
            ledger.execute('DELETE FROM results WHERE case_name=?', (completed[0],))
        elif change == 'result':
            ledger.execute('UPDATE results SET reference=? WHERE case_name=?',
                           (json.dumps(asdict(store.put('wrong result reference'))), completed[0]))
        else:
            if change == 'budget':
                record['case'] = tuple(36. if i == 3 else value for i, value in enumerate(record['case']))
            else:
                record['attempts'].pop()
            ledger.execute('UPDATE results SET audit=? WHERE case_name=?',
                           (json.dumps(asdict(store.put(record))), name))
    output = tmp_path/'second-recovery'
    with pytest.raises(ValueError, match=message):
        run(previous, output, **options)
    assert not output.exists()
    assert not (previous/'decoding-recovery.json').exists()


@pytest.mark.parametrize('change, message', [
    ('catalogue', 'differs from its saved catalogue'),
    ('thresholds', 'configuration differs'),
    ('event_slack', 'configuration differs'),
    ('identity', 'source identities differ'),
    ('guard', 'original memory guard'),
    ('budget', 'no remaining case budget'),
    ('unarchived', 'unarchived tasks or attempts remain'),
])
def test_decoding_recovery_rejects_changed_original_before_writing(tmp_path, change, message):
    from dataclasses import replace
    import pickle
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint

    source, selected = failed_decode_campaign(tmp_path)
    name, entry = next(iter(selected.items()))
    original_output = entry['output']
    store = BundleStore(original_output/'bundles')
    record = entry['record']
    if change in ('catalogue', 'thresholds', 'event_slack', 'budget'):
        if change == 'budget':
            record['case'] = tuple(record['case'][1] if i == 3 else item for i, item in enumerate(record['case']))
        else:
            task_row = list(record['tasks'][0])
            task = _loads_checkpoint(task_row[1])
            if change == 'catalogue':
                task.payload['catalogue'] = store.put('another catalogue')
            elif change == 'thresholds':
                task.payload['config'] = replace(task.payload['config'], threshold=.9)
            else:
                task.payload['config'] = replace(task.payload['config'], event_slack=3)
            task_row[1] = pickle.dumps(task, protocol=5)
            record['tasks'] = [tuple(task_row)]
        with sqlite3.connect(original_output/'ledger.sqlite') as db:
            db.execute("UPDATE results SET audit=? WHERE case_name=? AND stage='decoding'",
                       (json.dumps(asdict(store.put(record))), name))
    elif change == 'identity':
        path = original_output/'progress.json'
        path.write_text(json.dumps(dict(json.loads(path.read_text()), source='another-original-code')))
    elif change == 'unarchived':
        with sqlite3.connect(original_output/'bundles/queue.sqlite') as db:
            db.execute('INSERT INTO attempts VALUES (?,?,?,?,?,?,?)', record['attempts'][0])
    original_files = {path.relative_to(source): path.read_bytes() for path in source.rglob('*') if path.is_file()}
    output = tmp_path/'recovery'
    with pytest.raises(ValueError, match=message):
        run(source, output, local=True, memory_bytes=(4 if change == 'guard' else 2)*2**30,
            recover_stage='decoding')
    assert not output.exists()
    assert {path.relative_to(source): path.read_bytes() for path in source.rglob('*') if path.is_file()} == original_files
