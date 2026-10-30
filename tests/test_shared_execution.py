"""Scheduling changes must preserve full search graphs and decoded event sets."""
from dataclasses import replace
from pathlib import Path
import time

import pytest

from mappa.mappa_core.core.scheduler import Task, TaskQueue, BundleWriter, read_blob, run, connect, coordinator
from mappa.mappa_core.core.scheduler.search import plan_search, iter_cuts
from mappa.mappa_core.core.search.session import SearchSession
from mappa import AAMSearchConfig, AAMResult
from mappa.mappa_core.domain import AAMSearchMetrics
from mappa.mappa_core.core.search.graph import AAMSearchGraph
from mappa.mappa_core.core.search.branching.symmetry import finalize_graph_symmetry
from test_reused_aam_api import problem


IDENTITY = 'shared-execution-test-v1'
MEMORY = 1024*2**20


def task(key, **options):
    return Task(key, 'execution_helpers:value', options.pop('payload', {}), **options)


def test_dependency_identity_resume_and_single_acceptance(tmp_path):
    path = tmp_path/'queue.db'
    queue = TaskQueue(path, IDENTITY)
    a, b = task('a'), task('b', dependencies=('a',))
    queue.add([b, a])
    claim = queue.claim('first', MEMORY, IDENTITY)
    assert claim['task'] == a
    assert queue.claim('second', MEMORY, IDENTITY)['state'] == 'waiting'
    assert queue.finish('a', claim['token'], 'complete', {'value': 4})
    assert not queue.finish('a', claim['token'], 'complete', {'value': 5})
    queue.close()
    queue = TaskQueue(path, IDENTITY)
    queue.add([a, b])  # completed work is never requeued
    assert queue.claim('second', MEMORY, IDENTITY)['task'] == b
    with pytest.raises(ValueError, match='different work'):
        queue.add([replace(a, payload={'value': 20})])
    queue.close()
    with pytest.raises(ValueError, match='identity'):
        TaskQueue(path, 'different-source')


def test_budget_shared_across_tasks_and_explicit_lost_recovery(tmp_path, monkeypatch):
    import mappa.mappa_core.core.scheduler.queue as module
    clock = [100.]
    monkeypatch.setattr(module.time, 'time', lambda: clock[0])
    queue = TaskQueue(tmp_path/'queue.db', IDENTITY, lease_seconds=2)
    queue.add([task('a', case='same', case_seconds=10), task('b', case='same', case_seconds=10)])
    first = queue.claim('one', MEMORY, IDENTITY)
    clock[0] += 3
    assert queue.snapshot() == {'lost': 1, 'pending': 1}
    assert not queue.finish('a', first['token'], 'complete', {})
    queue.retry('a')
    resumed = queue.claim('two', MEMORY, IDENTITY)
    assert first['deadline'] == 110
    assert resumed['deadline'] == 111
    queue.finish('a', resumed['token'], 'complete', {})
    clock[0] = 111
    queue.lease_seconds = 20
    last = queue.claim('two', MEMORY, IDENTITY)
    assert last['deadline'] == 119  # idle time did not consume the remaining eight seconds
    clock[0] = 120
    assert queue.finish('b', last['token'], 'timeout', {})
    assert queue.snapshot() == {'complete': 1, 'timeout': 1}
    queue.close()


def test_resources_dependencies_and_backup(tmp_path):
    queue = TaskQueue(tmp_path/'queue.db', IDENTITY)
    queue.add([task('a', memory_bytes=2*MEMORY), task('b', dependencies=('a',))])
    assert queue.claim('tiny', MEMORY//2, IDENTITY)['state'] == 'insufficient_memory'
    with pytest.raises(ValueError, match='identity'):
        queue.claim('wrong', MEMORY, 'other')
    with pytest.raises(ValueError, match='Unknown dependency'):
        queue.add([task('missing', dependencies=('absent',))])
    from graphlib import CycleError
    with pytest.raises(CycleError):
        queue.add([task('c', dependencies=('d',)), task('d', dependencies=('c',))])
    claim = queue.claim('big', 2*MEMORY, IDENTITY)
    queue.finish('a', claim['token'], 'memory_limit', {})
    assert queue.snapshot() == {'memory_limit': 1, 'blocked': 1}
    queue.backup(tmp_path/'backup.db')
    restored = TaskQueue(tmp_path/'backup.db', IDENTITY)
    assert restored.results() == queue.results()
    restored.close()
    queue.close()


def test_checkpoint_does_not_block_worker_heartbeat(tmp_path, monkeypatch):
    import sqlite3
    import threading
    import mappa.mappa_core.core.scheduler.queue as module
    queue = TaskQueue(tmp_path/'queue.db', IDENTITY)
    queue.add([task('a')])
    claim = queue.claim('worker', MEMORY, IDENTITY)
    copying, release = threading.Event(), threading.Event()
    class SlowBackup(sqlite3.Connection):
        def backup(self, *args, **kwargs):
            copying.set()
            assert release.wait(5)
            return super().backup(*args, **kwargs)
    connect = sqlite3.connect
    def connection(*args, **kwargs):
        if kwargs.get('uri'):
            kwargs['factory'] = SlowBackup
        return connect(*args, **kwargs)
    monkeypatch.setattr(module.sqlite3, 'connect', connection)
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(2) as pool:
        backup = pool.submit(queue.backup, tmp_path/'snapshot.db')
        assert copying.wait(5)
        heartbeat = pool.submit(queue.heartbeat, 'a', claim['token'])
        try:
            assert heartbeat.result(timeout=2)
        finally:
            release.set()
        backup.result(timeout=5)
    restored = TaskQueue(tmp_path/'snapshot.db', IDENTITY)
    assert restored.snapshot() == {'running':1}
    restored.close();queue.close()


def test_memory_handoff_preserves_live_work_specs_and_deadline(tmp_path):
    queue = TaskQueue(tmp_path/'queue.db', IDENTITY)
    queue.add([task('a', case='same', memory_bytes=MEMORY),
               task('b', case='same', memory_bytes=MEMORY)])
    original = queue.db.execute('SELECT key,specification FROM tasks ORDER BY key').fetchall()
    active = queue.claim('old', MEMORY, IDENTITY)
    maintenance = TaskQueue(tmp_path/'queue.db', IDENTITY)
    assert maintenance.set_pending_memory(2*MEMORY, reason='Drain old workers') == 1
    assert queue.claim('old', MEMORY, IDENTITY)['state'] == 'waiting'
    queue.finish('a', active['token'], 'complete', {'value': 7})
    assert queue.claim('old', MEMORY, IDENTITY)['state'] == 'insufficient_memory'
    assert maintenance.set_pending_memory(MEMORY//2, reason='Lower-memory first pass') == 1
    next_task = queue.claim('new', MEMORY//2, IDENTITY)
    assert next_task['task'].key == 'b'
    assert next_task['deadline'] >= active['deadline']
    assert queue.db.execute('SELECT key,specification FROM tasks ORDER BY key').fetchall() == original
    assert queue.results(['a'])['a']['result'] == {'value': 7}
    maintenance.close()
    queue.close()


def test_larger_memory_stage_retries_only_authorized_limits(tmp_path):
    queue = TaskQueue(tmp_path/'queue.db', IDENTITY)
    queue.add([task('small', memory_bytes=MEMORY),
               task('dependent', dependencies=('small',), memory_bytes=MEMORY),
               task('broken', memory_bytes=MEMORY), task('already-large', memory_bytes=MEMORY)])
    first = queue.claim('first-pass', MEMORY, IDENTITY)
    queue.finish('small', first['token'], 'memory_limit', {'memory_guard_bytes': MEMORY})
    error = queue.claim('first-pass', MEMORY, IDENTITY)
    queue.finish('broken', error['token'], 'error', {'memory_guard_bytes': MEMORY})
    large = queue.claim('large', MEMORY, IDENTITY)
    queue.finish('already-large', large['token'], 'memory_limit', {'memory_guard_bytes': 4*MEMORY})
    assert queue.retry_memory_limits(MEMORY, 4*MEMORY, reason='Authorized memory second stage') == ['small']
    assert queue.claim('first-pass', MEMORY, IDENTITY)['state'] == 'insufficient_memory'
    second = queue.claim('second-pass', 4*MEMORY, IDENTITY)
    assert second['task'].key == 'small' and second['deadline'] >= first['deadline']
    queue.finish('small', second['token'], 'complete', {'value': 1})
    assert queue.claim('first-pass', MEMORY, IDENTITY)['task'].key == 'dependent'
    assert queue.results(['broken'])['broken']['state'] == 'error'
    assert queue.results(['already-large'])['already-large']['state'] == 'memory_limit'
    assert queue.db.execute("SELECT state FROM attempts WHERE task='small' ORDER BY started").fetchall() == [('memory_limit',), ('complete',)]
    assert queue.retry_memory_limits(MEMORY, 4*MEMORY, reason='No repeated escalation') == []
    queue.close()


def test_memory_stage_reuses_completed_outputs_and_refuses_live_source(tmp_path):
    import json
    from rebenchmark.shared import prepare_memory_stage
    source, output, scratch = (tmp_path/name for name in ('first', 'second', 'scratch'))
    for folder in (source, output, scratch):
        folder.mkdir()
    queue = TaskQueue(source/'queue.sqlite', IDENTITY)
    queue.add([task('saved'), task('limited')])
    saved = queue.claim('first', MEMORY, IDENTITY)
    queue.finish('saved', saved['token'], 'complete', {'published': 'unchanged'})
    limited = queue.claim('first', MEMORY, IDENTITY)
    queue.finish('limited', limited['token'], 'memory_limit', {'memory_guard_bytes': MEMORY})
    queue.backup(source/'final.sqlite')
    queue.close()
    (source/'final.sqlite').replace(source/'queue.sqlite')
    plan = dict(previous_stage=str(source), identity=IDENTITY,
                previous_memory_bytes=MEMORY, memory_bytes=4*MEMORY)
    (source/'progress.json').write_text(json.dumps(dict(finished=False, states={'running': 1})))
    with pytest.raises(ValueError, match='still owns'):
        prepare_memory_stage(plan, output, scratch)
    (source/'progress.json').write_text(json.dumps(dict(finished=True, states={'complete': 1, 'memory_limit': 1})))
    stage = prepare_memory_stage(plan, output, scratch)
    assert stage.snapshot() == {'complete': 1, 'pending': 1}
    assert stage.results(['saved'])['saved']['result'] == {'published': 'unchanged'}
    resumed = stage.claim('second', 4*MEMORY, IDENTITY)
    assert resumed['task'].key == 'limited' and resumed['deadline'] >= limited['deadline']
    assert json.loads((output/'selection.json').read_text())['keys'] == ['limited']
    stage.close()


def test_bundle_hash_and_file_count(tmp_path):
    writer = BundleWriter(tmp_path/'bundles', tmp_path)
    blobs = [writer.put({'index': n, 'data': bytes(range(255))}) for n in range(40)]
    assert len(list((tmp_path/'bundles').iterdir())) == 1
    for n, blob in enumerate(blobs):
        assert read_blob(blob, tmp_path)['index'] == n
    with open(blobs[3].path, 'r+b') as stream:
        stream.seek(blobs[3].offset+10)
        stream.write(b'bad')
    with pytest.raises(ValueError, match='checksum'):
        read_blob(blobs[3], tmp_path)
    assert read_blob(blobs[4], tmp_path)['index'] == 4


@pytest.mark.parametrize('workers,remote', [(1, False), (3, False), (3, True)])
def test_common_worker_dependencies_and_failure_isolation(tmp_path, workers, remote):
    queue = TaskQueue(tmp_path/'queue.db', IDENTITY)
    queue.add([task('bad', payload={'fail': True}), task('a', payload={'value': 3}),
        task('b', payload={'value': 7}),
        Task('sum', 'execution_helpers:total', ('a', 'b'), dependencies=('a', 'b'))])
    def execute(q):
        run(q, workers=workers, identity=IDENTITY, scratch=tmp_path,
            bundles=tmp_path/'bundles', memory_bytes=MEMORY)
    if remote:
        with coordinator(queue, ('127.0.0.1', 0), b'test-secret') as address:
            execute(connect(address, b'test-secret'))
    else:
        execute(queue)
    rows = queue.results()
    assert queue.snapshot() == {'complete': 3, 'error': 1}
    assert read_blob(rows['sum']['result']['blob'], tmp_path) == 10
    assert 'deliberate test failure' in rows['bad']['result']['error']
    queue.close()


def test_actual_deadline_does_not_block_other_cases(tmp_path):
    queue = TaskQueue(tmp_path/'queue.db', IDENTITY)
    queue.add([task('slow', payload={'sleep': 5}, case_seconds=.2), task('fast')])
    run(queue, identity=IDENTITY, scratch=tmp_path, bundles=tmp_path/'bundles', memory_bytes=MEMORY)
    assert queue.snapshot() == {'timeout': 1, 'complete': 1}
    queue.close()


def test_budget_counts_parallel_work_once_and_excludes_publication_and_idle(tmp_path, monkeypatch):
    import mappa.mappa_core.core.scheduler.queue as module
    clock = [100.]
    monkeypatch.setattr(module.time, 'time', lambda: clock[0])
    queue = TaskQueue(tmp_path/'queue.db', IDENTITY, lease_seconds=1000)
    queue.add([task(k, case='same', case_seconds=10) for k in ('a','b','c','d')])
    a = queue.claim('one', MEMORY, IDENTITY)
    clock[0] = 102
    b = queue.claim('two', MEMORY, IDENTITY)
    assert b['deadline'] == a['deadline'] == 110
    clock[0] = 104
    assert queue.computed('a', a['token'])
    assert not queue.computed('b', 'wrong-token')
    clock[0] = 106
    assert queue.computed('b', b['token'])
    assert queue.db.execute('SELECT spent,deadline FROM cases').fetchone() == (6, None)
    clock[0] = 200  # publication can finish much later without consuming computation time
    queue.finish('a', a['token'], 'complete', {})
    queue.finish('b', b['token'], 'complete', {})
    queue.backup(tmp_path/'backup.db')
    queue.close()
    clock[0] = 1000  # an allocation gap survives checkpoint/restart
    queue = TaskQueue(tmp_path/'backup.db', IDENTITY)
    c = queue.claim('three', MEMORY, IDENTITY)
    assert c['deadline'] == 1004
    clock[0] = 1005
    assert queue.snapshot() == {'complete': 2, 'running': 1, 'timeout': 1}
    queue.finish('c', c['token'], 'timeout', {})
    assert queue.snapshot() == {'complete': 2, 'timeout': 2}
    queue.close()


@pytest.mark.parametrize('batch', [1, 2, 5, 10])
@pytest.mark.parametrize('sizes', [(5, 7), (7, 5)])
def test_seed_batch_graph_and_prefix_equality(batch, sizes):
    case = problem(*sizes)
    config = AAMSearchConfig(seed_count=10, iso_tolerance=1., branch_limit=4)
    session = SearchSession(case, config)
    expected, counts = session.search_cut(())
    parts = [session.search_cut((), seed_indices=tuple(range(i, min(i+batch, 10))))[0]
             for i in range(0, 10, batch)]
    assert AAMSearchGraph.combine(parts) == expected
    from rebenchmark.search import select_contexts
    for n in range(1, 11):
        prefix = SearchSession(case, replace(config, seed_count=n)).search_cut(())[0]
        assert select_contexts(expected, range(n)) == prefix


@pytest.mark.parametrize('workers,remote,batch', [(1, False, 2), (3, False, 1), (3, True, 5)])
def test_multicase_tasks_equal_original_search_and_complete_decoding(tmp_path, workers, remote, batch):
    from mappa.mappa_core.post_processing.hierarchical_mapping import refine_mapping
    from mappa.mappa_core.post_processing.hierarchical_mapping import RefinementConfig
    queue = TaskQueue(tmp_path/'queue.db', IDENTITY)
    writer = BundleWriter(tmp_path/'inputs', tmp_path)
    expected = []
    for name, sizes in [('forward', (4, 5)), ('reverse', (5, 4)), ('balanced', (4, 4))]:
        case = problem(*sizes)
        config = AAMSearchConfig(seed_count=3, iso_tolerance=1., branch_limit=4)
        tasks, keys = plan_search(name, case, config, writer, seed_batch=batch, memory_bytes=MEMORY)
        queue.add(tasks)
        expected.append((case, config, keys))
    def execute(q):
        run(q, workers=workers, identity=IDENTITY, scratch=tmp_path,
            bundles=tmp_path/'bundles', memory_bytes=MEMORY)
    if remote:
        with coordinator(queue, ('127.0.0.1', 0), b'test-secret') as address:
            execute(connect(address, b'test-secret'))
    else:
        execute(queue)
    assert set(queue.snapshot()) == {'complete'}, queue.results()
    from mappa.mappa_core.core.search.planning import cut_sweep_items
    for case, config, keys in expected:
        session = SearchSession(case, config)
        for cut, observed in zip(cut_sweep_items(case.reactant.wbo, config.cut_floor), iter_cuts(queue, keys, tmp_path)):
            raw, _ = session.search_cut(cut)
            original, _ = finalize_graph_symmetry(raw, session.context[2], iso_tolerance=config.iso_tolerance)
            assert observed['graph'] == original
            if not cut and case.balanced:
                def decode(graph):
                    result = AAMResult(case, config, graph, AAMSearchMetrics.from_record({}, 0))
                    decoded = refine_mapping(result, RefinementConfig(max_events=8))
                    assert decoded.complete
                    return {str(sorted(c.events.items())) for c in decoded.candidates}
                assert decode(original) == decode(observed['graph'])
    queue.close()


def test_abrupt_exit_keeps_other_work_moving(tmp_path):
    queue = TaskQueue(tmp_path/'queue.db', IDENTITY)
    queue.add([task('dead', payload={'exit': True}), task('healthy')])
    run(queue, identity=IDENTITY, scratch=tmp_path, bundles=tmp_path/'bundles', memory_bytes=MEMORY)
    assert queue.snapshot() == {'error': 1, 'complete': 1}
    queue.close()


def test_parallel_reference_evaluation_equals_existing_evaluator(tmp_path):
    from golden_evaluation import prepare
    from mappa import AAMSearchPlan
    from rebenchmark.shared import plan_direction
    from rebenchmark.analysis import analyze_seed
    import json
    case, features, mapping = prepare('[CH3:1][OH:2]>>[CH3:1][OH:2]')
    config = AAMSearchConfig(seed_count=2, iso_tolerance=1, sweep_cuts=False)
    plan = AAMSearchPlan(case, case, config, False)
    reference = dict(features=features, mapping=sorted(mapping.items()))
    queue = TaskQueue(tmp_path/'queue.db', IDENTITY)
    writer = BundleWriter(tmp_path/'inputs', tmp_path)
    tasks, cuts, evaluation_keys = plan_direction('control', plan, reference, writer, seed_batch=1)
    queue.add(tasks)
    run(queue, workers=2, identity=IDENTITY, scratch=tmp_path,
        bundles=tmp_path/'bundles', memory_bytes=MEMORY)
    assert set(queue.snapshot()) == {'complete'}, queue.results()
    graph = next(iter_cuts(queue, cuts, tmp_path))['graph']
    rows = queue.results()
    for seed, key in enumerate(evaluation_keys):
        actual = read_blob(rows[key]['result']['blob'], tmp_path)
        expected = analyze_seed(plan, reference, graph, seed, tmp_path/f'direct-{seed}', 0)
        verdict = json.loads((tmp_path/f'direct-{seed}'/'evaluation.json').read_text())
        assert actual['summary'] == expected
        for field in ('reference_recovery', 'representative_covered'):
            assert actual['evaluation'].get(field) == verdict.get(field)
        assert actual['classes'] == json.loads((tmp_path/f'direct-{seed}'/'classes.json').read_text())
        assert actual['evidence'].size > 0
    queue.close()


def test_actual_memory_guard_releases_worker_for_other_cases(tmp_path):
    limit = 512*2**20
    queue = TaskQueue(tmp_path/'queue.db', IDENTITY)
    queue.add([task('large', payload={'allocate': 640*2**20, 'sleep': 5}, memory_bytes=limit),
               task('small', memory_bytes=limit)])
    run(queue, identity=IDENTITY, scratch=tmp_path, bundles=tmp_path/'bundles', memory_bytes=limit)
    assert queue.snapshot() == {'memory_limit': 1, 'complete': 1}
    assert queue.results()['large']['result']['peak_tree_rss_bytes'] > limit
    queue.close()


def test_scheduled_catalogue_conversion_equals_public_decoder(tmp_path):
    from golden_evaluation import prepare
    from mappa import search_aam
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig
    from mappa.mappa_core.post_processing.decoder.api import decode_events
    from mappa.mappa_core.core.scheduler.decoding import plan_event_decoding
    case, _, _ = prepare('[CH3:1][OH:2]>>[CH3:1][OH:2]')
    aam = search_aam(case, AAMSearchConfig(seed_count=2, sweep_cuts=False), execution='reused_native')
    catalogue = FinalBranchCatalogue(case).add_aam(aam)
    config = EventDecodeConfig()
    expected = decode_events(catalogue, config)
    assert expected.catalogue is catalogue
    writer = BundleWriter(tmp_path/'inputs', tmp_path)
    queue = TaskQueue(tmp_path/'queue.db', IDENTITY)
    tasks = plan_event_decoding('decode', case, catalogue, config, writer)
    queue.add(tasks)
    run(queue, workers=2, identity=IDENTITY, scratch=tmp_path,
        bundles=tmp_path/'bundles', memory_bytes=MEMORY)
    assert set(queue.snapshot()) == {'complete'}
    patterns = set()
    for row in queue.results().values():
        result = read_blob(row['result']['blob'], tmp_path)
        assert result['mapping_space_preserved']
        patterns.update(p['id'] for p in result['representatives'])
    assert patterns == {expected.evaluate(i)['id'] for i in range(len(expected.catalogue.families))}
    queue.close()


def test_incremental_queue_plan_keeps_history_and_fetches_selected_results(tmp_path):
    queue = TaskQueue(tmp_path/'queue.db', IDENTITY)
    queue.add([task('first')])
    claim = queue.claim('worker', MEMORY, IDENTITY)
    queue.finish('first', claim['token'], 'error', {'error': 'diagnosed'})
    queue.retry('first')
    claim = queue.claim('worker', MEMORY, IDENTITY)
    queue.finish('first', claim['token'], 'complete', {'ok': True})
    queue.add([task('second', dependencies=('first',))])
    assert list(queue.results(['second'])) == ['second']
    assert queue.claim('worker', MEMORY, IDENTITY)['task'].key == 'second'
    assert queue.db.execute("SELECT state FROM attempts WHERE task='first' ORDER BY started").fetchall() == [('error',), ('complete',)]
    queue.close()


def test_queue_backup_is_one_complete_copyable_file(tmp_path):
    import shutil
    queue=TaskQueue(tmp_path/'live.sqlite',IDENTITY)
    queue.add([task('one')])
    target=tmp_path/'snapshot.sqlite'
    queue.backup(target)
    claim=queue.claim('worker',MEMORY,IDENTITY)
    queue.finish('one',claim['token'],'complete',{'value':42})
    queue.backup(target)
    assert not target.with_name(target.name+'-wal').exists()
    copied=tmp_path/'copied.sqlite';shutil.copyfile(target,copied)
    restored=TaskQueue(copied,IDENTITY)
    assert restored.snapshot()=={'complete':1}
    assert restored.results()['one']['result']=={'value':42}
    restored.close();queue.close()


def test_queue_optimization_preserves_dependency_and_claim_results(tmp_path):
    queue=TaskQueue(tmp_path/'queue.sqlite',IDENTITY)
    queue.add([task('parent'),task('child',dependencies=('parent',))])
    queue.optimize()
    claim=queue.claim('worker',MEMORY,IDENTITY)
    assert claim['task'].key=='parent'
    queue.finish('parent',claim['token'],'error',{'reason':'test'})
    assert queue.claim('worker',MEMORY,IDENTITY)['state']=='finished'
    assert queue.snapshot()=={'error':1,'blocked':1}
    queue.close()


def test_independent_tasks_do_not_scan_payloads_during_dependency_refresh(tmp_path):
    queue = TaskQueue(tmp_path/'queue.sqlite', IDENTITY)
    queue.add([task(str(i), payload=b'x'*65536, case='batch') for i in range(600)])
    queue.optimize()
    calls = 0
    def instructions():
        nonlocal calls
        calls += 1
        return calls > 40
    queue.db.set_progress_handler(instructions, 100)
    try:
        # With no dependency edges, maintenance must not visit every payload.
        with queue.db:
            queue._refresh(time.time())
    finally:
        queue.db.set_progress_handler(None, 0)
        queue.close()


def test_checkpoint_finishes_under_continuous_heartbeat_writes(tmp_path, monkeypatch):
    import sqlite3
    import mappa.mappa_core.core.scheduler.queue as module
    queue = TaskQueue(tmp_path/'live.sqlite', IDENTITY)
    queue.add([task('large', payload=b'x'*(8*2**20))])
    claim = queue.claim('worker', MEMORY, IDENTITY)
    expires = queue.db.execute('SELECT expires FROM tasks').fetchone()[0]
    pages = queue.db.execute('PRAGMA page_count').fetchone()[0]
    writes = 0
    class UpdatingBackup(sqlite3.Connection):
        def backup(self, target, **options):
            def update(status, remaining, total):
                nonlocal writes
                writes += 1
                assert writes < pages // 256 + 5, 'Checkpoint restarted after heartbeat write'
                assert queue.heartbeat('large', claim['token'])
            return super().backup(target, progress=update, **options)
    connect = sqlite3.connect
    def connection(*args, **kwargs):
        if kwargs.get('uri'):
            kwargs['factory'] = UpdatingBackup
        return connect(*args, **kwargs)
    monkeypatch.setattr(module.sqlite3, 'connect', connection)
    target = tmp_path/'published.sqlite'
    try:
        queue.backup(target)
        with connect(target) as saved:
            assert saved.execute('SELECT expires FROM tasks').fetchone()[0] == expires
            assert saved.execute('PRAGMA quick_check').fetchone()[0] == 'ok'
        assert writes > 1
        assert queue.db.execute('SELECT expires FROM tasks').fetchone()[0] > expires
    finally:
        queue.close()


def test_stop_tolerates_child_exit_before_final_signal(monkeypatch):
    from mappa.mappa_core.core.scheduler.runtime import _stop
    from types import SimpleNamespace
    import mappa.mappa_core.core.scheduler.runtime as runtime
    child = SimpleNamespace(pid=12345, is_alive=lambda: True, join=lambda *args: None)
    monkeypatch.setattr(runtime.os, 'getpgid', lambda pid: pid)
    def gone(*args):
        raise ProcessLookupError()
    monkeypatch.setattr(runtime.os, 'killpg', gone)
    _stop(child)


def test_resource_logs_separate_computation_from_publication(tmp_path):
    queue = TaskQueue(tmp_path/'queue.sqlite', IDENTITY)
    queue.add([task('record', payload={'sleep': .3})])
    run(queue, identity=IDENTITY, scratch=tmp_path, bundles=tmp_path/'bundles', memory_bytes=MEMORY)
    result = queue.results()['record']['result']
    assert set(result['phase_peak_tree_rss_bytes']) == {'compute','publication'}
    assert result['peak_tree_rss_bytes'] == max(result['phase_peak_tree_rss_bytes'].values())
    assert result['phase_peak_tree_rss_bytes']['compute'] > 0
    assert result['compute_seconds'] >= .3
    assert result['total_cpu_seconds'] >= result['cpu_seconds'] + result['publish_cpu_seconds']
    assert result['process_lifetime_peak_rss_bytes'] > 0
    assert result['phase_at_exit'] == 'publication'
    queue.close()


def test_publication_memory_limit_preserves_compute_measurements(tmp_path):
    limit = 512*2**20
    queue = TaskQueue(tmp_path/'queue.sqlite', IDENTITY)
    queue.add([Task('publish', 'execution_helpers:publication_allocation', {}, memory_bytes=limit)])
    run(queue, identity=IDENTITY, scratch=tmp_path, bundles=tmp_path/'bundles', memory_bytes=limit)
    row = queue.results()['publish']
    assert row['state'] == 'memory_limit'
    result = row['result']
    assert result['phase_at_exit'] == 'publication'
    assert result['phase_peak_tree_rss_bytes']['publication'] > limit
    assert result['compute_seconds'] >= 0 and result['cpu_seconds'] >= 0
    queue.close()


def test_fresh_conversion_matches_direct_public_conversion_and_reuses_pilot(tmp_path):
    from dataclasses import asdict
    import hashlib
    import json
    from golden_evaluation import prepare
    from rebenchmark.shared import conversion_tasks, prepare_queue
    from rebenchmark.decoding import decode_graph
    case, _, _ = prepare('[CH3:1][OH:2]>>[CH3:1][OH:2]')
    config = AAMSearchConfig(seed_count=2, iso_tolerance=1, sweep_cuts=False)
    writer = BundleWriter(tmp_path/'inputs',tmp_path)
    tasks,cuts = plan_search('00000',case,config,writer,seed_batch=1,memory_bytes=MEMORY)
    tasks += list(conversion_tasks('00000',tasks,cuts))
    pilot = tmp_path/'pilot'; pilot.mkdir()
    queue = TaskQueue(pilot/'queue.sqlite',IDENTITY)
    queue.add(tasks)
    run(queue,workers=2,identity=IDENTITY,scratch=tmp_path,bundles=tmp_path/'bundles',memory_bytes=MEMORY)
    assert queue.snapshot() == {'complete':5}
    rows = queue.results()
    graph = read_blob(rows[cuts[0]]['result']['blob'],tmp_path)['graph']
    for task_spec in tasks[-2:]:
        result = read_blob(rows[task_spec.key]['result']['blob'],tmp_path)
        direct = decode_graph(case,graph,task_spec.payload)
        assert result['decoded'] == direct['decoded']
        assert result['decoded']['mapping_space_preserved']
        assert result['search_rerun'] is True
        assert result['finalized_sha256'] == rows[cuts[0]]['result']['blob'].sha256
    queue.backup(pilot/'published.sqlite')
    queue.close()
    # Move the consistent backup into the location used by production reuse.
    (pilot/'published.sqlite').replace(pilot/'queue.sqlite')
    index = tmp_path/'index.json'
    index.write_text(json.dumps(dict(directions=[dict(direction=0,blob=asdict(writer.put(tasks)))])))
    plan = dict(index=str(index),index_sha256=hashlib.sha256(index.read_bytes()).hexdigest(),
        identity=IDENTITY,commit='test-source',memory_bytes=MEMORY,pilots=[],pilot_commits=['test-source'])
    (pilot/'plan.json').write_text(json.dumps(plan))
    (pilot/'progress.json').write_text(json.dumps(dict(finished=True)))
    full = tmp_path/'full'; full.mkdir()
    scratch = tmp_path/'scratch'; scratch.mkdir()
    plan['pilots'] = [str(pilot)]
    reused = prepare_queue(plan,full,scratch)
    assert reused.snapshot() == {'complete':5}
    assert reused.db.execute('SELECT seconds FROM cases').fetchone()[0] < 3600
    assert len(json.loads((full/'reused.json').read_text())) == 5
    reused.close()
    plan['pilot_commits'] = []
    other = tmp_path/'other'; other.mkdir()
    with pytest.raises(ValueError,match='Pilot commit'):
        prepare_queue(plan,other,other)


def test_allocation_retirement_does_not_claim_or_reset_pending_work(tmp_path):
    queue = TaskQueue(tmp_path/'queue.sqlite', IDENTITY)
    queue.add([task('later', case='same',case_seconds=100)])
    run(queue,identity=IDENTITY,scratch=tmp_path,bundles=tmp_path/'bundles',
        memory_bytes=MEMORY,claim_until=time.time()-1)
    assert queue.snapshot() == {'pending':1}
    assert queue.db.execute('SELECT deadline FROM cases').fetchone()[0] is None
    assert queue.db.execute('SELECT count(*) FROM attempts').fetchone()[0] == 0
    run(queue,identity=IDENTITY,scratch=tmp_path,bundles=tmp_path/'bundles',memory_bytes=MEMORY)
    assert queue.snapshot() == {'complete':1}
    queue.close()


def test_dependency_refresh_skips_successful_parents(tmp_path):
    queue = TaskQueue(tmp_path/'queue.db', IDENTITY)
    queue.add([task('parent')])
    queue.db.execute("UPDATE tasks SET state='complete'")
    queue.db.executemany('INSERT INTO dependencies VALUES (?,?)',
        [(f'child/{i}', 'parent') for i in range(50000)])
    queue.db.commit()
    calls = [0]
    def budget():
        calls[0] += 1
        return calls[0] > 50
    queue.db.set_progress_handler(budget,1000)
    queue._refresh(time.time())
    queue.db.set_progress_handler(None,0)
    queue.close()


def test_enzyme_planner_keeps_input_independent_and_checks_both_directions(tmp_path, monkeypatch):
    from dataclasses import asdict
    from golden_evaluation import prepare
    from rebenchmark.shared import _initialize_planner, _prepare_record
    from mappa.mappa_core.core.scheduler import Blob
    import numpy as np
    monkeypatch.setenv('TMPDIR',str(tmp_path))
    _initialize_planner(str(tmp_path),str(tmp_path))
    config = AAMSearchConfig(seed_count=1,iso_tolerance=1,sweep_cuts=False)
    plan = dict(dataset='EnzymeMap',directions=[False,True],config=asdict(config),seed_batch=1,
        case_seconds=60,memory_bytes=MEMORY,slack=5,reference_checks=True)
    item = dict(index=0,rxn='CO>>CO',ground_truth='[CH3:1][OH:2]>>[CH3:1][OH:2]',original_id='example')
    rows = _prepare_record((plan,item,str(tmp_path)))
    queue = TaskQueue(tmp_path/'queue.db',IDENTITY)
    for row in rows:
        assert row['counts'] == dict(search_batch=1,finalize_cut=1,verify_seed=1,convert_seed=1)
        queue.add(read_blob(Blob(**row['blob']),tmp_path))
    run(queue,workers=2,identity=IDENTITY,scratch=tmp_path,bundles=tmp_path/'results',memory_bytes=MEMORY)
    assert queue.snapshot() == {'complete':8}
    checks = [read_blob(r['result']['blob'],tmp_path) for k,r in queue.results().items() if '/verify/' in k]
    assert all(r['reference_recovery']=='recovered' for r in checks)
    original,_,_ = prepare(item['rxn'])
    changed = _prepare_record((plan,dict(item,ground_truth='[CH4:1]>>[CH4:1]'),str(tmp_path)))
    tasks = read_blob(Blob(**changed[0]['blob']),tmp_path)
    record = read_blob(tasks[0].payload['input'],tmp_path)
    assert np.array_equal(record['reactant']['wbo'],original.reactant.wbo)
    reference = read_blob(Blob(**changed[0]['reference']),tmp_path)
    assert not reference['input_reference_index_agrees']
    queue.close()


def test_enzyme_preparation_freezes_complete_grid_in_worker_bundles(tmp_path, monkeypatch):
    import csv, gzip, hashlib, json
    from rebenchmark.shared import prepare_campaign
    monkeypatch.setenv('TMPDIR',str(tmp_path))
    dataset = tmp_path/'dataset.csv.gz'
    with gzip.open(dataset,'wt') as stream:
        writer = csv.DictWriter(stream,fieldnames=['r_id','rxn','ground_truth','original_id'])
        writer.writeheader()
        for i in range(3):
            writer.writerow(dict(r_id=i,rxn='CO>>CO',ground_truth='[CH3:1][OH:2]>>[CH3:1][OH:2]',original_id=i))
    manifest = tmp_path/'manifest.json'
    manifest.write_text(json.dumps(dict(cases=3,dataset_sha256=hashlib.sha256(dataset.read_bytes()).hexdigest())))
    plan = dict(dataset='EnzymeMap',dataset_path=str(dataset),input_manifest=str(manifest),
        input_manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest(),directions=[False,True],
        config=dict(seed_count=1,iso_tolerance=1,sweep_cuts=False),seed_batch=1,case_seconds=60,
        memory_bytes=MEMORY,slack=5,reference_checks=True,preparation_workers=2)
    (tmp_path/'plan.json').write_text(json.dumps(plan))
    prepare_campaign(tmp_path)
    index = json.loads((tmp_path/'index.json').read_text())
    assert index['cases']==3 and index['tasks']==24
    assert [r['direction'] for r in index['directions']]==list(range(6))
    assert len(list((tmp_path/'inputs').glob('*.bundle'))) <= 2


def test_repeat_stage_reuses_exact_parents_and_keeps_source_immutable(tmp_path):
    import hashlib, json
    from dataclasses import asdict
    from rebenchmark.shared import prepare_repeat_queue
    prior_path = tmp_path/'prior'; prior_path.mkdir()
    output = tmp_path/'repeat'; output.mkdir()
    scratch = tmp_path/'scratch'; scratch.mkdir()
    parent = task('00000/cut/00000/final', case='00000')
    child = replace(task('00000/cut/00000/seed/00', case='00000', dependencies=(parent.key,)),
                    operation='rebenchmark.shared:convert_seed')
    writer = BundleWriter(tmp_path/'bundles', scratch)
    blob = writer.put([parent, child])
    index = tmp_path/'index.json'
    index.write_text(json.dumps(dict(directions=[dict(direction=0,blob=asdict(blob))])))
    index_hash = hashlib.sha256(index.read_bytes()).hexdigest()
    (prior_path/'plan.json').write_text(json.dumps(dict(index_sha256=index_hash)))
    (prior_path/'progress.json').write_text(json.dumps(dict(finished=True,states=dict(complete=2))))
    prior = TaskQueue(prior_path/'queue.sqlite','prior')
    prior.add([parent, child])
    for expected in (parent,child):
        claim = prior.claim('old', MEMORY, 'prior')
        assert claim['task']==expected
        prior.finish(expected.key,claim['token'],'complete',dict(value=expected.key))
    prior.close()
    before = (prior_path/'queue.sqlite').read_bytes()
    plan = dict(index=str(index), index_sha256=index_hash,identity='new',pilots=[],
        repeat_from=str(prior_path),source_queue_sha256=hashlib.sha256(before).hexdigest(),
        repeat_operations=[child.operation])
    queue = prepare_repeat_queue(plan,output,scratch)
    assert queue.snapshot()==dict(complete=1,pending=1)
    assert queue.db.execute('SELECT count(*) FROM attempts').fetchone()[0]==0
    claim = queue.claim('new',MEMORY,'new')
    assert claim['task']==child
    assert claim['dependencies']=={parent.key:dict(value=parent.key)}
    assert queue.db.execute('SELECT spent FROM cases').fetchone()[0]==0
    queue.close()
    assert (prior_path/'queue.sqlite').read_bytes()==before
    plan['source_queue_sha256']='wrong'
    with pytest.raises(ValueError,match='hash differs'):
        prepare_repeat_queue(plan,output,scratch)
    (prior_path/'progress.json').write_text(json.dumps(dict(finished=False,states=dict(running=1))))
    with pytest.raises(ValueError,match='completed source'):
        prepare_repeat_queue(plan,output,scratch)
