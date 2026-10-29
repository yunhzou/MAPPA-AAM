"""Scheduling changes must preserve full search graphs and decoded event sets."""
from dataclasses import replace
from pathlib import Path
import time

import pytest

from mappa.mappa_core.core.scheduler import Task, TaskQueue, BundleWriter, read_blob, run, connect, coordinator
from mappa.mappa_core.core.scheduler.search import plan_search, iter_cuts
from mappa.mappa_core.core.search.session import SearchSession
from mappa import AAMSearchConfig, AAMResult
from mappa.domain import AAMSearchMetrics
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
    assert resumed['deadline'] == first['deadline'] == 110
    queue.finish('a', resumed['token'], 'complete', {})
    clock[0] = 111
    assert queue.claim('two', MEMORY, IDENTITY)['state'] == 'finished'
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
    from mappa.postprocessing import decode_events, EventDecodeConfig
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
                    decoded = decode_events(result, EventDecodeConfig(max_events=8))
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
    from mappa.final_branches import FinalBranchCatalogue
    from mappa.postprocessing import EventDecodeConfig, decode_events
    from mappa.mappa_core.core.scheduler.decoding import plan_event_decoding
    case, _, _ = prepare('[CH3:1][OH:2]>>[CH3:1][OH:2]')
    aam = search_aam(case, AAMSearchConfig(seed_count=2, sweep_cuts=False), execution='reused_native')
    catalogue = FinalBranchCatalogue(case).add_aam(aam)
    config = EventDecodeConfig(event_slack=5)
    expected = decode_events(catalogue, config)
    assert expected.complete
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
        assert result['complete']
        patterns.update(p['id'] for p in result['candidates'])
    assert patterns == {candidate.id for candidate in expected.candidates}
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
