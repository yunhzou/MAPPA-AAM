"""Stage boundaries, transport equivalence and bounded file creation."""
import json
from pathlib import Path

import pytest

from mappa import AAMSearchConfig, search_aam, decode_events
from mappa.execution import WorkPool, MemoryStore, BundleStore, StageResults
from mappa.mappa_core.core.search.graph import AAMSearchGraph
from mappa.mappa_core.utils.io.records import decoded_events_record
from test_reused_aam_api import problem


def normalize_provenance(record, cuts):
    """Translate per-archive local IDs to the combined graph's IDs for comparison."""
    record = json.loads(json.dumps(record))
    offsets = []
    contexts = states = transitions = 0
    for cut in cuts:
        offsets.append((contexts, states, transitions))
        contexts += len(cut.graph.contexts)
        states += len(cut.graph.states)
        transitions += len(cut.graph.transitions)
    for family in record['families']:
        for origin in family['provenance']:
            number = int(origin['archive'].rsplit('/', 1)[1])
            c, s, t = offsets[number]
            origin.update(archive='in_memory', context=origin['context']+c,
                          terminal=origin['terminal']+s,
                          transitions=[edge+t for edge in origin['transitions']])
    return record


@pytest.mark.parametrize('storage', ['ram', 'bundles'])
@pytest.mark.parametrize('workers', [1, 2])
def test_separate_stages_preserve_full_graph_catalogue_and_decoding(tmp_path, storage, workers):
    store = MemoryStore() if storage == 'ram' else BundleStore(tmp_path/'results')
    cases = {'forward': problem(), 'unbalanced': problem(5, 7)}
    config = AAMSearchConfig(seed_count=3, branch_limit=4)
    with WorkPool.local(workers, store=store) as pool:
        assert len(pool.queue.ready_workers()) == workers
        assert pool.queue.snapshot() == {}
        matched = pool.search(cases, config).wait_and_validate()
        matching_workers = {r['result']['worker'] for r in pool.queue.results().values()}
        catalogues = pool.deduplicate(matched).wait_and_validate()
        decoded = pool.decode(catalogues).wait_and_validate()
        assert {r['result']['worker'] for r in pool.queue.results().values()} == matching_workers
        for name, case in cases.items():
            reference = search_aam(case, config, execution='reused_native')
            saved = matched.load(name)
            cuts = [store.load(ref) for ref in saved.cuts]
            assert AAMSearchGraph.combine([cut.graph for cut in cuts]) == reference.graph
            catalogue = catalogues.load(name)
            assert catalogue.search_config == config
            assert catalogue.execution == 'reused_native'
            assert normalize_provenance(catalogue.catalogue.to_record(), cuts) == json.loads(json.dumps(reference.final_catalogue().to_record()))
            actual = dict(decoded.load(name).record)
            assert decoded.load(name).search_config == config
            actual['catalogue'] = normalize_provenance(actual['catalogue'], cuts)
            expected = decoded_events_record(decode_events(reference.final_catalogue()))
            assert json.loads(json.dumps(actual)) == json.loads(json.dumps(expected))
            assert len(catalogue.scope) == len(cuts)
            assert all(scope.seeds == (0, 1, 2) for scope in catalogue.scope)
        if storage == 'bundles':
            decoded.save_manifest(tmp_path/'decoded.json')
            restored = StageResults.open(tmp_path/'decoded.json', store=store)
            assert restored.load('forward').case == 'forward'
            # Inputs, manifests and many seed/cut outputs share worker bundles.
            assert len(list((tmp_path/'results').glob('*.bundle'))) <= workers+1
        else:
            with pytest.raises(ValueError, match='durable manifest'):
                decoded.save_manifest(tmp_path/'decoded.json')


def test_stage_barrier_and_explicit_decoding_input(tmp_path):
    from mappa.mappa_core.core.scheduler.queue import Task
    pool = WorkPool.local(store=MemoryStore())
    with pool:
        pool.queue.add([Task('pending', 'execution_helpers:value', {'sleep': .1})])
        pending = StageResults('matching', pool.store, pool=pool, tasks={'case': ('pending',)})
        with pytest.raises(RuntimeError, match='Finish and validate'):
            pool.deduplicate(pending)
        with pytest.raises(TypeError, match='deduplication'):
            pool.decode(pending)


def test_catalogue_manifest_decodes_in_a_fresh_pool(tmp_path):
    store = BundleStore(tmp_path/'results')
    with WorkPool.local(store=store) as pool:
        matched = pool.search({'case': problem()}, AAMSearchConfig(seed_count=2)).wait_and_validate()
        catalogues = pool.deduplicate(matched).wait_and_validate()
        catalogues.save_manifest(tmp_path/'catalogues.json')
    restored = StageResults.open(tmp_path/'catalogues.json', store=store)
    with WorkPool.local(store=store) as pool:
        decoded = pool.decode(restored).wait_and_validate()
        assert decoded.load('case').record['mapping_space_preserved']


def test_slurm_requires_an_allocation(monkeypatch):
    monkeypatch.delenv('SLURM_JOB_ID', raising=False)
    with pytest.raises(RuntimeError, match='existing Slurm allocation'):
        WorkPool.slurm()


@pytest.mark.parametrize('storage', ['ram', 'bundles'])
def test_one_slurm_startup_serves_all_stages_with_separate_remote_process(tmp_path, monkeypatch, storage):
    import subprocess
    import sys
    import mappa.execution.pool as module
    monkeypatch.setenv('SLURM_JOB_ID', 'test-allocation')
    monkeypatch.setenv('SLURM_JOB_NUM_NODES', '1')
    monkeypatch.setenv('SLURM_CPUS_PER_TASK', '2')
    monkeypatch.setenv('SLURM_MEM_PER_NODE', '4096')
    monkeypatch.setattr(module.socket, 'getfqdn', lambda: '127.0.0.1')
    launch = subprocess.Popen
    commands = []
    def start(command, **options):
        commands.append(command)
        assert command[0] == 'srun'
        # Exercise the actual network worker entry point locally; no real allocation is submitted.
        return launch(command[command.index(sys.executable):], **options)
    monkeypatch.setattr(module.subprocess, 'Popen', start)
    store = MemoryStore() if storage == 'ram' else BundleStore(tmp_path/'results')
    with WorkPool.slurm(store=store) as pool:
        assert len(pool.queue.ready_workers()) == 2
        assert pool.queue.snapshot() == {}
        matched = pool.search({'case': problem()}, AAMSearchConfig(seed_count=2)).wait_and_validate()
        pids = {r['result']['compute_pid'] for r in pool.queue.results().values()}
        catalogues = pool.deduplicate(matched).wait_and_validate()
        decoded = pool.decode(catalogues).wait_and_validate()
        assert decoded.load('case').record['mapping_space_preserved']
        assert {r['result']['compute_pid'] for r in pool.queue.results().values()} == pids
    assert len(commands) == 1


def test_failed_matching_blocks_the_next_stage():
    from mappa.mappa_core.core.scheduler.queue import Task
    with WorkPool.local() as pool:
        pool.queue.add([
            Task('failed', 'execution_helpers:value', {'fail': True}),
            Task('final', 'execution_helpers:value', {}, dependencies=('failed',)),
        ])
        matched = StageResults('matching', pool.store, pool=pool, tasks={'case': ('final',)})
        with pytest.raises(RuntimeError, match='did not complete'):
            matched.wait_and_validate()
        with pytest.raises(RuntimeError, match='Finish and validate'):
            pool.deduplicate(matched)
        assert len(pool.queue.results()) == 2


def test_archived_matching_retry_reuses_completed_cuts_and_attempts(tmp_path):
    from mappa.mappa_core.core.scheduler import TaskQueue
    from mappa.mappa_core.core.scheduler.runtime import Workspace
    from mappa.mappa_core.core.scheduler.search import plan_search, search_complete_cut
    store = BundleStore(tmp_path/'results')
    config = AAMSearchConfig(seed_count=1, branch_limit=10)
    tasks, keys = plan_search('original/case', problem(), config, store)
    assert len(tasks) > 1
    queue = TaskQueue(tmp_path/'old.sqlite', 'old-code')
    queue.add(tasks)
    workspace = Workspace(str(tmp_path), store=store)
    first = None
    for _ in tasks:
        claimed = queue.claim('old-worker', 2**30, 'old-code')
        task = claimed['task']
        if first is None:
            first = store.put(search_complete_cut(task.payload, workspace))
            queue.finish(task.key, claimed['token'], 'complete', dict(blob=first))
        else:
            queue.finish(task.key, claimed['token'], 'memory_limit', dict(memory_guard_bytes=2**30))
    audit = queue.archive_case(keys[0], store)
    before = store.load(audit)
    queue.close()
    with WorkPool.local(workers=2, store=store) as pool:
        matched = pool.retry_matching('case', audit).wait_and_validate()
        assert matched.load('case').cuts[0] == first
        assert pool.queue.db.execute('SELECT count(*) FROM attempts WHERE task=?', (keys[0],)).fetchone()[0] == 1
        assert pool.queue.db.execute('SELECT count(*) FROM attempts').fetchone()[0] == 2*len(tasks)-1
        assert pool.queue.db.execute('SELECT spent FROM cases').fetchone()[0] >= before['case'][3]
        catalogue = pool.deduplicate(matched).wait_and_validate()
        decoded = pool.decode(catalogue).wait_and_validate()
        assert decoded.load('case').record['mapping_space_preserved']
        again = store.load(pool.queue.archive_case(keys[0], store))
        old_tokens = {row[0] for row in before['attempts']}
        assert old_tokens <= {row[0] for row in again['attempts']}

    reference = search_aam(problem(), config, execution='reused_native')
    cuts = [store.load(ref) for ref in matched.load('case').cuts]
    assert AAMSearchGraph.combine([cut.graph for cut in cuts]) == reference.graph


def test_handoff_rejects_missing_cuts_and_different_memory_store(tmp_path):
    from dataclasses import replace
    from mappa.mappa_core.core.scheduler.decoding import reduce_matching
    from mappa.mappa_core.core.scheduler.runtime import Workspace
    store = MemoryStore()
    with WorkPool.local(store=store) as pool:
        matched = pool.search({'case': problem()}, AAMSearchConfig(seed_count=1)).wait_and_validate()
        saved = matched.load('case')
        malformed = store.put(replace(saved, cuts=saved.cuts[:-1]))
        with pytest.raises(ValueError, match='Missing matching cuts'):
            reduce_matching({'matching': malformed}, Workspace(tmp_path, store=store))
    other = WorkPool.local()
    with pytest.raises(ValueError, match='same MemoryStore'):
        other.deduplicate(matched)
    with pytest.raises(ValueError, match='same storage mode'):
        WorkPool.local(store=BundleStore(tmp_path/'bundles')).deduplicate(matched)


def test_pool_startup_timeout_closes_waiting_supervisor(monkeypatch):
    import time
    import mappa.execution.pool as module
    def never_ready(queue, **options):
        while queue.keep_open:
            time.sleep(.001)
        return {}
    monkeypatch.setattr(module, 'run', never_ready)
    pool = WorkPool.local(startup_seconds=.05)
    with pytest.raises(TimeoutError, match='0/1 workers'):
        with pool:
            pytest.fail('Unready pool returned')
    assert pool._runner.done()
