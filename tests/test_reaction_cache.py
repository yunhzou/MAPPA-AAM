"""Snapshots preserve scientific results through concurrency, restart and I/O failure."""
from pathlib import Path
import pytest

from mappa import AAMSearchConfig
from mappa.execution import BundleStore, WorkPool
from mappa.mappa_core.core.scheduler.cache import ReactionCaches, WorkerCache, cache_identity, read_section
from mappa.mappa_core.core.search.session import SearchSession
from test_reused_aam_api import problem


def session():
    return SearchSession(problem(), AAMSearchConfig(seed_count=2, branch_limit=4, growth_policy='existence'))


def execute(store, identity, search, worker, cut=()):
    worker.begin(store.describe(identity))
    worker.use(search.repair)
    graph, _ = search.search_cut(cut)
    proposal = worker.proposal()
    if proposal is not None:
        token = str(search.repair.stats()['calls']) + str(id(worker))
        plan = store.reserve(identity, proposal['version'], token)
        publication = worker.save(plan)
        if plan is not None:
            assert store.finish(token, publication)
    return graph


def test_latest_reload_and_immutable_baseline_after_failure(tmp_path, monkeypatch):
    store = ReactionCaches(tmp_path)
    identity = ('reaction', '0')
    first, worker = session(), WorkerCache()
    expected = execute(store, identity, first, worker)
    baseline = Path(store.describe(identity)['baseline'])
    first_bytes = baseline.read_bytes()
    # Its own successful publication needs neither a load nor another save.
    import mappa.mappa_core.core.scheduler.cache as module
    original = module.read_section
    monkeypatch.setattr(module, 'read_section', lambda *args: (_ for _ in ()).throw(AssertionError('redundant load')))
    assert execute(store, identity, first, worker) == expected
    assert worker.metrics['loaded_bytes'] == worker.metrics['saved_bytes'] == 0
    monkeypatch.setattr(module, 'read_section', original)
    execute(store, identity, first, worker, ((0,1),))
    assert baseline.read_bytes() == first_bytes
    latest = Path(store.describe(identity)['latest'])
    latest.write_bytes(b'truncated')
    store.close()
    restarted = ReactionCaches(tmp_path)
    recovered, loader = session(), WorkerCache()
    assert execute(restarted, identity, recovered, loader) == expected
    assert loader.metrics['baseline_loads'] == 1
    assert loader.metrics['errors']
    assert loader.metrics['fragment_hits'] == loader.metrics['fragment_calls']
    assert baseline.read_bytes() == first_bytes
    assert len(list(tmp_path.glob('*.cache'))) == 2
    assert not list(tmp_path.glob('*.pending'))


def test_competing_publishers_and_open_reader_keep_complete_snapshot(tmp_path):
    store = ReactionCaches(tmp_path)
    identity = ('reaction','0')
    search, worker = session(), WorkerCache()
    execute(store, identity, search, worker)
    descriptor = store.describe(identity)
    old_reader = open(descriptor['latest'], 'rb')
    old_bytes = old_reader.read()
    old_reader.seek(0)
    worker.begin(descriptor); worker.use(search.repair)
    search.search_cut(((0,1),))
    proposal = worker.proposal()
    plan = store.reserve(identity, proposal['version'], 'winner')
    assert store.reserve(identity, proposal['version'], 'loser') is None
    assert store.finish('winner', worker.save(plan))
    assert old_reader.read() == old_bytes
    old_reader.close()
    assert store.reserve(identity, proposal['version'], 'stale') is None
    assert len(list(tmp_path.glob('*.cache'))) == 2


def test_directions_share_two_files_but_never_matching_entries(tmp_path):
    store = ReactionCaches(tmp_path)
    forward, reverse = ('reaction','0'), ('reaction','1')
    execute(store, forward, session(), WorkerCache())
    baseline = store.describe(forward)['baseline']
    original, _ = read_section(baseline, *forward)
    execute(store, reverse, session(), WorkerCache())
    assert read_section(baseline, *forward)[0] == original
    assert read_section(baseline, *reverse)[0] is not None
    assert len(list(tmp_path.glob('*.cache'))) == 2


def test_input_config_and_direction_identity():
    from mappa import AAMProblem
    p = problem()
    config = AAMSearchConfig(growth_policy='existence')
    a = cache_identity(p, config, 'code1')
    b = cache_identity(AAMProblem(p.product,p.reactant), config, 'code1')
    assert a[0] == b[0] and a[1] != b[1]
    assert cache_identity(p, config, 'code2')[0] != a[0]
    assert cache_identity(p, AAMSearchConfig(), 'code1')[0] != a[0]


def test_fresh_pool_loads_previous_reaction_snapshot(tmp_path):
    store = BundleStore(tmp_path/'bundles')
    config = AAMSearchConfig(seed_count=2, branch_limit=4, growth_policy='existence')
    saved = []
    for _ in range(2):
        with WorkPool.local(2, store=store) as pool:
            result = pool.search({'reaction':problem()}, config).wait_and_validate().load('reaction')
            saved.append([store.load(ref).graph for ref in result.cuts])
            metrics = [r['result'].get('cache', {}) for r in pool.queue.results().values()]
            assert not [error for row in metrics for error in row.get('errors', [])]
    assert saved[0] == saved[1]
    assert sum(row.get('loaded_bytes', 0) for row in metrics) > 0
    assert sum(row.get('fragment_hits', 0) for row in metrics) > 0
    assert len(list((tmp_path/'bundles'/'cache').glob('*.cache'))) == 2


def test_missing_latest_uses_baseline_and_both_bad_are_explicit(tmp_path):
    store=ReactionCaches(tmp_path)
    identity=('reaction','0')
    expected=execute(store,identity,session(),WorkerCache())
    descriptor=store.describe(identity)
    Path(descriptor['latest']).unlink()
    loader=WorkerCache()
    assert execute(store,identity,session(),loader)==expected
    assert loader.metrics['baseline_loads']==1
    for name in ('latest','baseline'):
        Path(descriptor[name]).write_bytes(b'bad')
    store.close()
    restarted=ReactionCaches(tmp_path)
    with pytest.raises(RuntimeError,match='Neither reaction cache'):
        execute(restarted,identity,session(),WorkerCache())


def test_failed_publication_keeps_baseline_and_releases_writer(tmp_path):
    store=ReactionCaches(tmp_path)
    identity=('reaction','0')
    search,worker=session(),WorkerCache()
    execute(store,identity,search,worker)
    descriptor=store.describe(identity)
    original=Path(descriptor['baseline']).read_bytes()
    plan=store.reserve(identity,descriptor['version'],'failed')
    Path(plan['staging']).write_bytes(b'partial')
    assert not store.finish('failed')
    assert Path(descriptor['baseline']).read_bytes()==original
    assert store.reserve(identity,descriptor['version'],'next') is not None
    store.finish('next')
    assert not list(tmp_path.glob('*.pending'))
    with pytest.raises(RuntimeError,match='active coordinator'):
        ReactionCaches(tmp_path)
