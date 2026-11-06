"""Known fragments prove smaller queries without choosing growth trajectories."""
import numpy as np
import pytest

from growth_cases import graph, raw
from mappa import AAMProblem, AAMSearchConfig, MolecularEndpoint
from mappa.mappa_core.core.matching.reuse import FragmentPool
from mappa.mappa_core.core.search.session import SearchSession
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue


def fill(pool, seed=0, cut=()):
    return raw(pool.source, pool.target, seed, [-1]*len(pool.source), cut, pool=pool)


def test_completed_fragment_proves_different_starting_atom_and_subset():
    chain = graph(6, [(i, i+1) for i in range(5)])
    pool = FragmentPool(chain, chain)
    fill(pool)
    before = pool.stats()
    # Both queries use the same larger stored fragment, without growing it again.
    for atoms, edges in [([2, 3], [(2, 3)]), ([1, 2, 3], [(1, 2), (2, 3)])]:
        assert pool.engine.exists(atoms, edges, [-1]*6, 1.)
    assert pool.stats()['fragment_hits'] == before['fragment_hits'] + 2
    assert pool.stats()['entries'] == 1


def test_previous_other_fragment_lock_does_not_make_a_branch_cache_key():
    source = graph(6, [(0, 1), (2, 3), (3, 4), (4, 5)])
    pool = FragmentPool(source, source)
    raw(source, source, 2, [0, 1, -1, -1, -1, -1], (), pool=pool,
        islands=[(0, 1), (1, 1)])
    before = pool.stats()['fragment_hits']
    assert pool.engine.exists([3, 4], [(3, 4)], [1, 0, -1, -1, -1, -1], 1.)
    assert pool.stats()['fragment_hits'] == before+1


def test_available_targets_locks_and_new_required_edges_are_checked():
    source = graph(4, [(0, 1), (1, 2), (0, 2)])
    target = graph(4, [(0, 1), (1, 2)])
    pool = FragmentPool(source, target)
    fill(pool)
    before = pool.stats()['fragment_hits']
    # Same atom set, but the previous path does not prove a triangle in P.
    assert not pool.engine.exists([0, 1, 2], [(0, 1), (1, 2), (0, 2)], [-1]*4, 1.)
    # A target required by the known witness is now owned outside the fragment.
    assert not pool.engine.exists([0, 1, 2], [(0, 1), (1, 2)], [-1, -1, -1, 1], 1.)
    # A new anchor places the centre at an end, so the stored witness is invalid.
    assert not pool.engine.exists([0, 1, 2], [(0, 1), (1, 2)], [-1, 0, -1, -1], 1.)
    assert pool.stats()['fragment_hits'] == before


def test_tolerance_is_checked_on_the_proposed_edges():
    source = graph(3, [(0, 1, 1.), (1, 2, 1.)])
    target = graph(3, [(0, 1, 2.), (1, 2, 2.)])
    pool = FragmentPool(source, target)
    fill(pool)
    before = pool.stats()['fragment_hits']
    assert pool.engine.exists([0, 1], [(0, 1)], [-1]*3, 1.)
    assert not pool.engine.exists([0, 1], [(0, 1)], [-1]*3, .5)
    assert pool.stats()['fragment_hits'] == before+1


def test_snapshot_transfers_unique_fragments_and_preserves_local_knowledge():
    source = graph(6, [(0, 1), (1, 2), (3, 4), (4, 5)])
    first, second = FragmentPool(source, source), FragmentPool(source, source)
    fill(first, 0)
    fill(second, 3)
    assert second.loads(first.dumps())
    before = second.stats()['fragment_hits']
    assert second.engine.exists([0, 1], [(0, 1)], [-1]*6, 1.)
    assert second.engine.exists([3, 4], [(3, 4)], [-1]*6, 1.)
    assert second.stats()['fragment_hits'] == before+2
    state = second.dumps()
    with pytest.raises(ValueError):
        second.loads(state[:-2])
    assert second.dumps() == state


def test_budget_changes_reuse_not_answers():
    source = graph(6, [(i, i+1) for i in range(5)])
    pool = FragmentPool(source, source, cache_bytes=1)
    assert fill(pool) == raw(source, source, 0, [-1]*6, ())
    assert pool.stats()['entries'] == 0
    assert pool.engine.exists([1, 2], [(1, 2)], [-1]*6, 1.)
    assert pool.stats()['fragment_hits'] == 0


def test_eviction_removes_witness_index_and_snapshot_obeys_receiver_budget():
    source = graph(6, [(0, 1), (1, 2), (3, 4), (4, 5)])
    complete = FragmentPool(source, source)
    fill(complete, 0)
    budget = complete.stats()['resident_bytes']
    fill(complete, 3)
    limited = FragmentPool(source, source, cache_bytes=budget)
    limited.loads(complete.dumps())
    assert limited.stats()['entries'] == limited.stats()['evictions'] == 1
    assert limited.stats()['peak_resident_bytes'] <= budget
    assert limited.engine.exists([0, 1], [(0, 1)], [-1]*6, 1.)
    assert limited.stats()['fragment_hits'] == 0  # Eviction is not a negative answer.
    assert limited.engine.exists([3, 4], [(3, 4)], [-1]*6, 1.)
    assert limited.stats()['fragment_hits'] == 1


def test_cut_trial_streams_are_distinct_and_batching_does_not_change_search():
    weights = np.zeros((6, 6))
    for a, b in [(i, (i+1)%6) for i in range(6)]:
        weights[a, b] = weights[b, a] = 1.
    endpoint = MolecularEndpoint(('C',)*6, np.zeros((6, 3)), weights)
    problem = AAMProblem(endpoint, endpoint)
    config = AAMSearchConfig(seed_count=3, branch_limit=100)
    from mappa.mappa_core.core.search.graph import AAMSearchGraph
    session = SearchSession(problem, config)
    whole, _ = session.search_cut(())
    split = SearchSession(problem, config)
    parts = [split.search_cut((), seed_indices=(i,))[0] for i in range(3)]
    assert whole == AAMSearchGraph.combine(parts)
    cut, _ = session.search_cut(((0, 1),))
    assert len({c.random_seed for c in whole.contexts+cut.contexts}) == 6
    # Reduction uses the matching inventory's records, not reconstructed copies.
    catalogue = FinalBranchCatalogue(problem).add_graph(whole)
    assert all(fragment is whole.fragments[fragment.atom_sets]
               for fragment in catalogue._fragments.values())


def test_pool_never_replays_a_completed_growth_trajectory():
    source = graph(6, [(0, 1), (1, 2), (2, 3), (3, 0), (2, 4), (4, 5)])
    target = graph(6, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)])
    pool = FragmentPool(source, target)
    fill(pool)
    for seed in range(6):
        for cut in [(), ((0, 1),), ((2, 4),)]:
            assert raw(source, target, seed, [-1]*6, cut, pool=pool) == raw(
                source, target, seed, [-1]*6, cut)
