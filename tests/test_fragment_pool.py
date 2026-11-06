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


@pytest.mark.parametrize('elements,target_size', [(('C', 'C', 'C'), 2), (('C', 'Xe'), 1)])
def test_exhausted_element_supply_rejects_without_cache_or_matching(elements, target_size):
    source = graph(len(elements), [(r, r+1) for r in range(len(elements)-1)], elements)
    target = graph(target_size, [(r, r+1) for r in range(target_size-1)])
    pool = FragmentPool(source, target)
    result = pool.engine.grow([], 0, [-1]*len(elements), .2, 1., 1, 100, None, [], False)
    assert result['profile']['final_fragment_size'] == target_size
    assert result['profile']['deferred'] == 1
    assert result['profile']['feasibility_searches'] == pool.stats()['existence_queries'] == 0


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
    # Changing policy back must restore the indexed proof, including after load.
    restored = FragmentPool(source, target)
    restored.loads(pool.dumps())
    assert not restored.engine.exists([0, 1], [(0, 1)], [-1]*3, .5)
    assert restored.engine.exists([0, 1], [(0, 1)], [-1]*3, 1.)
    assert restored.stats()['fragment_hits'] == 1


def test_containment_across_multiple_bitset_blocks():
    source = graph(140, [(a, a+1) for a in range(0, 140, 2)])
    pool = FragmentPool(source, source)
    for a in range(0, 140, 2):
        fill(pool, a)
    assert pool.stats()['entries'] == 70
    before = pool.stats()['fragment_hits']
    for a in (0, 126, 128, 138):
        assert pool.engine.exists([a, a+1], [(a, a+1)], [-1]*140, 1.)
    assert pool.stats()['fragment_hits'] == before+4
    # Atom containment is across one stored fragment, not the union of records.
    before = pool.stats()['fragment_hits']
    assert pool.engine.exists([0, 138], [], [-1]*140, 1.)
    assert pool.stats()['fragment_hits'] == before
    assert pool.stats()['peak_resident_bytes'] <= pool.stats()['cache_budget_bytes']


def test_sparse_block_index_grows_past_one_word_and_restores():
    star = graph(14, [(0, r) for r in range(1, 14)], ['O']+['H']*13)
    pool = FragmentPool(star, star)
    # The last two source sets introduce atom 13 in block 64 (index word 1).
    for mask in range(1, 4098):
        cut = [(0, r) for r in range(1, 14) if not mask & (1 << (r-1))]
        fill(pool, cut=cut)
    assert pool.stats()['entries'] == 4097
    restored = FragmentPool(star, star)
    restored.loads(pool.dumps())
    assert restored.engine.exists([0, 1, 13], [(0, 1), (0, 13)], [-1]*14, 1.)
    assert restored.stats()['fragment_hits'] == 1
    # These atoms occur in separate records, which cannot be combined as a proof.
    assert restored.engine.exists([0, 2, 13], [(0, 2), (0, 13)], [-1]*14, 1.)
    assert restored.stats()['fragment_hits'] == 1
    assert restored.stats()['peak_resident_bytes'] <= restored.stats()['cache_budget_bytes']


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
