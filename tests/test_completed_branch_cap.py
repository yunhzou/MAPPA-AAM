"""Branch budgets apply to completed fragments, not temporary embeddings."""
import numpy as np
import pytest

from mappa import build_graph
from mappa.mappa_core.core.search.branching.construction import find_islands
from mappa.mappa_core.core.matching.growth.engine import grow_island
from mappa.mappa_core.core.matching.representation.records import IslandBranchLimitExceeded
from mappa.mappa_core.core.matching.growth import native


def prune_dag(edges, frontier, cap):
    from mappa.mappa_core.core.search.graph import (
        FragmentTransition, SearchState, prune_frontier)

    states = [SearchState(i, 0, (), (), ())
              for i in range(1 + max(max(edge) for edge in edges))]
    transitions = [FragmentTransition(i, source, target, None, (1, 0), None)
                   for i, (source, target) in enumerate(edges)]
    leaves, kept, cuts = prune_frontier(states, transitions, (0,), frontier, cap)
    return set(leaves), {edges[i] for i in kept}, list(cuts)


def test_recursive_cap_removes_the_largest_ancestor_subtree():
    # Root children B=(4+4) and C=(3+3) each fit cap10. The root must remove
    # all of B and keep C's6; dropping just one four-leaf group is not recursive.
    edges = [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)]
    edges += [(3, leaf) for leaf in range(7, 11)]
    edges += [(4, leaf) for leaf in range(11, 15)]
    edges += [(5, leaf) for leaf in range(15, 18)]
    edges += [(6, leaf) for leaf in range(18, 21)]
    leaves, kept, cuts = prune_dag(edges, tuple(range(7, 21)), 10)
    assert leaves == set(range(15, 21))
    assert kept == {(0, 2), (2, 5), (2, 6),
                    *((5, leaf) for leaf in range(15, 18)),
                    *((6, leaf) for leaf in range(18, 21))}
    assert cuts == [(0, 0, 14, 8)]


def test_recursive_cap_recounts_parent_after_pruning_below():
    # C's4+1+1 leaves shrink to2 before root A counts B's2 plus C's2.
    edges = [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6), (2, 7)]
    edges += [(5, leaf) for leaf in range(8, 12)] + [(6, 12), (7, 13)]
    leaves, kept, cuts = prune_dag(edges, (3, 4, 8, 9, 10, 11, 12, 13), 5)
    assert leaves == {3, 4, 12, 13}
    assert kept == {(0, 1), (0, 2), (1, 3), (1, 4),
                    (2, 6), (2, 7), (6, 12), (7, 13)}
    assert cuts == [(2, 4, 6, 4)]


def test_recursive_cap_repeats_removal_and_drops_later_equal_children():
    edges = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5),
             (2, 6), (2, 7), (3, 8), (3, 9)]
    leaves, kept, cuts = prune_dag(edges, (4, 5, 6, 7, 8, 9), 3)
    assert leaves == {4, 5}
    assert kept == {(0, 1), (1, 4), (1, 5)}
    assert cuts == [(0, 2, 6, 2), (0, 1, 4, 2)]


def test_recursive_cap_counts_a_shared_leaf_once_at_exact_limit():
    edges = [(0, 1), (0, 2), (1, 3), (1, 4), (2, 4), (2, 5)]
    leaves, kept, cuts = prune_dag(edges, (3, 4, 5), 3)
    assert leaves == {3, 4, 5}
    assert kept == set(edges)
    assert cuts == []


def test_recursive_cap_removes_rejected_history_but_keeps_shared_leaf():
    edges = [(0, 1), (0, 2), (1, 3), (1, 4), (2, 4), (2, 5), (2, 6)]
    leaves, kept, cuts = prune_dag(edges, (3, 4, 5, 6), 3)
    assert leaves == {3, 4}
    assert kept == {(0, 1), (1, 3), (1, 4)}
    # The removed child contains three leaves, but shared leaf4 stays via child1.
    assert cuts == [(0, 1, 4, 3)]


def test_recursive_cap_handles_shared_states_created_before_their_parents():
    # Hash-consing can reuse leaf2 after parent4 or parent3 was created. Node
    # numbers are not a traversal order; the root still removes the larger child.
    edges = [(0, 4), (0, 3), (4, 1), (4, 2), (3, 2), (3, 5), (3, 6)]
    leaves, kept, cuts = prune_dag(edges, (1, 2, 5, 6), 3)
    assert leaves == {1, 2}
    assert kept == {(0, 4), (4, 1), (4, 2)}
    assert cuts == [(0, 1, 4, 3)]


def graph(elements, edges):
    wbo = np.zeros((len(elements), len(elements)))
    for a, b in edges:
        wbo[a, b] = wbo[b, a] = 1.
    return build_graph(elements, wbo, bond_cut=.2)


@pytest.fixture(params=['python', 'native'])
def engine(request, monkeypatch):
    if request.param == 'native' and not native.built():
        pytest.skip('native growth engine not built')
    monkeypatch.setenv('MAPPA_NATIVE', '1' if request.param == 'native' else '0')
    return request.param


def target():
    # Four distinct N placements: only one can extend to O.
    return graph(['Si', 'N', 'N', 'N', 'N', 'O', 'F', 'Cl', 'Br'],
                 [(0, n) for n in range(1, 5)] + [(n, n+4) for n in range(1, 5)])


@pytest.mark.parametrize('initial_ambiguity', [True, False])
def test_temporary_ambiguity_can_resolve_below_cap(engine, initial_ambiguity):
    source = (graph(['N', 'O'], [(0, 1)]) if initial_ambiguity else
              graph(['Si', 'N', 'O'], [(0, 1), (1, 2)]))
    dest = target()
    kwargs = dict(iso_tol=.1)
    profile = []
    result = grow_island(source, dest, 0, {}, max_branches=1, profile=profile, **kwargs)
    unrestricted = grow_island(source, dest, 0, {}, max_branches=100, **kwargs)
    assert [r.key for r in result.matches] == [r.key for r in unrestricted.matches]
    assert result.deferred_edges == unrestricted.deferred_edges
    assert len(result.matches) == 1
    assert result.matches[0].source_atoms == tuple(source)
    candidates = max(p['max_cands_before'] for p in profile)
    assert candidates == 1
    assert dict(result.matches[0].witness) == ({0: 1, 1: 5} if initial_ambiguity else {0: 0, 1: 1, 2: 5})


def test_completed_fragment_still_obeys_cap(engine):
    source = graph(['Si', 'N'], [(0, 1)])
    dest = target()
    kwargs = dict(iso_tol=.1)
    assert len(grow_island(source, dest, 0, {}, max_branches=4, **kwargs).matches) == 4
    with pytest.raises(IslandBranchLimitExceeded) as exc:
        grow_island(source, dest, 0, {}, max_branches=3, **kwargs)
    assert (exc.value.count, exc.value.limit) == (4, 3)


@pytest.mark.smoke
def test_frontier_limit_combines_children_from_different_parents(engine):
    # Two distinct C placements, each with two distinct O placements. Each
    # completed fragment fits cap 3; their combined four leaves do not.
    source = graph(['C', 'O'], [])
    dest = graph(['C', 'C', 'O', 'O', 'N', 'F', 'Cl', 'Br'],
                 [(0, 4), (1, 5), (2, 6), (3, 7)])
    if engine == 'native':
        from mappa.mappa_core.core.search.branching.native import find_islands_native
        search = find_islands_native
    else:
        search = find_islands
    kwargs = dict(iso_tol=.1)
    full = search(source, dest, [0, 1], max_branches=4, **kwargs)
    capped = search(source, dest, [0, 1], max_branches=3, **kwargs)
    assert len(full.terminals) == 4
    assert [path.mapping for path in capped.paths()] == [
        {0: 0, 1: 2}, {0: 0, 1: 3}]
    cap, = [stop for stop in capped.stops if stop.reason == 'capped']
    assert (cap.state, capped.stop_context(cap), cap.seed, cap.step, cap.stage,
            cap.count, cap.limit, cap.occurrences) == (
                None, 0, 1, (1, 1), 'combined_live_leaves', 4, 3, 1)
    assert len(capped.transitions) == 3


@pytest.mark.parametrize('oversized_first', [False, True])
@pytest.mark.parametrize('limit', [3, 4])
def test_capped_child_preserves_its_smaller_sibling(engine, oversized_first, limit):
    # The first C–N fragment chooses one of two N atoms. The next O–N
    # fragment uses the remaining N: one choice has four placements, the
    # other has two. Discarding the four must leave both valid sibling leaves.
    source = graph(['C', 'N', 'O', 'N'], [(0, 1), (2, 3)])
    small, large = (1, 2) if oversized_first else (2, 1)
    # Different external attachments distinguish O placements globally without
    # changing the accepted C–N and O–N matching constraints.
    target = graph(['C', 'N', 'N'] + ['O'] * 6 + ['F', 'Cl', 'Br', 'I', 'P', 'S'],
                   [(0, 1), (0, 2), (small, 3), (small, 4)]
                   + [(large, oxygen) for oxygen in range(5, 9)]
                   + [(oxygen, oxygen+6) for oxygen in range(3, 9)])
    if engine == 'native':
        from mappa.mappa_core.core.search.branching.native import find_islands_native
        search = find_islands_native
    else:
        search = find_islands
    full = search(source, target, [0, 1, 2, 3], iso_tol=.1, max_branches=6)
    capped = search(source, target, [0, 1, 2, 3], iso_tol=.1, max_branches=limit)
    expected = [{0: 0, 1: large, 2: oxygen, 3: small} for oxygen in (3, 4)]
    assert len(full.terminals) == 6
    assert [path.mapping for path in capped.paths()] == expected
    assert all(mapping in [path.mapping for path in full.paths()] for mapping in expected)
    cap, = [stop for stop in capped.stops if stop.reason == 'capped']
    if limit == 3:
        assert (cap.stage, cap.seed, cap.count, cap.limit, cap.occurrences) == (
            'fragment_growth', 2, 4, 3, 1)
    else:
        assert (cap.stage, cap.seed, cap.count, cap.limit, cap.occurrences) == (
            'combined_live_leaves', 2, 6, 4, 1)


@pytest.mark.parametrize('larger_first', [False, True])
def test_three_fragment_search_prunes_the_whole_larger_ancestor(engine, larger_first):
    source = graph(['C', 'N', 'O', 'N', 'F'], [(0, 1), (2, 3)])
    small, large = (1, 2) if larger_first else (2, 1)
    target = graph(['C', 'N', 'N'] + ['O'] * 6 + ['F', 'Cl', 'Br', 'I', 'P', 'S', 'F'],
                   [(0, 1), (0, 2), (small, 3), (small, 4), (0, 15)]
                   + [(large, oxygen) for oxygen in range(5, 9)]
                   + [(oxygen, oxygen+6) for oxygen in range(3, 9)])
    if engine == 'native':
        from mappa.mappa_core.core.search.branching.native import find_islands_native
        search = find_islands_native
    else:
        search = find_islands
    # C–N gives two choices. O–N expands them to4 and2, then the independent F
    # doubles those subtrees to8 and4. Cap10 must remove the entire8-leaf ancestor.
    full = search(source, target, list(source), iso_tol=.1, max_branches=12)
    capped = search(source, target, list(source), iso_tol=.1, max_branches=10)
    expected = {tuple(sorted({0: 0, 1: large, 2: oxygen, 3: small, 4: fluorine}.items()))
                for oxygen in (3, 4) for fluorine in (9, 15)}
    assert len(full.terminals) == 12
    assert {tuple(sorted(path.mapping.items())) for path in capped.paths()} == expected
    assert all(len(path.fragments) == 3 for path in capped.paths())
    assert expected <= {tuple(sorted(path.mapping.items())) for path in full.paths()}
    cap, = [stop for stop in capped.stops if stop.reason == 'capped']
    assert (cap.stage, cap.seed, cap.count, cap.limit, cap.occurrences) == (
        'combined_live_leaves', 4, 12, 10, 1)


@pytest.mark.parametrize('checkpoint_only', [False, True])
def test_old_intermediate_cap_checkpoints_cannot_resume(tmp_path, checkpoint_only):
    import json
    from mappa import AAMProblem, AAMSearchConfig, search_aam
    from mappa.mappa_core.domain import MolecularEndpoint
    from mappa.mappa_core.core.search.planning import checkpoint_manifest
    from mappa.mappa_core.core.search.api import search_aam_checkpoints

    endpoint = MolecularEndpoint(('C',), np.zeros((1, 3)), np.zeros((1, 1)))
    problem = AAMProblem(endpoint, endpoint)
    config = AAMSearchConfig(sweep_cuts=False)
    old = checkpoint_manifest(problem, config)
    old.pop('branch_cap_policy')
    old['schema'] = 'rxn_core.aam_checkpoints/v2'
    (tmp_path/'manifest.json').write_text(json.dumps(old))
    run = search_aam_checkpoints if checkpoint_only else search_aam
    with pytest.raises(ValueError, match='configuration differs'):
        run(problem, config, intermediate_dir=tmp_path, resume=True)
