"""Branch budgets apply to completed fragments, not temporary embeddings."""
import numpy as np
import pytest

from mappa import build_graph
from mappa.mappa_core.core.search.branch import find_islands
from mappa.mappa_core.core.matching.growth import grow_island
from mappa.mappa_core.core.matching.records import IslandBranchLimitExceeded
from mappa.mappa_core.core.matching import native
from mappa.mappa_core.core.matching.orbits import _nauty_orbits


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
    kwargs = dict(iso_tol=.1, p_orbits=_nauty_orbits(dest, wbo_tol=.1))
    profile = []
    result = grow_island(source, dest, 0, {}, max_branches=1, profile=profile, **kwargs)
    unrestricted = grow_island(source, dest, 0, {}, max_branches=100, **kwargs)
    assert [dict(r) for r in result] == [dict(r) for r in unrestricted]
    assert len(result) == 1
    assert result[0].fragment == frozenset(source)
    assert max(p['max_cands_before'] for p in profile) > 1
    assert dict(result[0]) == ({0: 1, 1: 5} if initial_ambiguity else {0: 0, 1: 1, 2: 5})


def test_completed_fragment_still_obeys_cap(engine):
    source = graph(['Si', 'N'], [(0, 1)])
    dest = target()
    kwargs = dict(iso_tol=.1, p_orbits=_nauty_orbits(dest, wbo_tol=.1))
    assert len(grow_island(source, dest, 0, {}, max_branches=4, **kwargs)) == 4
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
        from mappa.mappa_core.core.search.native import find_islands_native
        search = find_islands_native
    else:
        search = find_islands
    kwargs = dict(iso_tol=.1, p_orbits=_nauty_orbits(dest, wbo_tol=.1))
    full = search(source, dest, [0, 1], max_branches=4, **kwargs)
    capped = search(source, dest, [0, 1], max_branches=3, **kwargs)
    assert len(full.terminals) == 4
    assert 0 < len(capped.terminals) <= 3
    assert any(s.reason == 'capped' and s.stage == 'combined_live_leaves'
               and s.count == 4 and s.limit == 3 for s in capped.stops)


@pytest.mark.parametrize('checkpoint_only', [False, True])
def test_old_intermediate_cap_checkpoints_cannot_resume(tmp_path, checkpoint_only):
    import json
    from mappa import AAMProblem, AAMSearchConfig, search_aam
    from mappa.domain import MolecularEndpoint
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
