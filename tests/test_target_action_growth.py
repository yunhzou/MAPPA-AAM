"""Exact target actions must transport future atoms, not discard their matches."""
from itertools import product

import numpy as np
import pytest

from mappa import FragmentMatchConfig, match_fragment
from mappa.matcher import _SymBlock, _SymCand, _dedup_sym_cands, _extend_sym_cands, _nauty_orbits
from mappa.matcher.state import candidate_from_record
from mappa.search_symmetry import SymmetryWorkspace
from test_growth_compression_exactness import graph, represented, assert_prefixes


def rings():
    cycle = [(i, (i + 1) % 6, 1) for i in range(6)]
    source = graph(['C'] * 6 + ['H'] * 6, cycle + [(i, 6 + i, 1) for i in range(6)])
    target = graph(['C'] * 6 + ['H'] * 12,
                   cycle + [(i, 6 + 2 * i + h, 1) for i in range(6) for h in range(2)])
    return source, target


def ring_mappings():
    # Independent construction: six rotations, two directions, two H per C.
    result = set()
    for offset, direction in product(range(6), (-1, 1)):
        carbons = [(offset + direction * i) % 6 for i in range(6)]
        for choices in product(range(2), repeat=6):
            images = carbons + [6 + 2 * p + h for p, h in zip(carbons, choices)]
            result.add(tuple(enumerate(images)))
    return result


@pytest.mark.parametrize('order', [list(range(12)), [a for i in range(6) for a in (i, 6 + i)]])
def test_ring_every_prefix_preserves_all_768_matches(order):
    source, target = rings()
    ro, po = _nauty_orbits(source, wbo_tol=1), _nauty_orbits(target, wbo_tol=1)
    candidates = [_SymCand({0: 0}, (_SymBlock((0,), tuple(range(6)), False),))]
    fragment = {0}
    for atom in order[1:]:
        candidates = _extend_sym_cands(candidates, fragment, atom, source, target,
                                       {}, 1, None, r_orbits=ro, p_orbits=po)
        fragment.add(atom)
        expected = {tuple((r, p) for r, p in mapping if r in fragment) for mapping in ring_mappings()}
        assert represented(candidates, target, 1) == expected
    assert len(candidates) == 1


@pytest.mark.parametrize('native', ['0', '1'])
def test_public_ring_match_is_one_branch_with_768_assignments(native, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', native)
    source, target = rings()
    result = match_fragment(source, target, seed=0,
                            config=FragmentMatchConfig(iso_tolerance=1, branch_limit=100))
    assert not result.capped
    assert len(result.matches) == 1
    candidate = candidate_from_record(result.matches[0].symmetry)
    assert set(candidate.mapping) == set(source)
    assert represented([candidate], target, 1) == ring_mappings()


def test_target_actions_fix_locked_images():
    _, target = rings()
    candidate = _SymCand({0: 0})
    generators = SymmetryWorkspace(target, 1).conditioned_generators(candidate, ((99, 2), (100, 10)))
    assert all(g[2] == 2 and g[10] == 10 for g in generators)


def test_weak_weight_breaks_apparent_ring_symmetry():
    source, target = rings()
    target.graph['wbo_matrix'][0, 2] = target.graph['wbo_matrix'][2, 0] = .05
    generators = SymmetryWorkspace(target, 1).conditioned_generators(_SymCand({0: 0}), ())
    weights = target.graph['wbo_matrix']
    assert all(np.array_equal(weights, weights[np.ix_(g, g)]) for g in generators)


def test_new_branch_is_not_constrained_by_old_automorph_domain_metadata():
    target = graph(['C'] * 4, [(i, (i + 1) % 4, 1) for i in range(4)])
    first = _SymCand({0: 0, 1: 1}, automorph_blocks=(_SymBlock((0,), (0, 2), False),))
    second = _SymCand({0: 1, 1: 2})
    candidates = _dedup_sym_cands([first, second], target, target)
    assert len(candidates) == 1


def test_asymmetric_future_neighbors_keep_non_equivalent_matches():
    source = graph(['C', 'C', 'O'], [(0, 1, 1), (1, 2, 1)])
    target = graph(['C', 'C', 'C', 'O'], [(0, 1, 1), (1, 2, 1), (2, 3, 1)])
    assert_prefixes(source, target, [0, 1, 2])
