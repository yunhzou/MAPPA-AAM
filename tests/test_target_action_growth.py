"""Local relaxed relations preserve correlated branch choices and target occupancy."""
from itertools import product, permutations

import numpy as np
import pytest

from mappa import FragmentMatchConfig, match_fragment
from test_growth_compression_exactness import graph, represented, discovered, assert_completed_constraints, complete


from growth_cases import rings, ring_mappings


@pytest.mark.parametrize('order', [list(range(12)), [a for i in range(6) for a in (i, 6 + i)]])
def test_ring_constraints_preserve_all_768_matches(order):
    source, target = rings()
    fragment = {0}
    bonds = []
    for atom in order[1:]:
        bonds.append((min(set(source.neighbors(atom)) & fragment), atom))
        fragment.add(atom)
        candidates = complete(source, target, fragment, bonds, 1)
        expected = {tuple((r, p) for r, p in mapping if r in fragment) for mapping in ring_mappings()}
        assert discovered(candidates) == expected
    assert len(candidates) == 64


@pytest.mark.parametrize('native', ['0', '1'])
@pytest.mark.smoke
def test_public_ring_match_has_64_atom_sets_with_768_assignments(native, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', native)
    source, target = rings()
    result = match_fragment(source, target, seed=0,
                            config=FragmentMatchConfig(iso_tolerance=1, branch_limit=100))
    assert not result.capped
    assert len(result.matches) == 64
    assert all(set(fragment.source_atoms) == set(source) for fragment in result.matches)
    assert represented(result.matches) == ring_mappings()


def test_target_local_actions_fix_fragment_anchors():
    source, target = rings()
    fragments = complete(source, target, {0, 1, 6, 7}, [(0, 1), (0, 6), (1, 7)],
                         locked={0: 0})
    assert fragments
    for fragment in fragments:
        relation = fragment.target_self_automorph
        anchor = relation.atoms.index(0)
        assert all(action[anchor] == anchor for action in relation.generators)


def test_asymmetric_future_neighbors_keep_non_equivalent_matches():
    source = graph(['C', 'C', 'O'], [(0, 1, 1), (1, 2, 1)])
    target = graph(['C', 'C', 'C', 'O'], [(0, 1, 1), (1, 2, 1), (2, 3, 1)])
    assert_completed_constraints(source, target, [0, 1, 2])


def test_decoding_keeps_fragment_atom_sets_and_later_occupancy():
    from mappa import AAMProblem, query_path
    from mappa.mappa_core.core.search.branching.construction import find_islands
    from mappa.mappa_workflow.mapping_analytics.family import project_path_assignments
    from test_family_query import endpoint

    source = graph(['C'] * 4, [(0, 1, 1), (1, 2, 1), (2, 0, 1)])
    target = graph(['C'] * 4, [(i, (i + 1) % 4, 1) for i in range(4)])
    problem = AAMProblem(endpoint('CCCC', source.edges), endpoint('CCCC', target.edges))
    found = find_islands(source, target, [0, 1, 2, 3], iso_tol=.2, max_branches=100)
    paths = list(found.paths())
    assert len(paths) == 4
    union = set()
    for path in paths:
        first, last = path.fragments
        expected = {tuple(sorted((*mapping, *last.witness))) for mapping in represented((first,))}
        union.update(expected)
        for images in permutations(range(4)):
            status, _ = query_path(path, problem, dict(enumerate(images)),
                                   source_atoms=tuple(range(4)), timeout_ms=2000)
            assert status == ('recovered' if tuple(enumerate(images)) in expected else 'not_recovered')
        projected = {tuple(mapping.items()) for mapping in project_path_assignments(
            path, problem, (0, 3), assignment_limit=100)}
        assert projected == {tuple((r, p) for r, p in mapping if r in (0, 3)) for mapping in expected}
    assert union == {tuple(enumerate(images)) for images in permutations(range(4))}


def test_locked_ring_then_free_hydrogen_preserves_all_joint_choices():
    import z3
    from mappa import AAMProblem
    from mappa.mappa_core.core.search.branching.construction import find_islands
    from mappa.mappa_workflow.mapping_analytics.family import compile_path
    from test_family_query import endpoint

    cycle = [(i, (i + 1) % 3, 1) for i in range(3)]
    source = graph(['C'] * 3 + ['H'] * 4, cycle + [(i, 3 + i, 1) for i in range(3)])
    target = graph(['C'] * 3 + ['H'] * 6,
                   cycle + [(i, 3 + 2 * i + h, 1) for i in range(3) for h in range(2)])
    problem = AAMProblem(endpoint('CCC' + 'H' * 4, source.edges),
                         endpoint('CCC' + 'H' * 6, target.edges))
    found = find_islands(source, target, list(range(7)), iso_tol=1, max_branches=100)
    paths = list(found.paths())
    assert len(paths) == 24  # Eight fixed ring/H sets, three remaining H images each.
    expected = set()
    for carbons in permutations(range(3)):
        for choices in product(range(2), repeat=3):
            images = tuple(carbons) + tuple(3 + 2 * p + h for p, h in zip(carbons, choices))
            for free in set(range(3, 9)) - set(images):
                expected.add(tuple(enumerate(images + (free,))))
    assert len(expected) == 144
    decoded = set()
    for path in paths:
        first, last = path.fragments
        assert set(first.source_atoms) == set(range(6)) and last.source_atoms == (6,)
        assert not set(first.target_atoms).intersection(last.target_atoms)
        compiled = compile_path(path, problem, {}, source_atoms=tuple(range(7)), complete_reference=False)
        compiled.solver.set(timeout=3000)
        while True:
            status = compiled.solver.check()
            assert status != z3.unknown, compiled.solver.reason_unknown()
            if status == z3.unsat:
                break
            mapping = tuple(map(tuple, compiled.realize(compiled.solver.model())['mapping']))
            assert mapping in expected and mapping not in decoded
            decoded.add(mapping)
            compiled.solver.add(z3.Or(*(value[0] != p
                for value, (_r, p) in zip(compiled.values, mapping))))
    assert decoded == expected
