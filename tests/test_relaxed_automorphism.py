"""Independent mapping oracles for constrained broad structural shuffles."""
from collections import deque
from itertools import permutations

import numpy as np
import pytest

from mappa.mappa_core.core.symmetry import relaxed_automorphism


def weights(size, edges):
    matrix = np.zeros((size, size))
    for a, b, value in edges:
        matrix[a, b] = matrix[b, a] = value
    return matrix


def permitted(elements, matrix, permutation, tolerance, threshold=.2, anchors=()):
    """Dense definition, independent of generators, responses and projections."""
    if any(elements[a] != elements[p] for a, p in enumerate(permutation)):
        return False
    if any(permutation[a] != a for a in anchors):
        return False
    for a in range(len(elements)):
        for b in range(a + 1, len(elements)):
            before = matrix[a, b]
            after = matrix[permutation[a], permutation[b]]
            if (before > threshold) != (after > threshold):
                return False
            if abs(before - after) > tolerance:
                return False
    return True


def broad_words(generators, size):
    """Expand only tiny test groups to audit the complete represented relation."""
    identity = tuple(range(size))
    words = {identity: ()}
    queue = deque([identity])
    while queue:
        previous = queue.popleft()
        for number, generator in enumerate(generators, 1):
            current = tuple(generator[a] for a in previous)
            if current not in words:
                words[current] = (*words[previous], number)
                queue.append(current)
                assert len(words) <= 720, 'test fixture exceeded its small-group scope'
    return words


def compare_complete_relation(relation, elements, matrix, tolerance, *, anchors=()):
    broad = broad_words(relation.generators, len(elements))
    expected_broad, expected_allowed = set(), set()
    for permutation in permutations(range(len(elements))):
        if permitted(elements, matrix, permutation, np.inf, anchors=anchors):
            expected_broad.add(permutation)
        expected = permitted(elements, matrix, permutation, tolerance, anchors=anchors)
        assert relation.contains(permutation) == expected, permutation
        if expected:
            expected_allowed.add(permutation)
    assert set(broad) == expected_broad
    for permutation, word in broad.items():
        expected = permutation in expected_allowed
        assert relation.allows(word) == expected
        inverse_word = tuple(-number for number in reversed(word))
        assert relation.allows(inverse_word) == expected
        if expected:
            assert relation.realize(word) == {
                relation.atoms[a]: relation.atoms[p] for a, p in enumerate(permutation)}
        else:
            with pytest.raises(ValueError):
                relation.realize(word)
    return broad, expected_allowed


def test_permitted_swaps_do_not_imply_permitted_composition():
    matrix = weights(4, [(0, 1, 1.), (0, 2, 1.6), (0, 3, 2.2)])
    relation = relaxed_automorphism('CHHH', matrix, tolerance=.7)
    words, allowed = compare_complete_relation(relation, 'CHHH', matrix, .7)
    left, right = (0, 2, 1, 3), (0, 1, 3, 2)
    assert allowed == {tuple(range(4)), left, right}
    assert len(words) == 6
    assert relation.allows(words[left]) and relation.allows(words[right])
    assert not relation.allows((*words[left], *words[right]))


@pytest.mark.parametrize('tolerance,expected', [(0., 6), (.2, 12)])
def test_ring_retains_correlated_rotations_and_reflections(tolerance, expected):
    matrix = weights(6, [(a, (a + 1) % 6, 1.1 if a % 2 else 1.2)
                         for a in range(6)])
    relation = relaxed_automorphism('CCCCCC', matrix, tolerance=tolerance)
    words, allowed = compare_complete_relation(relation, 'CCCCCC', matrix, tolerance)
    assert len(words) == 12
    assert len(allowed) == expected
    assert not relation.contains((1, 0, 2, 3, 4, 5))
    if tolerance == .2:
        assert relation.statistics['response_terms'] == 0
        assert relation.statistics['constrained_pairs'] == 0


def test_fragment_self_relationship_ignores_external_attachment_but_respects_anchors():
    # O-H-C-H has no global exchange. Inside H-C-H, the two H atoms can exchange.
    matrix = weights(4, [(0, 1, 1.), (1, 2, 1.1), (2, 3, 1.2)])
    global_relation = relaxed_automorphism('OHCH', matrix, tolerance=.2)
    _, global_allowed = compare_complete_relation(global_relation, 'OHCH', matrix, .2)
    assert global_allowed == {tuple(range(4))}

    for atoms in ((1, 2, 3), (3, 2, 1)):
        local = relaxed_automorphism('OHCH', matrix, tolerance=.2, atoms=atoms, anchors=(2,))
        local_elements = tuple('OHCH'[a] for a in local.atoms)
        local_weights = matrix[np.ix_(local.atoms, local.atoms)]
        center = local.atoms.index(2)
        words, allowed = compare_complete_relation(local, local_elements, local_weights,
                                                   .2, anchors=(center,))
        assert len(allowed) == 2
        swap = next(permutation for permutation in allowed if permutation != tuple(range(3)))
        assert local.realize(words[swap]) == {1: 3, 2: 2, 3: 1}

        anchored = relaxed_automorphism('OHCH', matrix, tolerance=.2,
                                       atoms=atoms, anchors=(1, 2))
        assert len(broad_words(anchored.generators, 3)) == 1
        assert anchored.realize() == {a: a for a in atoms}


@pytest.mark.parametrize('other,accepted', [(1.5, True), (np.nextafter(1.5, np.inf), False)])
def test_tolerance_includes_boundary_without_rounding(other, accepted):
    matrix = weights(3, [(0, 1, 1.), (0, 2, other)])
    relation = relaxed_automorphism('CHH', matrix, tolerance=.5)
    _, allowed = compare_complete_relation(relation, 'CHH', matrix, .5)
    assert ((0, 2, 1) in allowed) == accepted


def test_bond_presence_and_weight_tolerance_remain_distinct():
    matrix = weights(3, [(0, 1, .2), (0, 2, np.nextafter(.2, np.inf))])
    relation = relaxed_automorphism('CHH', matrix, tolerance=1.)
    _, allowed = compare_complete_relation(relation, 'CHH', matrix, 1.)
    assert allowed == {(0, 1, 2)}


def test_original_subthreshold_weights_still_constrain_relaxed_membership():
    matrix = weights(3, [(0, 1, .05), (0, 2, .1), (1, 2, .15)])
    relation = relaxed_automorphism('CCC', matrix, tolerance=.025)
    words, allowed = compare_complete_relation(relation, 'CCC', matrix, .025)
    assert len(words) == 6
    assert allowed == {(0, 1, 2)}


@pytest.mark.parametrize('seed', range(12))
def test_small_weighted_graphs_match_complete_independent_relation(seed):
    rng = np.random.default_rng(seed)
    matrix = weights(5, [(a, b, float(rng.choice((0., .1, .25, 1., 1.25))))
                         for a in range(5) for b in range(a + 1, 5)])
    tolerance = (0., .25, 1.)[seed % 3]
    relation = relaxed_automorphism('CCCCH', matrix, tolerance=tolerance)
    compare_complete_relation(relation, 'CCCCH', matrix, tolerance)


@pytest.mark.parametrize('atoms,anchors', [((), ()), ((2,), ()), ((2,), (2,))])
def test_empty_and_single_atom_fragments(atoms, anchors):
    matrix = weights(3, [(0, 1, 1.), (1, 2, 1.1)])
    relation = relaxed_automorphism('CCC', matrix, tolerance=.2, atoms=atoms, anchors=anchors)
    assert relation.atoms == atoms
    assert not relation.generators
    assert relation.contains(tuple(range(len(atoms))))
    assert relation.allows()
    assert relation.realize() == {atom: atom for atom in atoms}


def test_dendrimer_joint_words_match_dense_original_weight_checks():
    # A depth-three ternary tree represents 6**13 rearrangements. Test selected
    # complete choices without expanding that group in either test or production.
    matrix = weights(40, [((child - 1) // 3, child, 1.) for child in range(1, 40)])
    matrix[12, 39] = matrix[39, 12] = 1.8
    relation = relaxed_automorphism('C' * 40, matrix, tolerance=.5, anchors=(0,))
    assert relation.statistics['projected_pair_states'] <= 40 * 39
    assert relation.statistics['constrained_pairs'] <= 40 * 39 // 2
    rng = np.random.default_rng(83)
    words = [(), *((number,) for number in range(1, len(relation.generators) + 1))]
    words += [tuple(int(number) * int(sign) for number, sign in zip(
        rng.integers(1, len(relation.generators) + 1, size=8),
        rng.choice((-1, 1), size=8))) for _ in range(24)]
    outcomes = set()
    for word in words:
        permutation = tuple(range(40))
        for number in word:
            generator = relation.generators[abs(number) - 1]
            if number < 0:
                generator = tuple(np.argsort(generator))
            permutation = tuple(generator[a] for a in permutation)
        expected = permitted('C' * 40, matrix, permutation, .5, anchors=(0,))
        assert relation.contains(permutation) == expected
        assert relation.allows(word) == expected
        outcomes.add(expected)
    assert outcomes == {False, True}
