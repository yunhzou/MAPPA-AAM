"""Independent mapping oracles for binary structural shuffles."""
from collections import deque
from itertools import combinations, permutations

import numpy as np
import pytest

from mappa.mappa_core.core.symmetry import binary_symmetry


def weights(size, edges):
    matrix = np.zeros((size, size))
    for a, b, value in edges:
        matrix[a, b] = matrix[b, a] = value
    return matrix


def permitted(elements, matrix, permutation, threshold=.2, anchors=()):
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


def compare_complete_relation(relation, elements, matrix, *, anchors=()):
    words = broad_words(relation.generators, len(elements))
    expected = {p for p in permutations(range(len(elements)))
                if permitted(elements, matrix, p, anchors=anchors)}
    assert set(words) == expected
    for p in permutations(range(len(elements))):
        assert relation.contains(p) == (p in expected)
    for permutation, word in words.items():
        assert relation.realize(word) == {
            relation.atoms[a]: relation.atoms[p] for a, p in enumerate(permutation)}
    return words, expected


def test_binary_swaps_compose_regardless_of_original_weight_differences():
    matrix = weights(4, [(0, 1, 1.), (0, 2, 1.6), (0, 3, 2.2)])
    relation = binary_symmetry('CHHH', matrix)
    words, allowed = compare_complete_relation(relation, 'CHHH', matrix)
    assert len(allowed) == 6
    left, right = (0, 2, 1, 3), (0, 1, 3, 2)
    combined = relation.realize((*words[left], *words[right]))
    assert tuple(combined.values()) in allowed
    assert not hasattr(relation, 'tolerance') and not hasattr(relation, 'constraints')


def test_ring_retains_correlated_rotations_and_reflections():
    matrix = weights(6, [(a, (a + 1) % 6, 1.1 if a % 2 else 1.2) for a in range(6)])
    relation = binary_symmetry('CCCCCC', matrix)
    words, allowed = compare_complete_relation(relation, 'CCCCCC', matrix)
    assert len(words) == len(allowed) == 12
    assert not relation.contains((1, 0, 2, 3, 4, 5))


def test_fragment_self_relationship_ignores_external_attachment_but_respects_anchors():
    # O-H-C-H has no global exchange. Inside H-C-H, the two H atoms can exchange.
    matrix = weights(4, [(0, 1, 1.), (1, 2, 1.1), (2, 3, 1.2)])
    global_relation = binary_symmetry('OHCH', matrix)
    _, global_allowed = compare_complete_relation(global_relation, 'OHCH', matrix)
    assert global_allowed == {tuple(range(4))}

    for atoms in ((1, 2, 3), (3, 2, 1)):
        local = binary_symmetry('OHCH', matrix,  atoms=atoms, anchors=(2,))
        local_elements = tuple('OHCH'[a] for a in local.atoms)
        local_weights = matrix[np.ix_(local.atoms, local.atoms)]
        center = local.atoms.index(2)
        words, allowed = compare_complete_relation(local, local_elements, local_weights,
                                                   anchors=(center,))
        assert len(allowed) == 2
        swap = next(permutation for permutation in allowed if permutation != tuple(range(3)))
        assert local.realize(words[swap]) == {1: 3, 2: 2, 3: 1}

        anchored = binary_symmetry('OHCH', matrix,
                                       atoms=atoms, anchors=(1, 2))
        assert len(broad_words(anchored.generators, 3)) == 1
        assert anchored.realize() == {a: a for a in atoms}


@pytest.mark.parametrize('other', [1.5, np.nextafter(1.5, np.inf), 3.])
def test_original_weight_variation_does_not_filter_binary_symmetry(other):
    matrix = weights(3, [(0, 1, 1.), (0, 2, other)])
    relation = binary_symmetry('CHH', matrix)
    _, allowed = compare_complete_relation(relation, 'CHH', matrix)
    assert (0, 2, 1) in allowed


def test_bond_presence_uses_the_strict_threshold():
    matrix = weights(3, [(0, 1, .2), (0, 2, np.nextafter(.2, np.inf))])
    relation = binary_symmetry('CHH', matrix)
    _, allowed = compare_complete_relation(relation, 'CHH', matrix)
    assert allowed == {(0, 1, 2)}


def test_subthreshold_contacts_remain_available_without_restricting_symmetry():
    matrix = weights(3, [(0, 1, .05), (0, 2, .1), (1, 2, .15)])
    relation = binary_symmetry('CCC', matrix)
    words, allowed = compare_complete_relation(relation, 'CCC', matrix)
    assert len(words) == 6
    assert len(allowed) == 6
    assert np.array_equal(relation.weights, matrix)


@pytest.mark.parametrize('seed', range(12))
def test_small_weighted_graphs_match_complete_independent_relation(seed):
    rng = np.random.default_rng(seed)
    matrix = weights(5, [(a, b, float(rng.choice((0., .1, .25, 1., 1.25))))
                         for a in range(5) for b in range(a + 1, 5)])
    relation = binary_symmetry('CCCCH', matrix)
    compare_complete_relation(relation, 'CCCCH', matrix)


@pytest.mark.parametrize('atoms,anchors', [((), ()), ((2,), ()), ((2,), (2,))])
def test_empty_and_single_atom_fragments(atoms, anchors):
    matrix = weights(3, [(0, 1, 1.), (1, 2, 1.1)])
    relation = binary_symmetry('CCC', matrix,  atoms=atoms, anchors=anchors)
    assert relation.atoms == atoms
    assert not relation.generators
    assert relation.contains(tuple(range(len(atoms))))
    assert relation.realize() == {atom: atom for atom in atoms}


def test_dendrimer_joint_words_match_dense_original_weight_checks():
    # A depth-three ternary tree represents 6**13 rearrangements. Test selected
    # complete choices without expanding that group in either test or production.
    matrix = weights(40, [((child - 1) // 3, child, 1.) for child in range(1, 40)])
    matrix[12, 39] = matrix[39, 12] = 1.8
    relation = binary_symmetry('C' * 40, matrix,  anchors=(0,))
    assert len(relation.generators) < 40
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
        expected = permitted('C' * 40, matrix, permutation, anchors=(0,))
        assert relation.contains(permutation) == expected
        outcomes.add(expected)
    assert outcomes == {True}


@pytest.mark.parametrize('generators', [(), ((1, 2, 3, 4, 0),), ((1, 0, 3, 2, 4),)])
@pytest.mark.parametrize('selected', [None, ((1, 4), (2, 3))])
def test_pair_orbits_match_complete_small_group(generators, selected):
    from mappa.mappa_core.core.symmetry.groups import pair_orbits, permutation_group

    group = permutation_group(generators, 5)
    choices = broad_words(generators, 5)
    pairs = set(combinations(range(5), 2) if selected is None else selected)
    expected = {frozenset(tuple(sorted((choice[a], choice[b]))) for choice in choices)
                for a, b in pairs}
    ordered = sorted(expected, key=lambda orbit: min(orbit & pairs))
    assert list(pair_orbits(group, 5, selected)) == [tuple(sorted(orbit)) for orbit in ordered]


def test_placement_identity_avoids_group_search_and_nonidentity_remains_checked(monkeypatch):
    from mappa.mappa_core.core.symmetry import groups
    from mappa.mappa_core.core.symmetry.actions import apply_action

    action = ('group', ((1, 0, 3, 2),))
    images = np.array([0, 2])
    with monkeypatch.context() as patch:
        def unnecessary_search(*args):
            pytest.fail('the identity requires no group construction or membership search')
        patch.setattr(groups, 'permutation_group', unnecessary_search)
        assert np.array_equal(apply_action(images, action, {'permutation': (0, 1, 2, 3)}, 4), images)
    assert tuple(apply_action(images, action, {'permutation': (1, 0, 3, 2)}, 4)) == (1, 3)
    # Each swap is individually in an atom orbit, but only their joint action
    # belongs to this placement group.
    with pytest.raises(ValueError, match='outside the saved group'):
        apply_action(images, action, {'permutation': (1, 0, 2, 3)}, 4)
