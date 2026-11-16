"""Differential oracles for native bookkeeping, independent of chemistry cases."""
from itertools import permutations
import random
import subprocess
import sys

import pytest

from mappa._group_ops import (occupation_orbit, project_generators,
                                 conjugate_generators, OccupationLimitExceeded, SetOrbit,
                                 PermutationGroup)


def test_native_group_chain_matches_small_exhaustive_oracle():
    from growth_cases import group_actions
    rng = random.Random(1703)
    for degree in range(2, 8):
        for _ in range(12):
            generators = [rng.sample(range(degree), degree) for _ in range(rng.randrange(4))]
            group = PermutationGroup(degree, generators)
            expected = set(group_actions(generators, degree))
            assert group.order() == len(expected)
            for _ in range(12):
                selected = tuple(rng.sample(range(degree), degree))
                assert group.contains(selected) == (selected in expected)
            fixed = rng.sample(range(degree), rng.randrange(degree + 1))
            subgroup = group.stabilizer(fixed)
            assert set(group_actions(subgroup.generators, degree)) == {
                g for g in expected if all(g[p] == p for p in fixed)}
            point = rng.randrange(degree)
            assert set(group.orbit(point)) == {g[point] for g in expected}
            for image, action in group.transversal(point):
                assert action[point] == image and tuple(action) in expected
            # Compact factor choices must cover the group with correlated
            # permutations, including the orientation used by query compilation.
            represented = {tuple(range(degree))}
            for level in reversed(group.transversals()):
                represented = {tuple(g[p] for p in previous)
                               for previous in represented for _, g in level}
            assert represented == expected


def _check_large_group_chain():
    from math import factorial
    degree = 32
    generators = []
    for i in range(degree - 1):
        action = list(range(degree))
        action[i], action[i + 1] = action[i + 1], action[i]
        generators.append(action)
    group = PermutationGroup(degree, generators)
    assert group.order() == factorial(degree)
    assert group.contains(list(reversed(range(degree))))
    assert group.stabilizer(range(16)).order() == factorial(16)


def test_native_group_chain_keeps_factorial_symmetry_compact():
    subprocess.run([sys.executable, '-c',
        "import runpy,sys; runpy.run_path(sys.argv[1])['_check_large_group_chain']()",
        __file__], check=True, timeout=10)


def _check_independent_restrictions():
    # Fixing many independent pairs must not repeatedly rebuild one large chain.
    count = 128
    degree = 2 * count
    generators = []
    for i in range(count):
        action = list(range(degree))
        action[2*i], action[2*i+1] = action[2*i+1], action[2*i]
        generators.append(action)
    group = PermutationGroup(degree, generators)
    assert not group.stabilizer(range(degree), seconds=1.).generators
    remaining = group.stabilizer(range(degree-2), seconds=1.)
    assert remaining.generators == (tuple(generators[-1]),)


def test_native_restriction_does_not_rebuild_unrelated_pair_symmetries():
    subprocess.run([sys.executable, '-c',
        "import runpy,sys; runpy.run_path(sys.argv[1])['_check_independent_restrictions']()",
        __file__], check=True, timeout=10)


def test_native_restriction_preserves_coupled_and_composed_operations():
    from growth_cases import group_actions

    # A coupled pair exchange cannot become two independently usable swaps.
    coupled = PermutationGroup(4, [(1, 0, 3, 2)])
    assert not coupled.stabilizer((0,)).generators
    assert set(group_actions(coupled.stabilizer(()).generators, 4)) == {(0, 1, 2, 3), (1, 0, 3, 2)}
    # Keeping only individually eligible input generators would lose this swap.
    composed = PermutationGroup(3, [(1, 0, 2), (1, 2, 0)])
    assert set(group_actions(composed.stabilizer((0,)).generators, 3)) == {(0, 1, 2), (0, 2, 1)}


def test_native_restriction_deadline_interrupts_chain_construction():
    from mappa._group_ops import SearchDeadlineReached
    from time import perf_counter

    degree = 128
    generators = []
    for i in range(degree-1):
        action = list(range(degree))
        action[i], action[i+1] = action[i+1], action[i]
        generators.append(action)
    group = PermutationGroup(degree, generators)
    start = perf_counter()
    with pytest.raises(SearchDeadlineReached):
        group.stabilizer((0,), seconds=.001)
    assert perf_counter()-start < 1.
    # Interruption must not change the original group or publish partial output.
    assert group.generators == tuple(map(tuple, generators))
    assert not group.stabilizer(range(degree), seconds=1.).generators


def reference(witness, degree, stages, attachments, fragments, bonds, limit=-1):
    def key(images):
        return (tuple(sorted(images)), tuple(sorted(images[i] for i in attachments)),
                tuple(sorted((label, tuple(sorted(images[i] for i in positions)))
                             for label, positions in fragments)),
                tuple(sorted(tuple(sorted((images[a], images[b]))) for a, b in bonds)))
    states = {key(witness): (tuple(witness), tuple(range(degree)))}
    for stage in stages:
        queue = list(states.values())
        for images, action in queue:
            for generator in stage:
                moved = tuple(generator[p] for p in images)
                relation = key(moved)
                if relation in states:
                    continue
                if limit >= 0 and len(states) >= limit:
                    raise OccupationLimitExceeded()
                state = moved, tuple(generator[p] for p in action)
                states[relation] = state
                queue.append(state)
    return list(states.values())


def test_exact_orbit_matches_reference_including_first_witness_and_stage_order():
    rng = random.Random(416)
    for degree in range(1, 8):
        for _ in range(40):
            size = rng.randrange(1, degree + 1)
            witness = rng.sample(range(degree), size)
            stages = []
            for _ in range(3):
                stages.append([rng.sample(range(degree), degree) for _ in range(rng.randrange(3))])
            fragments = [(i % 2, list(range(i, size, 2))) for i in range(min(2, size))]
            attachments = list(range(0, size, 3))
            bonds = [(i, i + 1) for i in range(size - 1)]
            actual = occupation_orbit(witness, degree, stages, attachments, fragments, bonds, -1)
            assert [(tuple(images), tuple(action)) for images, action in actual] == reference(
                witness, degree, stages, attachments, fragments, bonds)


def test_limits_and_invalid_native_inputs_are_reported():
    with pytest.raises(OccupationLimitExceeded):
        occupation_orbit([0], 2, [[[1, 0]]], [], [(0, [0])], [], 1)
    with pytest.raises(ValueError):
        occupation_orbit([0], 2, [[[1, 1]]], [], [(0, [0])], [], -1)
    with pytest.raises(ValueError):
        occupation_orbit([0], 2, [], [], [(0, [1])], [], -1)
    with pytest.raises(ValueError):
        occupation_orbit([0], 2, [], [], [(0, [0])], [], -1, [2])


def _observed_relation(images, observed, attachments, fragments, bonds):
    values = tuple(p if p in observed else -1 for p in images)
    return (tuple(sorted(values)), tuple(sorted(values[i] for i in attachments)),
            tuple(sorted((label, tuple(sorted(values[i] for i in positions)))
                         for label, positions in fragments)),
            tuple(sorted(tuple(sorted((values[a], values[b]))) for a, b in bonds)))


def test_permanently_competitor_only_choices_stay_compressed():
    witness = [0, 2, 3]
    stages = [[[0, 1, 3, 2, 4, 5], [0, 1, 2, 4, 3, 5], [0, 1, 2, 3, 5, 4]]]
    fragments = [(0, [0]), (1, [1]), (1, [2])]
    full = occupation_orbit(witness, 6, stages, [], fragments, [], -1)
    compact = occupation_orbit(witness, 6, stages, [], fragments, [], -1, [0, 1])
    assert len(full) == 6
    assert len(compact) == 1
    assert compact[0] == full[0]


def test_future_actions_keep_current_competitors_distinct():
    # The first stage exchanges two competitors; the second sends them to
    # different target atoms. Projecting to the target prematurely loses atom 3.
    stages = [[[0, 2, 1, 3]], [[1, 0, 3, 2]]]
    compact = occupation_orbit([1], 4, stages, [], [(0, [0])], [], -1, [0, 3])
    full = occupation_orbit([1], 4, stages, [], [(0, [0])], [], -1)
    assert compact == full
    assert {images[0] for images, _ in compact} == {0, 1, 2, 3}


def test_observation_quotient_preserves_exhaustive_projected_relations():
    rng = random.Random(1976)
    for degree in range(2, 9):
        for _ in range(50):
            # Some actions mix target/competitor positions; others leave an
            # independent competitor component that can really be compressed.
            active = rng.randrange(1, degree + 1)
            stages = []
            for _ in range(3):
                generators = []
                for _ in range(rng.randrange(4)):
                    generators.append(rng.sample(range(active), active) +
                                      rng.sample(range(active, degree), degree - active))
                stages.append(generators)
            size = rng.randrange(1, degree + 1)
            witness = rng.sample(range(degree), size)
            observed = set(rng.sample(range(active), rng.randrange(active + 1)))
            fragments = [(i % 2, list(range(i, size, 2))) for i in range(min(2, size))]
            attachments = list(range(0, size, 3))
            bonds = [(i, i + 1) for i in range(size - 1)]
            full = reference(witness, degree, stages, attachments, fragments, bonds)
            compact = occupation_orbit(witness, degree, stages, attachments, fragments,
                                       bonds, -1, sorted(observed))
            key = lambda images: _observed_relation(images, observed, attachments, fragments, bonds)
            assert {key(images) for images, _ in compact} == {key(images) for images, _ in full}
            for images, action in compact:
                assert sorted(action) == list(range(degree))
                assert list(images) == [action[p] for p in witness]


def test_projection_preserves_generator_order_and_sparse_atom_labels():
    raw = [(1, 0, 3, 2), (0, 1, 2, 3), (1, 0, 3, 2)]
    assert project_generators(raw, [0, 3], [0, 1]) == ((3, 1, 2, 0),)
    with pytest.raises(ValueError):
        project_generators([(2, 1, 0)], [0, 1], [0, 1])


def test_conjugation_preserves_padded_frames_exhaustively():
    for action in permutations(range(4)):
        for generator in permutations(range(3)):
            expected = list(range(4))
            for i in range(4):
                expected[action[i]] = action[generator[i] if i < 3 else i]
            assert conjugate_generators([generator], action) == (tuple(expected),)


def test_point_lock_releases_independent_factors_from_a_correlated_exchange():
    from growth_cases import group_actions

    # Exchanging two central atoms carries both pairs attached to each center.
    # Fixing one center removes that coupled exchange but retains all four
    # independent pair swaps, including those on the other center's side.
    degree = 10
    generators = [(1, 0, 6, 7, 8, 9, 2, 3, 4, 5)]
    for a, b in ((2, 3), (4, 5), (6, 7), (8, 9)):
        action = list(range(degree))
        action[a], action[b] = b, a
        generators.append(tuple(action))
    original = group_actions(generators, degree)
    group = SetOrbit(degree, generators)
    for fixed, count in (((0,), 16), ((0, 2), 8), ((0, 2, 6), 4)):
        expected = {action for action in original if all(action[p] == p for p in fixed)}
        actual = set(group_actions(group.stabilizer(fixed), degree))
        assert actual == expected and len(actual) == count


def _check_correlated_set_transport(case):
    # A column permutation moves both rows together. Row counts and the number
    # of columns occupied in both rows are independent invariants of this group.
    width = 24
    generators = []
    for column in range(width - 1):
        action = list(range(2 * width))
        for row in (0, width):
            action[row + column], action[row + column + 1] = row + column + 1, row + column
        generators.append(action)
    group = SetOrbit(2 * width, generators)
    first = tuple(range(12)) + tuple(range(width, width + 12))
    if case == 'different_row_counts':
        second = tuple(range(11)) + tuple(range(width, width + 13))
        expected = False
    elif case == 'different_overlap':
        second = tuple(range(12)) + tuple(range(width + 12, 2 * width))
        expected = False
    else:
        second = tuple(range(12, width)) + tuple(range(width + 12, 2 * width))
        expected = True
    assert group.contains(first, second) is expected


@pytest.mark.parametrize('case', ['different_row_counts', 'different_overlap', 'joint_move'])
def test_correlated_set_transport_does_not_expand_equivalent_permutations(case):
    # A process deadline catches factorial recursion inside C++, where a Python
    # signal handler cannot interrupt the call. This is not a microbenchmark.
    subprocess.run([sys.executable, '-c',
        "import runpy,sys; runpy.run_path(sys.argv[1])['_check_correlated_set_transport'](sys.argv[2])",
        __file__, case], check=True, timeout=10)
