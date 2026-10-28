import numpy as np

from mappa.alignment.mapping_families import (
    compile_analytical_mapping_family,
)
from mappa.alignment.post_aam import (
    AtomBijection,
    AtomPermutation,
    PermutationGroup,
)


def test_bijection_group_action_is_target_o_mapping_o_source_inverse():
    mapping = AtomBijection((0, 1, 2, 3))
    source = AtomPermutation((1, 0, 2, 3))
    target = AtomPermutation((0, 1, 3, 2))

    transformed = mapping.act(source=source, target=target)

    assert transformed.images == (1, 0, 3, 2)


def test_permutation_group_orbits_come_from_generators_not_samples():
    group = PermutationGroup(6, (
        AtomPermutation((1, 0, 2, 3, 4, 5)),
        AtomPermutation((0, 2, 1, 3, 4, 5)),
        AtomPermutation((0, 1, 2, 4, 3, 5)),
    ))

    assert group.orbits() == ((0, 1, 2), (3, 4), (5,))


def test_analytical_mapping_family_dedupes_equal_cosets_not_distinct_cosets():
    elements = ["C"] * 4
    wbo = np.zeros((4, 4))
    hierarchy = {
        "fragments": [
            {"fragment_index": 0, "fragment": [0, 1]},
            {"fragment_index": 1, "fragment": [2, 3]},
        ],
    }
    identity = {0: 0, 1: 1, 2: 2, 3: 3}
    within_block = {0: 1, 1: 0, 2: 2, 3: 3}
    cross_block = {0: 0, 1: 2, 2: 1, 3: 3}

    first = compile_analytical_mapping_family(
        identity, hierarchy, elements, wbo, elements, wbo)
    same = compile_analytical_mapping_family(
        within_block, hierarchy, elements, wbo, elements, wbo)
    distinct = compile_analytical_mapping_family(
        cross_block, hierarchy, elements, wbo, elements, wbo)
    unrestricted = compile_analytical_mapping_family(
        identity, {"fragments": [{"fragment": [0, 1, 2, 3]}]},
        elements, wbo, elements, wbo)

    assert first.equivalent(same)
    assert same.equivalent(first)
    assert not first.equivalent(distinct)
    assert first.contains(identity)
    assert first.contains(within_block)
    assert not first.contains(cross_block)
    assert first.is_subset_of(unrestricted)
    assert distinct.is_subset_of(unrestricted)
    assert not unrestricted.is_subset_of(first)
