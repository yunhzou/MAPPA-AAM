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
