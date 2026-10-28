from itertools import permutations

from mappa.alignment.post_aam import (
    AAMHierarchy, AAMHierarchyView, AAMHierarchyChain, AtomPermutation, FragmentMatch)


def base():
    generator = AtomPermutation((1, 0, 2, 3))
    return AAMHierarchy(tuple(FragmentMatch(i, i, (i,),
        representative_assignments=((i, i),), target_generators=(generator,)) for i in (0, 1)))




def test_lazy_transport_equals_naive_conjugation_for_every_four_atom_action():
    hierarchy = base()
    for images in permutations(range(4)):
        view = hierarchy.relabel_target(dict(enumerate(images)))
        assert isinstance(view, AAMHierarchyView)
        assert "materialized" not in view.__dict__
        expected = [0] * 4
        original = hierarchy.fragments[0].target_generators[0].images
        for a in range(4):
            expected[images[a]] = images[original[a]]
        assert view.fragments[0].target_generators[0].images == tuple(expected)
        assert view.fragments[0].target_generators[0] is view.fragments[1].target_generators[0]
        assert view.fragments[0].representative_assignments == ((0, images[0]),)


def test_composed_views_flatten_and_pad_augmentation_frames():
    hierarchy = base()
    first = dict(enumerate((4, 1, 2, 3, 0)))
    second = dict(enumerate((0, 4, 2, 3, 1)))
    view = hierarchy.relabel_target(first).relabel_target(second)
    assert view.base is hierarchy
    assert "materialized" not in view.__dict__
    combined = {a: second[first[a]] for a in first}
    assert view.to_record() == hierarchy._materialize_target(combined).to_record()
    assert view.fragments[0].target_generators[0].degree == 5
