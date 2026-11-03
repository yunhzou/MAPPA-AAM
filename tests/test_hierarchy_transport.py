"""Coordinate transport preserves complete fragment relationships and sharing."""
from itertools import permutations
import io
import pickle

from mappa import match_fragment
from mappa.mappa_core.domain import AAMHierarchy, AAMHierarchyView, AAMHierarchyChain
from mappa.mappa_core.utils.io.checkpoints import _load_checkpoint
from growth_cases import graph, represented


def base():
    endpoint = graph(2, [(0, 1)])
    fragment, = match_fragment(endpoint, endpoint, seed=0).matches
    return AAMHierarchy((fragment, fragment))


def test_lazy_transport_preserves_every_mapping_in_each_atom_frame():
    hierarchy = base()
    original = hierarchy.fragments[0]
    expected = represented((original,))
    for images in permutations(range(4)):
        view = hierarchy.relabel_target(dict(enumerate(images)))
        if images == tuple(range(4)):
            assert view is hierarchy
        else:
            assert isinstance(view, AAMHierarchyView)
            assert 'materialized' not in view.__dict__
        transported = view.fragments[0]
        assert view.fragments[1] is transported
        assert transported.source_self_automorph is original.source_self_automorph
        assert transported.target_self_automorph.generators is original.target_self_automorph.generators
        assert transported.target_self_automorph.weights is original.target_self_automorph.weights
        assert represented((transported,)) == {
            tuple((r, images[p]) for r, p in pairs) for pairs in expected}
        assert transported.witness == tuple((r, images[p]) for r, p in original.witness)


def test_composed_views_flatten_and_keep_sparse_augmentation_frames():
    hierarchy = base()
    first = dict(enumerate((4, 1, 2, 3, 0)))
    second = dict(enumerate((0, 4, 2, 3, 1)))
    view = hierarchy.relabel_target(first).relabel_target(second)
    assert view.base is hierarchy
    assert 'materialized' not in view.__dict__
    combined = {a: second[first[a]] for a in first}
    assert view.to_record() == hierarchy._materialize_target(combined).to_record()
    assert set(view.fragments[0].target_atoms) == {1, 4}
    assert view.fragments[0].target_self_automorph.generators == ((1, 0),)


def test_saved_hierarchy_keeps_correlated_frames_after_module_relocation():
    hierarchy = base()
    chain = AAMHierarchyChain((hierarchy, hierarchy.relabel_target({0: 1, 1: 0})))
    saved = pickle.dumps(chain, protocol=0).replace(
        b'cmappa.mappa_core.domain\n', b'cmappa.alignment.post_aam\n')
    restored = _load_checkpoint(io.BytesIO(saved))
    assert isinstance(restored, AAMHierarchyChain)
    assert restored.parts[1].base is restored.parts[0]
    assert restored.parts[0].fragments[0] is restored.parts[0].fragments[1]
    assert restored.to_record() == chain.to_record()
