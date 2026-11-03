"""Declared fragment anchors constrain rematching; a stored witness does not."""
from dataclasses import replace

import numpy as np
import pytest

from mappa import BondTolerance, FragmentMatchConfig, FragmentMatchContext, match_fragment
from mappa.mappa_core.domain import AAMHierarchy
from mappa.mappa_core.core.matching.growth import native
from mappa.mappa_core.core.matching.representation.records import Fragment
from mappa.mappa_core.core.symmetry.relaxed import relaxed_automorphism
from mappa.mappa_core.core.symmetry.relations import RelaxedIsomorphism, DiscoveryEdge
from mappa.mappa_core.post_processing.final_branches import FinalFamily
from mappa.mappa_core.core.search.branching.construction import find_islands
from mappa.mappa_core.utils.bonds import bond_event_threshold
from growth_cases import graph as make_graph, represented, discovery_mappings


def graph(elements, edges):
    return make_graph(len(elements), edges, elements)


def event_threshold(source, a, b):
    return bond_event_threshold({n: source.nodes[n]['element'] for n in (a, b)}, a, b,
                               default_threshold=.5, metal_threshold=.3)


def saved_fragment(mapping, *, anchors=(), source=None, target=None):
    """Construct a complete saved relation; anchors are explicit scientific data."""
    atoms, images = tuple(sorted(mapping)), tuple(sorted(mapping.values()))
    anchors = tuple((r, mapping[r]) for r in anchors)
    def prepare(graph, selected, fixed):
        weights = (np.zeros((len(selected), len(selected))) if graph is None else
                   graph.graph['wbo_matrix'][np.ix_(selected, selected)])
        elements = tuple('C' if graph is None else graph.nodes[a]['element'] for a in selected)
        value = relaxed_automorphism(elements, weights, tolerance=1.5,
            anchors=tuple(selected.index(atom) for atom in fixed))
        return replace(value, atoms=selected, anchors=tuple(fixed))
    left = prepare(source, atoms, tuple(r for r, _ in anchors))
    right = prepare(target, images, tuple(p for _, p in anchors))
    domains = tuple((r, tuple(p for p in images if left.elements[atoms.index(r)] ==
                             right.elements[images.index(p)] and
                             (r not in dict(anchors) or dict(anchors)[r] == p))) for r in atoms)
    edges = (() if source is None else tuple(DiscoveryEdge((a,b),
             float(source.graph['wbo_matrix'][a,b]), 2.) for a,b in source.edges))
    target_edges = () if target is None else tuple(target.edges)
    return Fragment(RelaxedIsomorphism(left, right, tuple(mapping.items()), anchors,
                                      domains, edges, target_edges, .2))


def test_saved_witness_does_not_create_undeclared_anchors():
    fragment = saved_fragment({2: 8, 4: 6})
    context = FragmentMatchContext.from_saved_fragment(fragment)
    assert context.locked_mapping == context.islands == {}
    assert fragment.fixed_assignments == ()
    assert len(represented((fragment,))) == 2


def test_explicit_sparse_anchors_and_deferred_continuation_are_preserved():
    fragment = saved_fragment({2: 8, 4: 6}, anchors=(2, 4))
    context = FragmentMatchContext.from_saved_fragment(fragment, deferred_edges=((2, 4),))
    assert context.locked_mapping == {2: 8, 4: 6}
    assert context.islands == {2: 2, 4: 4}
    assert context.deferred_edges == ((2, 4),)
    assert FragmentMatchContext.from_saved_fragment(fragment).deferred_edges == ()
    assert represented((fragment,)) == {((2, 8), (4, 6))}


def test_anchor_keeps_other_assignments_free():
    fragment = saved_fragment({0: 0, 1: 1, 2: 2}, anchors=(0,))
    assert fragment.fixed_assignments == ((0, 0),)
    assert represented((fragment,)) == {
        ((0, 0), (1, 1), (2, 2)), ((0, 0), (1, 2), (2, 1))}


def test_anchors_survive_serialization_and_target_transport():
    fragment = saved_fragment({0: 0, 1: 1, 2: 2}, anchors=(0,))
    hierarchy = AAMHierarchy((fragment,))
    restored = AAMHierarchy.from_record(hierarchy.to_record())
    assert restored.fragments[0].fixed_assignments == ((0, 0),)
    transported = hierarchy.relabel_target({0: 2, 1: 0, 2: 1})
    assert transported.fragments[0].fixed_assignments == ((0, 2),)
    assert represented(transported.fragments) == {
        ((0, 2), (1, 0), (2, 1)), ((0, 2), (1, 1), (2, 0))}


def test_final_family_uses_saved_relations_without_freezing_the_witness():
    first = saved_fragment({1: 1, 2: 2})
    second = saved_fragment({3: 3, 4: 4})
    family = FinalFamily(tuple(enumerate(range(6))), (),
        (('relation', first), ('relation', second)), (.2, 1.))
    assert family.fixed_assignments == ((0, 0), (5, 5))


@pytest.mark.parametrize('backend', ['0', '1'])
def test_scheduler_splits_mixed_threshold_fragment(backend, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', backend)
    if backend == '1':
        assert native.built()
    source = graph(['Co', 'C', 'O'], [(0, 1, 1.), (1, 2, 1.)])
    target = graph(['Co', 'C', 'O'], [(0, 1, 1.4), (1, 2, 1.4)])
    result = find_islands(source, target, [1, 0, 2], iso_tol=.5,
        anchor_map={1: 1}, bond_tolerance=BondTolerance(event_threshold, inclusive=False))
    paths = list(result.paths())
    assert paths and not any(stop.reason == 'capped' for stop in result.stops)
    for path in paths:
        assert path.mapping == {0: 0, 1: 1, 2: 2}
        owners = dict(result.states[path.terminal].islands)
        assert owners[1] == owners[2] != owners[0]


@pytest.mark.parametrize('backend', ['0', '1'])
def test_tight_rematch_rediscover_mapping_instead_of_pinning_saved_witness(backend, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', backend)
    if backend == '1':
        assert native.built()
    source = graph(['O', 'C', 'C'], [(0, 1, 1.), (0, 2, 2.)])
    target = graph(['O', 'C', 'C'], [(0, 1, 1.), (0, 2, 2.)])
    fragment = saved_fragment({0: 0, 1: 2, 2: 1}, anchors=(0,), source=source, target=target)
    context = FragmentMatchContext.from_saved_fragment(fragment)
    assert context.locked_mapping == {0: 0}
    result = match_fragment(source, target, seed=0, context=context,
        config=FragmentMatchConfig(allow_mapped_seed=True,
                                  bond_tolerance=BondTolerance(event_threshold, inclusive=False)))
    assert not result.capped
    assert len(result.matches) == 1
    assert dict(result.matches[0].witness) == {0: 0, 1: 1, 2: 2}
    assert result.matches[0].source_atoms == (0, 1, 2)
    assert fragment.witness == ((0, 0), (1, 2), (2, 1))


@pytest.mark.parametrize('backend', ['0', '1'])
def test_rematch_rediscovers_compressed_symmetry(backend, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', backend)
    if backend == '1':
        assert native.built()
    source = graph(['C', 'H', 'H', 'H'], [(0, h, 1.) for h in (1, 2, 3)])
    target = graph(['C', 'H', 'H', 'H'], [(0, h, 1.) for h in (1, 2, 3)])
    fragment = saved_fragment({0: 0, 1: 2, 2: 3, 3: 1}, anchors=(0,), source=source, target=target)
    result = match_fragment(source, target, seed=0,
        context=FragmentMatchContext.from_saved_fragment(fragment),
        config=FragmentMatchConfig(allow_mapped_seed=True,
                                  bond_tolerance=BondTolerance(event_threshold, inclusive=False)))
    assert not result.capped and len(result.matches) == 1
    assert len(represented(result.matches)) == 6
