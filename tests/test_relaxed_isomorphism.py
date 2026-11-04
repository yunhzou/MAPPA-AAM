"""Independent finite oracles for the complete two-sided fragment relation."""
from dataclasses import FrozenInstanceError, replace
from itertools import permutations

import numpy as np
import pytest

from mappa.mappa_core.core.symmetry import (
    DiscoveryEdge, RelaxedIsomorphism, relaxed_automorphism)
from test_relaxed_automorphism import broad_words, permitted, weights


def relation(source, target, accepted, *, tolerance=0., self_tolerance=0.,
             anchors=(), domains=None, inclusive=True, target_edges=None, floor=.2):
    n = len(source)
    source_self = relaxed_automorphism('C' * n, source, tolerance=self_tolerance,
                                      anchors=tuple(r for r, _ in anchors))
    target_self = relaxed_automorphism('C' * n, target, tolerance=self_tolerance,
                                      anchors=tuple(p for _, p in anchors))
    edges = tuple(DiscoveryEdge((a, b), source[a, b], tolerance, inclusive)
                  for a, b in accepted)
    if target_edges is None:
        target_edges = tuple((a, b) for a in range(n) for b in range(a + 1, n)
                             if target[a, b] >= floor)
    return RelaxedIsomorphism(source_self, target_self, tuple(enumerate(range(n))),
        anchors, tuple((r, tuple(range(n))) for r in range(n)) if domains is None else domains,
        edges, target_edges, floor)


def discovery_oracle(value):
    """Literal definition of D; used only on tiny fixtures."""
    source, target = value.source_atoms, value.target_atoms
    positions = {p: i for i, p in enumerate(target)}
    domains = dict(value.domains)
    active = {frozenset(edge) for edge in value.target_edges}
    result = set()
    for images in permutations(target):
        mapping = dict(zip(source, images))
        if any(mapping[r] not in domains[r] for r in source):
            continue
        if any(mapping[r] != p for r, p in value.anchors):
            continue
        allowed = True
        for edge in value.accepted_edges:
            p, q = (mapping[r] for r in edge.atoms)
            weight = value.target.weights[positions[p], positions[q]]
            difference = abs(edge.weight - weight)
            within = difference <= edge.tolerance if edge.inclusive else difference < edge.tolerance
            allowed &= frozenset((p, q)) in active and weight >= value.graph_floor and within
        if allowed:
            result.add(images)
    return result


def self_oracle(value):
    positions = {atom: i for i, atom in enumerate(value.atoms)}
    anchors = tuple(positions[atom] for atom in value.anchors)
    return tuple(p for p in permutations(range(len(value.atoms)))
                 if permitted(value.elements, value.weights, p, value.tolerance,
                              value.bond_threshold, anchors))


def represented(value):
    """Enumerate only the independent test oracle, never production choices."""
    source_words = broad_words(value.source.generators, len(value.source_atoms))
    target_words = broad_words(value.target.generators, len(value.target_atoms))
    expected, actual = set(), set()
    for images in discovery_oracle(value):
        discovery = dict(zip(value.source_atoms, images))
        for alpha in self_oracle(value.source):
            inverse = tuple(alpha.index(i) for i in range(len(alpha)))
            for beta in self_oracle(value.target):
                target_position = {atom: i for i, atom in enumerate(value.target_atoms)}
                expected.add(tuple(value.target_atoms[beta[target_position[images[inverse[i]]]]]
                                   for i in range(len(images))))
                mapped = value.realize(discovery, source_word=source_words[alpha],
                                       target_word=target_words[beta])
                actual.add(tuple(mapped[r] for r in value.source_atoms))
    assert actual == expected
    return actual


def test_triangle_tail_retains_all_four_discoveries_and_six_two_sided_choices():
    matrix = weights(4, [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.), (0, 3, 1.)])
    value = relation(matrix, matrix, ((0, 1), (1, 2), (0, 3)))
    discoveries = discovery_oracle(value)
    assert len(discoveries) == 4
    for images in permutations(range(4)):
        assert value.allows_discovery(dict(enumerate(images))) == (images in discoveries)
    choices = represented(value)
    assert len(choices) == 6
    assert choices - discoveries == {(1, 3, 0, 2), (2, 3, 0, 1)}
    # One witness's self orbit contains only two choices and is incomplete.
    source_words = broad_words(value.source.generators, 4)
    target_words = broad_words(value.target.generators, 4)
    witness_orbit = {tuple(value.realize(source_word=a, target_word=b).values())
                     for a in source_words.values() for b in target_words.values()}
    assert len(witness_orbit) == 2


def test_source_inverse_precedes_discovery_and_target_action():
    matrix = weights(3, [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.)])
    value = relation(matrix, matrix, ((0, 1), (1, 2)))
    words = broad_words(value.source.generators, 3)
    discovery = {0: 2, 1: 0, 2: 1}
    source_word, target_word = words[(1, 2, 0)], words[(1, 0, 2)]
    assert value.realize(discovery, source_word=source_word) == {0: 1, 1: 2, 2: 0}
    assert value.realize(discovery, source_word=source_word, target_word=target_word) == {
        0: 0, 1: 2, 2: 1}


def test_self_tolerance_is_checked_jointly_and_relation_is_not_repeated_closure():
    matrix = weights(4, [(0, 1, 1.), (0, 2, 1.6), (0, 3, 2.2)])
    value = relation(matrix, matrix, ((0, 1), (0, 2), (0, 3)), self_tolerance=.7)
    assert discovery_oracle(value) == {(0, 1, 2, 3)}
    words = broad_words(value.source.generators, 4)
    left, right = words[(0, 2, 1, 3)], words[(0, 1, 3, 2)]
    moved = value.realize(source_word=left, target_word=right)
    assert not value.allows_discovery(moved)
    with pytest.raises(ValueError, match='accepted edge'):
        value.realize(moved)
    for options in (dict(source_word=(*left, *right)), dict(target_word=(*left, *right))):
        with pytest.raises(ValueError, match='self tolerance'):
            value.realize(**options)
    choices = represented(value)
    assert len(choices) == 5
    assert (0, 3, 2, 1) not in choices


def test_sparse_indices_anchors_and_cut_adjacency_use_their_own_frames():
    source_matrix = weights(8, [(2, 5, 1.), (2, 7, 1.1)])
    target_matrix = weights(9, [(1, 4, 1.1), (1, 8, 1.2)])
    source = relaxed_automorphism('C' * 8, source_matrix, tolerance=.2,
                                 atoms=(2, 5, 7), anchors=(2,))
    target = relaxed_automorphism('C' * 9, target_matrix, tolerance=.2,
                                 atoms=(1, 4, 8), anchors=(1,))
    value = RelaxedIsomorphism(source, target, ((2, 1), (5, 4), (7, 8)), ((2, 1),),
        ((2, (1,)), (5, (4, 8)), (7, (4, 8))),
        (DiscoveryEdge((2, 5), 1., .2),), ((1, 4),), .2)
    assert value.source is source and value.target is target
    assert value.source_atoms is source.atoms and value.target_atoms is target.atoms
    assert discovery_oracle(value) == {(1, 4, 8)}
    assert value.discovery_pairs(value.accepted_edges[0]) == ((1, 4),)
    assert represented(value) == {(1, 4, 8), (1, 8, 4)}
    # The removed discovery edge does not alter original fragment self symmetry.
    assert len(broad_words(target.generators, 3)) == 2


@pytest.mark.parametrize('inclusive', [False, True])
def test_discovery_predicate_preserves_exact_tolerance_boundary(inclusive):
    source = weights(3, [(0, 1, 1.), (0, 2, 1.5)])
    value = relation(source, source, ((0, 1),), tolerance=.5, inclusive=inclusive)
    swap = {0: 0, 1: 2, 2: 1}
    assert value.allows_discovery(swap) is inclusive
    pairs = value.discovery_pairs(value.accepted_edges[0])
    assert ((0, 2) in pairs) is inclusive
    raised = source.copy()
    raised[0, 2] = raised[2, 0] = np.nextafter(1.5, np.inf)
    other = relation(source, raised, ((0, 1),), tolerance=.5, inclusive=inclusive)
    assert not other.allows_discovery(swap)


def test_domains_anchors_injectivity_and_fixed_target_set_are_all_enforced():
    matrix = weights(3, [])
    value = relation(matrix, matrix, (), anchors=((0, 0),),
                     domains=((0, (0, 1)), (1, (1,)), (2, (2,))))
    assert value.realize() == {0: 0, 1: 1, 2: 2}
    for invalid in ({0: 0, 1: 2, 2: 1}, {0: 1, 1: 0, 2: 2},
                    {0: 0, 1: 1, 2: 1}, {0: 0, 1: 1, 2: 3}, {0: 0, 1: 1},
                    ((0, 0), (0, 1), (2, 2)), {0: 0., 1: 1, 2: 2}):
        assert not value.allows_discovery(invalid)
        with pytest.raises(ValueError):
            value.realize(invalid)
    with pytest.raises(ValueError, match='anchor projections'):
        replace(value, anchors=())
    with pytest.raises(FrozenInstanceError):
        value.graph_floor = 1.
    with pytest.raises(ValueError, match='atom domain'):
        replace(value, witness=((0, 0), (1, 2), (2, 1)))


def test_active_target_edge_and_floor_are_both_required():
    matrix = weights(3, [(0, 1, 1.), (0, 2, .2)])
    value = relation(matrix, matrix, ((0, 1),), tolerance=1.,
                     target_edges=((0, 1), (0, 2)), floor=.2)
    swap = {0: 0, 1: 2, 2: 1}
    assert value.allows_discovery(swap)
    assert not replace(value, graph_floor=np.nextafter(.2, np.inf)).allows_discovery(swap)
    assert not replace(value, target_edges=((0, 1),)).allows_discovery(swap)


def test_constructor_retains_constraints_without_rematching_or_expanding_self(monkeypatch):
    matrix = weights(15, [(0, atom, 1.) for atom in range(1, 15)])
    source = relaxed_automorphism('C' * 15, matrix, tolerance=0.)
    def forbidden(*args, **kwargs):
        raise AssertionError('construction must not expand self choices')
    monkeypatch.setattr(type(source), 'realize', forbidden)
    value = RelaxedIsomorphism(source, source, tuple(enumerate(range(15))), (),
        tuple((atom, tuple(range(15))) for atom in range(15)),
        tuple(DiscoveryEdge((0, atom), 1., 0.) for atom in range(1, 15)),
        tuple((0, atom) for atom in range(1, 15)), .2)
    assert value.source is value.target is source
    assert len(value.accepted_edges) == 14 and len(value.witness) == 15
    assert value.allows_discovery(dict(enumerate(range(15))))


def test_compiled_triangle_tail_membership_matches_all_twenty_four_bijections():
    import z3
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints, compile_relaxed_relation
    matrix = weights(4, [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.), (0, 3, 1.)])
    value = relation(matrix, matrix, ((0, 1), (1, 2), (0, 3)))
    expected = {(0, 1, 2, 3), (0, 2, 1, 3), (1, 0, 3, 2),
                (2, 0, 3, 1), (1, 3, 0, 2), (2, 3, 0, 1)}
    for images in permutations(range(4)):
        system = FiniteConstraints()
        compiled = compile_relaxed_relation(system, value, 'triangle-tail')
        for atom, image in enumerate(images):
            system.restrict(compiled.images[atom][0], {image})
        status, model, _ = system.solve(2000)
        assert (status == z3.sat) == (images in expected)
        assert status != z3.unknown
        if status == z3.sat:
            assert compiled.realize(model) == dict(enumerate(images))
            choice = compiled.choices(model)
            assert tuple(choice['discovery'].values()) in discovery_oracle(value)
            assert choice['source_permutation'] in self_oracle(value.source)
            assert choice['target_permutation'] in self_oracle(value.target)


def test_compiler_does_not_turn_nontransitive_self_tolerance_into_a_group():
    import z3
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints, compile_relaxed_relation
    matrix = weights(4, [(0, 1, 1.), (0, 2, 1.6), (0, 3, 2.2)])
    value = relation(matrix, matrix, ((0, 1), (0, 2), (0, 3)), self_tolerance=.7)
    for tail in permutations((1, 2, 3)):
        images = (0, *tail)
        system = FiniteConstraints()
        compiled = compile_relaxed_relation(system, value, 'weighted-star')
        for atom, image in enumerate(images):
            system.restrict(compiled.images[atom][0], {image})
        status, model, _ = system.solve(2000)
        assert (status == z3.sat) == (images != (0, 3, 2, 1))
        assert status != z3.unknown
        if status == z3.sat:
            assert compiled.realize(model) == dict(enumerate(images))
            choices = compiled.choices(model)
            assert choices['discovery'] == dict(enumerate(range(4)))
            assert choices['source_permutation'] in self_oracle(value.source)
            assert choices['target_permutation'] in self_oracle(value.target)


def test_compiler_reuses_source_factor_choices_for_their_inverse():
    import z3
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints, compile_relaxed_relation
    matrix = weights(3, [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.)])
    value = relation(matrix, matrix, ((0, 1), (1, 2)))
    system = FiniteConstraints()
    compiled = compile_relaxed_relation(system, value, 'cycle')
    for atom, image in enumerate((2, 0, 1)):
        system.restrict(compiled.discovery[atom][0], {image})
    for symbols, images in ((compiled.source_action, (1, 2, 0)),
                            (compiled.target_action, (1, 0, 2))):
        for symbol, image in zip(symbols, images):
            system.restrict(symbol[0], {image})
    status, model, _ = system.solve(2000)
    assert status == z3.sat
    assert compiled.realize(model) == {0: 0, 1: 2, 2: 1}
    assert compiled.choices(model) == dict(discovery={0: 2, 1: 0, 2: 1},
                                          source_permutation=(1, 2, 0),
                                          target_permutation=(1, 0, 2))


def test_compiler_preserves_sparse_original_indices_and_fixed_anchors():
    import z3
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints, compile_relaxed_relation
    source = relaxed_automorphism('C' * 8, weights(8, [(2, 5, 1.), (2, 7, 1.1)]),
                                 tolerance=.2, atoms=(2, 5, 7), anchors=(2,))
    target = relaxed_automorphism('C' * 9, weights(9, [(1, 4, 1.1), (1, 8, 1.2)]),
                                 tolerance=.2, atoms=(1, 4, 8), anchors=(1,))
    value = RelaxedIsomorphism(source, target, ((2, 1), (5, 4), (7, 8)), ((2, 1),),
        ((2, (1,)), (5, (4, 8)), (7, (4, 8))),
        (DiscoveryEdge((2, 5), 1., .2),), ((1, 4),), .2)
    for images in permutations((1, 4, 8)):
        system = FiniteConstraints()
        compiled = compile_relaxed_relation(system, value, 'sparse')
        for atom, image in zip(value.source_atoms, images):
            system.restrict(compiled.images[atom][0], {image})
        status, model, _ = system.solve(2000)
        assert (status == z3.sat) == (images[0] == 1)
        assert status != z3.unknown
        if status == z3.sat:
            assert compiled.realize(model) == dict(zip(value.source_atoms, images))


def test_large_symmetric_relation_uses_factors_without_generating_group_elements(monkeypatch):
    from math import factorial
    from sympy.combinatorics import PermutationGroup
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints, compile_relaxed_relation
    def forbidden(*args, **kwargs):
        raise AssertionError('compiler must not enumerate complete group elements')
    monkeypatch.setattr(PermutationGroup, 'generate_schreier_sims', forbidden)
    monkeypatch.setattr(PermutationGroup, 'generate_dimino', forbidden)
    n = 11
    matrix = weights(n, [(0, atom, 1.) for atom in range(1, n)])
    value = relation(matrix, matrix, tuple((0, atom) for atom in range(1, n)))
    compiled = compile_relaxed_relation(FiniteConstraints(), value, 'ten-leaf-star')
    factor_size = sum(len(factor) for _, _, _, factor in compiled.program)
    assert factor_size <= n * (n + 1)
    assert factor_size < factorial(n - 1)
    assert set(compiled.images) == set(range(n))


@pytest.mark.parametrize('fixture', ['unrestricted', 'partial_edges', 'weighted', 'domains'])
def test_image_only_compilation_absorbs_only_discovery_closed_self_actions(fixture):
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints, compile_relaxed_relation
    from mappa.mappa_core.core.symmetry.event_counts import EventCountProgram
    if fixture == 'weighted':
        matrix = weights(4, [(0, 1, 1.), (0, 2, 1.6), (0, 3, 2.2)])
        value = relation(matrix, matrix, ((0, 1), (0, 2), (0, 3)), self_tolerance=.7)
        assert value.absorbed_self_actions == (False, False)
    else:
        matrix = weights(4, [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.), (0, 3, 1.)])
        accepted = ((0, 1), (1, 2), (0, 3)) if fixture == 'partial_edges' else ()
        domains = ((0, (0,)), (1, (1, 2, 3)), (2, (1, 2, 3)), (3, (1, 2, 3))) if fixture == 'domains' else None
        value = relation(matrix, matrix, accepted, domains=domains)
        if fixture == 'unrestricted':
            assert value.absorbed_self_actions == (True, True)
        if fixture == 'partial_edges':
            assert value.absorbed_self_actions == (False, True)
    expected = represented(value)
    system = FiniteConstraints()
    compiled = compile_relaxed_relation(system, value, fixture, images_only=True)
    actual = EventCountProgram(system).project(tuple(compiled.images[a][0] for a in value.source_atoms))
    assert set(actual) == expected
    assert set(actual.values()) == {0}
