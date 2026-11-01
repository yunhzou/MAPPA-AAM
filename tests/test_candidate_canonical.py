import networkx as nx
from sympy.combinatorics import Permutation, PermutationGroup

from mappa.mappa_core.core.matching.symmetry.canonical import _CandidateAutomorphismCanonicalizer
from mappa.mappa_core.core.matching.symmetry.orbits import _nauty_orbits


def _carbon_path():
    graph = nx.Graph()
    graph.add_nodes_from((atom, {"element": "C"}) for atom in range(4))
    graph.add_edges_from(
        (left, right, {"wbo": 1.0})
        for left, right in ((0, 1), (1, 2), (2, 3)))
    return graph


def _mapping_orbit(canonicalizer, mapping):
    generators = canonicalizer.atom_generators({})
    group = PermutationGroup([Permutation(list(g)) for g in generators]
                             or [Permutation(list(range(len(canonicalizer.nodes))))])
    return {tuple((r, action(p)) for r, p in sorted(mapping.items()))
            for action in group.generate_schreier_sims()}


def test_global_action_preserves_semantic_roles():
    canonicalizer = _CandidateAutomorphismCanonicalizer(_carbon_path())
    assert _mapping_orbit(canonicalizer, {10: 0, 11: 1}) == {
        ((10, 0), (11, 1)), ((10, 3), (11, 2))}


def test_global_action_does_not_swap_distinct_source_roles():
    canonicalizer = _CandidateAutomorphismCanonicalizer(_carbon_path())
    assert ((10, 2), (11, 3)) not in _mapping_orbit(canonicalizer, {10: 0, 11: 1})
    coloring = canonicalizer._colored_vertices(
        None, candidate_roles={0: (('mapped', 10),), 1: (('mapped', 11),)})
    actions = canonicalizer.atom_generators(None, colored_vertices=coloring)
    assert all(action[0] == 0 and action[1] == 1 for action in actions)


def test_candidate_canonicalizer_reuses_product_base(monkeypatch):
    from mappa.mappa_core.core.matching.symmetry import canonical

    graph = _carbon_path()
    orbits = _nauty_orbits(graph)
    original = canonical._nauty_colored_wbo_graph
    calls = []

    def counted(*args, **kwargs):
        calls.append((args, kwargs))
        return original(*args, **kwargs)

    monkeypatch.setattr(canonical, '_nauty_colored_wbo_graph', counted)
    first = _CandidateAutomorphismCanonicalizer(
        graph, p_orbits=orbits, locked_mapping={10: 0})
    second = _CandidateAutomorphismCanonicalizer(
        graph, p_orbits=orbits, locked_mapping={10: 1})

    assert len(calls) == 1
    assert first.adjacency is second.adjacency
    assert first.locked_roles != second.locked_roles


def test_candidate_canonicalizer_reuses_explicit_operation_cache(monkeypatch):
    from mappa.mappa_core.core.matching.symmetry import canonical

    graph = _carbon_path()
    original = canonical._nauty_colored_wbo_graph
    calls = []

    def counted(*args, **kwargs):
        calls.append((args, kwargs))
        return original(*args, **kwargs)

    monkeypatch.setattr(canonical, '_nauty_colored_wbo_graph', counted)
    cache = {}
    first = _CandidateAutomorphismCanonicalizer(
        graph, locked_mapping={10: 0}, base_cache=cache)
    second = _CandidateAutomorphismCanonicalizer(
        graph, locked_mapping={10: 1}, base_cache=cache)

    assert len(calls) == 1
    assert first.adjacency is second.adjacency
