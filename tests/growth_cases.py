"""Small graph inputs and direct native comparison helpers."""
import numpy as np
from mappa.mappa_core.core.matching.graphs import build_graph
from mappa.mappa_core.core.matching.growth import native

def graph(n, edges, elements=None):
    w = np.zeros((n,n))
    for edge in edges:
        a,b = edge[:2]
        w[a,b] = w[b,a] = edge[2] if len(edge)==3 else 1.
    return build_graph(elements or ('C',)*n,w,bond_cut=.2)


def rings():
    cycle = [(i, (i + 1) % 6, 1) for i in range(6)]
    source = graph(12, cycle + [(i, 6 + i, 1) for i in range(6)], ['C'] * 6 + ['H'] * 6)
    target = graph(18, cycle + [(i, 6 + 2 * i + h, 1) for i in range(6) for h in range(2)],
                   ['C'] * 6 + ['H'] * 12)
    return source, target


def ring_mappings():
    from itertools import product
    result = set()
    for offset, direction in product(range(6), (-1, 1)):
        carbons = [(offset + direction * i) % 6 for i in range(6)]
        for choices in product(range(2), repeat=6):
            images = carbons + [6 + 2 * p + h for p, h in zip(carbons, choices)]
            result.add(tuple(enumerate(images)))
    return result


def discovery_mappings(fragment):
    """Independent VF2 enumeration of tiny saved discovery relations."""
    import networkx as nx
    relation = fragment.relation
    source, target = nx.Graph(), nx.Graph()
    domains = dict(relation.domains)
    anchors = dict(relation.anchors)
    source.add_nodes_from((atom, {'atom': atom}) for atom in fragment.source_atoms)
    target.add_nodes_from((atom, {'atom': atom}) for atom in fragment.target_atoms)
    for edge in relation.accepted_edges:
        source.add_edge(*edge.atoms, weight=edge.weight,
                        tolerance=edge.tolerance, inclusive=edge.inclusive)
    positions = {atom: i for i, atom in enumerate(fragment.target_atoms)}
    for a, b in relation.target_edges:
        target.add_edge(a, b, weight=float(relation.target.weights[positions[a], positions[b]]))
    def node_match(left, right):
        p, r = left['atom'], right['atom']
        return p in domains[r] and (r not in anchors or anchors[r] == p)
    def edge_match(left, right):
        difference = abs(left['weight'] - right['weight'])
        allowed = (difference <= right['tolerance'] if right['inclusive']
                   else difference < right['tolerance'])
        return left['weight'] >= relation.graph_floor and allowed
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        target, source, node_match=node_match, edge_match=edge_match)
    return {tuple(sorted((r, p) for p, r in mapping.items()))
            for mapping in matcher.subgraph_monomorphisms_iter()}


def self_actions(value):
    """Read original weights directly, independently of stored generators."""
    import networkx as nx
    structure = nx.Graph()
    for i, (atom, element) in enumerate(zip(value.atoms, value.elements)):
        structure.add_node(i, element=element, anchor=atom if atom in value.anchors else None)
    for i in range(len(value.atoms)):
        for j in range(i + 1, len(value.atoms)):
            if value.weights[i, j] > value.bond_threshold:
                structure.add_edge(i, j)
    matcher = nx.algorithms.isomorphism.GraphMatcher(structure, structure,
        node_match=lambda left, right: left == right)
    from sympy.combinatorics import Permutation, PermutationGroup
    structural = {tuple(mapping[i] for i in range(len(value.atoms)))
                  for mapping in matcher.isomorphisms_iter()}
    group = PermutationGroup([Permutation(list(g)) for g in value.generators]
        or [Permutation(list(range(len(value.atoms))))])
    assert {tuple(g(i) for i in range(len(value.atoms)))
            for g in group.generate_schreier_sims()} == structural
    actions = []
    for permutation in structural:
        moved = value.weights[np.ix_(permutation, permutation)]
        allowed = bool(np.all(np.abs(value.weights - moved) <= value.tolerance))
        assert value.contains(permutation) == allowed
        if allowed:
            actions.append({atom: value.atoms[permutation[i]] for i, atom in enumerate(value.atoms)})
    return actions


def group_actions(generators, degree):
    """Enumerate only the small test oracle, never production group elements."""
    from sympy.combinatorics import Permutation, PermutationGroup
    group = PermutationGroup([Permutation(list(g), size=degree) for g in generators]
        or [Permutation(list(range(degree)))])
    return tuple(tuple(action(i) for i in range(degree))
                 for action in group.generate_schreier_sims())


def placed(mappings, generators):
    """Move a whole mapping jointly under one saved placement group."""
    if not generators:
        return set(mappings)
    return {tuple((r, action[p]) for r, p in mapping)
            for action in group_actions(generators, len(generators[0]))
            for mapping in mappings}


def represented_local(fragments):
    """Expand tiny relations by independent discovery and two-sided self oracles."""
    result = set()
    for fragment in fragments:
        source = self_actions(fragment.source_self_automorph)
        target = self_actions(fragment.target_self_automorph)
        for pairs in discovery_mappings(fragment):
            discovery = dict(pairs)
            for alpha in source:
                inverse = {image: atom for atom, image in alpha.items()}
                for beta in target:
                    result.add(tuple((atom, beta[discovery[inverse[atom]]])
                                     for atom in fragment.source_atoms))
    return result


def represented(value):
    """Expand the complete operation or path without losing global placement."""
    from mappa.mappa_core.core.matching.representation.records import FragmentMatchResult
    from mappa.mappa_core.core.search.graph import AAMSearchGraph, SearchPath
    if isinstance(value, FragmentMatchResult):
        return placed(represented_local(value.matches), value.placement_generators)
    if isinstance(value, (AAMSearchGraph, SearchPath)):
        result = set()
        paths = value.paths() if isinstance(value, AAMSearchGraph) else (value,)
        for path in paths:
            mappings = {tuple(sorted(path.mapping.items()))}
            for index in reversed(path.transitions):
                edge = path.graph.transitions[index]
                if edge.match is None:
                    mappings = placed(mappings, edge.placement_generators)
                    continue
                fragment = edge.match
                choices = represented_local((fragment,))
                moved = set()
                for choice in choices:
                    images = dict(choice)
                    action = {p: images[r] for r, p in fragment.witness}
                    moved.update(tuple((r, action.get(p, p)) for r, p in mapping)
                                 for mapping in mappings)
                mappings = placed(moved, edge.placement_generators)
            result.update(mappings)
        return result
    return represented_local(value)


def raw(source,target,seed,mapping,cut,*,replay=None,cap=100,tol=1.,islands=None,deferred=()):
    if replay is None:
        cut_graph = source.copy()
        cut_graph.remove_edges_from(cut)
        view = native.source_graph(cut_graph)
        result = native._engine.grow_island(view.graph,native.target_graph(target).graph,
            seed,mapping,.2,tol,1,cap,islands,deferred,False)
    else:
        result = replay.engine.grow(cut,seed,mapping,.2,tol,1,cap,islands,deferred,False)
        result.pop('reused_extensions')
    result['profile'].pop('completion_seconds')  # Wall time is not a scientific result.
    for fragment in result['fragments']:
        for side in ('source', 'target'):
            fragment[side]['weights'] = np.asarray(fragment[side]['weights']).tolist()
    return result


def saved_fragment(problem, *, source_atoms=None, target_atoms=None, witness=None,
                   anchors=(), domains=None, accepted_edges=(),
                   iso_tolerance=1., self_tolerance=1.):
    """Construct an actual fixed-image relation with explicit scientific constraints."""
    from mappa.mappa_core.core.matching.representation.records import Fragment
    from mappa.mappa_core.core.symmetry import DiscoveryEdge, RelaxedIsomorphism, relaxed_automorphism
    source_atoms = tuple(range(problem.source_atom_count)) if source_atoms is None else tuple(source_atoms)
    target_atoms = source_atoms if target_atoms is None else tuple(target_atoms)
    witness = tuple(zip(source_atoms, target_atoms)) if witness is None else tuple(sorted(dict(witness).items()))
    source = relaxed_automorphism(problem.reactant.elements, problem.reactant.wbo,
        atoms=source_atoms, anchors=tuple(r for r, _ in anchors), tolerance=self_tolerance)
    target = relaxed_automorphism(problem.product.elements, problem.product.wbo,
        atoms=target_atoms, anchors=tuple(p for _, p in anchors), tolerance=self_tolerance)
    if domains is None:
        domains = tuple((r, tuple(p for p in target_atoms
            if problem.reactant.elements[r] == problem.product.elements[p])) for r in source_atoms)
    edges = tuple(DiscoveryEdge((a, b), problem.reactant.wbo[a, b], iso_tolerance)
                  for a, b in accepted_edges)
    target_edges = tuple((p, q) for p in target_atoms for q in target_atoms
                        if p < q and problem.product.wbo[p, q] >= .2)
    return Fragment(RelaxedIsomorphism(source, target, witness, tuple(anchors),
        tuple(domains), edges, target_edges, .2))


