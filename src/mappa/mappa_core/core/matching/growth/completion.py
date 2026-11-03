"""Discover target atom sets and retain their complete local relationships."""
from dataclasses import replace
from time import perf_counter

import numpy as np

from mappa.mappa_core.core.matching.growth.existence import FragmentQuery
from mappa.mappa_core.core.matching.policy import _edge_wbo, ElementNodeMatchPolicy
from mappa.mappa_core.core.matching.representation.records import Fragment, IslandBranchLimitExceeded
from mappa.mappa_core.core.symmetry.relaxed import relaxed_automorphism
from mappa.mappa_core.core.symmetry.relations import DiscoveryEdge, RelaxedIsomorphism


def _prepare_self(graph, atoms, anchors, tolerance, threshold, labels):
    """Original internal WBOs define self relationships, including cut bonds."""
    atoms = tuple(sorted(atoms))
    anchors = tuple(anchors)
    positions = {atom: index for index, atom in enumerate(atoms)}
    weights = np.asarray([[_edge_wbo(graph, a, b) for b in atoms] for a in atoms])
    prepared = relaxed_automorphism(
        tuple(labels[a] for a in atoms), weights,
        tolerance=tolerance, bond_threshold=threshold,
        anchors=tuple(positions[a] for a in anchors))
    return replace(prepared, atoms=atoms, anchors=tuple(sorted(anchors)))


def _self_labels(source, target, node_policy):
    if type(node_policy) is ElementNodeMatchPolicy:
        return tuple({a: graph.nodes[a].get('element') for a in graph}
                     for graph in (source, target))
    # Custom shuffles must preserve every compatibility row, including later
    # assignments transported by an earlier fragment choice.
    rows = {r: tuple(node_policy.compatible(source, r, target, p) for p in sorted(target))
            for r in sorted(source)}
    columns = {p: tuple(node_policy.compatible(source, r, target, p) for r in sorted(source))
               for p in sorted(target)}
    result = []
    for graph, signatures in ((source, rows), (target, columns)):
        colors = {}
        result.append({a: str(colors.setdefault((node_policy.key(graph, a), signatures[a]), len(colors)))
                       for a in sorted(graph)})
    return tuple(result)


def _global_self(target, labels, tolerance):
    """Reuse one original-WBO preparation within its owning endpoint graph."""
    threshold = target.graph.get('bond_cut', .2)
    key = threshold, tolerance, tuple(sorted(labels.items()))
    weights = target.graph.get('wbo_matrix')
    cached = target.graph.get('_relaxed_global')
    if cached is not None and cached[0] is target and cached[1] is weights and cached[2] == key:
        prepared = cached[3]
        # SearchSession owns read-only endpoint matrices. Explicitly mutable
        # graph inputs need a value check because they have no mutation stamp.
        if weights is not None and not np.asarray(weights).flags.writeable:
            return prepared
        if weights is None:
            current = np.asarray([[_edge_wbo(target,a,b) for b in prepared.atoms]
                                  for a in prepared.atoms])
        else:
            current = weights if prepared.atoms == tuple(range(len(weights))) else weights[
                np.ix_(prepared.atoms, prepared.atoms)]
        if np.array_equal(current, prepared.weights):
            return prepared
    prepared = _prepare_self(target, target, (), tolerance, threshold, labels)
    target.graph['_relaxed_global'] = target, weights, key, prepared
    return prepared


def complete_fragment(growth, source, target, locked, *, tolerance, self_tolerance,
                      node_policy, branch_limit, bond_tolerance=None, profile=None):
    """One Fragment per image; internal mappings remain a joint constraint set."""
    started = perf_counter()
    query = FragmentQuery(source, target, growth.source_atoms, growth.accepted_edges,
                          locked, node_policy, tolerance, bond_tolerance)
    source_labels, target_labels = _self_labels(source, target, node_policy)
    global_self = _global_self(target, target_labels, self_tolerance) if query.needs_global_symmetry() else None
    witnesses = query.images(branch_limit + 1, preferred=growth.witness, global_self=global_self)
    if profile is not None:
        profile.update(completion_seconds=perf_counter() - started,
                       completion_queries=query.nodes, completion_pruned=query.rejected,
                       completion_global_generators=query.global_generators,
                       completion_global_images=query.global_images,
                       completion_global_pruned=query.global_pruned,
                       completion_generated=len(witnesses), completion_peak_candidates=len(witnesses))
    if not witnesses:
        raise RuntimeError('completed fragment lost its proven matching witness')
    if len(witnesses) > branch_limit:
        raise IslandBranchLimitExceeded(len(witnesses), branch_limit)
    anchors = tuple(sorted((r, locked[r]) for r in query.atoms if r in locked))
    source_self = _prepare_self(source, query.atoms, dict(anchors), self_tolerance,
                                source.graph.get('bond_cut', .2), source_labels)
    edges = tuple(DiscoveryEdge((a, b), _edge_wbo(source, a, b),
                tolerance if bond_tolerance is None else bond_tolerance.limit(a, b),
                bond_tolerance is None or bond_tolerance.inclusive)
                  for a, b in sorted(growth.accepted_edges))
    fragments = []
    for witness in witnesses:
        image = frozenset(witness.values())
        domains = {r: tuple(p for p in query.domains[r] if p in image) for r in query.atoms}
        if not query._reduce(domains):
            raise RuntimeError('a proven target image lost its discovery support')
        target_self = _prepare_self(target, image, (p for _, p in anchors),
                                    self_tolerance, target.graph.get('bond_cut', .2), target_labels)
        relation = RelaxedIsomorphism(source_self, target_self, tuple(witness.items()), anchors,
            tuple(domains.items()),
            edges, tuple((a, b) for a, b in target.edges if a in image and b in image),
            target.graph.get('bond_cut', .2))
        fragments.append(Fragment(relation))
    if profile is not None:
        profile['completion_seconds'] = perf_counter() - started
    return tuple(fragments)
