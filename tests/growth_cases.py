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


def represented(candidates, target, tolerance):
    """Expand tiny test families independently of the production decoder."""
    from itertools import permutations, product
    from sympy.combinatorics import Permutation, PermutationGroup
    from mappa.mappa_core.core.matching.symmetry.workspace import SymmetryWorkspace

    workspace = SymmetryWorkspace(target, tolerance)
    result = set()
    for candidate in candidates:
        atoms = sorted(candidate.mapping)
        generators = (candidate.target_generators if candidate.target_generators is not None
                      else workspace.conditioned_generators(candidate, ()))
        group = PermutationGroup([Permutation(list(g)) for g in generators]
                                 or [Permutation(list(range(len(target))))])
        choices = [list(permutations(block.p_atoms, len(block.r_atoms)))
                   for block in candidate.blocks]
        for selection in product(*choices):
            mapping = dict(candidate.mapping)
            for block, images in zip(candidate.blocks, selection):
                mapping.update(zip(block.r_atoms, images))
            for action in group.generate_schreier_sims():
                result.add(tuple((r, action(mapping[r])) for r in atoms))
    return result


def raw(source,target,orbits,seed,mapping,cut,*,replay=None,cap=100,tol=1.,islands=None,deferred=()):
    if replay is None:
        cut_graph = source.copy()
        cut_graph.remove_edges_from(cut)
        view = native.source_graph(cut_graph)
        result = native._engine.grow_island(view.graph,native.target_graph(target,orbits).graph,
            seed,mapping,.2,tol,1,cap,islands,deferred,False)
    else:
        result = replay.engine.grow(cut,seed,mapping,.2,tol,1,cap,islands,deferred,False)
        result.pop('reused_extensions')
    result['profile'].pop('completion_seconds')  # Wall time is not a scientific result.
    return result


