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


