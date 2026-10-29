"""Shared graph preparation for fragment and subgraph matching."""
from __future__ import annotations

import networkx as nx
from mappa.frag import WeightedGraph


def _coerce_graph(graph, graph_floor):
    if isinstance(graph, WeightedGraph):
        return graph.to_networkx(bond_cut=graph_floor)
    if isinstance(graph, nx.Graph):
        if "wbo_matrix" not in graph.graph:
            raise ValueError("NetworkX graph must carry graph['wbo_matrix']")
        if "bond_cut" not in graph.graph:
            graph = graph.copy()
            graph.graph["bond_cut"] = float(graph_floor)
        return graph
    if hasattr(graph, "to_networkx"):
        return graph.to_networkx(bond_cut=graph_floor)
    raise TypeError("graph must be WeightedGraph or a NetworkX graph")
