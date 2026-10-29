"""Weighted graph representation and preparation shared by matching and search."""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any
import numpy as np
import networkx as nx

@dataclass
class WeightedNode:
    """Node record for generalized weighted-graph matching."""

    element: str | None = None
    label: str | None = None
    features: dict[str, Any] = field(default_factory=dict)
    attrs: dict[str, Any] = field(default_factory=dict)

    def as_attrs(self):
        out = dict(self.attrs)
        if self.element is not None:
            out["element"] = self.element
        if self.label is not None:
            out["label"] = self.label
        out["features"] = dict(self.features)
        for key, value in self.features.items():
            out.setdefault(key, value)
        return out


@dataclass
class WeightedGraph:
    """Typed weighted graph used by generalized subgraph matching."""

    nodes: list[WeightedNode | Mapping[str, Any] | str]
    weights: np.ndarray
    weight_name: str = "wbo"
    coords: np.ndarray | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_networkx(self, bond_cut=0.5):
        return build_weighted_graph(
            self.nodes, self.weights, bond_cut=bond_cut,
            weight_name=self.weight_name, coords=self.coords,
            metadata=self.metadata)


def _node_attrs(node):
    if isinstance(node, WeightedNode):
        return node.as_attrs()
    if isinstance(node, str):
        return {"element": node, "features": {}}
    if isinstance(node, Mapping):
        out = dict(node)
        features = dict(out.get("features") or {})
        out["features"] = features
        for key, value in features.items():
            out.setdefault(key, value)
        return out
    raise TypeError(f"unsupported weighted graph node: {node!r}")


def build_weighted_graph(nodes, weights, bond_cut=0.5, *,
                         weight_name="wbo", coords=None, metadata=None):
    """Build a NetworkX graph from arbitrary node descriptors and weights."""
    weights = np.asarray(weights, dtype=float)
    if weights.ndim != 2 or weights.shape[0] != weights.shape[1]:
        raise ValueError("weights must be a square matrix")
    if weights.shape[0] != len(nodes):
        raise ValueError("number of nodes must match weights shape")
    g = nx.Graph()
    g.graph["wbo_matrix"] = weights
    g.graph["bond_cut"] = float(bond_cut)
    g.graph["weight_name"] = str(weight_name)
    if coords is not None:
        g.graph["coords"] = np.asarray(coords, dtype=float)
    g.graph.update(dict(metadata or {}))
    for i, node in enumerate(nodes):
        g.add_node(i, **_node_attrs(node))
    n = len(nodes)
    for i in range(n):
        for j in range(i + 1, n):
            if weights[i, j] >= bond_cut:
                w = float(weights[i, j])
                g.add_edge(i, j, **{weight_name: w, "wbo": w})
    return g


def build_graph(elements, wbo, bond_cut=0.5):
    """Connectivity graph with element on each node and WBO weight on each
    edge. Bond exists iff WBO >= bond_cut; the full WBO matrix is retained on
    the graph so scoring and trace diagnostics can inspect exact WBO values."""
    g = nx.Graph()
    g.graph["wbo_matrix"] = np.asarray(wbo, dtype=float)
    g.graph["bond_cut"] = float(bond_cut)
    for i, e in enumerate(elements):
        g.add_node(i, element=e)
    n = len(elements)
    for i in range(n):
        for j in range(i + 1, n):
            if wbo[i, j] >= bond_cut:
                g.add_edge(i, j, wbo=float(wbo[i, j]))
    return g


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
