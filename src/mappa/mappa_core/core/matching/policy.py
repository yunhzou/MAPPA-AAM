"""Node compatibility and bond-tolerance rules for matching."""
from __future__ import annotations

from collections.abc import Callable, Hashable, Sequence
from dataclasses import dataclass
from typing import Any
from collections.abc import Callable
import math


def _node_attr(g, node, name, default=None):
    data = g.nodes[node]
    if name in data:
        return data.get(name, default)
    features = data.get("features")
    if isinstance(features, dict):
        return features.get(name, default)
    return default


class ElementNodeMatchPolicy:
    """Default chemistry policy: nodes are compatible iff elements match."""

    def compatible(self, g_query, query_node, g_target, target_node) -> bool:
        return (
            _node_attr(g_query, query_node, "element")
            == _node_attr(g_target, target_node, "element")
        )

    def key(self, g, node) -> Hashable:
        return _node_attr(g, node, "element")


@dataclass(frozen=True)
class AttributeNodeMatchPolicy:
    """Match nodes by one or more node attributes or feature entries."""

    fields: tuple[str, ...]

    def __init__(self, fields: str | Sequence[str]):
        if isinstance(fields, str):
            fields = (fields,)
        object.__setattr__(self, "fields", tuple(fields))
        if not self.fields:
            raise ValueError("AttributeNodeMatchPolicy requires at least one field")

    def key(self, g, node) -> Hashable:
        values = tuple(_node_attr(g, node, field) for field in self.fields)
        return values[0] if len(values) == 1 else values

    def compatible(self, g_query, query_node, g_target, target_node) -> bool:
        return self.key(g_query, query_node) == self.key(g_target, target_node)


@dataclass(frozen=True)
class CallableNodeMatchPolicy:
    """Policy wrapper for user-supplied node compatibility functions."""

    node_match: Callable[[Any, int, Any, int], bool]
    node_key: Callable[[Any, int], Hashable] | None = None

    def compatible(self, g_query, query_node, g_target, target_node) -> bool:
        return bool(self.node_match(g_query, query_node, g_target, target_node))

    def key(self, g, node) -> Hashable:
        if self.node_key is not None:
            return self.node_key(g, node)
        return DEFAULT_NODE_POLICY.key(g, node)


DEFAULT_NODE_POLICY = ElementNodeMatchPolicy()


def as_node_match_policy(policy=None, *, node_match=None, node_key=None):
    """Normalize policy/callback input to a policy object.

    Callers may pass:
    - ``None`` for the default same-element policy
    - an object with ``compatible(...)`` and ``key(...)``
    - a callable ``(g_query, i, g_target, j) -> bool``
    - explicit ``node_match=...`` and optional ``node_key=...``
    """
    if policy is not None and node_match is not None:
        raise ValueError("pass either policy or node_match, not both")
    if policy is None and node_match is None:
        if node_key is not None:
            raise ValueError("node_key requires a custom node_match or policy")
        return DEFAULT_NODE_POLICY
    if policy is None:
        return CallableNodeMatchPolicy(node_match, node_key)
    if isinstance(policy, str):
        if node_key is not None:
            raise ValueError("node_key is not used with string attribute policy")
        return AttributeNodeMatchPolicy(policy)
    if isinstance(policy, (tuple, list)):
        if node_key is not None:
            raise ValueError("node_key is not used with attribute policy")
        return AttributeNodeMatchPolicy(tuple(policy))
    if callable(policy) and not hasattr(policy, "compatible"):
        return CallableNodeMatchPolicy(policy, node_key)
    if not hasattr(policy, "compatible") or not hasattr(policy, "key"):
        raise TypeError("node policy must define compatible(...) and key(...)")
    if node_key is not None:
        raise ValueError("node_key is not used when policy object is supplied")
    return policy


@dataclass(frozen=True)
class BondTolerance:
    """A constant or ``value(source_graph, a, b)`` for each source bond.

    Conditions may use element, node features, or bond attributes. The rule
    must be symmetric in its endpoints and deterministic. It is evaluated
    once per active source bond, never per candidate mapping. ``inclusive``
    selects <= (default) or a strict < comparison.
    Target bonds must still satisfy the matcher's existing graph floor.
    """
    value: float | Callable
    inclusive: bool = True

    def resolve(self, source):
        values = {}
        for a, b in source.edges:
            a, b = sorted((a, b))
            value = float(self.value(source, a, b) if callable(self.value) else self.value)
            if not math.isfinite(value) or value < 0:
                raise ValueError('bond tolerances must be finite and nonnegative')
            values[a, b] = value
        return ResolvedBondTolerance(values, self.inclusive)


@dataclass(frozen=True)
class ResolvedBondTolerance:
    values: dict
    inclusive: bool

    def limit(self, a, b):
        return self.values[min(a, b), max(a, b)]

    def supports(self, a, b, source_weight, target_weight, floor):
        return _growth_edge_supported(source_weight, target_weight, self.limit(a, b),
                                      floor, inclusive=self.inclusive)

    def native_values(self, nodes):
        # Native indexing can differ from the graph's original atom indices.
        index = {atom: i for i, atom in enumerate(nodes)}
        return [(index[a], index[b], value) for (a, b), value in self.values.items()]


SYM_SUPPORT_MAX_STATES = 4096


def _edge_wbo(g, a, b):
    """WBO for any atom pair in the complete weighted graph."""
    if a == b:
        return 0.0
    mat = g.graph.get("wbo_matrix")
    if mat is not None:
        return float(mat[a, b])
    if g.has_edge(a, b):
        return float(g[a][b].get("wbo", 0.0))
    return 0.0


def _growth_edge_supported(w_R, w_P, iso_tol, graph_floor=0.2, *, inclusive=True):
    """Compatibility for one active R-side growth/vector edge.

    The popped edge is not special chemically; it is one active R-side pair in
    the weighted WBO vector checked during extension.  If an R pair is active,
    the target pair must also be an active graph edge before the WBO tolerance
    is considered.  This prevents a loose ``iso_tol`` from treating a real
    bond as matching a target nonbond.
    """
    w_P = float(w_P)
    difference = abs(float(w_R) - w_P)
    return w_P >= graph_floor and (difference <= iso_tol if inclusive else difference < iso_tol)


def _wbo_bucket(w):
    return int(round(float(w) * 5))


def _orbit_id(orbits, node):
    return orbits[node] if orbits is not None else node
