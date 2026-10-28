"""Conditional bond tolerances shared by Python and native fragment matching."""
from dataclasses import dataclass
from collections.abc import Callable
import math


@dataclass(frozen=True)
class BondTolerance:
    """A constant or ``value(source_graph, a, b)`` for each source bond.

    Conditions may use element, node features, or bond attributes. The rule
    must be symmetric in its endpoints and deterministic. It is evaluated
    once per active source bond, never per candidate mapping. ``inclusive``
    selects <= (ordinary growth) or < (unchanged under an event threshold).
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
        from .primitives import _growth_edge_supported
        return _growth_edge_supported(source_weight, target_weight, self.limit(a, b),
                                      floor, inclusive=self.inclusive)

    def native_values(self, nodes):
        # Native indexing can differ from the graph's original atom indices.
        index = {atom: i for i, atom in enumerate(nodes)}
        return [(index[a], index[b], value) for (a, b), value in self.values.items()]
