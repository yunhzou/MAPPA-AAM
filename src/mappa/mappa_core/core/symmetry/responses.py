"""Shared original-weight bond responses under correlated structural actions."""
from collections import OrderedDict
from dataclasses import dataclass
from functools import cached_property

import numpy as np

from .actions import ActionIndex


def response_bounds(values, weight, threshold):
    """Encode (-, zero, +) responses of sorted opposite weights with two cuts.

    Keep the original subtraction/comparisons, including floating boundaries.
    Algebraically moving the threshold to the other operand need not agree.
    """
    values, weight = np.asarray(values), np.asarray(weight)
    bounds = []
    # Search all distinct weights together, without a weights-by-weights table.
    for negative in (True, False):
        lo = np.zeros(weight.shape, dtype=int)
        hi = np.full(weight.shape, len(values), dtype=int)
        while np.any(lo < hi):
            mid = (lo + hi) // 2
            delta = values[np.minimum(mid, len(values) - 1)] - weight
            left = delta < -threshold if negative else delta <= threshold
            active = lo < hi
            lo = np.where(active & left, mid + 1, lo)
            hi = np.where(active & ~left, mid, hi)
        bounds.append(lo)
    return tuple(int(bound) if bound.ndim == 0 else bound for bound in bounds)


def bond_categories(source, target, threshold):
    """0 unchanged; 1..4 are exclusive categories using the existing bond cutoff."""
    source, target = np.asarray(source), np.asarray(target)
    decrease, increase = source - target > threshold, target - source > threshold
    return np.where(decrease, np.where(target > threshold, 4, 1),
                    np.where(increase, np.where(source > threshold, 3, 2), 0)).astype(np.uint8)


@dataclass(frozen=True, eq=False)
class BondResponse:
    """Shared responses to reachable bond weights, not independent choices."""
    source_weight: float
    threshold: float
    target_weights: np.ndarray
    categories: np.ndarray

    @cached_property
    def constant(self):
        return int(self.categories[0]) if np.all(self.categories == self.categories[0]) else None

    @cached_property
    def always_event(self):
        return bool(self.categories.size and self.categories.all())

    def lookup(self, weights):
        return self.categories[np.searchsorted(self.target_weights, weights)]


class BondResponseIndex(ActionIndex):
    """Reuse weight responses and pair projections for one fixed target matrix."""
    def __init__(self, weights):
        weights = np.asarray(weights)
        if (weights.ndim != 2 or weights.shape[0] != weights.shape[1]
                or not np.all(np.isfinite(weights))
                or not np.array_equal(weights, weights.T) or np.any(weights < 0)):
            raise ValueError('bond weights must be a finite nonnegative symmetric square matrix')
        super().__init__(len(weights))
        self.p = weights
        self.padded_product = np.pad(weights, ((0, 1), (0, 1)))
        self._event_tables = OrderedDict()
        self._domain_responses = OrderedDict()
        self._response_values = OrderedDict()

    def response_from_weights(self, weight, tau, weights):
        key = float(weight), float(tau), tuple(weights)
        if key not in self._response_values:
            codes = bond_categories(weight, weights, tau)
            weights.setflags(write=False)
            codes.setflags(write=False)
            self._response_values[key] = BondResponse(float(weight), float(tau), weights, codes)
            if len(self._response_values) > 4096:
                self._response_values.popitem(last=False)
        return self._response_values[key]

    def domain_response(self, weight, tau, xs, ys):
        """Certify constants on a superset before exploring any pair orbit."""
        xs, ys = sorted((xs, ys))
        key = float(weight), float(tau), xs, ys
        if key not in self._domain_responses:
            values = set()
            x, y = np.asarray(xs), np.asarray(ys)
            block = max(1, 65536 // len(y))
            for start in range(0, len(x), block):
                xx = x[start:start+block]
                valid = (xx[:,None] != y[None,:]) | (xx[:,None] == -1) | (y[None,:] == -1)
                values.update(np.unique(self.padded_product[np.ix_(xx,y)][valid]).tolist())
            self._domain_responses[key] = self.response_from_weights(weight, tau, np.array(sorted(values)))
            if len(self._domain_responses) > 4096:
                self._domain_responses.popitem(last=False)
        return self._domain_responses[key]

    def refine_response(self, weight, tau, pair, actions, response):
        """Refine structural pair orbits; relation responses keep bounded domains.

        Relation domains only prepare weight lookups. Their joint permissions
        remain in the saved constraints and are checked at choice evaluation.
        """
        if response.constant is not None or any(kind == 'relation' for kind, _ in actions):
            return response, 0
        key = float(weight), float(tau), pair, actions
        if key not in self._event_tables:
            pairs = self.endpoint_pairs(pair, actions)
            endpoints = np.asarray(pairs, dtype=int)
            weights = np.unique(self.padded_product[endpoints[:,0], endpoints[:,1]])
            self._event_tables[key] = self.response_from_weights(weight, tau, weights), len(pairs)
            if len(self._event_tables) > 16384:
                self._event_tables.popitem(last=False)
        return self._event_tables[key]
