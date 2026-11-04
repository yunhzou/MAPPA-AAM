"""Contract correlated event costs without expanding complete atom mappings.

Each factor describes a local constraint or bond response. Eliminating a choice
keeps the least cost for each remaining shared choice. Requested endpoint images
remain in the result, so sweep assembly uses the same internal correlation.
"""
from dataclasses import dataclass
from itertools import combinations, product
from math import log2

import numpy as np

from .constraints import Choice


@dataclass(frozen=True)
class _Factor:
    variables: tuple[int, ...]
    rows: np.ndarray
    costs: np.ndarray

    def reduce(self, variables):
        columns = [self.variables.index(v) for v in variables]
        return _reduce(variables, self.rows[:, columns], self.costs)


def _reduce(variables, rows, costs):
    """One minimum per shared assignment; discarded choices are minimized out."""
    variables = tuple(variables)
    rows = np.asarray(rows, dtype=np.int64).reshape(len(costs), len(variables))
    costs = np.asarray(costs, dtype=np.int64)
    if not len(costs):
        return _Factor(variables, rows, costs)
    if not variables:
        return _Factor((), np.empty((1, 0), dtype=np.int64), np.array([costs.min()]))
    unique, inverse = np.unique(rows, axis=0, return_inverse=True)
    minimum = np.full(len(unique), np.iinfo(np.int64).max, dtype=np.int64)
    np.minimum.at(minimum, inverse, costs)
    return _Factor(variables, unique, minimum)


def _keys(rows, columns):
    selected = np.ascontiguousarray(rows[:, columns])
    return selected.view(np.dtype((np.void, selected.dtype.itemsize * len(columns)))).ravel()


def _join(left, right):
    """Sparse natural join with addition of costs, performed on array indices."""
    variables = tuple(sorted(set(left.variables) | set(right.variables)))
    if not len(left.costs) or not len(right.costs):
        return _Factor(variables, np.empty((0, len(variables)), dtype=np.int64),
                       np.empty(0, dtype=np.int64))
    shared = tuple(v for v in left.variables if v in right.variables)
    if shared:
        a = _keys(left.rows, [left.variables.index(v) for v in shared])
        b = _keys(right.rows, [right.variables.index(v) for v in shared])
        order_a, order_b = np.argsort(a), np.argsort(b)
        keys_a, starts_a, lengths_a = np.unique(a[order_a], return_index=True, return_counts=True)
        keys_b, starts_b, lengths_b = np.unique(b[order_b], return_index=True, return_counts=True)
        _, ia, ib = np.intersect1d(keys_a, keys_b, assume_unique=True, return_indices=True)
        lengths = lengths_a[ia] * lengths_b[ib]
        blocks = np.repeat(np.arange(len(lengths)), lengths)
        offsets = np.arange(int(lengths.sum())) - np.repeat(np.cumsum(lengths) - lengths, lengths)
        rows_a = order_a[starts_a[ia][blocks] + offsets // lengths_b[ib][blocks]]
        rows_b = order_b[starts_b[ib][blocks] + offsets % lengths_b[ib][blocks]]
    else:
        rows_a = np.repeat(np.arange(len(left.costs)), len(right.costs))
        rows_b = np.tile(np.arange(len(right.costs)), len(left.costs))
    rows = np.empty((len(rows_a), len(variables)), dtype=np.int64)
    for column, variable in enumerate(variables):
        if variable in left.variables:
            rows[:, column] = left.rows[rows_a, left.variables.index(variable)]
        else:
            rows[:, column] = right.rows[rows_b, right.variables.index(variable)]
    return _Factor(variables, rows, left.costs[rows_a] + right.costs[rows_b])


class EventCountProgram:
    """Reusable min-sum response for one finite correlated fragment program.

    ``costs`` contains ``(left, right, pair_costs)`` records. Values are integers
    or ``Choice`` objects from ``system``; pair costs cover their domains.
    ``project`` retains requested endpoint values and minimizes only internal
    choices. Empty support means the original joint constraints are infeasible.
    """
    def __init__(self, system, costs=()):
        system.propagate()
        self.system = system
        self.factors = []
        self._projections = {}
        self._failed = system.failed
        if self._failed:
            return
        for kind, values, data in system.constraints:
            if kind in ('table', 'action'):
                rows = system.action_rows(values, data) if kind == 'action' else data
                self._add(values, rows)
            elif kind == 'distinct':
                for left, right in combinations(values, 2):
                    xs, ys = system.domain(left), system.domain(right)
                    if xs.isdisjoint(ys):
                        continue
                    self._add((left, right), ((x, y) for x in xs for y in ys if x != y))
            elif kind == 'lookup':
                selector, output = values
                for selected in system.domain(selector):
                    entry = data[selected]
                    # Conditional equalities avoid a factor containing every
                    # selectable entry, which would form a needlessly large clique.
                    self._add((selector, output, entry),
                        ((s, out, value) for s, out, value in product(
                            system.domain(selector), system.domain(output), system.domain(entry))
                         if s != selected or out == value))
            else:
                raise ValueError(f'Unknown finite constraint: {kind}')
        for left, right, entries in costs:
            rows = tuple(product(sorted(system.domain(left)), sorted(system.domain(right))))
            missing = set(rows).difference(entries)
            if missing:
                raise ValueError('Pair costs must cover both endpoint domains')
            self._add((left, right), rows, [entries[row] for row in rows])
        self.factors = tuple(self.factors)

    def _add(self, values, rows, costs=None):
        rows = tuple(rows)
        if costs is None:
            costs = np.zeros(len(rows), dtype=np.int64)
        rows = np.asarray(rows, dtype=np.int64).reshape(len(rows), len(values))
        costs = np.asarray(costs, dtype=np.int64)
        valid = np.ones(len(rows), dtype=bool)
        columns = {}
        for column, value in enumerate(values):
            domain = self.system.domain(value)
            valid &= np.isin(rows[:, column], tuple(domain))
            if isinstance(value, Choice) and len(domain) > 1:
                if value.index in columns:
                    valid &= rows[:, column] == rows[:, columns[value.index]]
                else:
                    columns[value.index] = column
        variables = tuple(sorted(columns))
        factor = _reduce(variables, rows[valid][:, [columns[v] for v in variables]], costs[valid])
        if not len(factor.costs):
            self._failed = True
        if len(factor.costs) and np.all(factor.costs == factor.costs[0]):
            size = 1
            for variable in variables:
                size *= len(self.system.domains[variable])
            if len(factor.costs) == size:
                factor = factor.reduce(())
                if factor.costs[0] == 0:
                    return
        self.factors.append(factor)

    def project(self, outputs=()):
        """Return ``{output_tuple: minimum_event_count}`` for requested images."""
        outputs = tuple(outputs)
        if outputs in self._projections:
            return self._projections[outputs]
        if self._failed:
            return {}
        retained = {v.index for v in outputs
                    if isinstance(v, Choice) and len(self.system.domain(v)) > 1}
        factors = list(self.factors)
        present = {v for factor in factors for v in factor.variables}
        for variable in retained - present:
            domain = sorted(self.system.domains[variable])
            factors.append(_Factor((variable,), np.asarray(domain, dtype=np.int64)[:, None],
                                   np.zeros(len(domain), dtype=np.int64)))
        remaining = {v for factor in factors for v in factor.variables} - retained
        while remaining:
            def order(variable):
                adjacent = {v for factor in factors if variable in factor.variables for v in factor.variables}
                width = sum(log2(len(self.system.domains[v])) for v in adjacent if v != variable)
                return width, len(adjacent), variable
            variable = min(remaining, key=order)
            bucket = [factor for factor in factors if variable in factor.variables]
            factors = [factor for factor in factors if variable not in factor.variables]
            bucket.sort(key=lambda factor: len(factor.costs))
            combined = bucket.pop(0)
            while bucket:
                selected = min(range(len(bucket)), key=lambda i: (
                    -len(set(combined.variables).intersection(bucket[i].variables)), len(bucket[i].costs)))
                combined = _join(combined, bucket.pop(selected))
                if not len(combined.costs):
                    self._projections[outputs] = {}
                    return {}
            factors.append(combined.reduce(tuple(v for v in combined.variables if v != variable)))
            remaining.remove(variable)
        combined = _Factor((), np.empty((1, 0), dtype=np.int64), np.array([0], dtype=np.int64))
        for factor in sorted(factors, key=lambda factor: len(factor.costs)):
            combined = _join(combined, factor)
        positions = {v: i for i, v in enumerate(combined.variables)}
        result = {}
        for row, cost in zip(combined.rows, combined.costs):
            key = tuple(int(row[positions[v.index]]) if isinstance(v, Choice) and v.index in positions
                        else int(next(iter(self.system.domain(v)))) for v in outputs)
            result[key] = int(cost)
        self._projections[outputs] = result
        return result
