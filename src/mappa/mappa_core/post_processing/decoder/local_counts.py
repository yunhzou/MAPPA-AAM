"""Lazy minima from shared fragment responses and their original joint choices.

Only unresolved local queries reach the solver. Integer action functions are
parsed in bulk, avoiding a Python Boolean object for every input/output pair.
"""
from collections import defaultdict

import numpy as np

from mappa.mappa_core.core.symmetry.constraints import (
    Choice, FiniteConstraints, compile_relaxed_relation)


def count_bounds(index, response):
    """Independent-term lower bound and the saved valid witness's count."""
    mapping = dict(response.mapping)
    lower = sum(bool(np.all(term.response.categories != 0)) for term in response.terms)
    upper = sum(abs(float(index.p[mapping[a], mapping[b]]) - term.response.source_weight)
                > term.response.threshold for term in response.terms for a, b in (term.atoms,))
    return int(lower), int(upper)


def _integer(value):
    value = int(value)
    return str(value) if value >= 0 else f'(- {-value})'


def _combine(operator, values):
    values = tuple(values)
    if not values:
        return {'and': 'true', 'or': 'false', '+': '0'}[operator]
    return values[0] if len(values) == 1 else f'({operator} {" ".join(values)})'


def _member(expression, domain):
    values = tuple(sorted(domain))
    if len(values) > 2 and values[-1] - values[0] + 1 == len(values):
        return f'(and (<= {_integer(values[0])} {expression}) (<= {expression} {_integer(values[-1])}))'
    return _combine('or', (f'(= {expression} {_integer(v)})' for v in values))


def _lookup(expression, entries):
    groups = defaultdict(list)
    for value, image in entries.items():
        groups[image].append(value)
    ordered = sorted(groups.items(), key=lambda row: len(row[1]), reverse=True)
    result = ordered[0][0]
    for image, values in ordered[1:]:
        result = f'(ite {_member(expression, values)} {image} {result})'
    return result


def _matrix(left, right, xs, ys, rows):
    return _lookup(left, {x: _lookup(right, {y: _integer(rows[i, j]) for j, y in enumerate(ys)})
                          for i, x in enumerate(xs)})


class _IntegerEncoding:
    def __init__(self, system):
        self.system = system
        self.lines = []
        self.factors = {}
        for i, domain in enumerate(system.domains):
            if len(domain) > 1:
                self.lines.append(f'(declare-const v{i} Int)')
                self.assertion(_member(f'v{i}', domain))

    def value(self, value):
        domain = self.system.domain(value)
        return _integer(next(iter(domain))) if len(domain) == 1 else f'v{value.index}'

    def assertion(self, expression):
        self.lines.append(f'(assert {expression})')

    def factor(self, factor):
        if factor not in self.factors:
            name = f'factor{len(self.factors)}'
            choices = {}
            for i, permutation in enumerate(factor):
                result = 'p'
                for atom, image in enumerate(permutation):
                    if atom != image:
                        result = f'(ite (= p {atom}) {image} {result})'
                choices[i] = result
            expression = _lookup('c', choices)
            self.lines.append(f'(define-fun {name} ((c Int) (p Int)) Int {expression})')
            self.factors[factor] = name
        return self.factors[factor]

    def table(self, values, rows):
        if not values:
            return 'true' if rows else 'false'
        groups = defaultdict(set)
        for row in rows:
            groups[row[0]].add(row[1:])
        shared = defaultdict(list)
        for item, suffixes in groups.items():
            shared[frozenset(suffixes)].append(item)
        return _combine('or', (_combine('and', (
            _member(self.value(values[0]), items), self.table(values[1:], suffixes)))
            for suffixes, items in shared.items()))

    def constraints(self):
        for kind, values, data in self.system.constraints:
            tokens = tuple(self.value(value) for value in values)
            if kind == 'action':
                source, choice, target = tokens
                self.assertion(f'(= {target} ({self.factor(data)} {choice} {source}))')
            elif kind == 'lookup':
                selector, output = tokens
                expression = _lookup(selector, {p: self.value(data[p])
                    for p in self.system.domain(values[0])})
                self.assertion(f'(= {output} {expression})')
            elif kind == 'table':
                self.assertion(self.table(values, data))
            elif kind == 'distinct':
                if len(tokens) > 1:
                    self.assertion(f'(distinct {" ".join(tokens)})')
            else:
                raise ValueError(f'Unknown finite constraint: {kind}')


class LocalCountQuery:
    """An ephemeral joint query; callers retain numeric minima, not solvers."""
    def __init__(self, solver, images, counts, domains):
        self.solver, self.images, self.counts = solver, images, counts
        self._domains = domains
        self._pair_count = 0

    def minimum(self, expression, lower, upper):
        """Find the exact minimum between a proven bound and attainable witness."""
        import z3
        if lower > upper:
            raise ValueError('Local minimum lower bound exceeds its witness')
        while lower < upper:
            middle = (lower + upper) // 2
            self.solver.push()
            try:
                self.solver.add(expression <= middle)
                status = self.solver.check()
                if status == z3.unknown:
                    raise RuntimeError(f'Local event minimum was inconclusive: {self.solver.reason_unknown()}')
                if status == z3.unsat:
                    lower = middle + 1
                else:
                    value = expression if isinstance(expression, int) else self.solver.model().eval(
                        expression, model_completion=True).as_long()
                    upper = int(value)
            finally:
                self.solver.pop()
        return int(lower)

    def pair_cost(self, left, right, rows):
        """Add one complete endpoint-pair cost table, including negative costs."""
        import z3
        xs = (left,) if isinstance(left, int) else self._domains[str(left)]
        ys = (right,) if isinstance(right, int) else self._domains[str(right)]
        if not {(x, y) for x in xs for y in ys}.issubset(rows):
            raise ValueError('Pair cost must cover both endpoint domains')
        values = np.array([[rows[x, y] for y in ys] for x in xs], dtype=int)
        if np.all(values == values[0, 0]):
            return int(values[0, 0])
        name = f'pair_cost{self._pair_count}'
        self._pair_count += 1
        output = z3.Int(name)
        declarations = {str(v): v for v in (left, right, output) if not isinstance(v, int)}
        expression = _matrix(_integer(left) if isinstance(left, int) else str(left),
            _integer(right) if isinstance(right, int) else str(right), xs, ys, values)
        self.solver.add(z3.parse_smt2_string(f'(assert (= {name} {expression}))', decls=declarations))
        return output


def compile_local(index, responses):
    """Compile independent fragment copies while preserving each joint program."""
    import z3
    system = FiniteConstraints()
    images = []
    for number, response in enumerate(responses):
        current = {a: system.constant(p) for a, p in response.mapping}
        for position, (kind, data) in enumerate(response.actions):
            label = number, position
            if kind == 'relation':
                compiled = compile_relaxed_relation(system, data.relation, label)
                transport = {p: compiled.images[r] for r, p in data.witness}
                current = {a: system.lookup(value, transport) for a, value in current.items()}
            elif kind == 'group':
                current = dict(zip(current, system.act(list(current.values()), data, label)))
            else:
                raise ValueError(f'Unknown fragment action: {kind}')
        images.append(current)
    system.propagate()
    if system.failed:
        raise ValueError('Saved fragment relationship has no valid witness')
    encoding = _IntegerEncoding(system)
    encoding.constraints()
    for number, (response, current) in enumerate(zip(responses, images)):
        terms = []
        for term in response.terms:
            if term.response.constant is not None:
                terms.append(str(int(term.response.constant != 0)))
                continue
            a, b = term.atoms
            left, right = current[a][0], current[b][0]
            xs, ys = tuple(sorted(system.domain(left))), tuple(sorted(system.domain(right)))
            # Projected response categories can omit unreachable Cartesian pairs.
            # Evaluate their original weights; the joint constraints remove them.
            changed = (np.abs(index.p[np.ix_(xs, ys)] - term.response.source_weight)
                       > term.response.threshold).astype(int)
            terms.append(_matrix(encoding.value(left), encoding.value(right), xs, ys, changed))
        encoding.lines.append(f'(declare-const events{number} Int)')
        encoding.assertion(f'(= events{number} {_combine("+", terms)})')
    solver = z3.Solver()
    solver.from_string('\n'.join(encoding.lines))
    domains = {}
    def physical(value):
        domain = tuple(sorted(system.domain(value[0])))
        if len(domain) == 1:
            return domain[0]
        name = encoding.value(value[0])
        domains[name] = domain
        return z3.Int(name)
    return LocalCountQuery(solver, tuple({a: physical(v) for a, v in current.items()} for current in images),
                           tuple(z3.Int(f'events{i}') for i in range(len(responses))), domains)
