"""Joint event counts from local endpoint-pair responses."""
import time
from dataclasses import dataclass
from mappa.mappa_core.utils.bonds import is_metal_element, bond_event_threshold

class _EventTerms(dict):
    """Signed terms plus exact Boolean event-presence indicators."""
    def __init__(self):
        super().__init__()
        self.constant_count = 0
        self.event_flags = []


def _event_model(compiled, index, deadline, *, internal_only=False):
    """Small local pair tables; never enumerate complete assignments."""
    import z3
    terms, objective, variable_pairs, table_entries, lower = _EventTerms(), [], 0, 0, 0
    constant_events = 0
    selectors = {}
    mapped = set(compiled.representative)
    terms.mapped_source = mapped
    terms.unmapped = {}
    system = compiled.constraints
    values = {a: (compiled.physical_values[a][0],
                  frozenset(system.domain(compiled.physical_values[a][0]))) for a in mapped}

    def select(value, subset):
        expr, support = value
        if len(subset) == len(support):
            return True
        if isinstance(expr, int):
            return expr in subset
        key = expr, subset
        if key not in selectors:
            selectors[key] = z3.Or(*(system.predicate(expr, x) for x in subset))
        return selectors[key]

    domains = {a: value[1] for a, value in values.items()}
    for a, b in sorted(index.neighboring_pairs(domains)):
        if internal_only and (a not in mapped or b not in mapped):
            continue
        if time.perf_counter() >= deadline:
            raise TimeoutError('event encoding budget')
        tau = bond_event_threshold(index.elements, a, b,
            default_threshold=index.threshold, metal_threshold=index.metal_threshold)
        if a not in mapped or b not in mapped:
            sign = -int((a in mapped or b in mapped) and index.r[a, b] > tau)
            terms[int(a), int(b)] = sign
            constant_events += int(sign != 0)
            lower += int(sign != 0)
            continue
        va, vb = values[a], values[b]
        default, rectangles, pair_lower, entries = index.pair_table(index.r[a, b], tau, va[1], vb[1])
        table_entries += entries
        term = default
        if not entries:
            compiled.solver.add(False)
        sign_conditions = {}
        for sign, blocks in rectangles.items():
            conditions = []
            for xs, ys in blocks:
                left, right = select(va, xs), select(vb, ys)
                conditions.append(right if left is True else left if right is True else z3.And(left, right))
            condition = z3.Or(*conditions)
            sign_conditions[sign] = condition
            term = z3.If(condition, sign, term)
        variable_pairs += bool(rectangles)
        lower += pair_lower
        terms[int(a), int(b)] = term
        if pair_lower:
            # This pair is always an event, even if its sign may change.
            constant_events += 1
        elif rectangles:
            # Count presence directly, without asking the arithmetic solver
            # to derive If(sign_condition, -1, +1) != 0 repeatedly.
            present = (z3.Or(*sign_conditions.values()) if default == 0
                       else z3.Not(sign_conditions[0]))
            terms.event_flags.append(present)
            objective.append(z3.If(present, 1, 0))
    if not internal_only and len(mapped) != index.m:
        # Each target boundary edge has one mapped end. Encode it once in that
        # source atom's frame, with the other end a native unmatched P index.
        for p in range(index.m):
            possible = [value[0] for value in values.values() if p in value[1]]
            terms.unmapped[p] = z3.simplify(z3.And(*(
                z3.Not(system.predicate(value, p)) for value in possible)))
        for a in sorted(mapped):
            value = values[a]
            for p in range(index.m):
                if time.perf_counter() >= deadline:
                    raise TimeoutError('event encoding budget')
                tau = (index.metal_threshold if index.metal_threshold is not None and
                       (is_metal_element(index.elements[a]) or
                        is_metal_element(index.problem.product.elements[p])) else index.threshold)
                support = tuple(x for x in value[1] if index.p[x, p] > tau)
                if not support:
                    continue
                condition = z3.simplify(z3.And(select(value, support), terms.unmapped[p]))
                if z3.is_false(condition):
                    continue
                if z3.is_true(condition):
                    terms[a, index.n+p] = 1
                    constant_events += 1
                    lower += 1
                else:
                    terms[a, index.n+p] = z3.If(condition, 1, 0)
                    terms.event_flags.append(condition)
                    objective.append(z3.If(condition, 1, 0))
                    variable_pairs += 1
    terms.constant_count = constant_events
    total = constant_events + z3.Sum(objective) if objective else z3.IntVal(constant_events)
    return terms, total, dict(
        variable_pairs=variable_pairs, pair_table_entries=table_entries, event_lower_bound=lower)


def _event_count_constraint(terms, count, *, at_most=False):
    """Exact Boolean cardinality, avoiding arithmetic search for an event sum."""
    import z3
    if isinstance(terms, _EventTerms):
        constant = terms.constant_count
        choices = [(flag, 1) for flag in terms.event_flags]
    else:
        constant = sum(int(t != 0) for t in terms.values() if isinstance(t, int))
        choices = [(t != 0, 1) for t in terms.values() if not isinstance(t, int)]
    remaining = count - constant
    if remaining == 0 and choices:
        return z3.And(*(z3.Not(flag) for flag, _ in choices))
    if not choices:
        return z3.BoolVal(remaining >= 0 if at_most else remaining == 0)
    return z3.PbLe(choices, remaining) if at_most else z3.PbEq(choices, remaining)


@dataclass(frozen=True)
class LocalEventSpectrum:
    minimum: int
    counts: tuple[int, ...]
    endpoint_counts: tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]

    def within(self, allowance):
        """Filter a proved support without rebuilding its local constraints."""
        if allowance is None or self.counts[-1] <= self.minimum + allowance:
            return self
        limit = self.minimum + allowance
        conditioned = []
        for images, counts in self.endpoint_counts:
            retained = tuple(count for count in counts if count <= limit)
            if retained:
                conditioned.append((images, retained))
        return LocalEventSpectrum(self.minimum,
            tuple(count for count in self.counts if count <= limit), tuple(conditioned))


class LocalEventCounts:
    """Joint internal counts, shared within one decoder without solver retention."""
    def __init__(self, index):
        self.index = index
        self.cache = {}
        self.minima = {}
        self.solver_queries = 0
        self.cache_hits = 0

    def spectrum(self, mapping, actions, *, endpoints=(), allowance=None):
        """Project only event totals and requested cut-endpoint images.

        A local relation retains its complete correlated discovery and self
        constraints. Other fragments' choices and boundary events are absent.
        Count gaps are preserved; this does not enumerate full bijections.
        """
        import z3
        from numbers import Integral
        from mappa.mappa_core.core.symmetry.family import compile_actions

        if allowance is not None and (isinstance(allowance, bool) or
                not isinstance(allowance, Integral) or allowance < 0):
            raise ValueError('local event allowance must be a nonnegative integer or None')
        mapping = tuple(sorted(dict(mapping).items()))
        actions, endpoints = tuple(actions), tuple(endpoints)
        if len(set(endpoints)) != len(endpoints) or not set(endpoints).issubset(dict(mapping)):
            raise ValueError('local endpoints must name distinct mapped source atoms')
        key = mapping, actions, endpoints
        if key in self.cache:
            covered, result = self.cache[key]
            if covered is None or (allowance is not None and allowance <= covered):
                self.cache_hits += 1
                return result.within(allowance)
        local_key = mapping, actions
        if allowance == 0 and not endpoints and local_key in self.minima:
            self.cache_hits += 1
            minimum = self.minima[local_key]
            result = LocalEventSpectrum(minimum, (minimum,), ())
            self.cache[key] = allowance, result
            return result
        compiled = compile_actions(actions, self.index.problem, mapping)
        terms, _, _ = _event_model(compiled, self.index, float('inf'), internal_only=True)
        solver = compiled.solver

        def check():
            status = solver.check()
            self.solver_queries += 1
            if status == z3.unknown:
                raise RuntimeError(f'Local event query was inconclusive: {solver.reason_unknown()}')
            return status

        def count(model):
            return terms.constant_count + sum(z3.is_true(model.eval(flag, model_completion=True))
                                               for flag in terms.event_flags)

        if local_key not in self.minima:
            if check() != z3.sat:
                raise ValueError('Saved local mapping relation has inconsistent constraints')
            lower, upper = terms.constant_count, count(solver.model())
            while lower < upper:
                bound = (lower + upper - 1) // 2
                solver.push()
                solver.add(_event_count_constraint(terms, bound, at_most=True))
                if check() == z3.unsat:
                    lower = bound + 1
                else:
                    upper = count(solver.model())
                solver.pop()
            self.minima[local_key] = upper
        else:
            self.cache_hits += 1
        minimum = self.minima[local_key]
        if not endpoints and (not terms.event_flags or allowance == 0):
            result = LocalEventSpectrum(minimum, (minimum,), ())
            # With no variable event flags this support is complete at every
            # allowance, even if the caller requested only the minimum.
            self.cache[key] = (None if not terms.event_flags else allowance), result
            return result
        if allowance is not None:
            solver.add(_event_count_constraint(terms, minimum + allowance, at_most=True))
        endpoint_values = tuple(compiled.physical_values[atom][0] for atom in endpoints)
        endpoint_expressions = tuple(compiled.constraints.expression(value) for value in endpoint_values)
        totals, conditioned = set(), {}
        count_conditions = {}
        while check() == z3.sat:
            model = solver.model()
            total = count(model)
            totals.add(total)
            images = tuple(value if isinstance(value, int) else
                           model.eval(value, model_completion=True).as_long()
                           for value in endpoint_expressions)
            if endpoints:
                conditioned.setdefault(images, set()).add(total)
            if total not in count_conditions:
                count_conditions[total] = _event_count_constraint(terms, total)
            solver.add(z3.Not(z3.And(count_conditions[total], *(
                compiled.constraints.predicate(value, image)
                for value, image in zip(endpoint_values, images)))))
        result = LocalEventSpectrum(minimum, tuple(sorted(totals)), tuple(
            (images, tuple(sorted(counts))) for images, counts in sorted(conditioned.items())))
        self.cache[key] = allowance, result
        return result
