"""Symbolic equality queries for signed event classes; never called by decoding."""
import time
from mappa.mappa_core.utils.bonds import is_metal_element

class _EventTerms(dict):
    """Signed terms plus exact Boolean event-presence indicators."""
    def __init__(self):
        super().__init__()
        self.constant_count = 0
        self.event_flags = []


def _event_model(compiled, index, deadline):
    """Small local pair tables; never enumerate complete assignments."""
    import z3
    terms, objective, variable_pairs, table_entries, lower = _EventTerms(), [], 0, 0, 0
    constant_events = 0
    selectors = {}
    mapped = set(compiled.representative)
    terms.mapped_source = mapped
    terms.unmapped = {}

    def select(value, subset):
        expr, support = value
        if len(subset) == len(support):
            return True
        if isinstance(expr, int):
            return expr in subset
        key = expr.get_id(), subset
        if key not in selectors:
            selectors[key] = z3.Or(*(expr == x for x in subset))
        return selectors[key]

    for a, b, tau in zip(index.a, index.b, index.tau):
        if time.perf_counter() >= deadline:
            raise TimeoutError('event encoding budget')
        if a not in mapped or b not in mapped:
            sign = -int((a in mapped or b in mapped) and index.r[a, b] > tau)
            terms[int(a), int(b)] = sign
            constant_events += int(sign != 0)
            lower += int(sign != 0)
            continue
        va, vb = compiled.values[a], compiled.values[b]
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
    if len(mapped) != index.m:
        # Each target boundary edge has one mapped end. Encode it once in that
        # source atom's frame, with the other end a native unmatched P index.
        for p in range(index.m):
            possible = [compiled.values[a][0] for a in mapped if p in compiled.values[a][1]]
            terms.unmapped[p] = z3.simplify(z3.And(*(v != p for v in possible)))
        for a in sorted(mapped):
            value = compiled.values[a]
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


def _same_event_pattern(terms, objective, pattern):
    """Exact sparse equality, also safe when different event counts are allowed."""
    import z3
    return z3.And(_event_count_constraint(terms, pattern['total']), *(terms[tuple(pair)] == sign
        for kind, sign in (('broken', -1), ('formed', 1)) for pair in pattern['events'][kind]))


def _same_event_orbit(terms, objective, pattern, index):
    """An entire signed-event class via factored source actions, not orbit lists.

    Equal count plus the transported signed edges implies equality, including
    every absent edge. Source actions only define output identity; they impose
    no extra assumptions about invariance of the original search family.
    """
    import z3
    from mappa.mappa_workflow.mapping_analytics.family import SymbolicActions
    extras = sorted({a for pairs in pattern['events'].values() for pair in pairs
                     for a in pair if a >= index.n})
    if pattern['total'] == 0 or (not index.source_generators and not extras):
        return _same_event_pattern(terms, objective, pattern)
    solver = z3.Solver()
    encoder = SymbolicActions(solver, fresh=True)
    atoms = {a for pairs in pattern['events'].values() for pair in pairs for a in pair if a < index.n}
    mapped = set(pattern.get('mapped_source', range(index.n)))
    if len(mapped) != index.n:
        atoms.update(mapped)
    atoms = sorted(atoms)
    moved = encoder.act([encoder.constant(a) for a in atoms], index.source_generators,
                        'signed_event_equivalence')
    images = dict(zip(atoms, moved))
    conditions = [_event_count_constraint(terms, pattern['total'])]
    if len(mapped) != index.n:
        conditions.extend(z3.Or(*(images[a][0] == b for b in terms.mapped_source)) for a in mapped)
    # Only boundary atoms used by the pattern need images. Distinct preserves
    # shared endpoints; it does not enumerate their permutations or whole AAMs.
    for a in extras:
        element = index.problem.product.elements[a-index.n]
        images[a] = encoder.symbol(index.n+p for p, e in enumerate(index.problem.product.elements) if e == element)
        conditions.append(z3.Or(*(z3.And(images[a][0] == index.n+p, terms.unmapped[p])
                                  for p, e in enumerate(index.problem.product.elements) if e == element)))
    if len(extras) > 1:
        values = [images[a][0] for a in extras]
        conditions.append(z3.Distinct(*(z3.IntVal(v) if isinstance(v, int) else v for v in values)))
    for kind, sign in (('broken', -1), ('formed', 1)):
        for a, b in pattern['events'][kind]:
            va, vb = images[a], images[b]
            rows = {x: encoder.lookup(vb, {
                y: (terms.get(tuple(sorted((x, y))), 0), frozenset((-1, 0, 1)))
                if x != y else encoder.constant(0) for y in vb[1]}) for x in va[1]}
            conditions.append(encoder.lookup(va, rows)[0] == sign)
    equality = z3.And(*solver.assertions(), *conditions)
    return z3.Exists(encoder.variables, equality) if encoder.variables else equality
