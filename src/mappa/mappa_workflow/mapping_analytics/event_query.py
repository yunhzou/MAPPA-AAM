"""Symbolic equality queries for signed event classes."""
from mappa.mappa_core.core.symmetry.event_counts import _event_count_constraint


def _same_event_pattern(terms, objective, pattern):
    """Exact sparse equality, also safe when different event counts are allowed."""
    import z3
    return z3.And(_event_count_constraint(terms, pattern['total']), *(terms.get(tuple(pair), 0) == sign
        for kind, sign in (('broken', -1), ('formed', 1)) for pair in pattern['events'][kind]))


def _same_event_orbit(terms, objective, pattern, index):
    """An entire signed-event class via factored source actions, not orbit lists.

    Equal count plus the transported signed edges implies equality, including
    every absent edge. Source actions only define output identity; they impose
    no extra assumptions about invariance of the original search family.
    """
    import z3
    from mappa.mappa_core.core.symmetry.constraints import SymbolicActions
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
