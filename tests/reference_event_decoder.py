"""Test-only signed-event enumeration oracle; never imported by MAPPA."""
import math
import time
from mappa.event_patterns import _event_model, _event_count_constraint, _same_event_orbit

def extract_path_events(path, search_problem, index, *, max_events=None, seconds=2.0,
                        max_patterns=32, on_pattern=None, query_log=None):
    """Bounded signed-event projection of one saved partial or complete path.

    Completeness concerns this path/window, not the archive or global search.
    Time is a soft budget; native canonicalization and path compilation require
    an external process watchdog for a hard deadline. Each exclusion removes
    an entire output class using quantified factored symmetry actions.
    ``max_patterns=None`` removes the class-count limit. ``on_pattern`` receives
    each newly discovered class so callers can persist witnesses before a
    process watchdog interrupts an expensive completeness proof.
    Optional query_log(text, status, cpu_seconds) persists exact solver evidence;
    its returned references are retained in the family report.
    """
    import z3
    from mappa.family_query import compile_path
    from mappa.search_validation import RepresentativeValidationWorkspace
    if (not math.isfinite(seconds) or seconds < 0 or
            (max_patterns is not None and
             (not isinstance(max_patterns, int) or max_patterns < 1))):
        raise ValueError('invalid extraction budget')
    if max_events is not None and (not isinstance(max_events, int) or max_events < 0):
        raise ValueError('max_events must be a nonnegative integer')
    if (search_problem.reactant.elements != index.problem.reactant.elements or
            search_problem.product.elements != index.problem.product.elements):
        raise ValueError('index and search path must use the same atom indexing')
    started = time.perf_counter()
    deadline = started + seconds
    validator = getattr(index, '_representative_validator', None)
    if validator is None or validator.problem is not search_problem:
        validator = RepresentativeValidationWorkspace(search_problem)
        index._representative_validator = validator
    validator.validate_path(path)
    patterns = {}
    queries = 0
    metrics = {}
    metrics['exclusion'] = 'signed_event_orbits'
    if query_log is not None:
        metrics['query_evidence'] = []

    def finish(complete, reason):
        return dict(patterns=sorted(patterns.values(), key=lambda p: (p['total'], p['id'])),
                    complete=complete, reason=reason, solver_queries=queries,
                    seconds=time.perf_counter() - started, max_events=max_events,
                    terminal=path.terminal, transitions=list(path.transitions), **metrics)

    def record(mapping, actions):
        pattern = index.describe(index.vector(mapping))
        pattern['actions'] = actions
        if pattern['id'] not in patterns:
            patterns[pattern['id']] = pattern
            if on_pattern is not None:
                on_pattern(pattern)
        return pattern

    representative = dict(path.mapping)
    score = int(index.counts([index.vector(representative)])[0].sum())
    seed = record(representative, []) if max_events is None or score <= max_events else None
    if time.perf_counter() >= deadline:
        return finish(False, 'time_budget')
    if index.invariant_path(path):
        return finish(True, 'signed_event_invariance')
    if time.perf_counter() >= deadline:
        return finish(False, 'time_budget')
    compiled = compile_path(path, search_problem, {}, source_atoms=(), complete_reference=False)
    metrics['constraint_encoding_seconds'] = compiled.encoding_seconds
    try:
        terms, objective, counts = _event_model(compiled, index, deadline)
    except TimeoutError:
        return finish(False, 'event_encoding_budget')
    metrics.update(counts)
    if max_events is not None:
        compiled.solver.add(_event_count_constraint(terms, max_events, at_most=True))
    if seed is not None:
        compiled.solver.add(z3.Not(_same_event_orbit(terms, objective, seed, index)))
    while time.perf_counter() < deadline:
        compiled.solver.set(timeout=max(1, int(1000 * (deadline - time.perf_counter()))))
        query_text = compiled.solver.to_smt2() if query_log is not None else None
        query_started = time.process_time()
        status = compiled.solver.check()
        query_seconds = time.process_time() - query_started
        if query_log is not None:
            metrics['query_evidence'].append(query_log(query_text, str(status), query_seconds))
        queries += 1
        if status == z3.unsat:
            return finish(True, 'all_event_orbits_excluded')
        if status == z3.unknown:
            return finish(False, compiled.solver.reason_unknown())
        if max_patterns is not None and len(patterns) >= max_patterns:
            return finish(False, 'pattern_budget')
        model = compiled.solver.model()
        witness = compiled.realize(model)
        pattern = record(dict(witness['mapping']), witness['actions'])
        if pattern['total'] != model.eval(objective).as_long():
            raise AssertionError('symbolic and scalar event counts disagree')
        compiled.solver.add(z3.Not(_same_event_orbit(terms, objective, pattern, index)))
    return finish(False, 'time_budget')
