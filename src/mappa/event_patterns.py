"""Opt-in signed-event decoding of complete, balanced saved AAM paths.

Search constraints stay at their recorded tolerance. Output identity uses
raw-WBO event responses, with no bijection or group-element expansion.
The existing Golden floor-based scorer and historical pattern IDs are separate.
"""
from collections import Counter, defaultdict, OrderedDict
from functools import cached_property
import hashlib
import json
import math
import time

import numpy as np
import pynauty

from .frag import classify_bonds, is_metal_element


def _response_bounds(values, weight, threshold):
    """Encode (-, zero, +) responses of sorted opposite weights with two cuts.

    Keep the original subtraction/comparisons, including floating boundaries.
    Algebraically moving the threshold to the other operand need not agree.
    """
    lo, hi = 0, len(values)
    while lo < hi:
        mid = (lo + hi) // 2
        if values[mid] - weight <= -threshold:
            lo = mid + 1
        else:
            hi = mid
    negative_end = lo
    lo, hi = 0, len(values)
    while lo < hi:
        mid = (lo + hi) // 2
        if values[mid] - weight >= threshold:
            hi = mid
        else:
            lo = mid + 1
    return negative_end, lo


class SignedEventIndex:
    """Reaction-local event identity; IDs are versioned separately from legacy.

    The raw problem must have equal element multisets and nonnegative, exactly
    symmetric WBOs. Partial mappings and historical floor-based events require
    their own output policy and are deliberately rejected here.
    """
    schema = 'rxn_core.signed_events/v1'

    def __init__(self, raw_problem, *, threshold=0.5, metal_threshold=0.3):
        for value in (threshold, metal_threshold):
            if value is not None and (not math.isfinite(value) or value <= 0):
                raise ValueError('event thresholds must be finite and positive')
        if threshold is None:
            raise ValueError('ordinary threshold is required')
        self.problem = raw_problem
        self.elements = raw_problem.reactant.elements
        if Counter(self.elements) != Counter(raw_problem.product.elements):
            raise ValueError('signed-event index requires balanced elements')
        self.n = len(self.elements)
        self.r, self.p = raw_problem.reactant.wbo, raw_problem.product.wbo
        for matrix in (self.r, self.p):
            if not np.array_equal(matrix, matrix.T) or np.any(matrix < 0):
                raise ValueError('raw WBOs must be nonnegative and exactly symmetric')
        self.threshold = float(threshold)
        self.metal_threshold = None if metal_threshold is None else float(metal_threshold)
        self.a, self.b = np.triu_indices(self.n, 1)
        metal = np.array([is_metal_element(e) for e in self.elements])
        self.tau = np.full(len(self.a), self.threshold)
        if self.metal_threshold is not None:
            self.tau[metal[self.a] | metal[self.b]] = self.metal_threshold
        self.cache = {}
        self._actions = {}
        self._pair_tables = OrderedDict()

    def pair_table(self, weight, tau, xs, ys):
        """Reuse exact local responses across pairs/families, never mappings.

        Raw weights remain in the key: no rounding or tolerance approximation.
        The bounded cache only saves work; eviction cannot change the result.
        """
        key = float(weight), float(tau), xs, ys
        if key in self._pair_tables:
            self._pair_tables.move_to_end(key)
            return self._pair_tables[key]
        xx, yy = sorted(xs), sorted(ys)
        delta = self.p[np.ix_(xx, yy)] - weight
        valid = np.asarray(xx)[:, None] != np.asarray(yy)[None, :]
        signs = np.where(delta <= -tau, -1, np.where(delta >= tau, 1, 0))
        counts = {int(k): int(v) for k, v in zip(*np.unique(signs[valid], return_counts=True))}
        rectangles = {}
        default = max(counts, key=counts.get) if counts else 0
        if len(counts) > 1:
            for sign in counts:
                if sign == default:
                    continue
                rows = defaultdict(list)
                for i, x in enumerate(xx):
                    selected = tuple(yy[j] for j in np.flatnonzero(valid[i] & (signs[i] == sign)))
                    if selected:
                        rows[selected].append(x)
                rectangles[sign] = tuple((tuple(row), col) for col, row in rows.items())
        result = default, rectangles, min((int(k != 0) for k in counts), default=0), sum(counts.values())
        self._pair_tables[key] = result
        if len(self._pair_tables) > 16384:
            self._pair_tables.popitem(last=False)
        return result

    @property
    def policy(self):
        return dict(schema=self.schema, threshold=self.threshold,
                    metal_threshold=self.metal_threshold, includes_hydrogens=True)

    @cached_property
    def family_certificates(self):
        from .event_certificates import EventFamilyCertificates
        return EventFamilyCertificates(self)

    def _vectors(self, vectors):
        values = np.asarray(vectors)
        if values.size == 0:
            if values.shape == (0,) or values.shape == (0, self.n):
                return np.empty((0, self.n), dtype=int)
            raise ValueError('expected complete integer mapping vectors')
        if values.ndim != 2 or values.shape[1] != self.n or values.dtype.kind not in 'iu':
            raise ValueError('expected complete integer mapping vectors')
        if not np.all(np.sort(values, axis=1) == np.arange(self.n)):
            raise ValueError('mappings must be complete bijections')
        if not np.all(np.asarray(self.problem.product.elements)[values] == np.asarray(self.elements)):
            raise ValueError('mappings must preserve elements')
        return values

    def counts(self, vectors):
        """Vectorized scoring; callers choose bounded batch sizes."""
        values = self._vectors(vectors)
        delta = self.p[values[:, self.a], values[:, self.b]] - self.r[self.a, self.b]
        return np.column_stack(((delta <= -self.tau).sum(axis=1),
                                (delta >= self.tau).sum(axis=1)))

    @cached_property
    def response_data(self):
        endpoints = (self.problem.reactant, self.problem.product)
        result = []
        for side in (0, 1):
            ep, other = endpoints[side], endpoints[1 - side]
            opposite_values = defaultdict(set)
            for a, b in zip(self.a, self.b):
                pair = tuple(sorted((other.elements[a], other.elements[b])))
                opposite_values[pair].add(float(other.wbo[a, b]))
            opposite_values = {k: tuple(sorted(v)) for k, v in opposite_values.items()}
            memo, pairs, counts = {}, [], defaultdict(Counter)
            for a, b in zip(self.a, self.b):
                pair = tuple(sorted((ep.elements[a], ep.elements[b])))
                tau = (self.metal_threshold if self.metal_threshold is not None
                       and any(is_metal_element(e) for e in pair) else self.threshold)
                key = pair, float(ep.wbo[a, b])
                if key not in memo:
                    memo[key] = _response_bounds(opposite_values[pair], key[1], tau)
                label = memo[key]
                pairs.append((int(a), int(b), pair, label))
                counts[pair][label] += 1
            defaults = {pair: min(c, key=lambda label: (-c[label], label))
                        for pair, c in counts.items()}
            labels = np.full((self.n, self.n), -1, dtype=int)
            intern = {}
            exceptional = []
            for a, b, pair, label in pairs:
                color = pair, label
                code = intern.setdefault(color, len(intern))
                labels[a, b] = labels[b, a] = code
                if label != defaults[pair]:
                    exceptional.append((a, b, color))
            labels.setflags(write=False)
            result.append(dict(elements=ep.elements, labels=labels,
                               pairs=tuple(exceptional), full_pair_count=len(pairs),
                               unique_weight_responses=len(memo)))
        return tuple(result)

    @cached_property
    def base(self):
        data = self.response_data[0]
        adjacency = {i: set() for i in range(self.n)}
        colors = defaultdict(set)
        for i, element in enumerate(data['elements']):
            colors[('atom', element)].add(i)
        for a, b, label in data['pairs']:
            vertex = len(adjacency)
            adjacency[vertex] = {a, b}
            adjacency[a].add(vertex)
            adjacency[b].add(vertex)
            colors[('pair', label)].add(vertex)
        return pynauty.Graph(len(adjacency), adjacency_dict={i: sorted(ns) for i, ns in adjacency.items()},
                            vertex_coloring=[colors[k] for k in sorted(colors, key=repr)])

    @cached_property
    def source_generators(self):
        return tuple(tuple(g[:self.n]) for g in pynauty.autgrp(self.base)[0])

    def preserves_target_action(self, images):
        images = tuple(images)
        if images not in self._actions:
            if len(images) != self.n or sorted(images) != list(range(self.n)):
                raise ValueError('invalid saved target permutation')
            elements = self.problem.product.elements
            labels = self.response_data[1]['labels']
            self._actions[images] = (all(elements[a] == elements[b] for a, b in enumerate(images))
                                    and np.array_equal(labels, labels[np.ix_(images, images)]))
        return self._actions[images]

    def invariant_path(self, path):
        """Sufficient full signed-event certificate, including boundary pairs."""
        for edge in path.transitions:
            placement = path.graph.fragment_placement(edge)
            if placement is None:
                continue
            if placement.target_generators is None:
                raise ValueError('finalize saved path symmetry first')
            for generator in placement.target_generators:
                if not self.preserves_target_action(generator.images):
                    return False
            for domain in placement.symmetry_domains:
                if domain.source == 'exact_automorph_group':
                    continue
                for a, b in zip(domain.p_atoms, domain.p_atoms[1:]):
                    images = list(range(self.n))
                    images[a], images[b] = images[b], images[a]
                    if not self.preserves_target_action(images):
                        return False
        return True

    def describe(self, vector):
        vector = self._vectors([vector])[0]
        delta = self.p[vector[self.a], vector[self.b]] - self.r[self.a, self.b]
        events = tuple(tuple((int(a), int(b)) for a, b in zip(self.a[mask], self.b[mask]))
                       for mask in (delta <= -self.tau, delta >= self.tau))
        if events not in self.cache:
            base = self.base
            adjacency = {i: set(ns) for i, ns in base.adjacency_dict.items()}
            colors = [set(cell) for cell in base.vertex_coloring]
            n = base.number_of_vertices
            for pairs in events:
                marker = n
                n += 1
                adjacency[marker] = set()
                colors.append({marker})
                gadgets = set()
                for a, b in pairs:
                    adjacency[n] = {a, b, marker}
                    for atom in (a, b, marker):
                        adjacency[atom].add(n)
                    gadgets.add(n)
                    n += 1
                if gadgets:
                    colors.append(gadgets)
            graph = pynauty.Graph(n, adjacency_dict={i: sorted(adjacency.get(i, ())) for i in range(n)},
                                  vertex_coloring=colors)
            payload = json.dumps(self.policy, sort_keys=True).encode() + pynauty.certificate(graph)
            self.cache[events] = self.schema + ':' + hashlib.sha256(payload).hexdigest()
        # Independently validate every emitted witness using the existing scalar rule.
        broken, formed, _, _ = classify_bonds(dict(enumerate(map(int, vector))), self.r, self.p,
            dwbo_threshold=self.threshold, elements_R=self.elements,
            elements_P=self.problem.product.elements, metal_dwbo_threshold=self.metal_threshold)
        inverse = {int(p): r for r, p in enumerate(vector)}
        if ({(a, b) for a, b, *_ in broken} != set(events[0]) or
                {tuple(sorted((inverse[a], inverse[b]))) for a, b, *_ in formed} != set(events[1])):
            raise AssertionError('vector and scalar event witnesses disagree')
        return dict(id=self.cache[events], mapping=vector.tolist(),
                    events=dict(zip(('broken', 'formed'), events)),
                    counts=dict(zip(('broken', 'formed'), map(len, events))), total=sum(map(len, events)),
                    policy=self.policy)


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
    from .family_query import SymbolicActions
    if not index.source_generators or pattern['total'] == 0:
        return _same_event_pattern(terms, objective, pattern)
    solver = z3.Solver()
    encoder = SymbolicActions(solver, fresh=True)
    atoms = sorted({a for pairs in pattern['events'].values() for pair in pairs for a in pair})
    moved = encoder.act([encoder.constant(a) for a in atoms], index.source_generators,
                        'signed_event_equivalence')
    images = dict(zip(atoms, moved))
    conditions = [_event_count_constraint(terms, pattern['total'])]
    for kind, sign in (('broken', -1), ('formed', 1)):
        for a, b in pattern['events'][kind]:
            va, vb = images[a], images[b]
            rows = {x: encoder.lookup(vb, {
                y: (terms[tuple(sorted((x, y)))], frozenset((-1, 0, 1)))
                if x != y else encoder.constant(0) for y in vb[1]}) for x in va[1]}
            conditions.append(encoder.lookup(va, rows)[0] == sign)
    equality = z3.And(*solver.assertions(), *conditions)
    return z3.Exists(encoder.variables, equality) if encoder.variables else equality


def extract_path_events(path, search_problem, index, *, max_events=None, seconds=2.0,
                        max_patterns=32, on_pattern=None):
    """Bounded signed-event projection of one complete saved path.

    Completeness concerns this path/window, not the archive or global search.
    Time is a soft budget; native canonicalization and path compilation require
    an external process watchdog for a hard deadline. Each exclusion removes
    an entire output class using quantified factored symmetry actions.
    ``max_patterns=None`` removes the class-count limit. ``on_pattern`` receives
    each newly discovered class so callers can persist witnesses before a
    process watchdog interrupts an expensive completeness proof.
    """
    import z3
    from .family_query import compile_path
    from .search_validation import RepresentativeValidationWorkspace
    if (not math.isfinite(seconds) or seconds < 0 or
            (max_patterns is not None and
             (not isinstance(max_patterns, int) or max_patterns < 1))):
        raise ValueError('invalid extraction budget')
    if max_events is not None and (not isinstance(max_events, int) or max_events < 0):
        raise ValueError('max_events must be a nonnegative integer')
    if (search_problem.reactant.elements != index.problem.reactant.elements or
            search_problem.product.elements != index.problem.product.elements or
            len(path.mapping) != index.n):
        raise ValueError('index and complete search path must use the same atom indexing')
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

    def finish(complete, reason):
        return dict(patterns=sorted(patterns.values(), key=lambda p: (p['total'], p['id'])),
                    complete=complete, reason=reason, solver_queries=queries,
                    seconds=time.perf_counter() - started, max_events=max_events,
                    terminal=path.terminal, transitions=list(path.transitions), **metrics)

    def record(mapping, actions):
        pattern = index.describe([mapping[a] for a in range(index.n)])
        pattern['actions'] = actions
        if pattern['id'] not in patterns:
            patterns[pattern['id']] = pattern
            if on_pattern is not None:
                on_pattern(pattern)
        return pattern

    representative = dict(path.mapping)
    score = int(index.counts([[representative[a] for a in range(index.n)]])[0].sum())
    seed = record(representative, []) if max_events is None or score <= max_events else None
    if time.perf_counter() >= deadline:
        return finish(False, 'time_budget')
    if index.invariant_path(path):
        return finish(True, 'signed_event_invariance')
    if time.perf_counter() >= deadline:
        return finish(False, 'time_budget')
    certificates = index.family_certificates
    actions = certificates.actions(path)
    vector = tuple(representative[a] for a in range(index.n))
    if max_events is not None:
        lower = certificates.lower_bound(vector, actions, max_events)
        metrics['reachable_image_event_lower_bound'] = lower
        if lower > max_events:
            return finish(True, 'reachable_image_event_lower_bound')
    if certificates.invariant_class(vector, actions):
        return finish(True, 'source_target_class_invariance')
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
        status = compiled.solver.check()
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
