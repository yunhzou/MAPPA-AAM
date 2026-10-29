"""Signed-event decoding on the mapped region and its boundary.

Search constraints stay at their recorded tolerance. Output identity uses
raw-WBO event responses, with no bijection or group-element expansion.
The existing Golden floor-based scorer and historical pattern IDs are separate.
"""
from collections import Counter, defaultdict, OrderedDict
from functools import cached_property
import hashlib
import json
import math

import numpy as np
import pynauty

from mappa.mappa_core.utils.bonds import is_metal_element, bond_event_threshold


def _response_bounds(values, weight, threshold):
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


class SignedEventIndex:
    """Reaction-local event identity; IDs are versioned separately from legacy.

    Endpoints have nonnegative, exactly symmetric WBOs. Mappings are injective
    and element-preserving; endpoint compositions may differ. Count differences
    touching mapped atoms and ignore bonds entirely within unmatched material.
    """
    schema = 'mappa.signed_events/v3'

    def __init__(self, raw_problem, *, threshold=0.5, metal_threshold=0.3):
        for value in (threshold, metal_threshold):
            if value is not None and (not math.isfinite(value) or value <= 0):
                raise ValueError('event thresholds must be finite and positive')
        if threshold is None:
            raise ValueError('ordinary threshold is required')
        self.problem = raw_problem
        self.elements = raw_problem.reactant.elements
        self.n = len(self.elements)
        self.m = len(raw_problem.product.elements)
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
        self.pa, self.pb = np.triu_indices(self.m, 1)
        target_metal = np.array([is_metal_element(e) for e in raw_problem.product.elements])
        self.ptau = np.full(len(self.pa), self.threshold)
        if self.metal_threshold is not None:
            self.ptau[target_metal[self.pa] | target_metal[self.pb]] = self.metal_threshold
        self.padded_product = np.pad(self.p, ((0, 1), (0, 1)))
        self.cache = {}
        self._actions = OrderedDict()
        self._pair_tables = OrderedDict()
        self._event_tables = OrderedDict()
        self._domain_weights = OrderedDict()
        self._response_values = OrderedDict()
        self._orbit_tables = OrderedDict()

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
        signs = np.where(delta < -tau, -1, np.where(delta > tau, 1, 0))
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

    def _event_orbits(self, action):
        if action not in self._orbit_tables:
            self._orbit_tables[action] = _action_orbits(action, self.m)
            if len(self._orbit_tables) > 1024:
                self._orbit_tables.popitem(last=False)
        return self._orbit_tables[action]

    def _event_response(self, weight, tau, xs, ys):
        from .records import BondEventResponse
        xs, ys = sorted((xs, ys))
        key = float(weight), float(tau), xs, ys
        if key not in self._event_tables:
            domains = xs, ys
            if domains not in self._domain_weights:
                values = set()
                x, y = np.array(xs), np.array(ys)
                block = max(1, 65536 // len(y))
                for start in range(0, len(x), block):
                    xx = x[start:start+block]
                    valid = (xx[:,None] != y[None,:]) | (xx[:,None] == -1) | (y[None,:] == -1)
                    values.update(np.unique(self.padded_product[np.ix_(xx,y)][valid]).tolist())
                weights = np.array(sorted(values))
                weights.setflags(write=False)
                self._domain_weights[domains] = weights
                if len(self._domain_weights) > 4096:
                    self._domain_weights.popitem(last=False)
            weights = self._domain_weights[domains]
            values_key = float(weight), float(tau), tuple(weights)
            if values_key not in self._response_values:
                codes = _event_categories(weight, weights, tau)
                codes.setflags(write=False)
                self._response_values[values_key] = BondEventResponse(float(weight), float(tau), weights, codes)
                if len(self._response_values) > 4096:
                    self._response_values.popitem(last=False)
            self._event_tables[key] = self._response_values[values_key]
            if len(self._event_tables) > 16384:
                self._event_tables.popitem(last=False)
        return self._event_tables[key]

    @property
    def policy(self):
        return dict(schema=self.schema, threshold=self.threshold,
                    metal_threshold=self.metal_threshold, includes_hydrogens=True,
                    event_comparison="absolute difference > threshold",
                    unmatched="ignore bonds with two unmapped endpoints",
                    event_frame="R indices; unmatched P atom p is encoded as source_atom_count+p",
                    source_atom_count=self.n, target_atom_count=self.m)

    def vector(self, mapping):
        """A missing source image is -1; this is bookkeeping, never an atom."""
        return [mapping.get(a, -1) for a in range(self.n)]

    @staticmethod
    def mapping(vector):
        return {a: int(p) for a, p in enumerate(vector) if p >= 0}

    def _vectors(self, vectors):
        values = np.asarray(vectors)
        if values.size == 0:
            if values.ndim == 2 and values.shape[1] == self.n == 0:
                return values.astype(int)
            if values.shape == (0,) or values.shape == (0, self.n):
                return np.empty((0, self.n), dtype=int)
        if values.ndim != 2 or values.shape[1] != self.n or values.dtype.kind not in 'iu':
            raise ValueError('expected integer mapping vectors, with -1 for unmatched atoms')
        if np.any(values < -1) or np.any(values >= self.m):
            raise ValueError('mapping image out of range')
        target_elements = np.asarray(self.problem.product.elements)
        source_elements = np.asarray(self.elements)
        for vector in values:
            mapped = vector >= 0
            images = vector[mapped]
            if len(set(images)) != len(images):
                raise ValueError('mappings must be injective')
            if not np.all(target_elements[images] == source_elements[mapped]):
                raise ValueError('mappings must preserve elements')
        return values

    def counts(self, vectors):
        """Vectorized bond differences; bounded batches also cover partial maps."""
        values = self._vectors(vectors)
        present = values >= 0
        both = present[:, self.a] & present[:, self.b]
        touching = present[:, self.a] | present[:, self.b]
        delta = self.padded_product[values[:, self.a], values[:, self.b]] - self.r[self.a, self.b]
        broken = ((delta < -self.tau) & touching).sum(axis=1)
        formed = ((delta > self.tau) & both).sum(axis=1)
        if self.n != self.m or not np.all(present):
            images = np.zeros((len(values), self.m), dtype=bool)
            rows, atoms = np.nonzero(present)
            images[rows, values[rows, atoms]] = True
            boundary = images[:, self.pa] ^ images[:, self.pb]
            formed += (boundary & (self.p[self.pa, self.pb] > self.ptau)).sum(axis=1)
        return np.column_stack((broken, formed))

    @cached_property
    def response_data(self):
        endpoints = (self.problem.reactant, self.problem.product)
        result = []
        for side in (0, 1):
            ep, other = endpoints[side], endpoints[1 - side]
            opposite_values = defaultdict(set)
            for a, b in zip(*np.triu_indices(other.atom_count, 1)):
                pair = tuple(sorted((other.elements[a], other.elements[b])))
                opposite_values[pair].add(float(other.wbo[a, b]))
            opposite_values = {k: tuple(sorted(v | {0.})) for k, v in opposite_values.items()}
            weights = defaultdict(set)
            for a, b in zip(*np.triu_indices(ep.atom_count, 1)):
                pair = tuple(sorted((ep.elements[a], ep.elements[b])))
                weights[pair].add(float(ep.wbo[a, b]))
            memo = {}
            for pair, values in weights.items():
                values = sorted(values)
                tau = (self.metal_threshold if self.metal_threshold is not None
                       and any(is_metal_element(e) for e in pair) else self.threshold)
                bounds = _response_bounds(opposite_values.get(pair, (0.,)), values, tau)
                memo.update(((pair, weight), (int(lo), int(hi)))
                            for weight, lo, hi in zip(values, *bounds))
            pairs, counts = [], defaultdict(Counter)
            for a, b in zip(*np.triu_indices(ep.atom_count, 1)):
                pair = tuple(sorted((ep.elements[a], ep.elements[b])))
                key = pair, float(ep.wbo[a, b])
                label = memo[key]
                pairs.append((int(a), int(b), pair, label))
                counts[pair][label] += 1
            defaults = {pair: min(c, key=lambda label: (-c[label], label))
                        for pair, c in counts.items()}
            labels = np.full((ep.atom_count, ep.atom_count), -1, dtype=int)
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
        return self.target_action_responses((tuple(images),))[0][0]

    def target_action_responses(self, generators, *, batch_entries=65536):
        """Certify generators over every bond response, never group elements.

        Return (preserves responses, changed pair count) per generator. Failure
        is conservative, not proof of different attainable event patterns.
        Temporary pair arrays have at most ``batch_entries`` entries each.
        """
        if batch_entries < 1:
            raise ValueError('batch_entries must be positive')
        generators = tuple(tuple(g) for g in generators)
        labels = self.response_data[1]['labels']
        elements = np.asarray(self.problem.product.elements)
        result = []
        # Bound both the generator matrix and the generator-by-pair scratch.
        batch_size = max(1, min(64, batch_entries // max(1, self.m)))
        for start in range(0, len(generators), batch_size):
            batch = generators[start:start + batch_size]
            missing = tuple(dict.fromkeys(g for g in batch if g not in self._actions))
            if missing:
                images = np.asarray(missing)
                if (images.shape != (len(missing), self.m) or (self.m and images.dtype.kind not in 'iu')
                        or not np.all(np.sort(images, axis=1) == np.arange(self.m))):
                    raise ValueError('invalid saved target permutation')
                images = images.astype(int, copy=False)
                changed = np.zeros(len(missing), dtype=int)
                moved = np.any(images != np.arange(self.m), axis=0)
                relevant = np.flatnonzero(moved[self.pa] | moved[self.pb])
                pair_batch = max(1, batch_entries // len(missing))
                for offset in range(0, len(relevant), pair_batch):
                    pairs = relevant[offset:offset + pair_batch]
                    a, b = self.pa[pairs], self.pb[pairs]
                    changed += np.count_nonzero(
                        labels[images[:, a], images[:, b]] != labels[a, b], axis=1)
                same_elements = np.all(elements[images] == elements, axis=1)
                for g, count, same in zip(missing, changed, same_elements):
                    self._actions[g] = (bool(same and count == 0), int(count))
            result.extend(self._actions[g] for g in batch)
            while len(self._actions) > 4096:
                self._actions.popitem(last=False)
        return tuple(result)

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
                    images = list(range(self.m))
                    images[a], images[b] = images[b], images[a]
                    if not self.preserves_target_action(images):
                        return False
        return True

    def describe(self, vector):
        vector = self._vectors([vector])[0]
        mapping = self.mapping(vector)
        broken, formed, _, _ = classify_bonds(mapping, self.r, self.p,
            dwbo_threshold=self.threshold, elements_R=self.elements,
            elements_P=self.problem.product.elements, metal_dwbo_threshold=self.metal_threshold)
        result = self._describe_events(vector, tuple((a, b) for a, b, *_ in broken),
                                       tuple((a, b) for a, b, *_ in formed))
        if tuple(self.counts([vector])[0]) != tuple(result['counts'].values()):
            raise AssertionError('vector and scalar event witnesses disagree')
        return result

    def _describe_events(self, vector, broken, formed):
        """Canonicalize already counted endpoint events without recounting bonds."""
        vector = self._vectors([vector])[0]
        mapping = self.mapping(vector)
        inverse = {p: r for r, p in mapping.items()}
        # Existing R coordinates remain valid for mapped pairs. New boundary
        # endpoints name actual P atoms explicitly, beyond the R index range.
        events = (tuple(sorted((a, b) for a, b, *_ in broken)),
                  tuple(sorted(tuple(sorted((inverse.get(a, self.n+a), inverse.get(b, self.n+b))))
                               for a, b, *_ in formed)))
        mapped = tuple(sorted(mapping))
        key = mapped, events
        if key not in self.cache:
            base = self.base
            adjacency = {i: set(ns) for i, ns in base.adjacency_dict.items()}
            colors = []
            # Event classes retain which source region is mapped, modulo exact
            # source symmetries. Unmatched target structure is outside the score.
            for cell in base.vertex_coloring:
                for matched in (False, True):
                    group = {a for a in cell if (a in mapping) == matched}
                    if group:
                        colors.append(group)
            n = base.number_of_vertices
            adjacency[n] = set(mapped)
            for a in mapped:
                adjacency[a].add(n)
            colors.append({n})
            n += 1
            extra = sorted({a for pairs in events for pair in pairs for a in pair if a >= self.n})
            vertices, extra_colors = {}, defaultdict(set)
            for a in extra:
                vertices[a] = n
                adjacency[n] = set()
                extra_colors[self.problem.product.elements[a-self.n]].add(n)
                n += 1
            colors.extend(extra_colors[k] for k in sorted(extra_colors))
            for pairs in events:
                marker = n
                n += 1
                adjacency[marker] = set()
                colors.append({marker})
                gadgets = set()
                for a, b in pairs:
                    a, b = vertices.get(a, a), vertices.get(b, b)
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
            self.cache[key] = self.schema + ':' + hashlib.sha256(payload).hexdigest()
        counts = dict(zip(('broken', 'formed'), map(len, events)))
        return dict(id=self.cache[key], mapping=vector.tolist(), mapped_source=mapped,
                    events=dict(zip(('broken', 'formed'), events)),
                    endpoint_events=dict(broken=tuple((a,b) for a,b,*_ in broken),
                                         formed=tuple((a,b) for a,b,*_ in formed)),
                    counts=counts, total=sum(counts.values()), policy=self.policy)

EVENT_CATEGORIES = ('breaking', 'forming', 'strengthening', 'weakening')


def _event_categories(source, target, threshold):
    """0 unchanged; 1..4 are exclusive categories using the existing bond cutoff."""
    source, target = np.asarray(source), np.asarray(target)
    decrease, increase = source - target > threshold, target - source > threshold
    return np.where(decrease, np.where(target > threshold, 4, 1),
                    np.where(increase, np.where(source > threshold, 3, 2), 0)).astype(np.uint8)


def _action_orbits(action, degree):
    """Atom orbits from generator edges; never construct permutation elements."""
    kind, data = action
    parent = list(range(degree))

    def root(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    edges = zip(data, data[1:]) if kind == 'pool' else (
        (a, b) for generator in data for a, b in enumerate(generator) if a != b)
    for a, b in edges:
        parent[root(a)] = root(b)
    cells = defaultdict(list)
    for a in range(degree):
        cells[root(a)].append(a)
    result = [None] * degree
    for cell in cells.values():
        orbit = frozenset(cell)
        for a in cell:
            result[a] = orbit
    return tuple(result)


class FragmentEventTables:
    """Sparse internal/boundary factors over ORIGINAL ordered family choices.

    Single-atom domains only screen impossible events. Their Cartesian product
    is an overapproximation, never a permission to combine independent images.
    Evaluation receives a realization checked by EventInterpretation.realize.
    """
    def __init__(self, index, family, fragments):
        from .records import BondEventTerm
        self.index, self.fragments = index, fragments
        owners = {a: i for i, atoms in enumerate(fragments) for a in atoms}
        mapping = dict(family.mapping)
        if len(owners) != sum(map(len, fragments)) or set(owners) != set(mapping):
            raise ValueError('saved fragments must partition the mapped source atoms')
        self.owners = np.array([owners.get(a, -1) for a in range(index.n)], dtype=int)
        domains = {a: frozenset((p,)) for a, p in family.mapping}
        dependencies = {a: set() for a in domains}
        for aid, action in enumerate(family.actions):
            orbits = index._event_orbits(action)
            for a, domain in domains.items():
                if any(len(orbits[p]) > 1 for p in domain):
                    dependencies[a].add(aid)
                    domains[a] = frozenset().union(*(orbits[p] for p in domain))
        self.atom_domains = {a: tuple(sorted(ds)) for a, ds in domains.items()}
        self.atom_actions = {a: tuple(sorted(ids)) for a, ids in dependencies.items()}
        self.internal, self.boundary = {}, {}
        groups = defaultdict(list)
        self._constant = Counter()
        # Source pairs cover both ordinary and manually cut boundaries. Product
        # bonds are included even when their source weight is zero.
        for a, b, tau in zip(index.a, index.b, index.tau):
            if a not in mapping and b not in mapping:
                continue
            response = index._event_response(index.r[a, b], tau,
                self.atom_domains.get(a, (-1,)), self.atom_domains.get(b, (-1,)))
            if response.constant == 0:
                continue
            left, right = int(self.owners[a]), int(self.owners[b])
            region = tuple(sorted((left, right)))
            term = BondEventTerm((int(a), int(b)),
                tuple(sorted(dependencies.get(a, set()) | dependencies.get(b, set()))), response)
            collection = self.internal if left == right else self.boundary
            key = left if left == right else region
            collection.setdefault(key, []).append(term)
            if response.constant is None:
                groups[response].append((int(a), int(b), region))
            else:
                self._constant[region, response.constant] += 1
        self.internal = {key: tuple(terms) for key, terms in self.internal.items()}
        self.boundary = {key: tuple(terms) for key, terms in self.boundary.items()}
        regions = sorted({region for terms in groups.values() for _,_,region in terms}
                         | {region for region, _ in self._constant})
        self._regions = tuple(regions)
        region_ids = {region:i for i,region in enumerate(regions)}
        self._fixed = np.zeros((len(regions), 5), dtype=np.int64)
        for (region, code), count in self._constant.items():
            self._fixed[region_ids[region], code] += count
        self._groups = tuple((response, np.array([(a,b) for a,b,_ in terms], dtype=int),
                              np.array([region_ids[region] for _,_,region in terms], dtype=int))
                             for response, terms in groups.items())
        # Product-to-unmatched boundaries depend on occupied target atoms, not
        # a fabricated correspondence between unmatched reactants and products.
        eligible = index.p[index.pa, index.pb] > index.ptau
        if len(mapping) == index.m:
            eligible[:] = False
        reachable = set().union(*domains.values()) if domains else set()
        eligible &= np.array([a in reachable or b in reachable for a,b in zip(index.pa,index.pb)])
        self.product_boundary_pairs = np.column_stack((index.pa[eligible], index.pb[eligible]))
        self.product_boundary_actions = tuple(range(len(family.actions)))

    def evaluate(self, vector, *, include_events=False):
        """Lookup only the requested joint choice; sum each pair once."""
        vector = np.asarray(vector, dtype=int)
        totals = self._fixed.copy()
        for response, pairs, regions in self._groups:
            codes = response.lookup(self.index.padded_product[vector[pairs[:,0]], vector[pairs[:,1]]])
            np.add.at(totals, (regions, codes), 1)
        inverse = np.full(self.index.m, -1, dtype=int)
        mapped = np.flatnonzero(vector >= 0)
        inverse[vector[mapped]] = mapped
        pairs = self.product_boundary_pairs
        source = inverse[pairs]
        boundary = (source[:,0] >= 0) ^ (source[:,1] >= 0)
        target_boundary = pairs[boundary]
        internal, boundary_result = {}, {}
        counts = totals[:,1:].sum(axis=0)
        for (left, right), row in zip(self._regions, totals[:,1:]):
            if np.any(row):
                collection = internal if left == right else boundary_result
                key = left if left == right else (left, right)
                collection[key] = row.copy()
        formed_owners = self.owners[np.maximum(source[boundary,0], source[boundary,1])]
        for owner, count in zip(*np.unique(formed_owners, return_counts=True)):
            boundary_result.setdefault((-1, int(owner)), np.zeros(4, dtype=int))[1] += count
            counts[1] += count
        def named(counts):
            return dict(zip(EVENT_CATEGORIES, map(int, counts)))
        result = dict(counts=named(counts), total=int(counts.sum()),
                      internal={k:named(v) for k,v in internal.items()},
                      boundary={k:named(v) for k,v in boundary_result.items()})
        if include_events:
            events = []
            for terms in (*self.internal.values(), *self.boundary.values()):
                for term in terms:
                    a,b = term.atoms
                    x,y = int(vector[a]), int(vector[b])
                    weight = self.index.padded_product[x,y]
                    code = int(term.response.lookup(weight))
                    if code:
                        events.append(dict(category=EVENT_CATEGORIES[code-1], reactant_pair=(a,b),
                            product_pair=tuple(sorted((x,y))) if x >= 0 and y >= 0 else None,
                            reactant_weight=float(self.index.r[a,b]), product_weight=float(weight)))
            events.extend(dict(category='forming', reactant_pair=None,
                product_pair=(int(a),int(b)), reactant_weight=0., product_weight=float(self.index.p[a,b]))
                for a,b in target_boundary)
            result['events'] = tuple(events)
        return result

    def counts(self, vectors):
        """Vectorized category lookup for caller-supplied, validated realizations."""
        vectors = np.asarray(vectors, dtype=int)
        result = np.broadcast_to(self._fixed[:,1:].sum(axis=0), (len(vectors), 4)).copy()
        for response, pairs, _ in self._groups:
            values = self.index.padded_product[vectors[:,pairs[:,0]], vectors[:,pairs[:,1]]]
            codes = response.lookup(values)
            for code in range(1, 5):
                result[:,code-1] += np.count_nonzero(codes == code, axis=1)
        pairs = self.product_boundary_pairs
        if len(pairs):
            present = np.zeros((len(vectors), self.index.m), dtype=bool)
            rows, atoms = np.nonzero(vectors >= 0)
            present[rows, vectors[rows,atoms]] = True
            result[:,1] += np.count_nonzero(present[:,pairs[:,0]] ^ present[:,pairs[:,1]], axis=1)
        return result

    @property
    def statistics(self):
        terms = [t for rows in (*self.internal.values(), *self.boundary.values()) for t in rows]
        return dict(fragments=len(self.fragments), internal_terms=sum(map(len,self.internal.values())),
                    boundary_terms=sum(map(len,self.boundary.values())),
                    variable_terms=sum(t.response.constant is None for t in terms),
                    constant_terms=sum(t.response.constant is not None for t in terms),
                    unique_responses=len({id(t.response) for t in terms}),
                    product_boundary_pairs=len(self.product_boundary_pairs))


class FragmentEventCounter:
    """Count each local alternative and interacting pair once per decoding call.

    Regions are converted fragment atom sets, including any internal splits.
    Unmapped material is a separate region: only its boundary is counted.
    No complete-mapping or Cartesian-combination cache is stored.
    """
    def __init__(self, index):
        self.index = index
        self.regions, self.region_ids, self.local = [], {}, []
        self.domains, self.domain_ids = [], {}
        self.topologies, self.boundaries = {}, {}
        self.bond_comparisons = self.branch_counts = 0
        self.edges = []
        for endpoint in (index.problem.reactant, index.problem.product):
            adjacency = [[] for _ in endpoint.elements]
            for a in range(endpoint.atom_count):
                for b in range(a + 1, endpoint.atom_count):
                    threshold = bond_event_threshold(endpoint.elements, a, b,
                        default_threshold=index.threshold, metal_threshold=index.metal_threshold)
                    weight = endpoint.wbo[a, b]
                    if weight > threshold:
                        adjacency[a].append((b, weight, threshold))
                        adjacency[b].append((a, weight, threshold))
            self.edges.append(adjacency)

    def _register(self, source, target, *, unmapped=False):
        key = tuple(sorted(source.items())), tuple(sorted(target.items()))
        if key not in self.region_ids:
            region = len(self.regions)
            self.region_ids[key] = region
            self.regions.append((source, target))
            domain = tuple(sorted(source)), tuple(sorted(target))
            self.domains.append(self.domain_ids.setdefault(domain, len(self.domain_ids)))
            self.local.append(((), ()) if unmapped else self._events(region, region))
        return self.region_ids[key]

    def _events(self, left, right):
        result = []
        for side, opposite in enumerate((self.index.p, self.index.r)):
            images, other = self.regions[left][side], self.regions[right][side]
            events = []
            for a, image in images.items():
                for b, weight, threshold in self.edges[side][a]:
                    if b not in other or (left == right and a >= b):
                        continue
                    partner = other[b]
                    if image is None and partner is None:
                        continue
                    self.bond_comparisons += 1
                    compared = 0. if image is None or partner is None else opposite[image, partner]
                    if weight - compared > threshold:
                        events.append((min(a, b), max(a, b)))
            result.append(tuple(sorted(events)))
        return tuple(result)

    def _partition(self, mapping, islands):
        owners = dict(islands)
        fragments = {}
        for a, b in mapping:
            fragments.setdefault(owners[a], {})[a] = b
        return tuple(self._register(part, {b: a for a, b in part.items()})
                     for part in fragments.values())

    def _topology(self, regions):
        key = tuple(self.domains[i] for i in regions)
        if key not in self.topologies:
            owners = [{a: i for i, region in enumerate(regions) for a in self.regions[region][side]}
                      for side in (0, 1)]
            absent = [{a: None for a in range(n) if a not in own}
                      for n, own in zip((self.index.n, self.index.m), owners)]
            unmatched = self._register(*absent, unmapped=True)
            pairs = set()
            for adjacency, own in zip(self.edges, owners):
                for a, edges in enumerate(adjacency):
                    for b, _, _ in edges:
                        left, right = own.get(a, len(regions)), own.get(b, len(regions))
                        if a < b and left != right:
                            pairs.add(tuple(sorted((left, right))))
            self.topologies[key] = unmatched, tuple(sorted(pairs))
        return self.topologies[key]

    def _contributions(self, mapping, islands):
        regions = self._partition(mapping, islands)
        unmatched, pairs = self._topology(regions)
        local = tuple(self.local[i] for i in regions)
        regions += (unmatched,)
        boundary = []
        for a, b in pairs:
            key = tuple(sorted((regions[a], regions[b])))
            if key not in self.boundaries:
                self.boundaries[key] = self._events(*key)
            boundary.append(self.boundaries[key])
        return local, tuple(boundary)

    def count(self, mapping, islands):
        local, boundary = self._contributions(mapping, islands)
        internal = sum(len(broken) + len(formed) for broken, formed in local)
        total = internal + sum(len(broken) + len(formed) for broken, formed in boundary)
        self.branch_counts += 1
        return total, internal

    def events(self, mapping, islands):
        local, boundary = self._contributions(mapping, islands)
        return tuple(tuple(sorted(pair for contribution in local + boundary
                                  for pair in contribution[side])) for side in (0, 1))

    @property
    def statistics(self):
        return dict(regions=len(self.regions), interacting_alternative_pairs=len(self.boundaries),
                    partitions=len(self.topologies), bond_comparisons=self.bond_comparisons,
                    branch_counts=self.branch_counts)


def classify_bonds(mapping, wbo_R, wbo_P, dwbo_threshold=0.5,
                   elements_R=None, elements_P=None,
                   metal_dwbo_threshold=None):
    """Bond classification by WBO change.

        broken iff (WBO_R - WBO_P) > pair_threshold
        formed iff (WBO_P - WBO_R) > pair_threshold

    `pair_threshold` is normally `dwbo_threshold`.  If element lists and
    `metal_dwbo_threshold` are supplied, any pair containing a metal uses the
    metal threshold instead.  Since WBO is non-negative, requiring
    (wR - wP) > pair_threshold automatically implies wR > pair_threshold
    (wP >= 0), so a single
    threshold both gates "is this a bond worth considering" and "did
    its order change enough." For pairs with exactly one endpoint unmapped, compare the bond with
    zero on the other side. Ignore pairs whose endpoints are both unmapped,
    symmetrically in R and P.

    Returns (broken_list, formed_list, core_R, core_P) where each bond
    record is (i, j, wbo_R_or_None, wbo_P_or_None)."""
    inv = {v: k for k, v in mapping.items()}
    default_threshold = float(dwbo_threshold)
    metal_threshold = (
        None if metal_dwbo_threshold is None
        else float(metal_dwbo_threshold))
    metal_R = (
        None if elements_R is None or metal_threshold is None
        else tuple(is_metal_element(element) for element in elements_R))
    metal_P = (
        None if elements_P is None or metal_threshold is None
        else tuple(is_metal_element(element) for element in elements_P))
    def eligible(matrix, metal):
        # Screen absent/subthreshold pairs in NumPy. The scalar subtraction
        # below remains the event oracle, including exact floating boundaries.
        threshold = default_threshold
        if metal is not None:
            flags = np.asarray(metal, dtype=bool)
            threshold = np.where(flags[:, None] | flags[None, :],
                                 metal_threshold, default_threshold)
        a, b = np.nonzero(np.triu(matrix > threshold, k=1))
        return zip(a.tolist(), b.tolist())

    broken, formed = [], []
    for i, j in eligible(wbo_R, metal_R):
        if i not in mapping and j not in mapping:
            continue
        threshold = (metal_threshold if metal_R is not None and (metal_R[i] or metal_R[j])
                     else default_threshold)
        wR = wbo_R[i, j]
        if i not in mapping or j not in mapping:
            broken.append((i, j, float(wR), None)); continue
        wP = wbo_P[mapping[i], mapping[j]]
        if wR - wP > threshold:
            broken.append((i, j, float(wR), float(wP)))
    for ip, jp in eligible(wbo_P, metal_P):
        if ip not in inv and jp not in inv:
            continue
        threshold = (metal_threshold if metal_P is not None and (metal_P[ip] or metal_P[jp])
                     else default_threshold)
        wP = wbo_P[ip, jp]
        if ip not in inv or jp not in inv:
            formed.append((ip, jp, None, float(wP))); continue
        wR = wbo_R[inv[ip], inv[jp]]
        if wP - wR > threshold:
            formed.append((ip, jp, float(wR), float(wP)))
    core_R = sorted({i for (i, j, _, _) in broken}
                    | {j for (i, j, _, _) in broken})
    core_P = sorted({i for (i, j, _, _) in formed}
                    | {j for (i, j, _, _) in formed})
    return broken, formed, core_R, core_P
