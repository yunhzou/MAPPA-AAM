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
import time

import numpy as np
import pynauty

from mappa.frag import classify_bonds, is_metal_element


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

    Endpoints have nonnegative, exactly symmetric WBOs. Mappings are injective
    and element-preserving; endpoint compositions may differ. Count differences
    touching mapped atoms and ignore bonds entirely within unmatched material.
    """
    schema = 'mappa.signed_events/v2'

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
                    metal_threshold=self.metal_threshold, includes_hydrogens=True,
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
        broken = ((delta <= -self.tau) & touching).sum(axis=1)
        formed = ((delta >= self.tau) & both).sum(axis=1)
        if self.n != self.m or not np.all(present):
            images = np.zeros((len(values), self.m), dtype=bool)
            rows, atoms = np.nonzero(present)
            images[rows, values[rows, atoms]] = True
            boundary = images[:, self.pa] ^ images[:, self.pb]
            formed += (boundary & (self.p[self.pa, self.pb] >= self.ptau)).sum(axis=1)
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
            memo, pairs, counts = {}, [], defaultdict(Counter)
            for a, b in zip(*np.triu_indices(ep.atom_count, 1)):
                pair = tuple(sorted((ep.elements[a], ep.elements[b])))
                tau = (self.metal_threshold if self.metal_threshold is not None
                       and any(is_metal_element(e) for e in pair) else self.threshold)
                key = pair, float(ep.wbo[a, b])
                if key not in memo:
                    memo[key] = _response_bounds(opposite_values.get(pair, (0.,)), key[1], tau)
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
        images = tuple(images)
        if images not in self._actions:
            if len(images) != self.m or sorted(images) != list(range(self.m)):
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
                    images = list(range(self.m))
                    images[a], images[b] = images[b], images[a]
                    if not self.preserves_target_action(images):
                        return False
        return True

    def describe(self, vector):
        vector = self._vectors([vector])[0]
        mapping = self.mapping(vector)
        inverse = {p: r for r, p in mapping.items()}
        broken, formed, _, _ = classify_bonds(mapping, self.r, self.p,
            dwbo_threshold=self.threshold, elements_R=self.elements,
            elements_P=self.problem.product.elements, metal_dwbo_threshold=self.metal_threshold)
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
        if tuple(self.counts([vector])[0]) != tuple(counts.values()):
            raise AssertionError('vector and scalar event witnesses disagree')
        return dict(id=self.cache[key], mapping=vector.tolist(), mapped_source=mapped,
                    events=dict(zip(('broken', 'formed'), events)),
                    endpoint_events=dict(broken=tuple((a,b) for a,b,*_ in broken),
                                         formed=tuple((a,b) for a,b,*_ in formed)),
                    counts=counts, total=sum(counts.values()), policy=self.policy)
