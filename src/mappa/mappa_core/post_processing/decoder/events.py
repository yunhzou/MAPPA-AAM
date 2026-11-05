"""Signed-event decoding on the mapped region and its boundary.

Saved fragment relations retain their joint constraints. Output identity uses
raw-WBO event responses, with no bijection or group-element expansion.
The existing Golden floor-based scorer and historical pattern IDs are separate.
"""
from collections import Counter, defaultdict, OrderedDict
from dataclasses import dataclass
from functools import cached_property
import hashlib
import json
import math

import numpy as np
import pynauty

from mappa.mappa_core.utils.bonds import is_metal_element, bond_event_threshold
from mappa.mappa_core.core.symmetry.responses import BondResponseIndex, response_bounds


class SignedEventIndex(BondResponseIndex):
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
        super().__init__(raw_problem.product.wbo)
        self.r = raw_problem.reactant.wbo
        if not np.array_equal(self.r, self.r.T) or np.any(self.r < 0):
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
        self.reactant_bonds = np.column_stack((self.a, self.b))[
            self.r[self.a, self.b] > self.tau]
        self.product_bonds = np.column_stack((self.pa, self.pb))[
            self.p[self.pa, self.pb] > self.ptau]
        self.cache = OrderedDict()
        self._actions = OrderedDict()
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

    def neighboring_pairs(self, domains, *, owners=None):
        """Source bonds plus source pairs that can reach a product bond.

        Nonnegative weights can differ by more than the event threshold only
        if at least one endpoint bond exceeds it. Inverse atom domains retain
        possible forming bonds without visiting unrelated source nonbonds.
        Bitset unions visit each reachable source pair once. When owners are
        supplied, only cross-fragment pairs are assembled; unique fragment
        preparation already owns the internal pairs. Joint action constraints
        remain attached to every retained response.
        """
        pairs = {(int(a), int(b)) for a, b in self.reactant_bonds
                 if (a in domains or b in domains)
                 and (owners is None or owners[a] != owners[b])}
        sources = [0] * self.m
        for atom, images in domains.items():
            for image in images:
                sources[image] |= 1 << int(atom)
        neighbors = [0] * self.m
        for x, y in self.product_bonds:
            neighbors[x] |= sources[y]
            neighbors[y] |= sources[x]
        owner_masks = {}
        if owners is not None:
            for atom in domains:
                owner = owners[atom]
                owner_masks[owner] = owner_masks.get(owner, 0) | (1 << int(atom))
        domain_neighbors = {}
        for a, images in domains.items():
            a = int(a)
            images = tuple(images)
            if images not in domain_neighbors:
                reached = 0
                for image in images:
                    reached |= neighbors[image]
                domain_neighbors[images] = reached
            reached = domain_neighbors[images] & ~((1 << (a + 1)) - 1)
            if owners is not None:
                reached &= ~owner_masks[owners[a]]
            while reached:
                bit = reached & -reached
                pairs.add((a, bit.bit_length() - 1))
                reached ^= bit
        return pairs

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
                bounds = response_bounds(opposite_values.get(pair, (0.,)), values, tau)
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
        elements = np.asarray(self.problem.product.elements)
        result = []

        def changed_pairs(matrix, images):
            changed = np.zeros(len(images), dtype=int)
            moved = np.any(images != np.arange(self.m), axis=0)
            relevant = np.flatnonzero(moved[self.pa] | moved[self.pb])
            pair_batch = max(1, batch_entries // len(images))
            for offset in range(0, len(relevant), pair_batch):
                pairs = relevant[offset:offset + pair_batch]
                a, b = self.pa[pairs], self.pb[pairs]
                changed += np.count_nonzero(
                    matrix[images[:, a], images[:, b]] != matrix[a, b], axis=1)
            return changed

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
                # Exact original weights certify every event response directly;
                # build response labels only for actions that change a weight.
                same_elements = np.all(elements[images] == elements, axis=1)
                varying = np.flatnonzero((changed_pairs(self.p, images) > 0) | ~same_elements)
                changed = np.zeros(len(missing), dtype=int)
                if len(varying):
                    changed[varying] = changed_pairs(self.response_data[1]['labels'], images[varying])
                for g, count, same in zip(missing, changed, same_elements):
                    self._actions[g] = (bool(same and count == 0), int(count))
            result.extend(self._actions[g] for g in batch)
            while len(self._actions) > 4096:
                self._actions.popitem(last=False)
        return tuple(result)

    def invariant_path(self, path):
        """Certify a fixed path; return None when relation responses are unknown.

        A self generator certificate cannot cover all accepted discovery maps.
        Conservative domains certify identity only when no atom can move.
        """
        for edge in path.transitions:
            transition = path.graph.transitions[edge]
            fragment = transition.match
            if fragment is not None:
                domains = self.action_domains(('relation', fragment))
                if any(len(domains[p]) > 1 for p in fragment.target_atoms):
                    return None
            if any(not preserves for preserves, _ in
                   self.target_action_responses(transition.placement_generators)):
                return None
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
            if len(self.cache) > 4096:
                self.cache.popitem(last=False)
        counts = dict(zip(('broken', 'formed'), map(len, events)))
        return dict(id=self.cache[key], mapping=vector.tolist(), mapped_source=mapped,
                    events=dict(zip(('broken', 'formed'), events)),
                    endpoint_events=dict(broken=tuple((a,b) for a,b,*_ in broken),
                                         formed=tuple((a,b) for a,b,*_ in formed)),
                    counts=counts, total=sum(counts.values()), policy=self.policy)

EVENT_CATEGORIES = ('breaking', 'forming', 'strengthening', 'weakening')


@dataclass(frozen=True)
class BoundaryEvents:
    """Fixed boundary ledger and the reservation for one swept correspondence."""
    source_owners: np.ndarray
    target_owners: np.ndarray
    source_pairs: np.ndarray
    target_pairs: np.ndarray
    fixed_counts: dict
    cut_bonds: frozenset
    sweep_cut: tuple | None
    source_cut_counted: bool
    sweep_removed: int
    sweep_width: int

    @property
    def count(self):
        return len(self.source_pairs) + len(self.target_pairs)

    @property
    def lower(self):
        return self.count - self.sweep_removed

    @cached_property
    def target_pair_set(self):
        return frozenset(map(tuple, self.target_pairs))


def prepare_boundary_events(index, mapping, pieces, cut_bonds=()):
    """Prepare fixed ordinary events without preparing any local symmetry."""
    mapping = dict(mapping)
    owners = {a: i for i, piece in enumerate(pieces) for a in piece.source_atoms}
    if len(owners) != sum(len(piece.source_atoms) for piece in pieces) or set(owners) != set(mapping):
        raise ValueError('saved fragments must partition the mapped source atoms')
    cuts = frozenset(tuple(sorted(pair)) for pair in cut_bonds)
    if len(cuts) > 1:
        raise ValueError('each decoder occurrence requires at most one sweep cut')
    source_owners = np.array([owners.get(a, -1) for a in range(index.n)], dtype=int)
    target_owners = np.full(index.m, -1, dtype=int)
    for a, p in mapping.items():
        target_owners[p] = owners[a]
    source_pairs = index.reactant_bonds[
        source_owners[index.reactant_bonds[:, 0]] != source_owners[index.reactant_bonds[:, 1]]]
    target_pairs = index.product_bonds[
        target_owners[index.product_bonds[:, 0]] != target_owners[index.product_bonds[:, 1]]]
    fixed = {}
    for pairs, labels, category in ((source_pairs, source_owners, 0), (target_pairs, target_owners, 1)):
        for a, b in pairs:
            region = tuple(sorted((int(labels[a]), int(labels[b]))))
            fixed.setdefault(region, np.zeros(4, dtype=int))[category] += 1
    cut = next(((a, b) for a, b in cuts
                if a in mapping and b in mapping and owners[a] != owners[b]), None)
    counted = False
    removed = width = 0
    if cut is not None:
        a, b = cut
        tau = bond_event_threshold(index.elements, a, b,
            default_threshold=index.threshold, metal_threshold=index.metal_threshold)
        counted = bool(index.r[a, b] > tau)
        reserve = int(bool(len(target_pairs)))
        removed, width = int(counted) + reserve, 1 + reserve
    return BoundaryEvents(source_owners, target_owners, source_pairs, target_pairs,
                          fixed, cuts, cut, counted, removed, width)


class FragmentEventTables:
    """Shared local responses, fixed ordinary boundaries and one sweep response.

    Ordinary boundary counts are an explicit approximation. Local responses and
    the selected sweep still evaluate the original ordered, validated choices.
    """
    model = 'fixed_ordinary_boundaries_with_sweep'
    approximate = True

    def __init__(self, local, family, fragments, cut_bonds=(), *, boundary=None):
        index = local.index
        self.index = index
        self.local, self.family = local, family
        self.fragments = tuple(piece.source_atoms for piece in fragments)
        mapping = dict(family.mapping)
        boundary = boundary or prepare_boundary_events(index, mapping, fragments, cut_bonds)
        self.cut_bonds = boundary.cut_bonds
        self.boundary_events = boundary
        self.owners, self.target_owners = boundary.source_owners, boundary.target_owners
        self.reactant_boundary_pairs, self.product_boundary_pairs = boundary.source_pairs, boundary.target_pairs
        self.fixed_boundary = boundary.fixed_counts
        self.product_boundary_set = boundary.target_pair_set
        self.boundary_count = boundary.count
        self.sweep_cut = boundary.sweep_cut
        self.sweep_removed, self.sweep_width = boundary.sweep_removed, boundary.sweep_width
        program = local.programs.program(family.actions)
        self.atom_domains, self.atom_actions = local.programs.bind(program, family.mapping)
        self.internal, self.boundary = {}, {}
        self.fragment_responses = tuple(local.fragment(fragment, mapping, program,
            self.atom_domains, self.atom_actions, self.cut_bonds) for fragment in fragments)
        # A branch binds shared local action positions to its original program.
        # Only these bindings and boundary relationships belong to the branch.
        for owner, (response, binding) in enumerate(self.fragment_responses):
            self.atom_domains.update(response.atom_domains)
            if response.terms:
                self.internal[owner] = local.bind(response.terms, binding)
        self.projected_pair_states = sum(response.projected_pair_states
                                         for response, _ in self.fragment_responses)
        self.neighboring_pair_count = int(self.sweep_cut is not None) + sum(
            response.neighboring_pair_count for response, _ in self.fragment_responses)
        if self.sweep_cut is not None:
            a, b = self.sweep_cut
            left, right = int(self.owners[a]), int(self.owners[b])
            pair = mapping[a], mapping[b]
            action_ids = tuple(sorted(set(self.atom_actions[a]) | set(self.atom_actions[b])))
            actions = local.programs.select(program, action_ids)
            term, projected = local.bond((a, b), pair, actions,
                (self.atom_domains[a], self.atom_domains[b]), action_ids, True)
            self.projected_pair_states += projected
            self.boundary[tuple(sorted((left, right)))] = (term,)
        # Preserve the public record order of the original sorted bond walk.
        self.internal = dict(sorted(self.internal.items(), key=lambda item: item[1][0].atoms))
        self.product_boundary_actions = ()
        self.components = self._dependency_components()
        lower, upper = boundary.lower, boundary.count
        for response, _ in self.fragment_responses:
            _, local_upper = local.bounds(response)
            lower += response.event_counts.lower
            upper += local_upper
        if self.sweep_cut is not None:
            term = next(iter(self.boundary.values()))[0]
            a, b = self.sweep_cut
            upper += int(term.response.lookup(index.p[mapping[a], mapping[b]]) != 0)
            upper -= int(boundary.source_cut_counted)
            upper -= int(tuple(sorted((mapping[a], mapping[b]))) in self.product_boundary_set)
        self.count_bounds = int(lower), int(upper)

    @cached_property
    def program(self):
        """Bind shared fragment costs; only the sweep adds a boundary choice."""
        from mappa.mappa_core.core.symmetry.program import OrderedSymmetryEvents
        local, family = self.local, self.family
        actions = tuple((action, local.prepare(action[1]).response if action[0] == 'relation' else None)
                        for action in family.actions)
        fixed = self.boundary_events.lower
        projections = tuple((response.event_counts.program, binding)
                            for response, binding in self.fragment_responses)
        sweep = None
        if self.sweep_cut is not None:
            a, b = self.sweep_cut
            tau = bond_event_threshold(self.index.elements, a, b,
                default_threshold=self.index.threshold, metal_threshold=self.index.metal_threshold)
            reserve = self.sweep_removed - int(self.boundary_events.source_cut_counted)
            sweep = self.sweep_cut, float(self.index.r[a,b]), tau, reserve, self.product_boundary_set
        return OrderedSymmetryEvents(family.mapping, actions, self.index.p, (),
            fixed_events=fixed, sweep=sweep, projections=projections)

    @cached_property
    def minimum_result(self):
        return self.program.solve(seconds=self.local.fragment_seconds)

    def composed_minimum(self):
        """Compose saved local witnesses, checking shared choices agree."""
        from mappa.mappa_core.core.symmetry.actions import apply_action
        choices = {}
        for response, binding in self.fragment_responses:
            for position, choice in zip(binding, response.event_counts.minimum_choices):
                if position in choices and choices[position] != choice:
                    return None
                choices[position] = choice
        selections = tuple(choices.get(i, dict(discovery=dict(self.local.prepare(action[1]).optimization.discovery_witness))
            if action[0] == 'relation' else None) for i, action in enumerate(self.family.actions))
        images = np.asarray([p for _, p in self.family.mapping])
        for action, choice in zip(self.family.actions, selections):
            images = apply_action(images, action, choice, self.index.m)
        mapping = tuple((a, int(p)) for (a, _), p in zip(self.family.mapping, images))
        count = self.boundary_events.lower + sum(
            response.event_counts.minimum for response, _ in self.fragment_responses)
        if self.sweep_cut is not None:
            a, b = self.sweep_cut
            images = dict(mapping)
            pair = tuple(sorted((images[a], images[b])))
            term = next(iter(self.boundary.values()))[0]
            count += int(term.response.lookup(self.index.p[pair]) != 0)
            count += self.sweep_removed - int(self.boundary_events.source_cut_counted)
            count -= int(pair in self.product_boundary_set)
        return count, mapping, selections

    @cached_property
    def _evaluation(self):
        """Compile array lookups only when a caller evaluates joint choices."""
        groups, constant = defaultdict(list), Counter()
        for terms in (*self.internal.values(), *self.boundary.values()):
            for term in terms:
                a, b = term.atoms
                region = tuple(sorted((int(self.owners[a]), int(self.owners[b]))))
                if term.response.constant is None:
                    groups[term.response].append((a, b, region))
                else:
                    constant[region, term.response.constant] += 1
        regions = tuple(sorted({region for terms in groups.values() for _,_,region in terms}
                               | {region for region, _ in constant}))
        region_ids = {region:i for i,region in enumerate(regions)}
        fixed = np.zeros((len(regions), 5), dtype=np.int64)
        for (region, code), count in constant.items():
            fixed[region_ids[region], code] += count
        variable = tuple((response, np.array([(a,b) for a,b,_ in terms], dtype=int),
                          np.array([region_ids[region] for _,_,region in terms], dtype=int))
                         for response, terms in groups.items())
        return fixed, variable, regions

    def _dependency_components(self):
        """Separate contribution groups without separating shared choices.

        Local relations and the sweep keep their original action bindings.
        Fixed ordinary boundary counts introduce no additional choice dependency.
        """
        parent = {('fragment', i): ('fragment', i) for i in range(len(self.fragments))}
        def root(node):
            parent.setdefault(node, node)
            while parent[node] != node:
                parent[node] = parent[parent[node]]
                node = parent[node]
            return node
        def join(owners, actions):
            nodes = [('fragment', i) for i in owners if i >= 0] + [('action', i) for i in actions]
            if nodes:
                first = root(nodes[0])
                for node in nodes[1:]:
                    parent[root(node)] = first
        # Atom rows join every owner to its actions. Boundary terms only add
        # the connection between their two owners, once per variable region.
        for owners, terms in self.boundary.items():
            if any(term.response.constant is None for term in terms):
                join(owners, ())
        for atom, actions in self.atom_actions.items():
            if actions:
                join((int(self.owners[atom]),), actions)
        groups = defaultdict(lambda: ([], []))
        for kind,i in tuple(parent):
            groups[root((kind,i))][kind == 'action'].append(i)
        return tuple(sorted((tuple(sorted(fs)), tuple(sorted(actions)))
                            for fs,actions in groups.values() if fs))

    def evaluate(self, vector, *, include_events=False):
        """Compose the approximate boundary model for one valid joint choice."""
        vector = np.asarray(vector, dtype=int)
        fixed, groups, regions = self._evaluation
        totals = fixed.copy()
        for response, pairs, region_ids in groups:
            codes = response.lookup(self.index.padded_product[vector[pairs[:,0]], vector[pairs[:,1]]])
            np.add.at(totals, (region_ids, codes), 1)
        internal = {}
        boundary_result = {region: row.copy() for region, row in self.fixed_boundary.items()}
        for (left, right), row in zip(regions, totals[:,1:]):
            if np.any(row):
                collection = internal if left == right else boundary_result
                key = left if left == right else (left, right)
                collection.setdefault(key, np.zeros(4, dtype=int))[:] += row
        sweep_target = None
        if self.sweep_cut is not None:
            a, b = self.sweep_cut
            region = tuple(sorted((int(self.owners[a]), int(self.owners[b]))))
            boundary_result.setdefault(region, np.zeros(4, dtype=int))[0] -= int(
                self.boundary_events.source_cut_counted)
            sweep_target = tuple(sorted((int(vector[a]), int(vector[b]))))
            if sweep_target in self.product_boundary_set:
                region = tuple(sorted(int(self.target_owners[p]) for p in sweep_target))
                boundary_result[region][1] -= 1
        counts = sum(internal.values(), np.zeros(4, dtype=int)) + sum(
            boundary_result.values(), np.zeros(4, dtype=int))
        def named(counts):
            return dict(zip(EVENT_CATEGORIES, map(int, counts)))
        result = dict(counts=named(counts), total=int(counts.sum()),
                      internal={k:named(v) for k,v in internal.items()},
                      boundary={k:named(v) for k,v in boundary_result.items() if np.any(v)},
                      model=self.model, approximate=self.approximate)
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
            events.extend(dict(category='breaking', reactant_pair=(int(a), int(b)),
                product_pair=None, reactant_weight=float(self.index.r[a,b]), product_weight=0.)
                for a,b in self.reactant_boundary_pairs if (a,b) != self.sweep_cut)
            events.extend(dict(category='forming', reactant_pair=None,
                product_pair=(int(a),int(b)), reactant_weight=0., product_weight=float(self.index.p[a,b]))
                for a,b in self.product_boundary_pairs if (a,b) != sweep_target)
            result['events'] = tuple(events)
        return result

    def counts(self, vectors):
        """Vectorized approximate counts for caller-supplied valid choices."""
        vectors = np.asarray(vectors, dtype=int)
        fixed, groups, _ = self._evaluation
        base = fixed[:,1:].sum(axis=0) + sum(self.fixed_boundary.values(), np.zeros(4, dtype=int))
        result = np.broadcast_to(base, (len(vectors), 4)).copy()
        for response, pairs, _ in groups:
            values = self.index.padded_product[vectors[:,pairs[:,0]], vectors[:,pairs[:,1]]]
            codes = response.lookup(values)
            for code in range(1, 5):
                result[:,code-1] += np.count_nonzero(codes == code, axis=1)
        if self.sweep_cut is not None:
            result[:,0] -= int(self.boundary_events.source_cut_counted)
            pairs = np.sort(vectors[:, self.sweep_cut], axis=1)
            keys = pairs[:,0] * self.index.m + pairs[:,1]
            counted = self.product_boundary_pairs
            result[:,1] -= np.isin(keys, counted[:,0] * self.index.m + counted[:,1])
        return result

    @property
    def statistics(self):
        terms = [t for rows in (*self.internal.values(), *self.boundary.values()) for t in rows]
        return dict(fragments=len(self.fragments), internal_terms=sum(map(len,self.internal.values())),
                    neighboring_pairs=self.neighboring_pair_count,
                    boundary_terms=sum(map(len,self.boundary.values())),
                    fixed_boundary_events=self.boundary_count,
                    variable_terms=sum(t.response.constant is None for t in terms),
                    constant_terms=sum(t.response.constant is not None for t in terms),
                    unique_responses=len({id(t.response) for t in terms}),
                    product_boundary_pairs=len(self.product_boundary_pairs),
                    cut_terms=sum(t.manual_cut for t in terms),
                    boundary_cases=dict(Counter(t.endpoint_case for ts in self.boundary.values() for t in ts)),
                    variable_boundary_terms=sum(t.response.constant is None
                                                for ts in self.boundary.values() for t in ts),
                    projected_pair_states=self.projected_pair_states,
                    contribution_components=len(self.components))


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
