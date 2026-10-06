"""Share exact representative-constraint checks without unfolding valid DAG paths."""

from collections import Counter
import numpy as np


class RepresentativeValidationWorkspace:
    """Reaction-local reuse; no event thresholds or symmetry equivalence assumptions.

    All histories ending at a state share its representative mapping. Checking
    that mapping against the union of ancestor preservation constraints proves
    every incoming path valid. A failed union is checked pathwise so rejection
    counts and diagnostics agree with the original validator.
    """

    def __init__(self, problem, *, cache_entries=4096):
        self.problem = problem
        self.n = problem.source_atom_count
        self.a, self.b = np.triu_indices(self.n, 1)
        self.weights = problem.reactant.wbo[self.a, self.b]
        self.pair_bits = {
            (int(a), int(b)): 1 << i for i, (a, b) in enumerate(zip(self.a, self.b))
        }
        self.r_elements = np.asarray(problem.reactant.elements)
        self.p_elements = np.asarray(problem.product.elements)
        self.fragment_masks = {}
        self.source_bonds = {}
        self.mask_arrays = {}
        self.cache_entries = cache_entries
        self.stats = Counter()

    def _remember(self, cache, key, value):
        if len(cache) >= self.cache_entries:
            cache.clear()
        cache[key] = value
        return value

    def _edge_mask(self, edge, context):
        if edge.match is None:
            return 0
        fragment = tuple(sorted(edge.match["fragment"]))
        key = context.graph_floor, fragment
        if key in self.fragment_masks:
            self.stats["fragment_mask_hits"] += 1
            mask = self.fragment_masks[key]
        else:
            self.stats["fragment_mask_builds"] += 1
            if context.graph_floor not in self.source_bonds:
                self.source_bonds[context.graph_floor] = tuple(
                    (int(a), int(b), self.pair_bits[int(a), int(b)])
                    for a, b, weight in zip(self.a, self.b, self.weights)
                    if weight >= context.graph_floor
                )
            atoms = set(fragment)
            mask = 0
            for a, b, bit in self.source_bonds[context.graph_floor]:
                if a in atoms and b in atoms:
                    mask |= bit
            self._remember(self.fragment_masks, key, mask)
        excluded = 0
        for a, b in (*context.cuts, *edge.match.get("deferred_edges", ())):
            if a != b:
                excluded |= self.pair_bits[min(a, b), max(a, b)]
        return mask & ~excluded

    def validate_path(self, path):
        """Validate a single complete path using cached preservation masks."""
        mask = 0
        context = path.context
        for edge_id in path.transitions:
            mask |= self._edge_mask(path.graph.transitions[edge_id], context)
        self.stats["path_checks"] += 1
        assert self._valid(path.graph.states[path.terminal], mask, context)

    def _arrays(self, mask):
        if mask in self.mask_arrays:
            self.stats["constraint_array_hits"] += 1
            return self.mask_arrays[mask]
        self.stats["constraint_array_builds"] += 1
        bits, ids = mask, []
        while bits:
            bit = bits & -bits
            ids.append(bit.bit_length() - 1)
            bits -= bit
        ids = np.asarray(ids, dtype=int)
        return self._remember(
            self.mask_arrays, mask, (self.a[ids], self.b[ids], self.weights[ids])
        )

    def _valid(self, state, mask, context):
        pairs = dict(state.mapping)
        if set(pairs) != set(range(self.n)):
            return False
        values = np.asarray([pairs[i] for i in range(self.n)], dtype=int)
        if len(set(values)) != self.n:
            return False
        if not np.all(self.r_elements == self.p_elements[values]):
            return False
        a, b, weights = self._arrays(mask)
        target = self.problem.product.wbo[values[a], values[b]]
        return bool(
            np.all(target >= context.graph_floor)
            and np.all(np.abs(weights - target) <= context.iso_tolerance + 1e-9)
        )

    def _pathwise(self, graph, terminal=None):
        from .family_scoring import validate_representative

        valid, invalid = 0, []
        for path in graph.paths(terminal):
            if len(path.mapping) != self.n:
                continue
            self.stats["fallback_path_checks"] += 1
            try:
                validate_representative(path, self.problem)
                valid += 1
            except AssertionError:
                invalid.append(path.terminal)
        return valid, invalid

    def validate_graph(self, graph):
        """Return (valid full-path count, invalid terminal IDs with multiplicity)."""
        # Native scheduling numbers states topologically. Other graph layouts
        # use the existing oracle; this optimization never repairs their order.
        if len(graph.contexts) != 1 or any(
            e.source >= e.target for e in graph.transitions
        ):
            return self._pathwise(graph)
        roots = set(graph.roots)
        incoming = [[] for _ in graph.states]
        for edge in graph.transitions:
            incoming[edge.target].append(edge)
        masks, paths = [0] * len(graph.states), [0] * len(graph.states)
        context = graph.contexts[0]
        for state in graph.states:
            if state.id in roots:
                paths[state.id] = 1
                continue
            for edge in incoming[state.id]:
                if not paths[edge.source]:
                    continue
                paths[state.id] += paths[edge.source]
                masks[state.id] |= masks[edge.source] | self._edge_mask(edge, context)
        valid, invalid = 0, []
        for terminal in graph.terminals:
            state = graph.states[terminal]
            if len(state.mapping) != self.n or not paths[terminal]:
                continue
            self.stats["terminal_checks"] += 1
            if self._valid(state, masks[terminal], context):
                valid += paths[terminal]
                self.stats["certified_paths"] += paths[terminal]
            else:
                good, bad = self._pathwise(graph, terminal)
                valid += good
                invalid.extend(bad)
        return valid, invalid
