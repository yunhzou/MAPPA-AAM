"""Opt-in native conditioning over one immutable target symmetry graph.

The Python fragment and AAM graph abstractions are unchanged. The native object
reuses role colors and complete ordered partitions, returning the same ordered
generators as the reference finalizer, without permutation enumeration.
"""
from ._engine import ConditionedAutGraph
from .matcher.canonical import _CandidateAutomorphismCanonicalizer
from .search_symmetry import SymmetryWorkspace


class ConditionedSymmetryWorkspace(SymmetryWorkspace):
    def __init__(self, target, iso_tolerance, *, cache_bytes=64 * 1024 * 1024):
        super().__init__(target, iso_tolerance)
        self.canonicalizer = _CandidateAutomorphismCanonicalizer(target, wbo_tol=iso_tolerance)
        self.cache_bytes = cache_bytes
        self.engines = {}
        self.identity_shortcuts = 0

    def _engine(self, continuation):
        if continuation in self.engines:
            return self.engines[continuation]
        c = self.canonicalizer
        edge_colors = c.exact_edge_color_classes if continuation else c.edge_color_classes
        engine = ConditionedAutGraph(
            c.n_vertices,
            [(a, b) for a, neighbors in c.adjacency.items() for b in neighbors if a <= b],
            list(c.nodes), [c.atom_index[a] for a in c.nodes],
            [repr(c.atom_base_color[c.atom_index[a]]) for a in c.nodes],
            [(repr(('edge', color)), sorted(vertices)) for color, vertices in edge_colors],
            self.cache_bytes // 2,
        )
        self.engines[continuation] = engine
        return engine

    def conditioned_generators(self, candidate, locked, intern=tuple):
        continuation = candidate.compression_version >= 2
        # Version 2 fixes every outside/pool atom pointwise. Without an
        # automorphism block, each remaining atom has its unique source role.
        # Thus the conditioned group is identity; local pool permutations
        # remain represented by candidate.blocks, independently of this group.
        if continuation and not candidate.automorph_blocks:
            self.identity_shortcuts += 1
            return intern(())
        c = self.canonicalizer
        role_builder = c._growth_roles if continuation else c._candidate_roles
        roles = role_builder(candidate, group_domains=True)
        return intern(self._engine(continuation).generators(
            locked, {p: repr(role) for p, role in roles.items()}))

    def coloring_count(self):
        return sum(engine.stats()['group_solves'] for engine in self.engines.values())

    def stats(self):
        totals = dict.fromkeys(('calls', 'partition_hits', 'group_solves', 'color_hits',
                              'cache_clears', 'resident_bytes', 'cache_budget_bytes'), 0)
        for engine in self.engines.values():
            for key, value in engine.stats().items():
                totals[key] += value
        totals['identity_shortcuts'] = self.identity_shortcuts
        return totals
