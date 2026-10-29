"""Exact conditioned symmetry workspaces for fragment matching."""
from __future__ import annotations

from collections import defaultdict
from mappa.mappa_core.core.matching.symmetry.canonical import _CandidateAutomorphismCanonicalizer


class SymmetryWorkspace:
    """Exact finalization reuse for one fixed target/tolerance, not a search cap.

    Owned by an explicit caller (typically one reaction), never process-global.
    Keys include locked atom roles and the complete candidate symmetry record.
    """
    def __init__(self, target, iso_tolerance):
        self.target = target
        self.iso_tolerance = iso_tolerance
        self.cache = {}
        self.coloring_cache = {}
        self.generators = {}
        self.groups = {}
        self.canonicalizers = {}

    def canonicalizer_for(self, version):
        if version not in self.canonicalizers:
            self.canonicalizers[version] = _CandidateAutomorphismCanonicalizer(
                self.target, wbo_tol=self.iso_tolerance, compression_version=version)
        return self.canonicalizers[version]

    def conditioned_generators(self, candidate, locked, intern=tuple):
        canonicalizer = self.canonicalizer_for(candidate.compression_version)
        locked_roles = defaultdict(list)
        for r, p in sorted(locked):
            if p in canonicalizer.atom_index:
                locked_roles[p].append(('locked', r))
        coloring = canonicalizer._colored_vertices(candidate, group_domains=True,
            locked_roles={p: tuple(roles) for p, roles in locked_roles.items()})
        key = candidate.compression_version, tuple(vertices for _label, vertices in coloring)
        if key not in self.coloring_cache:
            self.coloring_cache[key] = intern(canonicalizer.atom_generators(candidate, colored_vertices=coloring))
        return self.coloring_cache[key]

    def coloring_count(self):
        return len(self.coloring_cache)


class ConditionedSymmetryWorkspace(SymmetryWorkspace):
    def __init__(self, target, iso_tolerance, *, cache_bytes=64 * 1024 * 1024):
        super().__init__(target, iso_tolerance)
        self.cache_bytes = cache_bytes
        self.engines = {}
        self.identity_shortcuts = 0

    def _engine(self, version):
        from mappa._engine import ConditionedAutGraph
        if version in self.engines:
            return self.engines[version]
        c = self.canonicalizer_for(version)
        edge_colors = c.exact_edge_color_classes if version >= 2 else c.edge_color_classes
        engine = ConditionedAutGraph(
            c.n_vertices,
            [(a, b) for a, neighbors in c.adjacency.items() for b in neighbors if a <= b],
            list(c.nodes), [c.atom_index[a] for a in c.nodes],
            [repr(c.atom_base_color[c.atom_index[a]]) for a in c.nodes],
            [(repr(('edge', color)), sorted(vertices)) for color, vertices in edge_colors],
            self.cache_bytes // 3,
        )
        self.engines[version] = engine
        return engine

    def conditioned_generators(self, candidate, locked, intern=tuple):
        version = candidate.compression_version
        # Version 2 fixes every outside/pool atom pointwise. Without an
        # automorphism block, each remaining atom has its unique source role.
        # Thus the conditioned group is identity; local pool permutations
        # remain represented by candidate.blocks, independently of this group.
        if version == 2 and not candidate.automorph_blocks:
            self.identity_shortcuts += 1
            return intern(())
        c = self.canonicalizer_for(version)
        role_builder = c._growth_roles if version >= 2 else c._candidate_roles
        roles = role_builder(candidate, group_domains=True)
        return intern(self._engine(version).generators(
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
