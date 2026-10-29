"""Prepared search state and per-cut scientific operations."""

from __future__ import annotations
from mappa.mappa_core.core.search.planning import cut_seed, _generate_seed_orders
import gc, time
from mappa.mappa_core.core.search.branch import find_islands
from mappa.frag import build_graph
from mappa.mappa_core.core.matching.symmetry.orbits import _nauty_orbits
from mappa.mappa_core.core.search.graph import AAMSearchGraph


class _GrowthCounts:
    """Reduce worker profile events online instead of retaining every row."""
    def __init__(self):
        self.maximum = 0

    def append(self, row):
        self.maximum = max(self.maximum, row.get('max_cands_before', 0))


class SearchSession:
    """Prepared matching state reused by local or remote workers.

    Seed batches select from the original full seed plan. They never reseed a
    separate search. Sessions own their caches and are not shared across threads.
    """
    def __init__(self, problem, config, execution='reused_native'):
        if execution not in ('reference', 'reused_native', 'shared_policies'):
            raise ValueError('unknown AAM execution backend')
        if execution != 'reference':
            from mappa.mappa_core.core.matching.growth.native import available
            if not available():
                raise ValueError('native execution requires the built native engine')
        self.execution = execution
        target = build_graph(problem.product.elements, problem.product.wbo,
                             bond_cut=config.graph_floor)
        self.context = (problem, config, target,
                        _nauty_orbits(target, wbo_tol=config.iso_tolerance))
        self.repair = None
        if execution in ('reused_native', 'shared_policies'):
            from mappa.mappa_core.core.search.reuse import FragmentRepair
            source = build_graph(problem.reactant.elements, problem.reactant.wbo,
                                 bond_cut=config.graph_floor)
            self.repair = FragmentRepair(source, target, self.context[3])

    def search_cut(self, cut, *, seed_indices=None):
        started = time.perf_counter()
        problem, config, target, target_orbits = self.context
        view = (self.repair.for_cut(tuple(edge for edge in cut if self.repair.source.has_edge(*edge)))
                if self.repair is not None else None)
        if view is None:
            source = build_graph(problem.reactant.elements, problem.reactant.wbo,
                                 bond_cut=config.graph_floor)
            source.remove_edges_from(cut)
            matcher = find_islands
        else:
            from mappa.mappa_core.core.search.native import find_islands_native
            source = view.source
            matcher = find_islands_native
        if self.execution == 'shared_policies':
            if seed_indices is not None:
                raise ValueError('Shared-policy search is an indivisible cut task')
            from mappa.mappa_core.core.search.shared import AdaptiveSearchCondition
            from mappa.mappa_core.core.search.shared import SharedSeedSearch
            condition = AdaptiveSearchCondition(source, target, (), target_orbits, None, cut, view)
            gc_enabled = gc.isenabled()
            try:
                # These DAG/cache records are acyclic; avoid repeatedly scanning
                # the growing cut heap. Restore the caller's setting at the boundary.
                gc.disable()
                profile = _GrowthCounts()
                session = SharedSeedSearch(problem, config, condition=condition, profile=profile)
                while session.advance():
                    pass
                graph = session.builder.finish()
                return graph, {'search_seconds':time.perf_counter()-started,
                    'max_live_branches':session.max_live_branches,
                    'max_growth_candidates':profile.maximum}
            finally:
                if gc_enabled:
                    gc.enable()
        source_orbits = _nauty_orbits(source, wbo_tol=config.iso_tolerance)
        graphs, profile = [], _GrowthCounts()
        orders = tuple(_generate_seed_orders(source, n_trials=config.seed_count,
            rng_seed=cut_seed(cut, config.random_seed), seed_selection=config.seed_selection))
        indices = tuple(range(len(orders))) if seed_indices is None else tuple(seed_indices)
        if not indices or tuple(sorted(set(indices))) != indices or any(i < 0 or i >= len(orders) for i in indices):
            raise ValueError('Seed indices must be a nonempty ordered subset of the original seed plan')
        for index in indices:
            order = orders[index]
            graphs.append(matcher(source, target, order,
                graph_floor=config.graph_floor, iso_tol=config.iso_tolerance,
                max_branches=config.branch_limit, p_orbits=target_orbits,
                r_orbits=source_orbits, anchor_map=dict(config.anchors),
                profile=profile, cuts=cut, growth_replay=view, random_seed=config.random_seed))
        graph = AAMSearchGraph.combine(graphs)
        return graph, {
            'search_seconds': time.perf_counter() - started,
            'max_live_branches': max((len(g.terminals) for g in graphs), default=0),
            'max_growth_candidates': profile.maximum,
        }
