"""Prepared search state and per-cut scientific operations."""

from __future__ import annotations
from mappa.mappa_core.core.search.planning import cut_seed, growth_seed, _generate_seed_orders
import gc, time
from mappa.mappa_core.core.search.branching.construction import find_islands
from mappa.mappa_core.core.matching.graphs import build_graph
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
        self.context = (problem, config, target)
        self.fragments = None
        if execution in ('reused_native', 'shared_policies'):
            from mappa.mappa_core.core.matching.reuse import FragmentPool
            source = build_graph(problem.reactant.elements, problem.reactant.wbo,
                                 bond_cut=config.graph_floor)
            self.fragments = FragmentPool(source, target)

    def search_cut(self, cut, *, seed_indices=None):
        started = time.perf_counter()
        problem, config, target = self.context
        view = (self.fragments.for_cut(tuple(edge for edge in cut if self.fragments.source.has_edge(*edge)))
                if self.fragments is not None else None)
        if view is None:
            source = build_graph(problem.reactant.elements, problem.reactant.wbo,
                                 bond_cut=config.graph_floor)
            source.remove_edges_from(cut)
            matcher = find_islands
        else:
            from mappa.mappa_core.core.search.branching.native import find_islands_native
            source = view.source
            matcher = find_islands_native
        if self.execution == 'shared_policies':
            if seed_indices is not None:
                raise ValueError('Shared-policy search is an indivisible cut task')
            from mappa.mappa_core.core.search.branching.shared import AdaptiveSearchCondition
            from mappa.mappa_core.core.search.branching.shared import SharedSeedSearch
            condition = AdaptiveSearchCondition(source, target, (), cut, view)
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
        graphs, seed_metrics, profile = [], [], _GrowthCounts()
        orders = tuple(_generate_seed_orders(source, n_trials=config.seed_count,
            rng_seed=cut_seed(cut, config.random_seed), seed_selection=config.seed_selection))
        indices = tuple(range(len(orders))) if seed_indices is None else tuple(seed_indices)
        if not indices or tuple(sorted(set(indices))) != indices or any(i < 0 or i >= len(orders) for i in indices):
            raise ValueError('Seed indices must be a nonempty ordered subset of the original seed plan')
        cache_before = self.fragments.stats() if self.fragments is not None else None
        for index in indices:
            order = orders[index]
            seed_started, seed_cpu = time.perf_counter(), time.process_time()
            graphs.append(matcher(source, target, order,
                graph_floor=config.graph_floor, iso_tol=config.iso_tolerance,
                max_branches=config.branch_limit,
                anchor_map=dict(config.anchors),
                profile=profile, cuts=cut, fragment_cache=view, random_seed=growth_seed(cut, index, config.random_seed)))
            measured = dict(index=index, wall_seconds=time.perf_counter()-seed_started,
                            cpu_seconds=time.process_time()-seed_cpu)
            if self.fragments is not None:
                cache_after = self.fragments.stats()
                measured['cache'] = {key: cache_after[key]-cache_before[key]
                                     for key in ('existence_queries', 'fragment_hits', 'set_blocks_checked')}
                measured['cache'].update(retained_bytes=cache_after['resident_bytes'],
                                         entries=cache_after['entries'])
                cache_before = cache_after
            seed_metrics.append(measured)
        graph = AAMSearchGraph.combine(graphs)
        return graph, {
            'search_seconds': time.perf_counter() - started,
            'max_live_branches': max((len(g.terminals) for g in graphs), default=0),
            'max_growth_candidates': profile.maximum,
            'seeds': tuple(seed_metrics),
        }
