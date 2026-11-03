"""Reusable fragment growth interface."""
from __future__ import annotations

from dataclasses import replace
from mappa.mappa_core.core.matching.representation.records import FragmentMatchContext, FragmentMatchConfig, FragmentMatchResult
import time
from mappa.mappa_core.core.matching.representation.records import IslandBranchLimitExceeded
from mappa.mappa_core.core.matching.growth.engine import grow_island


def match_fragment(source, target, *, seed, context=None, config=None,
                   events=None, profile=None, profile_context=None):
    """Grow one fragment, preserving compressed placements and cap evidence.

    The caller chooses seeds and commits results. This function neither groups
    mechanisms nor requires balanced endpoint composition. Inputs are prepared
    weighted NetworkX graphs, as used by the underlying growth kernel.
    """
    context = context or FragmentMatchContext()
    config = config or FragmentMatchConfig()
    started = time.perf_counter()
    from mappa.mappa_core.core.matching.graphs import _coerce_graph
    source = _coerce_graph(source, config.graph_floor)
    target = _coerce_graph(target, config.graph_floor)
    if not set(context.placement_fixed_targets).issubset(target):
        raise ValueError('fixed placement target is outside the target graph')
    try:
        result = grow_island(
            source, target, seed, context.locked_mapping,
            graph_floor=config.graph_floor, iso_tol=config.iso_tolerance,
            min_lock_size=config.minimum_size, max_branches=config.branch_limit,
            node_policy=config.node_policy, islands_R=context.islands,
            prior_deferred_edges=context.deferred_edges,
            allow_mapped_seed=config.allow_mapped_seed,
            fragment_cache=context.fragment_cache, random_seed=config.random_seed,
            events=events, profile=profile, profile_context=profile_context,
            bond_tolerance=config.bond_tolerance, self_tolerance=config.self_tolerance,
            placement_fixed_targets=context.placement_fixed_targets)
    except IslandBranchLimitExceeded as exc:
        return FragmentMatchResult((), True, exc.count, exc.limit,
                                   time.perf_counter() - started)
    return replace(result, elapsed_seconds=time.perf_counter() - started)
