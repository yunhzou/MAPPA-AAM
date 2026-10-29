"""Reusable fragment growth interface."""
from __future__ import annotations

from mappa.mappa_core.core.matching.representation.records import FragmentPlacement, FragmentMatchContext, FragmentMatchConfig, FragmentMatchResult
import time
from mappa.mappa_core.core.matching.representation.records import IslandBranchLimitExceeded
from mappa.mappa_core.core.matching.growth.engine import grow_island
from mappa.mappa_core.core.matching.symmetry.orbits import _nauty_orbits


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
    source_orbits, target_orbits = context.source_orbits, context.target_orbits
    if config.orbit_dedup:
        if source_orbits is None:
            source_orbits = _nauty_orbits(source, wbo_tol=config.iso_tolerance,
                                          node_policy=config.node_policy)
        if target_orbits is None:
            target_orbits = _nauty_orbits(target, wbo_tol=config.iso_tolerance,
                                          node_policy=config.node_policy)
    try:
        matches = grow_island(
            source, target, seed, context.locked_mapping,
            graph_floor=config.graph_floor, iso_tol=config.iso_tolerance,
            min_lock_size=config.minimum_size, max_branches=config.branch_limit,
            node_policy=config.node_policy, islands_R=context.islands,
            p_orbits=target_orbits, r_orbits=source_orbits,
            prior_deferred_edges=context.deferred_edges,
            allow_mapped_seed=config.allow_mapped_seed,
            replay=context.growth_replay, random_seed=config.random_seed,
            events=events, profile=profile, profile_context=profile_context,
            bond_tolerance=config.bond_tolerance)
    except IslandBranchLimitExceeded as exc:
        return FragmentMatchResult((), True, exc.count, exc.limit,
                                   time.perf_counter() - started)
    placements = tuple(FragmentPlacement.from_match(match, source) for match in matches)
    return FragmentMatchResult(placements, False, len(matches),
                               config.branch_limit, time.perf_counter() - started)
