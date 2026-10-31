"""Complete a proven fragment using the established compressed representation."""
from dataclasses import dataclass
from time import perf_counter

from mappa.mappa_core.core.matching.growth.existence import FragmentQuery
from mappa.mappa_core.core.matching.growth.extend import _extend_sym_cands
from mappa.mappa_core.core.matching.growth.frontier import _frontier_boundary_edges


@dataclass(frozen=True)
class AcceptedGrowthStep:
    anchor: int
    atom: int
    weight: float
    deferred: tuple


def complete_fragment(candidates, seed, steps, fragment, bonds, source, target,
                      locked, islands, source_orbits, target_orbits, graph_floor,
                      tolerance, node_policy, bond_tolerance, profile=None):
    """Preserve the old result order, witnesses, pools and symmetry metadata.

    Future constraints may discard an impossible candidate, but never repack
    surviving candidates. The full query is only needed for ambiguous prefixes.
    """
    started = perf_counter() if profile is not None else None
    query = None

    def retain_completable(candidates, prefix):
        nonlocal query
        if len(candidates) <= 1 or prefix == fragment:
            return candidates
        if query is None:
            query = FragmentQuery(source, target, fragment, bonds, locked,
                                  node_policy, tolerance, bond_tolerance)
        return [c for c in candidates if query.find(candidate=c) is not None]

    prefix = {seed}
    candidates = retain_completable(candidates, prefix)
    for step in steps:
        enlarged = prefix | {step.atom}
        if step.atom in locked and islands is not None and step.atom in islands:
            enlarged.update(r for r, label in islands.items() if label == islands[step.atom])
        candidates = _extend_sym_cands(candidates, prefix, step.atom, source, target, locked,
            tolerance, islands, p_orbits=target_orbits, r_orbits=source_orbits,
            deferred_edges=step.deferred, anchor_u=step.anchor, anchor_wbo=step.weight,
            dedupe_edges=set(step.deferred) | _frontier_boundary_edges(source, enlarged, graph_floor),
            node_policy=node_policy, bond_tolerance=bond_tolerance)
        prefix = enlarged
        candidates = retain_completable(candidates, prefix)
    if not candidates:
        raise RuntimeError('completed fragment lost its proven matching witness')
    if profile is not None:
        profile['completion_seconds'] = perf_counter() - started
        profile['completion_steps'] = len(steps)
        profile['completion_queries'] = 0 if query is None else query.calls
        profile['completion_pruned'] = 0 if query is None else query.rejected
    return candidates
