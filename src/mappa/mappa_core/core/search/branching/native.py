"""Native branch-search adapter."""

from __future__ import annotations
from dataclasses import asdict
import gc, time
from mappa.mappa_core.core.search.branching.construction import _normalize_anchor_map
from mappa.mappa_core.core.matching.growth import native
from mappa.mappa_core.core.matching.symmetry.orbits import _nauty_orbits
from mappa.mappa_core.core.search.graph import AAMSearchGraph, SearchContext


def find_islands_native(source, target, seed_order, *, graph_floor=.2, iso_tol=1.,
                        max_branches=100, p_orbits=None, r_orbits=None, cuts=(),
                        anchor_map=None, core_R=None, stop_when_core_mapped=False,
                        growth_replay=None, profile=None, random_seed=42):
    order=list(dict.fromkeys(seed_order))
    anchors=_normalize_anchor_map(anchor_map,source,target)
    p_orbits=p_orbits if p_orbits is not None else _nauty_orbits(target,wbo_tol=iso_tol)
    if not native.available():
        raise ValueError('native scheduling requires the built native engine')
    repair=None;repair_cuts=()
    if growth_replay is not None:
        if growth_replay.source is not source or growth_replay.target is not target:
            raise ValueError('fragment repair view belongs to another graph pair')
        repair=growth_replay.engine
        if not isinstance(repair,native._engine.FragmentRepair):
            raise ValueError('native scheduling accepts a FragmentRepair session, not prefix replay')
        repair_cuts=growth_replay.cuts
    r=native.source_graph(source)
    p=growth_replay.target_view if growth_replay is not None else native.target_graph(target,p_orbits)
    if r is None or p is None:
        raise ValueError('native scheduling requires WBO graphs and exact target orbits')
    core=tuple(sorted(set(core_R or ())))
    context=SearchContext(tuple(sorted(source)),tuple(sorted(target)),tuple(order),tuple(cuts),
        core,tuple(sorted(anchors.items())),float(graph_floor),float(iso_tol),int(max_branches),
        objective='requested_core' if stop_when_core_mapped and core else 'full_source',
        random_seed=random_seed)
    cpu=time.process_time() if profile is not None else 0.
    wall=time.perf_counter() if profile is not None else 0.
    enabled = gc.isenabled()
    try:
        # Native records and the resulting decision DAG contain no cycles.
        # Avoid rescanning retained cuts while constructing this one record;
        # reference counting still releases temporary acyclic objects normally.
        gc.disable()
        record=native._engine.search_fragments(r.graph,p.graph,[r.index[a] for a in order],
            float(graph_floor),float(iso_tol),int(max_branches),
            [(r.index[a],p.index[b]) for a,b in sorted(anchors.items())],
            [r.index[a] for a in core],bool(stop_when_core_mapped),r.nodes,p.nodes,repair,repair_cuts,int(random_seed) & ((1<<64)-1))
        native_cpu=time.process_time()-cpu if profile is not None else 0.
        native_wall=time.perf_counter()-wall if profile is not None else 0.
        growth_calls=record.pop('growth_calls')
        max_growth_candidates=record.pop('max_growth_candidates')
        record.update(schema='rxn_core.aam_search_graph/v1',contexts=(asdict(context),))
        cpu=time.process_time() if profile is not None else 0.
        graph=AAMSearchGraph.from_record(record,copy=False)
        del record
    finally:
        if enabled:
            gc.enable()
    if profile is not None:
        profile.append(dict(growth_calls=growth_calls,max_cands_before=max_growth_candidates,native_schedule_and_export_cpu=native_cpu,
            python_graph_cpu=time.process_time()-cpu,extend_elapsed_sec=native_wall))
    return graph
