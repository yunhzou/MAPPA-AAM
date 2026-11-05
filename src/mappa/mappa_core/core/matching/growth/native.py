"""Fragment matching and growth with compressed symmetry."""
from __future__ import annotations

import os
from importlib import import_module
import time
from dataclasses import dataclass
import numpy as np
from mappa.mappa_core.core.matching.policy import ElementNodeMatchPolicy
from mappa.mappa_core.core.matching.representation.records import IslandBranchLimitExceeded, FragmentMatchResult, fragment_from_native


try:  # pragma: no cover - depends on the build
    _engine = import_module('mappa._engine')
except ModuleNotFoundError as exc:  # pragma: no cover
    if exc.name != 'mappa._engine':
        raise
    _engine = None


def built():
    return _engine is not None


def available():
    selected = os.environ.get("MAPPA_NATIVE", "1")
    if selected not in ('0', '1'):
        raise ValueError('MAPPA_NATIVE must be 0 (Python) or 1 (native)')
    if selected == '0':
        return False
    if not built():
        raise ImportError('Native growth is selected but mappa._engine is not built; '
                          'build it or explicitly select MAPPA_NATIVE=0')
    return True


@dataclass(frozen=True)
class _NativeGraphView:
    graph: object
    nodes: tuple[int, ...]
    index: dict[int, int]
    bond_tolerance: object = None


def _graph_edges(g, index):
    return sorted(
        (min(index[a], index[b]), max(index[a], index[b]))
        for a, b in g.edges()
    )


def _graph_inputs(graph):
    matrix = graph.graph.get('wbo_matrix')
    if matrix is None:
        return None
    nodes = tuple(sorted(map(int, graph)))
    index = {atom: position for position, atom in enumerate(nodes)}
    rows = np.asarray(matrix, dtype=float)
    if nodes != tuple(range(len(nodes))) or rows.shape != (len(nodes), len(nodes)):
        rows = rows[np.ix_(nodes, nodes)]
    signature = (nodes, tuple(str(graph.nodes[v].get('element')) for v in nodes),
                 tuple(_graph_edges(graph, index)), float(graph.graph.get('bond_cut', .2)),
                 rows.tobytes())
    return signature, index, rows


def source_graph(g_R, bond_tolerance=None):
    """Native view of the (possibly cut) reactant graph, cached on the graph."""
    inputs = _graph_inputs(g_R)
    if inputs is None:
        return None
    signature, index, rows = inputs
    nodes, elements, edges, bond_cut, _ = signature
    tolerances = () if bond_tolerance is None else tuple(bond_tolerance.native_values(nodes))
    inclusive = True if bond_tolerance is None else bond_tolerance.inclusive
    signature = signature, tolerances, inclusive
    cached = g_R.graph.get("_native_source")
    if cached is not None and cached[0] is g_R and cached[1] == signature:
        return cached[2]
    native = _engine.SourceGraph(elements, rows.tolist(), bond_cut, edges, tolerances, inclusive)
    view = _NativeGraphView(native, nodes, index, bond_tolerance)
    g_R.graph["_native_source"] = g_R, signature, view
    return view


def target_graph(g_P):
    """One native target view of the active graph and original WBO matrix."""
    inputs = _graph_inputs(g_P)
    if inputs is None:
        return None
    signature, index, rows = inputs
    entry = g_P.graph.get('_native_target')
    if entry is not None and entry[0] is g_P and entry[1] == signature:
        return entry[2]
    nodes, elements, edges, bond_cut, _ = signature
    native = _engine.TargetGraph(elements, rows.tolist(), bond_cut, edges)
    view = _NativeGraphView(native, nodes, index)
    g_P.graph['_native_target'] = g_P, signature, view
    return view


def validate_inputs(node_policy, events):
    """Reject unsupported native requests; never select another implementation."""
    if events is not None:
        raise ValueError('Native growth does not support trace events; explicitly select MAPPA_NATIVE=0')
    if type(node_policy) is not ElementNodeMatchPolicy:
        raise ValueError('Native growth requires the element node policy; explicitly select MAPPA_NATIVE=0 for custom policies')


def grow_island(g_R, g_P, seed, mapping, *, graph_floor, iso_tol, min_lock_size,
                max_branches, islands_R, prior_deferred_edges,
                allow_mapped_seed, profile, profile_context, fragment_cache=None, random_seed=42, bond_tolerance=None,
                placement_fixed_targets=()):
    """Run the selected native engine or report unsupported inputs explicitly."""
    if fragment_cache is not None and not isinstance(getattr(fragment_cache, 'engine', None), _engine.FragmentRepair):
        raise TypeError('fragment cache must be a completed-fragment cache view')
    if fragment_cache is not None and (fragment_cache.source is not g_R or fragment_cache.target is not g_P):
        raise ValueError("fragment cache view belongs to different source/target graphs")
    source_view = fragment_cache.source_view if fragment_cache is not None else source_graph(g_R, bond_tolerance)
    target_view = fragment_cache.target_view if fragment_cache is not None else target_graph(g_P)
    if fragment_cache is not None:
        recorded = source_view.bond_tolerance
        if (recorded is None) != (bond_tolerance is None) or (
                recorded is not None and (recorded.inclusive != bond_tolerance.inclusive or
                any(recorded.limit(a, b) != value for (a, b), value in bond_tolerance.values.items()))):
            raise ValueError("fragment cache bond tolerances differ from its source policy")
    if source_view is None or target_view is None:
        raise ValueError('Native growth requires original WBO matrices')
    source_index = source_view.index
    target_index = target_view.index
    if any(p not in target_index for p in placement_fixed_targets):
        raise ValueError('Placement fixed targets contain atoms outside the target graph')
    fixed_targets = sorted({target_index[int(p)] for p in placement_fixed_targets})
    if seed not in source_index:
        raise ValueError('Growth seed is not in the source graph')
    n_r = len(source_view.nodes)
    image = [-1] * n_r
    for r, p in mapping.items():
        if r not in source_index or p not in target_index:
            raise ValueError('Locked mapping contains atoms outside the endpoint graphs')
        image[source_index[int(r)]] = target_index[int(p)]
    islands = None
    if islands_R is not None:
        if any(r not in source_index for r in islands_R):
            raise ValueError('Locked islands contain atoms outside the source graph')
        islands = [(source_index[int(r)], int(k))
                   for r, k in islands_R.items()]
    deferred = []
    for a, b in prior_deferred_edges or ():
        if a not in source_index or b not in source_index:
            raise ValueError('Deferred edges contain atoms outside the source graph')
        deferred.append((source_index[int(a)], source_index[int(b)]))
    started = time.perf_counter()
    arguments = (source_index[int(seed)], image, float(graph_floor), float(iso_tol),
                 int(min_lock_size), int(max_branches), islands, deferred, bool(allow_mapped_seed), int(random_seed) & ((1<<64)-1))
    if fragment_cache is None:
        out = _engine.grow_island(source_view.graph, target_view.graph, *arguments,
                                  placement_fixed_targets=fixed_targets)
    else:
        out = fragment_cache.engine.grow(fragment_cache.cuts, *arguments,
                                        placement_fixed_targets=fixed_targets)
    elapsed = time.perf_counter() - started
    if profile is not None:
        prof = {'seed': int(seed)}
        prof.update(out["profile"])
        prof['extend_elapsed_sec'] = elapsed
        prof['max_extend_elapsed_sec'] = 0.0
        prof['slowest_extend'] = None
        prof['elapsed_sec'] = elapsed
        if profile_context:
            prof.update(profile_context)
        profile.append(prof)
    if out["capped"]:
        raise IslandBranchLimitExceeded(
            out["cap_count"], out["cap_limit"], seed=int(seed))
    self_pool = {}
    results = tuple(fragment_from_native(record, source_view.nodes, target_view.nodes, self_pool)
                    for record in out['fragments'])
    deferred = frozenset(tuple(sorted((source_view.nodes[int(a)], source_view.nodes[int(b)])))
                         for a, b in out["deferred_edges"])
    generators = out['placement_generators']
    if target_view.nodes != tuple(range(len(target_view.nodes))):
        translated = []
        for generator in generators:
            values = list(range(max(target_view.nodes) + 1))
            for atom, image in zip(target_view.nodes, generator):
                values[atom] = target_view.nodes[image]
            translated.append(tuple(values))
        generators = tuple(translated)
    return FragmentMatchResult(tuple(results), False, len(results), max_branches,
                               elapsed, deferred, generators)
