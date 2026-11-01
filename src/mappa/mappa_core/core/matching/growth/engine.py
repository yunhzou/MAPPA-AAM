"""Fragment matching and growth with compressed symmetry."""
from __future__ import annotations

import time
from mappa.mappa_core.core.matching.representation.state import _cand_map, _cand_possible_p_atoms, _symmetry_state
from mappa.mappa_core.core.matching.policy import as_node_match_policy
from mappa.mappa_core.core.matching.growth import native
from mappa.mappa_core.core.matching.growth.frontier import _push_edges_from, _set_unique, _pop_random_atom, GrowthRandom
from mappa.mappa_core.core.matching.representation.records import IslandBranchLimitExceeded, _IsoResult
from mappa.mappa_core.core.matching.growth.trace import cand_possible_values, cands_pattern_sample, cands_sample, heap_snapshot, pool_by_frag_atom, represented_assignment_expr
from mappa.mappa_core.core.matching.growth.existence import extend_witness
from mappa.mappa_core.core.matching.growth.completion import complete_fragment


def grow_island(g_R, g_P, seed, mapping,
                graph_floor=0.2,
                iso_tol=0.5,
                min_lock_size=1,
                max_branches=1_000_000,
                events=None,
                islands_R=None,
                p_orbits=None,
                prior_deferred_edges=None,
                node_policy=None,
                allow_mapped_seed=False,
                profile=None,
                profile_context=None, fragment_cache=None, random_seed=42, bond_tolerance=None):
    """
    Grow a fragment from `seed` by uniformly sampling the live boundary atoms.

    Bond order is used in compatibility checks, never in selection priority.
    ``random_seed`` makes each conditional fragment call reproducible.

    ``max_branches`` bounds distinct completed fragment placements.
    Growth carries one witness; completion constructs the compressed choices.

    Returns a list of isos:
      []            -- failed (no initial cands, or fragment too small)
      [single_iso]  -- locked successfully (set-unique or single cand)
      [iso_a, ...]  -- non-set-unique saturation; caller branches

    Optional `events` list receives diagnostic events (seed_start /
    commit / deferred / merge / seed_end) compatible with the
    existing trace_run.HTML viewer.
    """
    if native.available():
        native.validate_inputs(p_orbits, as_node_match_policy(node_policy), events)
        bond_tolerance = None if bond_tolerance is None else bond_tolerance.resolve(g_R)
        return native.grow_island(
            g_R, g_P, seed, mapping, graph_floor=graph_floor, iso_tol=iso_tol,
            min_lock_size=min_lock_size, max_branches=max_branches,
            islands_R=islands_R, p_orbits=p_orbits,
            prior_deferred_edges=prior_deferred_edges,
            allow_mapped_seed=allow_mapped_seed, profile=profile,
            profile_context=profile_context, fragment_cache=fragment_cache, random_seed=random_seed,
            bond_tolerance=bond_tolerance)
    if fragment_cache is not None:
        raise ValueError('Fragment cache requires explicit native execution')
    return grow_island_python(g_R, g_P, seed, mapping, graph_floor, iso_tol,
        min_lock_size, max_branches, events, islands_R, p_orbits,
        prior_deferred_edges, node_policy, allow_mapped_seed, profile,
        profile_context, random_seed, bond_tolerance)


def grow_island_python(g_R, g_P, seed, mapping, graph_floor=.2, iso_tol=.5,
                      min_lock_size=1, max_branches=1_000_000, events=None,
                      islands_R=None, p_orbits=None,
                      prior_deferred_edges=None, node_policy=None,
                      allow_mapped_seed=False, profile=None, profile_context=None,
                      random_seed=42, bond_tolerance=None):
    """Explicit Python growth, also used by diagnostic trajectory replay."""
    node_policy = as_node_match_policy(node_policy)
    bond_tolerance = None if bond_tolerance is None else bond_tolerance.resolve(g_R)
    record = events is not None
    prof = None
    profile_t0 = None
    if profile is not None:
        profile_t0 = time.perf_counter()
        prof = {
            'seed': int(seed),
            'seed_targets': 0,
            'init_cands': 0,
            'initial_heap': 0,
            'heap_pops': 0,
            'stale_pops': 0,
            'fragment_skip_pops': 0,
            'extend_calls': 0,
            'extend_elapsed_sec': 0.0,
            'max_extend_elapsed_sec': 0.0,
            'max_cands_before': 0,
            'max_cands_after': 0,
            'max_fragment_size': 0,
            'max_heap_len': 0,
            'commits': 0,
            'deferred': 0,
            'merge_calls': 0,
            'free_extend_calls': 0,
            'slowest_extend': None,
            'feasibility_calls': 0,
            'feasibility_nodes': 0,
            'feasibility_searches': 0,
            'completion_seconds': 0.,
            'completion_queries': 0,
            'completion_pruned': 0,
        }
        if profile_context:
            prof.update(profile_context)

    def _finish_profile(result, cands_count=0, fragment_set=None, branches=0):
        if prof is None:
            return
        prof['result'] = result
        prof['final_cands'] = int(cands_count)
        prof['final_fragment_size'] = int(len(fragment_set or ()))
        prof['branches'] = int(branches)
        prof['elapsed_sec'] = time.perf_counter() - profile_t0
        profile.append(prof)

    locked_p_atoms = set(mapping.values())
    if seed in mapping:
        if not allow_mapped_seed:
            _finish_profile('already_mapped')
            return []
        seed_targets = [mapping[seed]]
    else:
        seed_targets = [v for v in g_P.nodes()
                        if v not in locked_p_atoms
                        and node_policy.compatible(g_R, seed, g_P, v)]
    if prof is not None:
        prof['seed_targets'] = int(len(seed_targets))
    cands = [{seed: seed_targets[0]}] if seed_targets else []
    if not cands:
        _finish_profile('no_initial_cands')
        if record:
            events.append({'type': 'seed_start', 'seed': int(seed),
                           'init_cands': 0, 'fragment': [int(seed)],
                           'p_atoms': []})
            events.append({'type': 'seed_end', 'result': 'no_initial_cands',
                           'final_cands': 0, 'fragment': [int(seed)],
                           'iso': None})
        return []
    fragment = {seed}
    preserved_bonds = set()
    distance = {seed: 0}
    used_edges = set()
    deferred_edges = {tuple(sorted(e)) for e in (prior_deferred_edges or ())}
    heap = []
    _push_edges_from(heap, used_edges, g_R, seed, fragment, graph_floor)
    if prof is not None:
        prof['init_cands'] = int(len(cands))
        prof['initial_heap'] = int(len(heap))
        prof['max_cands_before'] = int(len(cands))
        prof['max_cands_after'] = int(len(cands))
        prof['max_fragment_size'] = int(len(fragment))
        prof['max_heap_len'] = int(len(heap))

    if record:
        events.append({
            'type': 'seed_start',
            'candidate_kind': 'witness',
            'seed': int(seed),
            'init_cands': len(cands),
            'represented_assignments': represented_assignment_expr(cands),
            'fragment': [int(seed)],
            'p_atoms': sorted({int(v) for c in cands
                               for v in _cand_possible_p_atoms(c)}),
            'cand_patterns': cands_pattern_sample(cands, 5),
        })

    rng = GrowthRandom(random_seed)
    while heap:
        selected = _pop_random_atom(heap, used_edges, fragment, rng)
        if selected is None:
            break
        neg_w, u, n = selected
        if prof is not None:
            prof['heap_pops'] += 1
        wbo = -neg_w
        edge = frozenset({u, n})
        if edge in used_edges:
            if prof is not None:
                prof['stale_pops'] += 1
            continue
        used_edges.add(edge)
        if n in fragment:
            if prof is not None:
                prof['fragment_skip_pops'] += 1
            if record:
                events.append({
                    'type': 'pop_skip',
                    'edge': {'frag_atom': int(u), 'ext_atom': int(n), 'wbo': round(wbo, 3)},
                    'reason': 'ext_atom already in fragment',
                })
            continue

        n_in_mapping = n in mapping
        if prof is not None:
            prof['extend_calls'] += 1
            if n_in_mapping:
                prof['merge_calls'] += 1
            else:
                prof['free_extend_calls'] += 1
        if record:
            events.append({
                'type': 'pop',
                'edge': {'frag_atom': int(u), 'ext_atom': int(n),
                         'wbo': round(wbo, 3),
                         'ext_element': g_R.nodes[n].get('element')},
                'scenario': ('merge_island' if n_in_mapping else 'extend_free'),
                'island_id_at_ext': (int(islands_R[n]) if islands_R and n in islands_R else None),
                'island_image': (int(mapping[n]) if n_in_mapping else None),
                'pre_state': {
                    'fragment_size': len(fragment),
                    'fragment': sorted(int(x) for x in fragment),
                    'cands_count': len(cands),
                    'represented_assignments': represented_assignment_expr(cands),
                    'cands_sample': cands_sample(cands, 5),
                    'cands_pattern_sample': cands_pattern_sample(cands, 5),
                    'p_atoms_in_cands': sorted({int(v) for c in cands
                                                for v in _cand_possible_p_atoms(c)}),
                },
                'heap_top_after_pop': heap_snapshot(heap, used_edges, fragment, mapping, 8),
                'pool_by_frag_atom': pool_by_frag_atom(heap, used_edges, fragment, mapping, g_R),
            })

        # Tentatively add n to the shared fragment. If any weighted extension
        # survives, commit it; otherwise record this edge as deferred boundary.
        old_count = len(cands)
        if prof is not None:
            prof['max_cands_before'] = max(
                int(prof['max_cands_before']), int(old_count))
        old_fragment = set(fragment)
        candidate_fragment = fragment | {n}
        enlarged_fragment = set(candidate_fragment)
        if n_in_mapping and islands_R is not None and n in islands_R:
            target_iid = islands_R[n]
            enlarged_fragment |= {
                r for r, k in islands_R.items() if k == target_iid
            }
        extend_t0 = time.perf_counter() if prof is not None else None
        witness, nodes, searched = extend_witness(g_R, g_P, cands[0],
            u, n, enlarged_fragment, preserved_bonds, mapping, node_policy, iso_tol, bond_tolerance)
        new_cands = [] if witness is None else [witness]
        if prof is not None:
            prof['feasibility_calls'] += 1
            prof['feasibility_nodes'] += nodes
            prof['feasibility_searches'] += searched
        if prof is not None:
            extend_elapsed = time.perf_counter() - extend_t0
            prof['extend_elapsed_sec'] += extend_elapsed
            if extend_elapsed > prof['max_extend_elapsed_sec']:
                prof['max_extend_elapsed_sec'] = extend_elapsed
                prof['slowest_extend'] = {
                    'frag_atom': int(u),
                    'ext_atom': int(n),
                    'wbo': float(wbo),
                    'scenario': (
                        'merge_island' if n_in_mapping else 'extend_free'
                    ),
                    'elapsed_sec': extend_elapsed,
                    'cands_before': int(old_count),
                    'cands_after': int(len(new_cands or ())),
                    'fragment_size_before': int(len(fragment)),
                    'candidate_fragment_size': int(len(candidate_fragment)),
                }
        if new_cands:
            preserved_bonds.add(tuple(sorted((u, n))))
            if prof is not None:
                prof['commits'] += 1
                prof['max_cands_after'] = max(
                    int(prof['max_cands_after']), int(len(new_cands)))
            cands = new_cands
            ref_dist = distance.get(u, 0) + 1
            if n_in_mapping and islands_R is not None and n in islands_R:
                target_iid = islands_R[n]
                whole_island = [r for r, k in islands_R.items() if k == target_iid]
                candidate_fragment = candidate_fragment | set(whole_island)
                added_extra = [r for r in whole_island if r not in distance]
                for r in added_extra:
                    distance[r] = ref_dist
            fragment = candidate_fragment
            for r in fragment - old_fragment:
                if r not in distance:
                    distance[r] = ref_dist
                _push_edges_from(heap, used_edges, g_R, r, fragment, graph_floor)
            if prof is not None:
                prof['max_fragment_size'] = max(
                    int(prof['max_fragment_size']), int(len(fragment)))
                prof['max_heap_len'] = max(
                    int(prof['max_heap_len']), int(len(heap)))
            if record:
                added_atoms = sorted(int(r) for r in fragment - old_fragment)
                events.append({
                    'type': 'commit',
                    'scenario': ('merge_island' if n_in_mapping else 'extend_free'),
                    'edge': {'frag_atom': int(u), 'ext_atom': int(n), 'wbo': round(wbo, 3)},
                    'added': int(n), 'element': g_R.nodes[n].get('element'),
                    'atoms_added': added_atoms,
                    'island_size_absorbed': len(added_atoms) if n_in_mapping else None,
                    'island_image': int(mapping[n]) if n_in_mapping else None,
                    'cand_n_value_set': sorted({int(v) for c in cands
                                                for v in cand_possible_values(c, n)}),
                    'cands_before': old_count,
                    'cands_after': len(cands),
                    'represented_assignments_after': represented_assignment_expr(cands),
                    'cands_sample_after': cands_sample(cands, 5),
                    'cands_pattern_after': cands_pattern_sample(cands, 5),
                    'fragment': sorted(int(x) for x in fragment),
                    'p_atoms': sorted({int(v) for c in cands
                                       for v in _cand_possible_p_atoms(c)}),
                    'distance_from_seed': distance[n],
                    'bonds_to_fragment': [(int(u), round(wbo, 3))],
                    'heap_remaining': len(heap),
                    'heap_top': heap_snapshot(heap, used_edges, fragment, mapping, 8),
                    'pool_by_frag_atom': pool_by_frag_atom(heap, used_edges, fragment, mapping, g_R),
                })
        else:
            if prof is not None:
                prof['deferred'] += 1
            deferred_edges.add(tuple(sorted((u, n))))
            if record:
                events.append({
                    'type': 'consumed',
                    'deferred': True,
                    'scenario': ('merge_island' if n_in_mapping else 'extend_free'),
                    'edge': {'frag_atom': int(u), 'ext_atom': int(n), 'wbo': round(wbo, 3),
                             'ext_element': g_R.nodes[n].get('element')},
                    'reason': ('merge_failed' if n_in_mapping else 'cut_all_cands'),
                    'island_image': int(mapping[n]) if n_in_mapping else None,
                    'island_id': int(islands_R[n]) if islands_R and n in islands_R else None,
                    'cands_count': len(cands),
                    'represented_assignments': represented_assignment_expr(cands),
                    'cands_sample': cands_sample(cands, 5),
                    'cands_pattern_sample': cands_pattern_sample(cands, 5),
                    'fragment': sorted(int(x) for x in fragment),
                    'heap_remaining': len(heap),
                    'heap_top': heap_snapshot(heap, used_edges, fragment, mapping, 8),
                    'pool_by_frag_atom': pool_by_frag_atom(heap, used_edges, fragment, mapping, g_R),
                })

    cands = complete_fragment(fragment, preserved_bonds, g_R, g_P, mapping,
                              p_orbits, iso_tol, node_policy, bond_tolerance, prof)

    # heap empty
    if not cands or len(fragment) < min_lock_size:
        _finish_profile(
            'no_cands' if not cands else 'too_small',
            len(cands), fragment)
        if record:
            events.append({
                'type': 'seed_end',
                'result': ('no_cands' if not cands else 'too_small'),
                'final_cands': len(cands),
                'fragment': sorted(int(x) for x in fragment),
                'iso': None,
            })
        return []
    if _set_unique(cands):
        _finish_profile('success', len(cands), fragment, 1)
        if record:
            iso_map = _cand_map(cands[0])
            iso = {int(k): int(v) for k, v in iso_map.items()}
            events.append({
                'type': 'seed_end', 'result': 'success',
                'final_cands': len(cands),
                'fragment': sorted(int(x) for x in fragment),
                'iso': iso,
                'all_isos': [iso],
                'cand_patterns': cands_pattern_sample(cands, 5),
            })
        return [_IsoResult(_cand_map(cands[0]),
                           deferred_edges=deferred_edges,
                           fragment=fragment, preserved_bonds=preserved_bonds,
                           symmetry=_symmetry_state(cands[0]))]
    # Reaction search also bounds the combined live leaves across parents.
    if len(cands) > max_branches:
        _finish_profile('subtree_branch_cap', len(cands), fragment, len(cands))
        if record:
            events.append({
                'type': 'seed_end',
                'result': 'subtree_branch_cap',
                'final_cands': len(cands),
                'max_branches': int(max_branches),
                'fragment': sorted(int(x) for x in fragment),
                'iso': None,
            })
        raise IslandBranchLimitExceeded(
            len(cands), max_branches, seed=int(seed))
    branches = [
        _IsoResult(
            _cand_map(c), deferred_edges=deferred_edges,
            fragment=fragment, preserved_bonds=preserved_bonds, symmetry=_symmetry_state(c))
        for c in cands
    ]
    _finish_profile('branched', len(cands), fragment, len(branches))
    if record:
        events.append({
            'type': 'seed_end', 'result': 'branched',
            'final_cands': len(cands), 'n_branches': len(branches),
            'fragment': sorted(int(x) for x in fragment),
            'iso': {int(k): int(v) for k, v in branches[0].items()},
            'all_isos': [{int(k): int(v) for k, v in c.items()}
                         for c in branches],
            'cand_patterns': cands_pattern_sample(cands, 10),
        })
    return branches
