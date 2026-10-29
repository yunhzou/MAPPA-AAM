"""Cut planning and saved-branch helpers shared by AAM and explicit analysis.

Search execution belongs to mappa.aam; event conversion to mappa.postprocessing.
"""
from __future__ import annotations

import copy
from collections import deque

from mappa.mappa_core.core.matching.state import _SymBlock, _sym_block_assignment_expr
from mappa.mappa_core.core.matching.orbits import _atom_tuple_orbit, _nauty_atom_generators, _nauty_orbits
from mappa.mappa_core.core.matching.canonical import _CandidateAutomorphismCanonicalizer
from mappa.mappa_core.core.matching.state import candidate_from_record
from .branch import BranchLimitExceeded


def _canon_pair(a, b):
    return (a, b) if a <= b else (b, a)


def _strong_edges(wboR, cut_floor):
    n = int(wboR.shape[0])
    return [(i, j) for i in range(n) for j in range(i + 1, n)
            if float(wboR[i, j]) >= cut_floor]


def _core_mapping_variants(branch, core_R, max_variants, *,
                           g_P, p_orbits=None):
    """Expand final candidate-carried core ambiguity with strict pynauty."""
    core_R = tuple(sorted(set(int(r) for r in core_R)))
    base = {
        int(r): int(branch.mapping[r])
        for r in core_R
        if r in branch.mapping
    }
    if len(base) != len(core_R):
        return []

    blocks = []
    covered = set()
    final_symmetry = _branch_symmetry_record(branch)
    for block in final_symmetry.get('blocks', ()):
        source = block.get('source') or 'sym_block'
        if source not in {
                'sym_block', 'exact_automorph_group'}:
            continue
        r_atoms = tuple(sorted(
            int(r) for r in block.get('r_atoms', ())
            if int(r) in set(core_R) and int(r) not in covered
        ))
        if not r_atoms:
            continue
        p_atoms = tuple(sorted(int(p) for p in block.get('p_atoms', ())))
        if len(p_atoms) < len(r_atoms):
            continue
        blocks.append((r_atoms, frozenset(p_atoms)))
        covered.update(r_atoms)
    if not blocks:
        return [base]

    fixed = {r: base[r] for r in core_R if r not in covered}
    tag_parts = {}

    def add_tag(atom, tag):
        tag_parts.setdefault(int(atom), []).append(tag)

    for r, p in fixed.items():
        add_tag(p, ('fixed', int(r), int(p)))
    for block_index, (_r_atoms, p_atoms) in enumerate(blocks):
        for p in p_atoms:
            add_tag(p, ('block', int(block_index)))

    generators = _nauty_atom_generators(
        g_P,
        wbo_tol=float(getattr(p_orbits, 'wbo_tol', 0.2)),
        atom_color_tags={
            atom: tuple(tags) for atom, tags in tag_parts.items()
        },
    )
    seed = tuple(base[r] for r in core_R)
    states = _atom_tuple_orbit(seed, generators)
    out = []
    for state in states:
        mapping = {int(r): int(p) for r, p in zip(core_R, state)}
        if any(mapping[r] != p for r, p in fixed.items()):
            continue
        if any(mapping[r] not in p_atoms
               for r_atoms, p_atoms in blocks for r in r_atoms):
            continue
        out.append(mapping)
    if len(out) > int(max_variants):
        raise BranchLimitExceeded(
            max_variants,
            branch_count=len(out),
            stage='core_symmetry_expansion',
        )
    return out


def _branch_symmetry_record(branch, symmetry_fragments=None):
    fragments = []
    blocks = []
    p_to_r = {
        int(p): int(r)
        for r, p in getattr(branch, 'mapping', {}).items()
    }
    if symmetry_fragments is None:
        paths = branch.symmetry_paths
        if len(paths) != 1:
            raise ValueError("select an explicit search path before projecting its hierarchy")
        symmetry_fragments = paths[0]
    for frag_index, fragment in enumerate(symmetry_fragments):
        record = {
            'fragment_index': int(frag_index),
            'island_idx': int(fragment.get('island_idx', frag_index)),
            'fragment': [int(r) for r in fragment.get('fragment', ())],
            'deferred_edges': [
                [int(a), int(b)]
                for a, b in fragment.get('deferred_edges', ())
            ],
            'symmetry': fragment.get('symmetry') or {},
        }
        fragments.append(record)
        symmetry = record['symmetry']
        for block_index, block in enumerate(symmetry.get('blocks') or ()):
            r_atoms = {
                int(r) for r in block.get('r_atoms', ())
            }
            p_atoms = [int(p) for p in block.get('p_atoms', ())]
            for p in p_atoms:
                if p in p_to_r:
                    r_atoms.add(p_to_r[p])
            r_atoms = sorted(r_atoms)
            if len(p_atoms) <= 1 and len(r_atoms) <= 1:
                continue
            normalized_block = _SymBlock(tuple(r_atoms), tuple(p_atoms),
                                         extendable=False)
            source = block.get('source')
            assignments = (
                block.get('assignments')
                if source in {'island_automorph', 'sym_block', 'interbranch',
                              'exact_automorph_group'}
                else _sym_block_assignment_expr(normalized_block)
            )
            item = {
                'fragment_index': int(frag_index),
                'block_index': int(block_index),
                'island_idx': record['island_idx'],
                'r_atoms': r_atoms,
                'p_atoms': p_atoms,
                'extendable': bool(block.get('extendable', False)),
                'open': bool(normalized_block.open),
                'assignments': assignments,
            }
            if source:
                item['source'] = source
            blocks.append(item)
    return {
        'rule': 'branch_symmetry_blocks',
        'fragments': fragments,
        'blocks': blocks,
    }


def _color_groups_from_blocks(blocks):
    groups = []
    by_key = {}
    for block in blocks:
        r_atoms = sorted(int(r) for r in block.get('r_atoms', ()))
        p_atoms = sorted(int(p) for p in block.get('p_atoms', ()))
        if len(r_atoms) <= 1 and len(p_atoms) <= 1:
            continue
        source = block.get('source') or 'sym_block'
        key = (tuple(r_atoms), tuple(p_atoms))
        if key in by_key:
            by_key[key]['sources'].add(source)
            continue
        group = {
            'r_atoms': r_atoms,
            'p_atoms': p_atoms,
            'sources': {source},
        }
        by_key[key] = group
        groups.append(group)
    for group in groups:
        group['sources'] = sorted(group['sources'])
    return groups


def _generator_orbit_map(g, generators):
    """Return exact atom orbits generated by pynauty permutations."""
    parent = {int(atom): int(atom) for atom in g.nodes()}

    def find(atom):
        atom = int(atom)
        while parent[atom] != atom:
            parent[atom] = parent[parent[atom]]
            atom = parent[atom]
        return atom

    def union(left, right):
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for generator in generators:
        for atom, image in generator.items():
            union(atom, image)
    return {atom: find(atom) for atom in parent}


def _perm_compose(first, second):
    """Return ``second o first`` for dense image permutations."""
    return tuple(second[first[atom]] for atom in range(len(first)))


def _perm_inverse(permutation):
    inverse = [0] * len(permutation)
    for atom, image in enumerate(permutation):
        inverse[image] = atom
    return tuple(inverse)


def _point_stabilizer_generators(generators, point, degree):
    """Exact Schreier generators for the subgroup fixing ``point``."""
    identity = tuple(range(int(degree)))
    generators = tuple(dict.fromkeys(
        tuple(map(int, generator)) for generator in generators
        if tuple(map(int, generator)) != identity
    ))
    transversals = {int(point): identity}
    queue = deque([int(point)])
    while queue:
        current = queue.popleft()
        transversal = transversals[current]
        for generator in generators:
            image = generator[current]
            if image in transversals:
                continue
            transversals[image] = _perm_compose(transversal, generator)
            queue.append(image)
    stabilizers = []
    seen = set()
    for current, transversal in transversals.items():
        for generator in generators:
            image = generator[current]
            stabilizer = _perm_compose(
                _perm_compose(transversal, generator),
                _perm_inverse(transversals[image]))
            if stabilizer != identity and stabilizer not in seen:
                seen.add(stabilizer)
                stabilizers.append(stabilizer)
    return tuple(stabilizers)


def _stored_generator_orbit_map(atoms, generators):
    parent = {int(atom): int(atom) for atom in atoms}

    def find(atom):
        while parent[atom] != atom:
            parent[atom] = parent[parent[atom]]
            atom = parent[atom]
        return atom

    def union(left, right):
        left, right = find(left), find(right)
        if left != right:
            parent[right] = left

    for generator in generators:
        for atom in tuple(parent):
            image = int(generator[atom])
            if image in parent:
                union(atom, image)
    return {atom: find(atom) for atom in parent}


def complete_chosen_automorphism_groups(branch_symmetry, mapping, g_R, g_P,
                                        wbo_tol, *,
                                        exact_target_generators=None,
                                        exact_branch_mappings=()):
    """Complete chosen-candidate display blocks using exact stabilizers.

    Growth blocks are a compressed search trace and can omit an atom assigned
    earlier as a singleton.  For each mapped center, intersect the exact R and
    P neighbor orbits under automorphisms that fix that center.  The resulting
    sibling groups are therefore mutable at both endpoints and belong to the
    chosen mapping, rather than to discarded branches or orbit history.
    """
    mapping = {int(r): int(p) for r, p in dict(mapping).items()}
    result = dict(branch_symmetry or {})
    complete = []
    if exact_target_generators is not None:
        degree = len(mapping)
        exact_target_generators = tuple(dict.fromkeys(
            tuple(map(int, generator))
            for generator in exact_target_generators))
        if any(len(generator) != degree
               or set(generator) != set(range(degree))
               for generator in exact_target_generators):
            raise ValueError("stored AAM target generator is not bijective")
    exact_branch_mappings = tuple(
        {int(r): int(p) for r, p in dict(candidate).items()}
        for candidate in exact_branch_mappings or ())

    def stabilizer_orbits(g, center, cache):
        if center not in cache:
            generators = _nauty_atom_generators(
                g, wbo_tol=float(wbo_tol),
                atom_color_tags={int(center): ('fixed_center',)})
            cache[center] = _generator_orbit_map(g, generators)
        return cache[center]

    fragments = result.get('fragments') or [{
        'fragment_index': 0,
        'island_idx': 0,
        'fragment': sorted(mapping),
    }]
    for fragment_position, fragment in enumerate(fragments):
        fragment_R = {
            int(atom) for atom in fragment.get('fragment', ())
            if int(atom) in mapping and int(atom) in g_R
        }
        fragment_P = {mapping[atom] for atom in fragment_R}
        if not fragment_R or not fragment_P:
            continue
        # Final mutability belongs to this chosen fragment.  Edges crossing
        # its boundary are masked, just as they were for the candidate.
        sub_R = g_R.subgraph(fragment_R).copy()
        sub_P = g_P.subgraph(fragment_P).copy()
        if exact_target_generators is not None:
            fragment_group_start = len(complete)
            for center_R in sorted(fragment_R):
                center_P = mapping[center_R]
                neighbors_R = [
                    int(atom) for atom in sub_R.neighbors(center_R)
                    if mapping.get(int(atom)) in fragment_P
                    and sub_P.has_edge(center_P, mapping[int(atom)])
                ]
                if len(neighbors_R) < 2:
                    continue
                stabilizers = _point_stabilizer_generators(
                    exact_target_generators, center_P, degree)
                orbit_map = _stored_generator_orbit_map(
                    fragment_P, stabilizers)
                groups = {}
                for atom_R in neighbors_R:
                    atom_P = mapping[atom_R]
                    groups.setdefault(orbit_map[atom_P], []).append(
                        (atom_R, atom_P))
                for pairs in groups.values():
                    if len(pairs) <= 1:
                        continue
                    complete.append({
                        'fragment_index': int(fragment.get(
                            'fragment_index', fragment_position)),
                        'island_idx': int(fragment.get(
                            'island_idx', fragment_position)),
                        'center_R': int(center_R),
                        'center_P': int(center_P),
                        'r_atoms': sorted(r for r, _ in pairs),
                        'p_atoms': sorted(p for _, p in pairs),
                        'extendable': False,
                        'open': False,
                        'assignments': f"{len(pairs)}!",
                        'source': 'stored_AAM_branch_mapping_group',
                    })
            if (len(complete) == fragment_group_start
                    and len(fragment_R) > 1):
                orbit_map = _stored_generator_orbit_map(
                    fragment_P, exact_target_generators)
                groups = {}
                for atom_R in sorted(fragment_R):
                    atom_P = mapping[atom_R]
                    groups.setdefault(orbit_map[atom_P], []).append(
                        (atom_R, atom_P))
                for pairs in groups.values():
                    if len(pairs) <= 1:
                        continue
                    complete.append({
                        'fragment_index': int(fragment.get(
                            'fragment_index', fragment_position)),
                        'island_idx': int(fragment.get(
                            'island_idx', fragment_position)),
                        'r_atoms': sorted(r for r, _ in pairs),
                        'p_atoms': sorted(p for _, p in pairs),
                        'extendable': False,
                        'open': False,
                        'assignments': f"{len(pairs)}!",
                        'source': 'stored_AAM_branch_mapping_group',
                    })
            # Separate exact branch families can carry correlated assignment
            # changes absent from every individual branch stabilizer.  Record
            # only permutations directly witnessed between deduplicated AAM
            # branch representatives with the same center and neighbor set.
            for center_R in sorted(fragment_R):
                center_P = mapping[center_R]
                neighbors_R = tuple(sorted(sub_R.neighbors(center_R)))
                selected_images = tuple(mapping[r] for r in neighbors_R)
                selected_by_image = {
                    mapping[r]: int(r) for r in neighbors_R}
                local_parent = {int(r): int(r) for r in neighbors_R}

                def local_find(atom):
                    while local_parent[atom] != atom:
                        local_parent[atom] = local_parent[local_parent[atom]]
                        atom = local_parent[atom]
                    return atom

                def local_union(left, right):
                    left, right = local_find(left), local_find(right)
                    if left != right:
                        local_parent[right] = left

                for alternative in exact_branch_mappings:
                    if alternative.get(center_R) != center_P:
                        continue
                    if any(r not in alternative for r in neighbors_R):
                        continue
                    alternative_images = tuple(
                        alternative[r] for r in neighbors_R)
                    if set(alternative_images) != set(selected_images):
                        continue
                    for r, image in zip(neighbors_R, alternative_images):
                        local_union(int(r), selected_by_image[image])
                groups = {}
                for atom_R in neighbors_R:
                    groups.setdefault(local_find(int(atom_R)), []).append(
                        int(atom_R))
                for r_atoms in groups.values():
                    if len(r_atoms) <= 1:
                        continue
                    complete.append({
                        'fragment_index': int(fragment.get(
                            'fragment_index', fragment_position)),
                        'island_idx': int(fragment.get(
                            'island_idx', fragment_position)),
                        'center_R': int(center_R),
                        'center_P': int(center_P),
                        'r_atoms': sorted(r_atoms),
                        'p_atoms': sorted(mapping[r] for r in r_atoms),
                        'extendable': False,
                        'open': False,
                        'assignments': 'exact_branch_relation',
                        'source': 'AAM_cross_branch_assignment',
                    })
            continue
        r_cache = {}
        p_cache = {}
        fragment_group_start = len(complete)
        for center_R in sorted(fragment_R):
            center_P = mapping[center_R]
            if sub_R.degree(center_R) < 2 or sub_P.degree(center_P) < 2:
                continue
            r_orbits = stabilizer_orbits(sub_R, center_R, r_cache)
            p_orbits = stabilizer_orbits(sub_P, center_P, p_cache)
            groups = {}
            for atom_R in sub_R.neighbors(center_R):
                atom_P = mapping.get(int(atom_R))
                if atom_P is None or not sub_P.has_edge(center_P, atom_P):
                    continue
                key = (int(r_orbits[atom_R]), int(p_orbits[atom_P]))
                groups.setdefault(key, []).append((int(atom_R), int(atom_P)))
            for pairs in groups.values():
                if len(pairs) <= 1:
                    continue
                r_atoms = sorted(r for r, _ in pairs)
                p_atoms = sorted(p for _, p in pairs)
                complete.append({
                    'fragment_index': int(fragment.get(
                        'fragment_index', fragment_position)),
                    'island_idx': int(fragment.get(
                        'island_idx', fragment_position)),
                    'center_R': int(center_R),
                    'center_P': int(center_P),
                    'r_atoms': r_atoms,
                    'p_atoms': p_atoms,
                    'extendable': False,
                    'open': False,
                    'assignments': f"{len(r_atoms)}!",
                    'source': 'chosen_candidate_automorph',
                })

        # Some fragments have no degree-2 center (a symmetric diatomic or a
        # pair of equivalent disconnected atoms).  In that case the exact
        # whole-fragment orbit is the only meaningful fragment-level group.
        if len(complete) == fragment_group_start and len(fragment_R) > 1:
            r_orbits = _nauty_orbits(sub_R, wbo_tol=float(wbo_tol))
            p_orbits = _nauty_orbits(sub_P, wbo_tol=float(wbo_tol))
            groups = {}
            for atom_R in sorted(fragment_R):
                atom_P = mapping[atom_R]
                key = (int(r_orbits[atom_R]), int(p_orbits[atom_P]))
                groups.setdefault(key, []).append((atom_R, atom_P))
            for pairs in groups.values():
                if len(pairs) <= 1:
                    continue
                r_atoms = sorted(r for r, _ in pairs)
                p_atoms = sorted(p for _, p in pairs)
                complete.append({
                    'fragment_index': int(fragment.get(
                        'fragment_index', fragment_position)),
                    'island_idx': int(fragment.get(
                        'island_idx', fragment_position)),
                    'r_atoms': r_atoms,
                    'p_atoms': p_atoms,
                    'extendable': False,
                    'open': False,
                    'assignments': f"{len(r_atoms)}!",
                    'source': 'chosen_fragment_automorph',
                })

    unique_complete = {}
    for block in complete:
        key = (tuple(block['r_atoms']), tuple(block['p_atoms']))
        unique_complete.setdefault(key, block)
    complete = list(unique_complete.values())

    result['rule'] = 'chosen_analytical_branch_fragment_automorphisms'
    result.pop('witnesses', None)
    result.pop('matching_generators', None)
    result.pop('matching_blocks', None)
    result.pop('dedup_witness_count', None)
    result.pop('selected_witness_index', None)
    result['blocks'] = complete
    result['color_groups'] = _color_groups_from_blocks(complete)
    return result


def _freeze_analytical(value):
    if isinstance(value, dict):
        return tuple(sorted((str(key), _freeze_analytical(item))
                            for key, item in value.items()))
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_analytical(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return tuple(sorted((_freeze_analytical(item) for item in value),
                            key=repr))
    return value


def attach_completed_candidate_groups(branches, g_P, *, wbo_tol,
                                      node_policy=None, return_metrics=False):
    """Attach exact groups after completed branch-family reduction.

    The cache key is the complete locked prefix plus candidate state.  Live
    candidates and discarded growth paths never enter this bounded stage.
    """
    cache = {}
    canonical_base_cache = {}
    completed = []
    metrics = {
        'completed_candidate_group_requests': 0,
        'completed_candidate_group_calculations': 0,
        'completed_candidate_group_cache_hits': 0,
    }
    for raw_branch in branches:
        branch = copy.deepcopy(raw_branch)
        hierarchy = branch.get('hierarchy') or {}
        locked = {}
        for fragment in hierarchy.get('fragments') or ():
            metrics['completed_candidate_group_requests'] += 1
            state = fragment.get('symmetry') or {}
            candidate = candidate_from_record(state)
            key = (tuple(sorted(locked.items())),
                   _freeze_analytical(state))
            generators = cache.get(key)
            if generators is None:
                metrics['completed_candidate_group_calculations'] += 1
                canonicalizer = _CandidateAutomorphismCanonicalizer(
                    g_P, locked_mapping=locked, node_policy=node_policy,
                    wbo_tol=float(wbo_tol),
                    base_cache=canonical_base_cache)
                generators = canonicalizer.atom_generators(candidate)
                cache[key] = generators
            else:
                metrics['completed_candidate_group_cache_hits'] += 1
            state = dict(state)
            state['automorph_generators'] = [
                list(generator) for generator in generators]
            state['automorph_group_source'] = (
                'completed_candidate_after_branch_family_reduction')
            fragment['symmetry'] = state
            for r, p in candidate.mapping.items():
                prior = locked.get(int(r))
                if prior is not None and prior != int(p):
                    raise RuntimeError(
                        "completed AAM fragment conflicts with locked prefix")
                locked[int(r)] = int(p)
        completed.append(branch)
    return (completed, metrics) if return_metrics else completed


def cut_sweep_items(wboR, cut_floor=0.2, *, elements=None,
                    heavy_only=False):
    """Return independent no-cut plus one-edge cut work items.

    ``heavy_only`` keeps explicit hydrogens in AAM mapping and final scoring;
    it only avoids using X-H bonds as artificial graph-search cuts.
    """
    edges = _strong_edges(wboR, float(cut_floor))
    if heavy_only:
        if elements is None:
            raise ValueError("elements are required for heavy-only cut sweep")
        edges = [
            (i, j) for i, j in edges
            if str(elements[i]) != 'H' and str(elements[j]) != 'H'
        ]
    return [()] + [((int(i), int(j)),) for i, j in edges]
