"""Partial core assignments and their saved fragment hierarchy."""

from __future__ import annotations
import time
from collections import Counter
from mappa.mappa_core.core.search.branching.construction import BranchLimitExceeded, find_islands
from mappa.mappa_core.core.search.planning import _generate_seed_orders
from mappa.mappa_core.domain import AAMHierarchy
from mappa.domain import (
    AAMSearchConfig,
    AtomAssignment,
    CoreAAMBranch,
    CoreAAMResult,
    MolecularEndpoint,
)
from mappa.frag import build_graph
from mappa.mappa_core.core.matching.symmetry.orbits import _nauty_orbits, _atom_tuple_orbit, _nauty_atom_generators
from mappa.mappa_core.core.matching.representation.state import _SymBlock, _sym_block_assignment_expr


def search_core_assignments(
        source: MolecularEndpoint, target: MolecularEndpoint, core_atoms, *,
        config: AAMSearchConfig | None = None,
        assignment_limit: int = 20_000) -> CoreAAMResult:
    """Map a source core into a target without losing tuple correlations.

    Search stops once the requested source core has been assigned.  Every
    final branch is expanded only by its candidate-carried exact target
    automorphism action.  Deduplication is then performed on complete core
    tuples, never on independent vertex orbits.
    """
    config = config or AAMSearchConfig()
    core = tuple(sorted(set(map(int, core_atoms))))
    if not core:
        raise ValueError("partial AAM requires a non-empty source core")
    if source.atom_count != target.atom_count:
        raise ValueError("partial AAM currently requires equal atom counts")
    if Counter(source.elements) != Counter(target.elements):
        raise ValueError("partial AAM requires equal endpoint compositions")
    started = time.perf_counter()
    graph_source = build_graph(
        source.elements, source.wbo, bond_cut=config.graph_floor)
    graph_target = build_graph(
        target.elements, target.wbo, bond_cut=config.graph_floor)
    source_orbits = _nauty_orbits(
        graph_source, wbo_tol=config.iso_tolerance)
    target_orbits = _nauty_orbits(
        graph_target, wbo_tol=config.iso_tolerance)

    typed_branches = []
    unique = {}
    capped = 0
    for seed_index, order in enumerate(_generate_seed_orders(
            graph_source, n_trials=config.seed_count, seed_selection=config.seed_selection)):
        graph = find_islands(
            graph_source, graph_target, list(order),
            iso_tol=config.iso_tolerance,
            max_branches=config.branch_limit,
            dwbo_threshold=config.event_threshold,
            metal_dwbo_threshold=config.metal_event_threshold,
            symmetry_wbo_tol=config.iso_tolerance,
            core_R=core,
            stop_when_core_mapped=True,
            p_orbits=target_orbits,
            r_orbits=source_orbits,
        )
        capped += int(graph.capped)
        for raw in graph.paths():
            if not all(atom in raw.mapping for atom in core):
                continue
            variants = _core_mapping_variants(
                raw, core, assignment_limit,
                g_P=graph_target, p_orbits=target_orbits)
            exact = tuple(AtomAssignment.from_mapping(mapping)
                          for mapping in variants)
            representative = AtomAssignment.from_mapping({
                atom: raw.mapping[atom] for atom in core})
            hierarchy = AAMHierarchy.from_record(
                _branch_symmetry_record(raw))
            typed_branches.append(CoreAAMBranch(
                representative=representative,
                hierarchy=hierarchy,
                exact_assignments=exact,
                seed_index=seed_index,
                search_path=raw,
            ))
            for assignment in exact:
                unique.setdefault(assignment.pairs, assignment)
                if len(unique) > int(assignment_limit):
                    raise BranchLimitExceeded(
                        assignment_limit, branch_count=len(unique),
                        stage="exact_core_assignment_union")

    if capped and not typed_branches:
        raise BranchLimitExceeded(
            config.branch_limit, branch_count=config.branch_limit + 1,
            stage="all_core_seed_paths_capped")
    return CoreAAMResult(
        source=source,
        target=target,
        core_atoms=core,
        branches=tuple(typed_branches),
        assignments=tuple(unique[key] for key in sorted(unique)),
        elapsed_seconds=time.perf_counter() - started,
        capped_seed_count=capped,
    )


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
