"""Partial core assignments and their saved fragment hierarchy."""

from __future__ import annotations
import time
from collections import Counter
from mappa.mappa_core.core.search.branching.construction import BranchLimitExceeded, find_islands
from mappa.mappa_core.core.search.branching.symmetry import finalize_graph_symmetry
from mappa.mappa_core.core.search.planning import _generate_seed_orders
from mappa.mappa_core.domain import AAMHierarchy, AAMProblem
from mappa.mappa_core.domain import AAMSearchConfig, AtomAssignment, CoreAAMBranch, CoreAAMResult, MolecularEndpoint
from mappa.mappa_core.core.matching.graphs import build_graph
from mappa.mappa_core.core.matching.symmetry.orbits import _nauty_orbits
from mappa.mappa_core.core.matching.symmetry.workspace import SymmetryWorkspace
from mappa.mappa_core.core.matching.representation.state import _SymBlock, _sym_block_assignment_expr
from mappa.mappa_workflow.mapping_analytics.family import project_path_assignments


def search_core_assignments(
        source: MolecularEndpoint, target: MolecularEndpoint, core_atoms, *,
        config: AAMSearchConfig | None = None,
        assignment_limit: int = 20_000) -> CoreAAMResult:
    """Map a source core into a target without losing tuple correlations.

    Search stops once the requested source core has been assigned.  Every
    final branch is projected through its ordered fragment pools and exact
    conditioned target actions. Only requested core tuples are enumerated;
    internal full mappings are never expanded.
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
    target_orbits = _nauty_orbits(
        graph_target, wbo_tol=config.iso_tolerance)
    symmetry_workspace = SymmetryWorkspace(graph_target, config.iso_tolerance)
    problem = AAMProblem(source, target)

    typed_branches = []
    unique = {}
    capped = 0
    for seed_index, order in enumerate(_generate_seed_orders(
            graph_source, n_trials=config.seed_count, seed_selection=config.seed_selection)):
        graph = find_islands(
            graph_source, graph_target, list(order),
            iso_tol=config.iso_tolerance,
            max_branches=config.branch_limit,
            core_R=core,
            stop_when_core_mapped=True,
            p_orbits=target_orbits,
        )
        capped += int(graph.capped)
        graph, _ = finalize_graph_symmetry(
            graph, graph_target, iso_tolerance=config.iso_tolerance,
            workspace=symmetry_workspace)
        for raw in graph.paths():
            if not all(atom in raw.mapping for atom in core):
                continue
            variants = project_path_assignments(
                raw, problem, core, assignment_limit=assignment_limit)
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
