"""Convert saved fragment pairs to the bond-event matching thresholds."""
from dataclasses import dataclass, field
from itertools import product
import time

from .alignment.branch import find_islands
from .final_branches import FinalBranchCatalogue, FinalFamily, canonical_actions
from .frag import bond_event_threshold, build_graph
from .matcher.bond_policy import BondTolerance
from .search_symmetry import finalize_graph_symmetry


@dataclass
class FragmentRematch:
    graph_floor: float
    source_atoms: tuple
    target_atoms: tuple
    anchors: tuple
    graph: object
    paths: tuple
    families: tuple
    complete: bool
    seconds: float

    def to_record(self):
        return dict(graph_floor=self.graph_floor, source_atoms=self.source_atoms,
                    target_atoms=self.target_atoms, anchors=self.anchors,
                    complete=self.complete, capped=self.graph.capped,
                    seconds=self.seconds, search=self.graph.to_record())


@dataclass(frozen=True)
class ConvertedBranch:
    source_family_id: int
    family_id: int
    fragment_paths: tuple
    islands: tuple


@dataclass
class FragmentConversion:
    source: FinalBranchCatalogue
    catalogue: FinalBranchCatalogue
    rematches: list = field(default_factory=list)
    branches: list = field(default_factory=list)
    reused_pairs: int = 0
    seconds: float = 0.

    @property
    def complete(self):
        return (self.source.incomplete_paths == 0
                and all(item.complete for item in self.rematches))


def convert_fragments(aam, *, threshold=.5, metal_threshold=.3,
                      branch_limit=100, random_seed=42):
    """Rematch each distinct saved pair/anchor context with the existing matcher.

    Free choices are rediscovered inside the recorded atom sets. Recombination
    retains compressed actions and fragment partitions, not all bijections.
    The branch limit applies to each rematch's existing growth scheduler.
    """
    from .postprocessing import EventDecodeConfig
    EventDecodeConfig(threshold=threshold, metal_threshold=metal_threshold,
                      branch_limit=branch_limit, random_seed=random_seed)
    source = aam if isinstance(aam, FinalBranchCatalogue) else FinalBranchCatalogue(aam.problem).add_aam(aam)
    if source.families and not (source.branches or source.coupled_branches):
        raise ValueError('Conversion requires saved fragment atom sets, not flattened families alone')
    problem = source.problem
    converted = FinalBranchCatalogue(problem)
    converted.mode = 'anchored_fragment_conversion'
    result = FragmentConversion(source, converted)
    graphs, cache = {}, {}
    started = time.perf_counter()

    def tolerance(graph, a, b):
        return bond_event_threshold(problem.reactant.elements, a, b,
            default_threshold=threshold, metal_threshold=metal_threshold)

    rule = BondTolerance(tolerance, inclusive=False)
    entries = [(policy, pairs, fid) for (policy, pairs), ids in source.branches.items()
               for fid in sorted(ids)]
    entries += [(source.families[fid].policy, pairs, fid)
                for fid, pairs_list in source.coupled_branches.items() for pairs in sorted(pairs_list)]
    for policy, pairs, fid in entries:
        anchors = dict(source.families[fid].fixed_assignments)
        if policy[0] not in graphs:
            graphs[policy[0]] = tuple(build_graph(ep.elements, ep.wbo, bond_cut=policy[0])
                                    for ep in (problem.reactant, problem.product))
        pair_ids = []
        for rs, ps in pairs:
            fixed = tuple((r, anchors[r]) for r in rs if r in anchors)
            key = policy[0], rs, ps, fixed
            if key not in cache:
                reactant, target = (g.subgraph(atoms).copy()
                                   for g, atoms in zip(graphs[policy[0]], (rs, ps)))
                pair_start = time.perf_counter()
                graph = find_islands(reactant, target, sorted(rs), graph_floor=policy[0],
                    iso_tol=threshold, max_branches=branch_limit, anchor_map=dict(fixed),
                    bond_tolerance=rule, random_seed=random_seed)
                graph, _ = finalize_graph_symmetry(graph, target, iso_tolerance=threshold)
                paths = tuple(graph.paths())
                # Flatten each path separately so its compressed choices stay
                # attached to the corresponding recombined branch.
                families = []
                for path in paths:
                    local = FinalBranchCatalogue(problem).add_graph(graph, paths=(path,))
                    families.append(local.families[0])
                complete = (not graph.capped and bool(paths)
                            and all(len(path.mapping) == len(rs) for path in paths))
                cache[key] = len(result.rematches)
                result.rematches.append(FragmentRematch(policy[0], rs, ps, fixed,
                    graph, paths, tuple(families), complete, time.perf_counter()-pair_start))
            else:
                result.reused_pairs += 1
            pair_ids.append(cache[key])
        for selected in product(*(range(len(result.rematches[i].paths)) for i in pair_ids)):
            pieces = [result.rematches[i] for i in pair_ids]
            families = [piece.families[j] for piece, j in zip(pieces, selected)]
            states = [piece.graph.states[piece.paths[j].terminal] for piece, j in zip(pieces, selected)]
            mapping = tuple(sorted(pair for state in states for pair in state.mapping))
            if len(dict(mapping)) != len(mapping) or len({p for _, p in mapping}) != len(mapping):
                raise ValueError('Saved fragment atom sets overlap during recombination')
            if any(dict(mapping).get(r) != p for r, p in anchors.items()):
                raise ValueError('Rematching changed a saved fixed assignment')
            owners = tuple(sorted((r, (i, island)) for i, state in enumerate(states)
                                  for r, island in state.islands))
            required = tuple(sorted({edge for family in families for edge in family.required_edges}))
            actions = canonical_actions(tuple(action for family in families for action in family.actions))
            # Local target graphs can omit the highest indexed endpoint atoms.
            actions = tuple((kind, tuple(tuple(g) + tuple(range(len(g), problem.target_atom_count))
                                         for g in data) if kind == 'group' else data)
                            for kind, data in actions)
            family_id = len(converted.families)
            refs = tuple(zip(pair_ids, selected))
            provenance = [dict(source_family_id=fid, fragment_paths=refs)]
            family = FinalFamily(mapping, required, actions, (policy[0], threshold), provenance)
            converted.families.append(family)
            groups = {}
            for r, owner in owners:
                groups.setdefault(owner, []).append(r)
            assignment = dict(mapping)
            new_pairs = tuple(sorted((tuple(atoms), tuple(sorted(assignment[r] for r in atoms)))
                                     for atoms in groups.values()))
            converted.branches.setdefault((family.policy, new_pairs), set()).add(family_id)
            result.branches.append(ConvertedBranch(fid, family_id, refs, owners))
    converted.path_count = len(result.branches)
    converted.incomplete_paths = source.incomplete_paths + sum(not item.complete for item in result.rematches)
    result.seconds = time.perf_counter()-started
    return result
