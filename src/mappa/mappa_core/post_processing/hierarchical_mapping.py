"""Experimental hierarchical mapping: rematch fragments at tighter tolerances.

This changes the retained mapping space. It is explicitly requested and never
used by event decoding of an existing AAM result.
"""
from __future__ import annotations
import math
from dataclasses import dataclass, field
from itertools import chain, product
import time
from .decoder.events import SignedEventIndex
from mappa.mappa_core.core.search.branching.construction import find_islands
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue, FinalFamily, canonical_actions
from mappa.mappa_core.utils.bonds import bond_event_threshold
from mappa.mappa_core.core.matching.graphs import build_graph
from mappa.mappa_core.core.matching.policy import BondTolerance
from mappa.mappa_core.core.search.branching.symmetry import finalize_graph_symmetry


@dataclass(frozen=True)
class RefinementConfig:
    threshold: float = 0.5
    metal_threshold: float | None = 0.3
    max_events: int | None = None
    branch_limit: int = 100
    random_seed: int = 42
    event_slack: int | None = None

    def __post_init__(self):
        for value in (self.threshold, self.metal_threshold):
            if value is not None and (not math.isfinite(value) or value <= 0):
                raise ValueError("event thresholds must be finite and positive")
        if self.threshold is None:
            raise ValueError("ordinary event threshold is required")
        if not isinstance(self.branch_limit, int) or self.branch_limit < 1:
            raise ValueError("branch_limit must be a positive integer")
        if not isinstance(self.random_seed, int):
            raise ValueError("random_seed must be an integer")
        for name, minimum in [("max_events", 0), ("event_slack", 0)]:
            value = getattr(self, name)
            if value is not None and (not isinstance(value, int) or value < minimum):
                raise ValueError(f"invalid {name}")
        if self.max_events is not None and self.event_slack is not None:
            raise ValueError("choose an absolute max_events or relative event_slack")


@dataclass
class RefinementCandidate:
    id: str
    mapping: dict[int, int]
    events: dict
    total: int
    family_ids: list[int] = field(default_factory=list)
    endpoint_events: dict = field(default_factory=dict)


@dataclass
class RefinedMapping:
    candidates: list[RefinementCandidate]
    catalogue: FinalBranchCatalogue
    index: SignedEventIndex
    family_reports: list[dict]
    config: RefinementConfig
    conversion: FragmentConversion
    counting_statistics: dict = field(default_factory=dict)

    @property
    def observed_minimum(self):
        return min((c.total for c in self.candidates), default=None)

    @property
    def minimum_proven(self):
        """Minimum across converted branches, not exhaustive old-family decoding."""
        return self.observed_minimum is not None and self.complete

    @property
    def event_window(self):
        if self.config.event_slack is not None and self.observed_minimum is not None:
            return self.observed_minimum + self.config.event_slack
        return self.config.max_events

    @property
    def complete(self):
        return self.conversion.complete and all(
            r["complete"] for r in self.family_reports
        )

    @property
    def minimum_candidates(self):
        """Lowest direct event counts among retained converted branches."""
        return (
            [c for c in self.candidates if c.total == self.candidates[0].total]
            if self.candidates
            else []
        )

    def symmetry(self, candidate):
        """Compressed ordered actions and constraints of supporting families.

        Pool/group actions may change events. They are not independent allowed
        atom swaps; use query_assignments() to certify a proposed joint assignment.
        """
        self._check_candidate(candidate)
        return tuple(
            dict(
                family_id=i,
                actions=self.catalogue.families[i].actions,
                required_edges=self.catalogue.families[i].required_edges,
                policy=self.catalogue.families[i].policy,
                provenance=tuple(self.catalogue.families[i].provenance),
            )
            for i in candidate.family_ids
        )

    def _check_candidate(self, candidate):
        if not any(candidate is c for c in self.candidates):
            raise ValueError("candidate must belong to this decoded result")


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
    source_family_ids: tuple[int, ...]
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
    RefinementConfig(threshold=threshold, metal_threshold=metal_threshold,
                      branch_limit=branch_limit, random_seed=random_seed)
    source = aam if isinstance(aam, FinalBranchCatalogue) else FinalBranchCatalogue(aam.problem).add_aam(aam)
    if source.families and not (source.branches or source.coupled_branches):
        raise ValueError('Conversion requires saved fragment atom sets, not flattened families alone')
    problem = source.problem
    converted = FinalBranchCatalogue(problem)
    converted.mode = 'anchored_fragment_conversion'
    result = FragmentConversion(source, converted)
    graphs, cache, programs = {}, {}, {}
    started = time.perf_counter()

    def tolerance(graph, a, b):
        return bond_event_threshold(problem.reactant.elements, a, b,
            default_threshold=threshold, metal_threshold=metal_threshold)

    rule = BondTolerance(tolerance)
    entries = ((policy, pairs, fid) for (policy, pairs), ids in source.branches.items()
               for fid in sorted(ids))
    coupled = ((source.families[fid].policy, pairs, fid)
               for fid, pairs_list in source.coupled_branches.items() for pairs in sorted(pairs_list))
    entries = tuple(chain(entries, coupled))
    # A saved uniform matching policy already satisfies the requested rule.
    # Preserve its compressed choices and partitions without another search.
    if (source.mode == 'lossless_saved_family_union'
            and all(policy[1] == threshold == (threshold if metal_threshold is None else metal_threshold)
                    for policy, _, _ in entries)):
        result.catalogue = source
        seen = set()
        for _, pairs, fid in entries:
            owners = tuple(sorted((r, i) for i, (rs, _) in enumerate(sorted(pairs)) for r in rs))
            if (fid, owners) not in seen:
                result.branches.append(ConvertedBranch((fid,), fid, (), owners))
                seen.add((fid, owners))
        result.seconds = time.perf_counter() - started
        return result
    groups = {}
    for policy, pairs, fid in entries:
        pairs = tuple(sorted((tuple(sorted(rs)), tuple(sorted(ps))) for rs, ps in pairs))
        fixed = tuple(sorted(source.families[fid].fixed_assignments))
        groups.setdefault((policy, pairs, fixed), set()).add(fid)
    for (policy, pairs, fixed), ids in groups.items():
        source_ids = tuple(sorted(ids))
        anchors = dict(fixed)
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
                    family = local.families[0]
                    # Pad local actions once per rematch, then share identical
                    # immutable factors across recombined branches.
                    factors = []
                    for kind, data in family.actions:
                        if kind == 'group':
                            data = tuple(tuple(g) + tuple(range(len(g), problem.target_atom_count))
                                         for g in data)
                        factors.append(converted._intern(converted._factors, (kind, data)))
                    family.actions = tuple(factors)
                    families.append(family)
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
            required = converted._intern(converted._required,
                tuple(sorted({edge for family in families for edge in family.required_edges})))
            program = tuple(action for family in families for action in family.actions)
            if program not in programs:
                programs[program] = converted._intern(converted._programs, canonical_actions(program))
            actions = programs[program]
            family_id = len(converted.families)
            refs = tuple(zip(pair_ids, selected))
            provenance = [dict(source_family_ids=source_ids, fragment_paths=refs)]
            family = FinalFamily(mapping, required, actions, (policy[0], threshold), provenance)
            converted.families.append(family)
            islands = {}
            for r, owner in owners:
                islands.setdefault(owner, []).append(r)
            assignment = dict(mapping)
            new_pairs = tuple(sorted((tuple(atoms), tuple(sorted(assignment[r] for r in atoms)))
                                     for atoms in islands.values()))
            converted.branches.setdefault((family.policy, new_pairs), set()).add(family_id)
            result.branches.append(ConvertedBranch(source_ids, family_id, refs, owners))
    converted.path_count = len(result.branches)
    converted.incomplete_paths = source.incomplete_paths + sum(not item.complete for item in result.rematches)
    result.seconds = time.perf_counter()-started
    return result

class FragmentEventCounter:
    """Count each local alternative and interacting pair once per decoding call.

    Regions are converted fragment atom sets, including any internal splits.
    Unmapped material is a separate region: only its boundary is counted.
    No complete-mapping or Cartesian-combination cache is stored.
    """
    def __init__(self, index):
        self.index = index
        self.regions, self.region_ids, self.local = [], {}, []
        self.domains, self.domain_ids = [], {}
        self.topologies, self.boundaries = {}, {}
        self.bond_comparisons = self.branch_counts = 0
        self.edges = []
        for endpoint in (index.problem.reactant, index.problem.product):
            adjacency = [[] for _ in endpoint.elements]
            for a in range(endpoint.atom_count):
                for b in range(a + 1, endpoint.atom_count):
                    threshold = bond_event_threshold(endpoint.elements, a, b,
                        default_threshold=index.threshold, metal_threshold=index.metal_threshold)
                    weight = endpoint.wbo[a, b]
                    if weight > threshold:
                        adjacency[a].append((b, weight, threshold))
                        adjacency[b].append((a, weight, threshold))
            self.edges.append(adjacency)

    def _register(self, source, target, *, unmapped=False):
        key = tuple(sorted(source.items())), tuple(sorted(target.items()))
        if key not in self.region_ids:
            region = len(self.regions)
            self.region_ids[key] = region
            self.regions.append((source, target))
            domain = tuple(sorted(source)), tuple(sorted(target))
            self.domains.append(self.domain_ids.setdefault(domain, len(self.domain_ids)))
            self.local.append(((), ()) if unmapped else self._events(region, region))
        return self.region_ids[key]

    def _events(self, left, right):
        result = []
        for side, opposite in enumerate((self.index.p, self.index.r)):
            images, other = self.regions[left][side], self.regions[right][side]
            events = []
            for a, image in images.items():
                for b, weight, threshold in self.edges[side][a]:
                    if b not in other or (left == right and a >= b):
                        continue
                    partner = other[b]
                    if image is None and partner is None:
                        continue
                    self.bond_comparisons += 1
                    compared = 0. if image is None or partner is None else opposite[image, partner]
                    if weight - compared > threshold:
                        events.append((min(a, b), max(a, b)))
            result.append(tuple(sorted(events)))
        return tuple(result)

    def _partition(self, mapping, islands):
        owners = dict(islands)
        fragments = {}
        for a, b in mapping:
            fragments.setdefault(owners[a], {})[a] = b
        return tuple(self._register(part, {b: a for a, b in part.items()})
                     for part in fragments.values())

    def _topology(self, regions):
        key = tuple(self.domains[i] for i in regions)
        if key not in self.topologies:
            owners = [{a: i for i, region in enumerate(regions) for a in self.regions[region][side]}
                      for side in (0, 1)]
            absent = [{a: None for a in range(n) if a not in own}
                      for n, own in zip((self.index.n, self.index.m), owners)]
            unmatched = self._register(*absent, unmapped=True)
            pairs = set()
            for adjacency, own in zip(self.edges, owners):
                for a, edges in enumerate(adjacency):
                    for b, _, _ in edges:
                        left, right = own.get(a, len(regions)), own.get(b, len(regions))
                        if a < b and left != right:
                            pairs.add(tuple(sorted((left, right))))
            self.topologies[key] = unmatched, tuple(sorted(pairs))
        return self.topologies[key]

    def _contributions(self, mapping, islands):
        regions = self._partition(mapping, islands)
        unmatched, pairs = self._topology(regions)
        local = tuple(self.local[i] for i in regions)
        regions += (unmatched,)
        boundary = []
        for a, b in pairs:
            key = tuple(sorted((regions[a], regions[b])))
            if key not in self.boundaries:
                self.boundaries[key] = self._events(*key)
            boundary.append(self.boundaries[key])
        return local, tuple(boundary)

    def count(self, mapping, islands):
        local, boundary = self._contributions(mapping, islands)
        internal = sum(len(broken) + len(formed) for broken, formed in local)
        total = internal + sum(len(broken) + len(formed) for broken, formed in boundary)
        self.branch_counts += 1
        return total, internal

    def events(self, mapping, islands):
        local, boundary = self._contributions(mapping, islands)
        return tuple(tuple(sorted(pair for contribution in local + boundary
                                  for pair in contribution[side])) for side in (0, 1))

    @property
    def statistics(self):
        return dict(regions=len(self.regions), interacting_alternative_pairs=len(self.boundaries),
                    partitions=len(self.topologies), bond_comparisons=self.bond_comparisons,
                    branch_counts=self.branch_counts)



def refine_mapping(aam, config=None):
    """Convert saved fragment branches and count their bond differences.

    Accept an AAMResult or a FinalBranchCatalogue with saved fragment atom sets.
    Distinct pair/anchor contexts are rematched once at the event thresholds.
    Counts and optional minimum+n filtering apply to the new retained branches;
    no hidden-permutation optimization or solver enumeration is performed.
    """
    config = config or RefinementConfig()
    conversion = convert_fragments(aam, threshold=config.threshold,
        metal_threshold=config.metal_threshold, branch_limit=config.branch_limit,
        random_seed=config.random_seed)
    catalogue = conversion.catalogue
    index = SignedEventIndex(catalogue.problem, threshold=config.threshold,
                             metal_threshold=config.metal_threshold)
    candidates, reports = {}, []
    counter = FragmentEventCounter(index)
    for branch in conversion.branches:
        family = catalogue.families[branch.family_id]
        total, internal = counter.count(family.mapping, branch.islands)
        complete = internal == 0 and all(conversion.rematches[i].complete
                                        for i, _ in branch.fragment_paths)
        reports.append(dict(family_id=branch.family_id, source_family_ids=branch.source_family_ids,
            complete=complete, internal_events=internal, total=total, solver_queries=0,
            reason='converted_branch' if complete else 'incomplete_fragment_conversion'))
    best = min((r['total'] for r in reports), default=None)
    bound = (best + config.event_slack if best is not None and config.event_slack is not None
             else config.max_events)
    descriptions = {}
    for branch, report in zip(conversion.branches, reports):
        if bound is not None and report['total'] > bound:
            continue
        family = catalogue.families[report['family_id']]
        if family.mapping not in descriptions:
            broken, formed = counter.events(family.mapping, branch.islands)
            descriptions[family.mapping] = index._describe_events(
                index.vector(dict(family.mapping)), broken, formed)
        pattern = descriptions[family.mapping]
        candidate = candidates.setdefault(pattern['id'], RefinementCandidate(pattern['id'],
            dict(family.mapping), pattern['events'], pattern['total'],
            endpoint_events=pattern['endpoint_events']))
        if report['family_id'] not in candidate.family_ids:
            candidate.family_ids.append(report['family_id'])
    selected = list(candidates.values())
    return RefinedMapping(sorted(selected, key=lambda c: (c.total, c.id)),
                         catalogue, index, reports, config, conversion, counter.statistics)
