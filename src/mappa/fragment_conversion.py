"""Refine saved fragment boundaries without replacing their mapping relation.

The saved families remain the authority for allowed assignments. Rematching
chooses a finer partition in a representative coordinate frame; the original
correlated actions travel with that partition. A boundary count is certified
only when those actions preserve the event responses and no event is internal.
Unresolved families are retained explicitly, never sent to another decoder.
"""
from dataclasses import dataclass, field
from functools import partial
import time

from .fragment import FragmentMatchConfig, FragmentMatchContext, match_fragment
from .matcher.bond_policy import BondTolerance
from .matcher import _nauty_orbits
from .frag import bond_event_threshold, build_graph, classify_bonds
from .event_patterns import SignedEventIndex
from .final_branches import FinalBranchCatalogue


@dataclass(frozen=True)
class ConvertedFragmentBranch:
    family_id: int
    original_pairs: tuple
    fragment_pairs: tuple
    broken: tuple
    formed: tuple
    complete: bool
    reason: str = ''

    @property
    def events(self):
        return len(self.broken) + len(self.formed) if self.complete else None


@dataclass
class ConvertedFragments:
    catalogue: FinalBranchCatalogue
    branches: list
    statistics: dict
    threshold: float
    metal_threshold: float | None

    @property
    def complete(self):
        return (bool(self.branches) and not self.catalogue.incomplete_paths
                and {b.family_id for b in self.branches} == set(range(len(self.catalogue.families)))
                and all(b.complete for b in self.branches))

    @property
    def minimum(self):
        return min(b.events for b in self.branches) if self.complete else None

    def to_record(self):
        from dataclasses import asdict
        return dict(schema='mappa.fragment_conversion/v1',
                    original_catalogue=self.catalogue.to_record(),
                    branches=[asdict(b) | {'events': b.events} for b in self.branches],
                    statistics=self.statistics, complete=self.complete,
                    minimum=self.minimum, threshold=self.threshold,
                    metal_threshold=self.metal_threshold)


@dataclass
class _FragmentPair:
    source: object
    target: object
    source_orbits: object
    target_orbits: object
    partitions: dict = field(default_factory=dict)


class FragmentRematcher:
    """Share intrinsic pair preparation across saved branch constraints."""

    def __init__(self, problem, threshold=.5, metal_threshold=.3):
        self.problem = problem
        self.threshold = threshold
        self.metal_threshold = metal_threshold
        rule = partial(bond_event_threshold, problem.reactant.elements,
                       default_threshold=threshold, metal_threshold=metal_threshold)
        self.policy = BondTolerance(lambda graph, a, b: rule(a, b), inclusive=False)
        self.pairs = {}
        self.graphs = {}
        self.calls = 0
        self.cache_hits = 0

    def refine(self, policy, pair, mapping):
        key = (policy, pair)
        if key not in self.pairs:
            floor, iso = policy
            if floor not in self.graphs:
                self.graphs[floor] = tuple(build_graph(ep.elements, ep.wbo, bond_cut=floor)
                    for ep in (self.problem.reactant, self.problem.product))
            source = self.graphs[floor][0].subgraph(pair[0]).copy()
            target = self.graphs[floor][1].subgraph(pair[1]).copy()
            self.pairs[key] = _FragmentPair(source, target,
                _nauty_orbits(source, wbo_tol=iso), _nauty_orbits(target, wbo_tol=iso))
        prepared = self.pairs[key]
        assignments = tuple((a, mapping[a]) for a in pair[0])
        if assignments in prepared.partitions:
            self.cache_hits += 1
            return prepared.partitions[assignments]
        # These anchors fix a coordinate frame, not the final allowed relation.
        # The result retains the original family's actions and constraints.
        context = FragmentMatchContext(locked_mapping=dict(assignments),
            islands={a: a for a in pair[0]}, source_orbits=prepared.source_orbits,
            target_orbits=prepared.target_orbits)
        config = FragmentMatchConfig(graph_floor=policy[0], iso_tolerance=policy[1],
                                     allow_mapped_seed=True, bond_tolerance=self.policy)
        remaining, pieces = set(pair[0]), []
        while remaining:
            result = match_fragment(prepared.source, prepared.target,
                                    seed=min(remaining), context=context, config=config)
            self.calls += 1
            if result.capped or len(result.matches) != 1:
                raise ValueError('fixed-frame rematching did not return one uncapped placement')
            atoms = tuple(sorted(result.matches[0].fragment & remaining))
            if not atoms:
                raise ValueError('fixed-frame rematching made no progress')
            pieces.append((atoms, tuple(sorted(mapping[a] for a in atoms))))
            remaining.difference_update(atoms)
        result = tuple(sorted(pieces))
        prepared.partitions[assignments] = result
        return result


def convert_aam_fragments(aam, *, threshold=.5, metal_threshold=.3):
    """Convert saved families with certified event-invariant actions.

    Families whose event responses still vary stay explicit and unresolved.
    This initial conversion does not enumerate permutations or use SMT. It
    never treats a representative-only count as a proven family minimum.
    """
    started, cpu = time.perf_counter(), time.process_time()
    catalogue = aam if isinstance(aam, FinalBranchCatalogue) else FinalBranchCatalogue(aam.problem).add_aam(aam)
    problem = catalogue.problem
    index = SignedEventIndex(problem, threshold=threshold, metal_threshold=metal_threshold)
    rematcher = FragmentRematcher(problem, threshold, metal_threshold)
    branches, invariance = [], {}
    entries = [(policy, pairs, family_id) for (policy, pairs), ids in catalogue.branches.items()
               for family_id in sorted(ids)]
    entries += [(catalogue.families[i].policy, pairs, i)
                for i, sources in catalogue.coupled_branches.items() for pairs in sorted(sources)]
    for policy, pairs, family_id in entries:
        family = catalogue.families[family_id]
        if family_id not in invariance:
            family.validate_representative(problem)
            invariance[family_id] = index.invariant_path(family.as_path(problem))
        if not invariance[family_id]:
            branches.append(ConvertedFragmentBranch(family_id, pairs, (), (), (), False,
                'saved symmetry still has event-dependent alternatives'))
            continue
        mapping = dict(family.mapping)
        refined = tuple(sorted(piece for pair in pairs
                               for piece in rematcher.refine(policy, pair, mapping)))
        broken, formed, _, _ = classify_bonds(mapping, problem.reactant.wbo, problem.product.wbo,
            dwbo_threshold=threshold, metal_dwbo_threshold=metal_threshold,
            elements_R=problem.reactant.elements, elements_P=problem.product.elements)
        r_owner = {a: i for i, (rs, _) in enumerate(refined) for a in rs}
        p_owner = {a: i for i, (_, ps) in enumerate(refined) for a in ps}
        internal = any(a in owners and b in owners and owners[a] == owners[b]
                       for events, owners in ((broken, r_owner), (formed, p_owner))
                       for a, b, _, _ in events)
        if internal:
            branches.append(ConvertedFragmentBranch(family_id, pairs, refined, (), (), False,
                'rematched fragment contains a target-only or weak-pair event'))
        else:
            branches.append(ConvertedFragmentBranch(family_id, pairs, refined,
                tuple(broken), tuple(formed), True))
    statistics = dict(families=len(catalogue.families), original_branches=len(entries),
        converted_branches=sum(b.complete for b in branches),
        distinct_fragment_pairs=len(rematcher.pairs), matching_calls=rematcher.calls,
        reused_partitions=rematcher.cache_hits, solver_queries=0,
        expanded_permutations=0, wall_seconds=time.perf_counter()-started,
        cpu_seconds=time.process_time()-cpu)
    return ConvertedFragments(catalogue, branches, statistics, threshold, metal_threshold)
