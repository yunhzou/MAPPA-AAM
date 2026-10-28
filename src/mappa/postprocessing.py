"""Anchored fragment conversion, direct bond counts, and explicit shuffle queries.

Minimum and completeness concern retained converted branches. Decoding does not
optimize over hidden bijections or call a constraint solver.
"""

from dataclasses import dataclass, field
import math
import time

from .event_patterns import SignedEventIndex
from .frag import classify_bonds
from .fragment_conversion import FragmentConversion, convert_fragments
from .final_branches import FinalBranchCatalogue


@dataclass(frozen=True)
class EventDecodeConfig:
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
class EventCandidate:
    id: str
    mapping: dict[int, int]
    events: dict
    total: int
    family_ids: list[int] = field(default_factory=list)
    endpoint_events: dict = field(default_factory=dict)


@dataclass(frozen=True)
class ShuffleQuery:
    status: str  # allowed / forbidden / unknown
    mapping: dict | None = None
    actions: tuple = ()
    family_id: int | None = None
    reason: str = ""


@dataclass
class DecodedEvents:
    candidates: list[EventCandidate]
    catalogue: FinalBranchCatalogue
    index: SignedEventIndex
    family_reports: list[dict]
    config: EventDecodeConfig
    conversion: FragmentConversion

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
        atom swaps; use query() to certify a proposed joint assignment.
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

    def query(self, candidate, assignments, *, same_events=True, seconds=10.0):
        """Can these R->P assignments coexist in a retained mapping family?

        Unspecified atoms may move jointly. To test exactly one swap, supply
        the entire modified witness. By default retain the candidate's event
        class. Search the converted catalogue, never the old free-choice program. Negative answers concern this saved catalogue only. Budgets
        are soft (encoding and native calls are not preemptible).
        """
        import z3
        from .family_query import compile_path
        from .event_patterns import _event_model, _same_event_orbit

        self._check_candidate(candidate)
        if not math.isfinite(seconds) or seconds < 0:
            raise ValueError("seconds must be finite and nonnegative")
        fixed = dict(assignments)
        problem = self.catalogue.problem
        if any(
            not isinstance(a, int)
            or not isinstance(b, int)
            or a not in range(problem.source_atom_count)
            or b not in range(problem.target_atom_count)
            for a, b in fixed.items()
        ):
            raise ValueError("assignments require valid integer endpoint indices")
        if len(set(fixed.values())) != len(fixed) or any(
            problem.reactant.elements[a] != problem.product.elements[b]
            for a, b in fixed.items()
        ):
            return ShuffleQuery(
                "forbidden", reason="noninjective or element-incompatible assignments"
            )
        deadline = time.perf_counter() + seconds
        pattern = self.index.describe(
            self.index.vector(candidate.mapping)
        )
        unknown = False
        for i, family in enumerate(self.catalogue.families):
            if time.perf_counter() >= deadline:
                return ShuffleQuery("unknown", reason="time budget")
            compiled = compile_path(
                family.as_path(problem),
                problem,
                fixed,
                source_atoms=tuple(fixed),
                complete_reference=False,
            )
            if same_events:
                try:
                    terms, total, _ = _event_model(compiled, self.index, deadline)
                except TimeoutError:
                    return ShuffleQuery("unknown", reason="event encoding budget")
                compiled.solver.add(
                    _same_event_orbit(terms, total, pattern, self.index)
                )
            remaining = deadline - time.perf_counter()
            if remaining <= 0:
                return ShuffleQuery("unknown", reason="encoding budget")
            compiled.solver.set(timeout=max(1, int(remaining * 1000)))
            status = compiled.solver.check()
            if status == z3.sat:
                witness = compiled.realize(compiled.solver.model())
                mapping = dict(witness["mapping"])
                assert all(mapping[a] == b for a, b in fixed.items())
                if same_events:
                    assert (
                        self.index.describe(self.index.vector(mapping))[
                            "id"
                        ]
                        == candidate.id
                    )
                return ShuffleQuery("allowed", mapping, tuple(witness["actions"]), i)
            unknown |= status == z3.unknown
        return ShuffleQuery(
            "unknown" if unknown else "forbidden", reason="saved-family query"
        )


def decode_events(aam, config=None):
    """Convert saved fragment branches and count their bond differences.

    Accept an AAMResult or a FinalBranchCatalogue with saved fragment atom sets.
    Distinct pair/anchor contexts are rematched once at the event thresholds.
    Counts and optional minimum+n filtering apply to the new retained branches;
    no hidden-permutation optimization or solver enumeration is performed.
    """
    config = config or EventDecodeConfig()
    conversion = convert_fragments(aam, threshold=config.threshold,
        metal_threshold=config.metal_threshold, branch_limit=config.branch_limit,
        random_seed=config.random_seed)
    catalogue = conversion.catalogue
    index = SignedEventIndex(catalogue.problem, threshold=config.threshold,
                             metal_threshold=config.metal_threshold)
    candidates, reports, counted = {}, [], {}
    problem = catalogue.problem
    for branch in conversion.branches:
        family = catalogue.families[branch.family_id]
        if family.mapping not in counted:
            broken, formed, _, _ = classify_bonds(dict(family.mapping),
                problem.reactant.wbo, problem.product.wbo,
                dwbo_threshold=config.threshold, metal_dwbo_threshold=config.metal_threshold,
                elements_R=problem.reactant.elements, elements_P=problem.product.elements)
            counted[family.mapping] = dict(broken=tuple((a,b) for a,b,*_ in broken),
                                          formed=tuple((a,b) for a,b,*_ in formed))
        events = counted[family.mapping]
        total = sum(map(len, events.values()))
        owners = dict(branch.islands)
        target_owners = {p: owners[r] for r, p in family.mapping}
        internal = sum(a in own and b in own and own[a] == own[b]
            for kind, own in (('broken', owners), ('formed', target_owners))
            for a, b in events[kind])
        complete = internal == 0 and all(conversion.rematches[i].complete
                                        for i, _ in branch.fragment_paths)
        reports.append(dict(family_id=branch.family_id, source_family_ids=branch.source_family_ids,
            complete=complete, internal_events=internal, total=total, solver_queries=0,
            reason='converted_branch' if complete else 'incomplete_fragment_conversion'))
    best = min((r['total'] for r in reports), default=None)
    bound = (best + config.event_slack if best is not None and config.event_slack is not None
             else config.max_events)
    descriptions = {}
    for report in reports:
        if bound is not None and report['total'] > bound:
            continue
        family = catalogue.families[report['family_id']]
        if family.mapping not in descriptions:
            descriptions[family.mapping] = index.describe(index.vector(dict(family.mapping)))
        pattern = descriptions[family.mapping]
        candidate = candidates.setdefault(pattern['id'], EventCandidate(pattern['id'],
            dict(family.mapping), pattern['events'], pattern['total'],
            endpoint_events=pattern['endpoint_events']))
        candidate.family_ids.append(report['family_id'])
    selected = list(candidates.values())
    return DecodedEvents(sorted(selected, key=lambda c: (c.total, c.id)),
                         catalogue, index, reports, config, conversion)
