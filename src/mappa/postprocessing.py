"""Separate signed bond-event decoding and conditional symmetry inspection.

Candidates are event classes, not all atom bijections. A completeness flag only
covers the retained complete families within the requested event window.
"""

from dataclasses import dataclass, field
import math
import time

from .event_patterns import SignedEventIndex, extract_path_events
from .final_branches import FinalBranchCatalogue


@dataclass(frozen=True)
class EventDecodeConfig:
    threshold: float = 0.5
    metal_threshold: float | None = 0.3
    max_events: int | None = None
    seconds_per_family: float = 10.0
    max_patterns_per_family: int | None = None

    def __post_init__(self):
        for value in (self.threshold, self.metal_threshold):
            if value is not None and (not math.isfinite(value) or value <= 0):
                raise ValueError("event thresholds must be finite and positive")
        if self.threshold is None:
            raise ValueError("ordinary event threshold is required")
        if not math.isfinite(self.seconds_per_family) or self.seconds_per_family < 0:
            raise ValueError("seconds_per_family must be finite and nonnegative")
        for name, minimum in [("max_events", 0), ("max_patterns_per_family", 1)]:
            value = getattr(self, name)
            if value is not None and (not isinstance(value, int) or value < minimum):
                raise ValueError(f"invalid {name}")


@dataclass
class EventCandidate:
    id: str
    mapping: dict[int, int]
    events: dict
    total: int
    family_ids: list[int] = field(default_factory=list)


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

    @property
    def complete(self):
        return self.catalogue.incomplete_paths == 0 and all(
            r["complete"] for r in self.family_reports
        )

    @property
    def minimum_candidates(self):
        """Observed minima; proven within the saved window only if complete."""
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
        class. Search all catalogue families, including families not yet fully
        decoded. Negative answers concern this saved catalogue only. Budgets
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
            [candidate.mapping[a] for a in range(self.index.n)]
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
                        self.index.describe([mapping[a] for a in range(self.index.n)])[
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
    """Deduplicate final branches and return one witness per signed-event class.

    Accept one AAMResult or a FinalBranchCatalogue containing multiple searches
    of the *same directed problem*. Search settings are never modified. This
    balanced, complete-mapping policy includes explicit H; partial mappings
    remain available through the raw AAM graph, not this decoder.
    """
    config = config or EventDecodeConfig()
    catalogue = (
        aam
        if isinstance(aam, FinalBranchCatalogue)
        else FinalBranchCatalogue(aam.problem).add_aam(aam)
    )
    index = SignedEventIndex(
        catalogue.problem,
        threshold=config.threshold,
        metal_threshold=config.metal_threshold,
    )
    candidates, reports = {}, []
    source_edges = {}
    for i, family in enumerate(catalogue.families):
        floor = family.policy[0]
        if floor not in source_edges:
            source_edges[floor] = tuple((int(a), int(b)) for a, b in zip(index.a, index.b)
                                       if index.r[a, b] >= floor)
        report = extract_path_events(
            family.as_path(catalogue.problem, source_edges[floor]),
            catalogue.problem,
            index,
            max_events=config.max_events,
            seconds=config.seconds_per_family,
            max_patterns=config.max_patterns_per_family,
        )
        for pattern in report["patterns"]:
            candidate = candidates.setdefault(
                pattern["id"],
                EventCandidate(
                    pattern["id"],
                    dict(enumerate(pattern["mapping"])),
                    pattern["events"],
                    pattern["total"],
                ),
            )
            candidate.family_ids.append(i)
        reports.append(
            dict(family_id=i, **{k: v for k, v in report.items() if k != "patterns"})
        )
    return DecodedEvents(
        sorted(candidates.values(), key=lambda c: (c.total, c.id)),
        catalogue,
        index,
        reports,
        config,
    )
