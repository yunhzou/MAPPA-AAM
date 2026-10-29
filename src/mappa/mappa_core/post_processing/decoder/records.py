"""Decoder configuration, witnesses and scoped completion records."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import TYPE_CHECKING
import math
if TYPE_CHECKING:
    from .conversion import FragmentConversion
    from .events import SignedEventIndex
    from ..final_branches import FinalBranchCatalogue

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


@dataclass
class DecodedEvents:
    candidates: list[EventCandidate]
    catalogue: FinalBranchCatalogue
    index: SignedEventIndex
    family_reports: list[dict]
    config: EventDecodeConfig
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
