"""Decoder configuration, witnesses and scoped completion records."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import TYPE_CHECKING
import math
from numbers import Integral
import numpy as np
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


@dataclass(frozen=True)
class SymmetryEventResponse:
    """One original action, annotated without expanding or removing it."""
    action: tuple
    generators: tuple  # (response-preserving certificate, changed pair count)

    @property
    def response_invariant(self):
        return all(preserves for preserves, _ in self.generators)


@dataclass
class EventInterpretation:
    """Bond responses attached to the unchanged, ordered saved mapping space.

    This is a compressed event function, not an enumerated pattern catalogue.
    Invariance concerns bond responses/counts; labelled unmatched target atoms
    can still move. No minimum or exhaustive event-window claim is made.
    """
    catalogue: FinalBranchCatalogue
    index: SignedEventIndex
    family_actions: tuple[tuple[SymmetryEventResponse, ...], ...]

    def realize(self, family_id, choices=None):
        """Apply explicitly requested choices and enforce original constraints.

        One choice per saved action, in saved order: a pool choice lists the
        images of its stored atoms; a group choice is a word of signed 1-based
        generator numbers (negative means inverse). None selects identity.
        No group closure or Cartesian product is constructed.
        """
        family = self.catalogue.families[family_id]
        choices = (None,) * len(family.actions) if choices is None else tuple(choices)
        if len(choices) != len(family.actions):
            raise ValueError('provide one choice per saved action')
        atoms = np.asarray([r for r, _ in family.mapping], dtype=int)
        images = np.asarray([p for _, p in family.mapping], dtype=int)
        for (kind, data), choice in zip(family.actions, choices):
            if choice is None:
                continue
            if kind == 'pool':
                selected = tuple(choice)
                if (any(not isinstance(p, Integral) for p in selected)
                        or sorted(selected) != sorted(data)):
                    raise ValueError('pool choice must permute exactly its saved atoms')
                permutation = np.arange(self.index.m)
                permutation[list(data)] = selected
                images = permutation[images]
            elif kind == 'group':
                for number in choice:
                    if not isinstance(number, Integral) or not 1 <= abs(number) <= len(data):
                        raise ValueError('group word must use signed 1-based generator numbers')
                    permutation = np.asarray(data[abs(number) - 1], dtype=int)
                    images = (permutation if number > 0 else np.argsort(permutation))[images]
            else:
                raise ValueError(f'unknown saved action: {kind}')
        mapping = {int(r): int(p) for r, p in zip(atoms, images)}
        vector = self.index.vector(mapping)
        self.index._vectors([vector])
        # Actions preserve correlation; required edges additionally restrict
        # their allowed compositions at the ORIGINAL matching tolerance.
        if family.required_edges:
            a, b = np.asarray(family.required_edges, dtype=int).T
            vector = np.asarray(vector)
            weights = self.index.p[vector[a], vector[b]]
            floor, iso = family.policy
            if np.any(weights < floor) or np.any(np.abs(self.index.r[a, b] - weights) > iso + 1e-9):
                raise ValueError('choices violate original fragment constraints')
        return mapping

    def evaluate(self, family_id, choices=None):
        """Return events for one requested realization, without rematching."""
        mapping = self.realize(family_id, choices)
        return self.index.describe(self.index.vector(mapping))

    def counts(self, family_id, choices, *, batch_size=64):
        """Stream bounded vector batches for caller-supplied choices only."""
        if batch_size < 1:
            raise ValueError('batch_size must be positive')
        batch = []
        for choice in choices:
            batch.append(self.index.vector(self.realize(family_id, choice)))
            if len(batch) == batch_size:
                yield from self.index.counts(batch)
                batch.clear()
        if batch:
            yield from self.index.counts(batch)
