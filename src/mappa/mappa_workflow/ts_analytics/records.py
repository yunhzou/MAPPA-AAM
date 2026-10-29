"""Typed transition-state inputs, scoring settings and results."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Mapping
import numpy as np
from mappa.mappa_core.domain import MolecularEndpoint, VibrationalModes, AtomAssignment, CoreAAMResult, ResolvedMechanism, ReactionContext

@dataclass(frozen=True)
class TransitionStateTarget:
    """A molecular target and its vibrational analysis."""

    molecule: MolecularEndpoint
    vibrations: VibrationalModes
    kind: str = "candidate"
    reactant_atom_mapping: tuple[int, ...] | None = None

    def __post_init__(self):
        if self.vibrations.displacements.shape[1] != self.molecule.atom_count:
            raise ValueError("normal modes and target atom counts differ")
        object.__setattr__(self, "kind", str(self.kind))
        if self.reactant_atom_mapping is not None:
            images = tuple(self.reactant_atom_mapping)
            if any(not isinstance(i, (int, np.integer)) or isinstance(i, bool) for i in images):
                raise ValueError("reactant_atom_mapping requires integer indices")
            object.__setattr__(self, "reactant_atom_mapping", images)


@dataclass(frozen=True)
class TSScore:
    """Mode score for one exact mechanism-local core assignment."""

    assignment: AtomAssignment
    sources: frozenset[str]
    score: float
    wbo_progress: float
    mode_index: int
    frequency: float
    event_terms: tuple[Mapping[str, Any], ...]
    diagnostics: Mapping[str, Any]


@dataclass(frozen=True)
class TSScoringConfig:
    prefer_endpoint_consensus: bool = True
    core_assignment_limit: int = 20_000

    def __post_init__(self):
        if self.core_assignment_limit < 1:
            raise ValueError("core assignment limit must be positive")


@dataclass(frozen=True)
class TSMechanismResult:
    mechanism: ResolvedMechanism
    target: TransitionStateTarget
    reactant_core_aam: CoreAAMResult | None
    product_core_aam: CoreAAMResult | None
    candidates: tuple[TSScore, ...]
    selected: TSScore | None
    status: str = "scored"
    reason: str | None = None


@dataclass(frozen=True)
class TSResult:
    reaction: ReactionContext
    mechanisms: tuple[TSMechanismResult, ...]
    elapsed_seconds: float

    def __post_init__(self):
        object.__setattr__(self, "mechanisms", tuple(self.mechanisms))
        object.__setattr__(self, "elapsed_seconds", float(self.elapsed_seconds))
