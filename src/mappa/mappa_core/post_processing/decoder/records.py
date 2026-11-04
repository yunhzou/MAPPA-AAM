"""Decoder configuration, witnesses and scoped completion records."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import TYPE_CHECKING
import numpy as np
from mappa.mappa_core.core.symmetry.responses import BondResponse
from mappa.mappa_core.core.symmetry.actions import apply_action
if TYPE_CHECKING:
    from .events import SignedEventIndex
    from .fragments import FragmentInventory, FragmentDecoder
    from ..final_branches import FinalBranchCatalogue

@dataclass(frozen=True)
class EventDecodeConfig:
    threshold: float = 0.5
    metal_threshold: float | None = 0.3

    def __post_init__(self):
        import math
        if self.threshold is None or any(
            v is not None and (not math.isfinite(v) or v <= 0)
            for v in (self.threshold, self.metal_threshold)
        ):
            raise ValueError('event thresholds must be finite and positive')

    def __setstate__(self, state):
        validated = type(self)(**state)
        for name, value in vars(validated).items():
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class SymmetryEventResponse:
    """One original action, annotated without expanding or removing it."""
    action: tuple
    generators: tuple | None  # None means no response-invariance certificate.

    @property
    def response_invariant(self):
        return self.generators is not None and all(preserves for preserves, _ in self.generators)


@dataclass(frozen=True, slots=True)
class BondEventTerm:
    atoms: tuple[int, int]
    action_ids: tuple[int, ...]
    response: BondResponse
    initial_pair: tuple[int, int]
    projected_actions: tuple
    index: SignedEventIndex = field(repr=False, compare=False)
    endpoint_case: str
    manual_cut: bool

    @property
    def target_pairs(self):
        """Inspect pair bounds lazily for relations; joint choices remain required."""
        return self.index.endpoint_pairs(self.initial_pair, self.projected_actions)


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
    inventory: FragmentInventory
    local: FragmentDecoder
    event_problem_sha256: str
    _tables: dict = field(default_factory=dict, repr=False)

    @property
    def matching_problem_sha256(self):
        return self.catalogue.problem_sha256

    @property
    def event_problem(self):
        return self.index.problem

    @property
    def fragment_partitions(self):
        return self.inventory.partitions

    def _table_key(self, family_id, partition):
        family = self.catalogue.families[family_id]
        fragments = self.fragment_partitions[family_id][partition]
        cut_bonds = tuple(sorted({tuple(pair) for origin in family.provenance
                                  for pair in origin.get('cuts', ())}))
        actions = self.local.programs.program(family.actions)
        return family.mapping, actions, fragments, cut_bonds

    def _prepare_table(self, family_id, partition, key):
        from .events import FragmentEventTables
        if key in self._tables:
            return self._tables[key]
        return FragmentEventTables(self.local, self.catalogue.families[family_id],
                                   self.inventory.uses[family_id][partition], key[-1])

    def tables(self, family_id, partition=0):
        """Prepare one saved partition lazily; share bond responses across branches."""
        key = self._table_key(family_id, partition)
        if key not in self._tables:
            self._tables[key] = self._prepare_table(family_id, partition, key)
        return self._tables[key]

    def _export_tables(self):
        """Yield stable table references without retaining export-only working state.

        A shared table is supplied only at its first occurrence. Tables explicitly
        requested through tables() remain owned by the interpretation.
        """
        references = {}
        for family_id, partitions in enumerate(self.fragment_partitions):
            for partition in range(len(partitions)):
                key = self._table_key(family_id, partition)
                if key in references:
                    yield family_id, references[key], None
                    continue
                number = references[key] = len(references)
                table = self._prepare_table(family_id, partition, key)
                yield family_id, number, table
                del table

    def contributions(self, family_id, choices=None, *, partition=0, include_events=False):
        """Compose internal and boundary categories for one permitted choice."""
        vector = self.index.vector(self.realize(family_id, choices))
        return self.tables(family_id, partition).evaluate(vector, include_events=include_events)

    def realize(self, family_id, choices=None):
        """Apply ordered local choices and joint global placement words.

        A relation choice names discovery and source/target self rearrangements;
        a placement choice names signed generator numbers. None retains the
        witness. Discovery predicates constrain D; they are not reapplied to the
        two-sided relaxed alternative.
        """
        family = self.catalogue.families[family_id]
        choices = (None,) * len(family.actions) if choices is None else tuple(choices)
        if len(choices) != len(family.actions):
            raise ValueError('provide one choice per saved action')
        atoms = np.asarray([r for r, _ in family.mapping], dtype=int)
        images = np.asarray([p for _, p in family.mapping], dtype=int)
        for action, choice in zip(family.actions, choices):
            images = apply_action(images, action, choice, self.index.m)
        mapping = {int(r): int(p) for r, p in zip(atoms, images)}
        vector = self.index.vector(mapping)
        self.index._vectors([vector])
        return mapping

    def evaluate(self, family_id, choices=None):
        """Explicit whole-reaction report with a canonical ID for one realization.

        Use contributions() to compose local events without canonicalization.
        """
        mapping = self.realize(family_id, choices)
        return self.index.describe(self.index.vector(mapping))

    def counts(self, family_id, choices, *, batch_size=64, categories=False):
        """Stream requested counts; categories=True returns breaking/forming/strengthening/weakening."""
        if batch_size < 1:
            raise ValueError('batch_size must be positive')
        pair_count = len(self.index.a) + len(self.index.pa)
        batch_size = min(batch_size, max(1, 65536 // max(1, pair_count)))
        count = self.tables(family_id).counts if categories else self.index.counts
        batch = []
        for choice in choices:
            batch.append(self.index.vector(self.realize(family_id, choice)))
            if len(batch) == batch_size:
                yield from count(batch)
                batch.clear()
        if batch:
            yield from count(batch)
