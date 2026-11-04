"""Decoder configuration and approximate fragment event results."""
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
    from .window import EventWindowSelection

@dataclass(frozen=True)
class EventDecodeConfig:
    threshold: float = 0.5
    metal_threshold: float | None = 0.3
    event_slack: int | None = 6

    def __post_init__(self):
        import math
        from numbers import Integral
        if self.threshold is None or any(
            v is not None and (not math.isfinite(v) or v <= 0)
            for v in (self.threshold, self.metal_threshold)
        ):
            raise ValueError('event thresholds must be finite and positive')
        if self.event_slack is not None and (
                isinstance(self.event_slack, bool) or
                not isinstance(self.event_slack, Integral) or self.event_slack < 0):
            raise ValueError('event_slack must be a nonnegative integer or None')

    def __setstate__(self, state):
        if 'event_slack' not in state:
            raise ValueError('historical decoder configuration requires its producing revision')
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

    @property
    def target_pairs(self):
        """Inspect pair bounds lazily for relations; joint choices remain required."""
        return self.index.endpoint_pairs(self.initial_pair, self.projected_actions)


@dataclass
class EventInterpretation:
    """Independent-fragment event counts on the unchanged saved mapping space.

    Ordinary boundaries contribute fixed events; only the sweep correction
    couples local event contributions. Minima and bins describe this model.
    Requested atom mappings still satisfy the original joint relation.
    """
    catalogue: FinalBranchCatalogue
    index: SignedEventIndex
    family_actions: tuple[tuple[SymmetryEventResponse, ...], ...]
    inventory: FragmentInventory
    local: FragmentDecoder
    event_problem_sha256: str
    window: EventWindowSelection
    _tables: dict = field(default_factory=dict, repr=False)

    @property
    def retained_family_ids(self):
        return self.window.retained_family_ids

    @property
    def minimum_events(self):
        return self.window.minimum

    @property
    def event_window(self):
        return self.window.limit

    @property
    def bin_witnesses_guaranteed(self):
        return False

    def _check_partition(self, family_id, partition):
        if family_id not in self.window.retained_family_ids:
            raise ValueError('family lies outside the requested event window')
        if partition not in self.window.retained_partitions[family_id]:
            raise ValueError('partition lies outside the requested event window')

    def attainable_counts(self, family_id, *, partition=0):
        """Return model bins, which need not have a complete joint mapping witness."""
        self._check_partition(family_id, partition)
        return self.window.event_counts[family_id][partition]

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
        cut_bonds = self.inventory.cuts[family_id][partition]
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
        self._check_partition(family_id, partition)
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
        for family_id in self.retained_family_ids:
            for partition in self.window.retained_partitions[family_id]:
                key = self._table_key(family_id, partition)
                if key in references:
                    yield family_id, references[key], None
                    continue
                number = references[key] = len(references)
                table = self._prepare_table(family_id, partition, key)
                yield family_id, number, table
                del table

    def contributions(self, family_id, choices=None, *, partition=0, include_events=False):
        """Compose approximate categories for one permitted joint choice."""
        vector = self.index.vector(self.realize(family_id, choices, partition=partition))
        return self.tables(family_id, partition).evaluate(vector, include_events=include_events)

    def realize(self, family_id, choices=None, *, partition=0):
        """Apply ordered local choices and joint global placement words.

        A relation choice names discovery and source/target self rearrangements;
        a placement choice names signed generator numbers. None retains the
        witness. Discovery predicates constrain D; they are not reapplied to the
        two-sided relaxed alternative.
        """
        self._check_partition(family_id, partition)
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
        if (self.window.limit is not None and
                int(self.tables(family_id, partition).counts([vector])[0].sum()) > self.window.limit):
            raise ValueError('joint choice lies outside the requested event window')
        return mapping

    def evaluate(self, family_id, choices=None, *, partition=0):
        """Report the permitted mapping and its approximate event contributions."""
        mapping = self.realize(family_id, choices, partition=partition)
        result = self.tables(family_id, partition).evaluate(
            self.index.vector(mapping), include_events=True)
        return dict(result, mapping=mapping, model='independent_fragments_with_sweep',
                    approximate=True, bin_witnesses_guaranteed=False)

    def counts(self, family_id, choices, *, partition=0, batch_size=64, categories=False):
        """Stream requested counts; categories=True returns breaking/forming/strengthening/weakening."""
        if batch_size < 1:
            raise ValueError('batch_size must be positive')
        pair_count = len(self.index.a) + len(self.index.pa)
        batch_size = min(batch_size, max(1, 65536 // max(1, pair_count)))
        table = self.tables(family_id, partition)
        def count(vectors):
            rows = table.counts(vectors)
            return rows if categories else rows[:, (0, 1)] + rows[:, (3, 2)]
        batch = []
        for choice in choices:
            batch.append(self.index.vector(self.realize(family_id, choice, partition=partition)))
            if len(batch) == batch_size:
                yield from count(batch)
                batch.clear()
        if batch:
            yield from count(batch)
