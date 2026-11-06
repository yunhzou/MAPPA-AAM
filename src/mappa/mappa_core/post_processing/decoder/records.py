"""Decoder configuration, witnesses and scoped completion records."""
from __future__ import annotations
from dataclasses import dataclass, field, replace
from collections.abc import Sequence
from functools import cached_property
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
    event_slack: int = 6
    fragment_seconds: float = 5.

    def __post_init__(self):
        import math
        if self.threshold is None or any(
            v is not None and (not math.isfinite(v) or v <= 0)
            for v in (self.threshold, self.metal_threshold)
        ):
            raise ValueError('event thresholds must be finite and positive')
        if isinstance(self.event_slack, bool) or not isinstance(self.event_slack, int) or self.event_slack < 0:
            raise ValueError('event_slack must be a nonnegative integer')
        if not math.isfinite(self.fragment_seconds) or self.fragment_seconds <= 0:
            raise ValueError('fragment_seconds must be finite and positive')

    def __setstate__(self, state):
        validated = type(self)(**state)
        for name, value in vars(validated).items():
            object.__setattr__(self, name, value)


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


@dataclass(frozen=True, slots=True)
class BoundEventTerms(Sequence):
    """One shared term sequence and its branch action positions.

    Explicit inspection resolves a term's positions lazily. Computation and
    persistence use the shared terms directly, without copying every bond.
    """
    terms: tuple[BondEventTerm, ...]
    binding: tuple[int, ...]

    def __len__(self):
        return len(self.terms)

    def __getitem__(self, position):
        if isinstance(position, slice):
            return type(self)(self.terms[position], self.binding)
        term = self.terms[position]
        positions = tuple(self.binding[i] for i in term.action_ids)
        return term if positions == term.action_ids else replace(term, action_ids=positions)


@dataclass
class EventInterpretation:
    """Selected fragment responses beside the unchanged matching catalogue.

    Local conditional responses compose with fixed ordinary boundaries and one
    sweep correction. The minimum and window belong to this approximate model,
    not to physical whole-reaction counts or unsearched mappings.
    """
    catalogue: FinalBranchCatalogue
    index: SignedEventIndex
    inventory: FragmentInventory
    local: FragmentDecoder
    event_problem_sha256: str
    event_slack: int = 6
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

    @cached_property
    def cut_contexts(self):
        """Keep alternative sweep histories separate from each other."""
        result = []
        for family, partitions in zip(self.catalogue.families, self.fragment_partitions):
            origins = tuple(family.provenance.metadata) or ({'cuts': ()},)
            choices = [set() for _ in partitions]
            distinct_cuts = {tuple(sorted(tuple(sorted(pair)) for pair in origin.get('cuts', ())))
                             for origin in origins}
            for origin in origins:
                cuts = tuple(sorted(tuple(sorted(pair)) for pair in origin.get('cuts', ())))
                if len(cuts) > 1:
                    raise ValueError('each decoding context must identify at most one sweep bond')
                saved = origin.get('fragment_partition')
                if saved is not None:
                    partition = tuple(sorted(tuple(sorted(atoms)) for atoms in saved))
                    choices[partitions.index(partition)].add(cuts)
                elif len(partitions) == 1 or len(distinct_cuts) == 1:
                    # Older records are unambiguous when only one side varies.
                    for alternatives in choices:
                        alternatives.add(cuts)
                else:
                    raise ValueError('catalogue lacks the association between partitions and sweep histories')
            result.append(tuple(tuple(sorted(alternatives)) for alternatives in choices))
        return tuple(result)

    @cached_property
    def window(self):
        from .window import EventWindow
        return EventWindow(self)

    @property
    def minimum_events(self):
        return self.window.minimum

    @property
    def event_window(self):
        return self.window.maximum

    @property
    def minimum_proven(self):
        return self.window.minimum_proven

    @property
    def minimum_lower_bound(self):
        return self.window.lower_bound

    @property
    def fragments(self):
        """Prepared fragment responses retained by this decoding window."""
        _ = self.window
        return tuple(self.local.prepared[fragment.key] for fragment in self.inventory.fragments
                     if fragment.key in self.local.prepared)

    def select(self, family_id, lower=None, upper=None, *, partition=0, cut=0):
        """Retain a branch's choices inside a requested count or inclusive window."""
        if lower is None:
            lower = self.minimum_lower_bound
            upper = self.event_window if upper is None else upper
        return self.tables(family_id, partition, cut=cut).program.select(lower, upper)

    def minimum(self, family_id, *, partition=0, cut=0):
        """Return a best actual branch witness and its scoped proof status."""
        return self.tables(family_id, partition, cut=cut).minimum_result

    def _table_key(self, family_id, partition, cut=0):
        family = self.catalogue.families[family_id]
        fragments = self.fragment_partitions[family_id][partition]
        cut_bonds = self.cut_contexts[family_id][partition][cut]
        actions = self.local.programs.program(family.actions)
        return family.mapping, actions, fragments, cut_bonds

    def _prepare_table(self, family_id, partition, key, *, boundary=None):
        from .events import FragmentEventTables
        if key in self._tables:
            return self._tables[key]
        return FragmentEventTables(self.local, self.catalogue.families[family_id],
                                   self.inventory.uses[family_id][partition], key[-1], boundary=boundary)

    def tables(self, family_id, partition=0, *, cut=0):
        """Prepare one saved partition lazily; share bond responses across branches."""
        key = self._table_key(family_id, partition, cut)
        if key not in self._tables:
            self._tables[key] = self._prepare_table(family_id, partition, key)
        return self._tables[key]

    def _export_tables(self):
        """Yield stable table references without retaining export-only working state.

        A shared table is supplied only at its first occurrence. Tables explicitly
        requested through tables() remain owned by the interpretation.
        """
        contexts = self.window.contexts
        # Every retained context is prepared by the window. Its responses own
        # the terms; do not keep their construction indexes beside saved copies.
        self.local.release_bond_lookups()
        references = {}
        for family_id, partition, cut in contexts:
            key = self._table_key(family_id, partition, cut)
            if key in references:
                yield family_id, references[key], None
                continue
            number = references[key] = len(references)
            table = self._prepare_table(family_id, partition, key)
            yield family_id, number, table
            del table

    def contributions(self, family_id, choices=None, *, partition=0, cut=0, include_events=False):
        """Compose internal and boundary categories for one permitted choice."""
        vector = self.index.vector(self.realize(family_id, choices))
        result = self.tables(family_id, partition, cut=cut).evaluate(vector, include_events=include_events)
        result['within_window'] = self.minimum_lower_bound <= result['total'] <= self.event_window
        return result

    def realize(self, family_id, choices=None):
        """Apply ordered local choices and joint global placement words.

        Explicit proofs are validated against the original matching catalogue.
        This also permits inspecting discoveries outside the selected decoded
        orbit; it does not assert membership in an EventSelection. An omitted
        discovery uses the selected start. Placement choices remain joint.
        """
        family = self.catalogue.families[family_id]
        choices = (None,) * len(family.actions) if choices is None else tuple(choices)
        if len(choices) != len(family.actions):
            raise ValueError('provide one choice per saved action')
        atoms = np.asarray([r for r, _ in family.mapping], dtype=int)
        images = np.asarray([p for _, p in family.mapping], dtype=int)
        for action, choice in zip(family.actions, choices):
            if action[0] == 'relation':
                selected = self.local.prepare(action[1]).optimization.discovery_witness
                choice = {} if choice is None else dict(choice)
                choice.setdefault('discovery', dict(selected))
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

    def counts(self, family_id, choices, *, partition=0, cut=0, batch_size=64, categories=False):
        """Stream requested counts; categories=True returns breaking/forming/strengthening/weakening."""
        if batch_size < 1:
            raise ValueError('batch_size must be positive')
        pair_count = len(self.index.a) + len(self.index.pa)
        batch_size = min(batch_size, max(1, 65536 // max(1, pair_count)))
        def count(vectors):
            values = self.tables(family_id, partition, cut=cut).counts(vectors)
            return values if categories else np.column_stack((values[:,0] + values[:,3],
                                                              values[:,1] + values[:,2]))
        batch = []
        for choice in choices:
            batch.append(self.index.vector(self.realize(family_id, choice)))
            if len(batch) == batch_size:
                yield from count(batch)
                batch.clear()
        if batch:
            yield from count(batch)
