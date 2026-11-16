"""Direction-aware stage calls over the existing single-direction operations."""
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Generic, TypeVar

from mappa.mappa_core.core.search.api import DirectedAAMResult
from mappa.mappa_core.domain import AAMProblem
from mappa.mappa_core.post_processing.decoder.api import decode_events as _decode_events
from mappa.mappa_core.post_processing.final_branches import (
    FinalBranchCatalogue, deduplicate_final_branches,
)
from mappa.mappa_core.utils.io.records import (
    decoded_events_record as _decoded_events_record,
    decoded_events_from_record as _decoded_events_from_record,
)


T = TypeVar('T')


@dataclass(frozen=True)
class DirectionalResults(Mapping[str, T], Generic[T]):
    """Stage results in their own base/target frames, accessible without loops.

    Forward uses the input endpoint order; reverse exchanges the endpoints.
    Only the requested directions are present. Compressed relations and event
    IDs remain direction-local; concrete mappings can be returned to input order.
    """

    forward: T | None = None
    reverse: T | None = None

    def __post_init__(self):
        if self.forward is None and self.reverse is None:
            raise ValueError('at least one direction is required')

    def __getitem__(self, direction):
        if direction not in ('forward', 'reverse') or getattr(self, direction) is None:
            raise KeyError(direction)
        return getattr(self, direction)

    def __iter__(self):
        return (direction for direction in ('forward', 'reverse')
                if getattr(self, direction) is not None)

    def __len__(self):
        return int(self.forward is not None) + int(self.reverse is not None)

    def _map(self, operation):
        return DirectionalResults(**{direction: operation(value, direction)
                                     for direction, value in self.items()})

    def to_input_mapping(self, mapping, *, direction):
        """Convert one concrete mapping to the original endpoint order."""
        self[direction]
        pairs = dict(mapping)
        return {b: a for a, b in pairs.items()} if direction == 'reverse' else pairs


def _directional_search(value):
    """Attach direction names to the existing search_aam_directions result."""
    if not isinstance(value, (tuple, list)):
        return value
    if not value or any(not isinstance(run, DirectedAAMResult) for run in value):
        raise TypeError('expected nonempty search_aam_directions results')
    identities = {FinalBranchCatalogue(run.plan.input_problem).problem_sha256 for run in value}
    if len(identities) != 1:
        raise ValueError('directions must belong to the same endpoint problem')
    results = {}
    for run in value:
        direction = 'reverse' if run.plan.reversed else 'forward'
        if direction in results:
            raise ValueError(f'duplicate {direction} result')
        results[direction] = run.aam
    return DirectionalResults(**results)


def deduplicate(aam):
    """Finish matching reduction for one result or all search directions."""
    aam = _directional_search(aam)
    if isinstance(aam, DirectionalResults):
        return aam._map(lambda result, _: deduplicate_final_branches(result))
    return deduplicate_final_branches(aam)


def decode_events(aam, config=None, *, threshold=.5, metal_threshold=.3,
                  event_slack=6, fragment_seconds=5., event_problem=None):
    """Decode one result or every supplied direction, preserving separate windows.

    For directional input, event_problem uses the original input atom order;
    its endpoints are exchanged automatically for the reverse direction.
    """
    aam = _directional_search(aam)
    settings = dict(threshold=threshold, metal_threshold=metal_threshold,
                    event_slack=event_slack, fragment_seconds=fragment_seconds)
    if isinstance(aam, DirectionalResults):
        def decode(result, direction):
            scoring = event_problem
            if scoring is not None and direction == 'reverse':
                scoring = AAMProblem(scoring.product, scoring.reactant, scoring.name)
            return _decode_events(result, config, event_problem=scoring, **settings)
        return aam._map(decode)
    return _decode_events(aam, config, event_problem=event_problem, **settings)


def decoded_events_record(result, *, array_lists=True):
    """Save one decoded result or all directions in one record."""
    if isinstance(result, DirectionalResults):
        return dict(schema='mappa.directional_decoded_events/v1', directions={
            direction: _decoded_events_record(value, array_lists=array_lists)
            for direction, value in result.items()})
    return _decoded_events_record(result, array_lists=array_lists)


def decoded_events_from_record(record):
    """Restore a single or directional decoded result without new searches."""
    if record['schema'] == 'mappa.directional_decoded_events/v1':
        return DirectionalResults(**{direction: _decoded_events_from_record(value)
                                     for direction, value in record['directions'].items()})
    return _decoded_events_from_record(record)
