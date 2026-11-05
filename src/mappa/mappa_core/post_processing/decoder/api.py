"""Queryable binary-symmetry responses attached to retained matching catalogues."""
import numpy as np
from .records import EventDecodeConfig, EventInterpretation
from .events import SignedEventIndex
from .fragments import FragmentInventory, FragmentDecoder
from ..final_branches import FinalBranchCatalogue


def decode_events(aam, config=None, *, threshold=.5, metal_threshold=.3,
                  event_slack=6, fragment_seconds=5., event_problem=None):
    """Prepare the independent-fragment event model and its minimum+n window.

    Unique fragments share bounded starting-map improvement and binary-symmetry
    responses. The original matching catalogue remains intact; event selections
    cover the chosen starts' symmetry orbits. Ordinary boundaries are fixed;
    the selected sweep correction retains its endpoint dependence. Original
    scoring weights may be supplied for the same indexed atoms in event_problem.
    """
    config = config or EventDecodeConfig(threshold, metal_threshold, event_slack, fragment_seconds)
    catalogue = (aam if isinstance(aam, FinalBranchCatalogue)
                 else FinalBranchCatalogue(aam.problem).add_aam(aam))
    return EventDecodeSession(catalogue.problem, config,
                              event_problem=event_problem).decode(catalogue)


def _validate_atom_identity(problem, event_problem):
    for side in ('reactant', 'product'):
        matching, scoring = getattr(problem, side), getattr(event_problem, side)
        if (matching.elements != scoring.elements or
                not np.array_equal(matching.coordinates, scoring.coordinates)):
            raise ValueError(f'event_problem must preserve {side} atom indexing and coordinates')


class EventDecodeSession:
    """Reuse responses for fixed matching endpoints, scoring endpoints and policy."""
    def __init__(self, problem, config=None, *, event_problem=None):
        config = config or EventDecodeConfig()
        self.config = config
        self.problem_sha256 = FinalBranchCatalogue(problem).problem_sha256
        self.event_problem_sha256 = self.problem_sha256
        if event_problem is None:
            event_problem = problem
        else:
            _validate_atom_identity(problem, event_problem)
            self.event_problem_sha256 = FinalBranchCatalogue(event_problem).problem_sha256
        self.index = SignedEventIndex(event_problem, threshold=config.threshold,
                                     metal_threshold=config.metal_threshold)
        self.fragments = {}

    def decode(self, catalogue):
        if catalogue.problem_sha256 != self.problem_sha256:
            raise ValueError('Catalogue belongs to different endpoints')
        if self.event_problem_sha256 != self.problem_sha256:
            _validate_atom_identity(catalogue.problem, self.index.problem)
        inventory = FragmentInventory(catalogue)
        local = FragmentDecoder(self.index, inventory, fragment_seconds=self.config.fragment_seconds,
                                prepared=self.fragments)
        return EventInterpretation(catalogue, self.index, inventory, local, self.event_problem_sha256,
                                   self.config.event_slack)
