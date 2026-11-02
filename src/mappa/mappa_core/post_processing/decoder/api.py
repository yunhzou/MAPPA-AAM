"""Event interpretation on the original compressed atom mappings."""
from collections import OrderedDict
import numpy as np
from .records import EventDecodeConfig, EventInterpretation, SymmetryEventResponse
from .events import SignedEventIndex
from .fragments import FragmentInventory, FragmentDecoder
from ..final_branches import FinalBranchCatalogue


def decode_events(aam, config=None, *, threshold=.5, metal_threshold=.3,
                  event_problem=None):
    """Attach tolerance-dependent responses to original compressed families.

    Unique saved actions share checks. Pool generating swaps and group
    generators are checked directly; no group elements, fragment rematches,
    complete mappings or Cartesian combinations are enumerated. Invariant
    actions remain in place because later actions may depend on their motion.
    event_problem may supply different scoring weights for the same indexed
    atoms; the catalogue's matching constraints remain unchanged.
    """
    config = config or EventDecodeConfig(threshold, metal_threshold)
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
        self.actions = OrderedDict()

    def decode(self, catalogue):
        if catalogue.problem_sha256 != self.problem_sha256:
            raise ValueError('Catalogue belongs to different endpoints')
        if self.event_problem_sha256 != self.problem_sha256:
            _validate_atom_identity(catalogue.problem, self.index.problem)
        inventory = FragmentInventory(catalogue)
        local = FragmentDecoder(self.index, inventory)
        families = []
        for family in catalogue.families:
            responses = []
            for original, action in zip(family.actions, local.programs.program(family.actions)):
                if action not in self.actions:
                    kind, data = action
                    if kind == 'group':
                        generators = data
                    elif kind == 'pool':
                        generators = []
                        for a, b in zip(data, data[1:]):
                            images = list(range(self.index.m))
                            images[a], images[b] = images[b], images[a]
                            generators.append(tuple(images))
                    else:
                        raise ValueError(f'unknown saved action: {kind}')
                    self.actions[action] = SymmetryEventResponse(
                        original, self.index.target_action_responses(generators))
                    if len(self.actions) > 1024:
                        self.actions.popitem(last=False)
                self.actions.move_to_end(action)
                responses.append(self.actions[action])
            families.append(tuple(responses))
        return EventInterpretation(catalogue, self.index, tuple(families),
                                   inventory, local, self.event_problem_sha256)
