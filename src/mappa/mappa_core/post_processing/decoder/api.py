"""Approximate local event interpretation on saved atom mappings."""
from collections import OrderedDict
import numpy as np
from .records import EventDecodeConfig, EventInterpretation, SymmetryEventResponse
from .events import SignedEventIndex
from .fragments import FragmentInventory, FragmentDecoder
from ..final_branches import FinalBranchCatalogue

_UNSET = object()


def decode_events(aam, config=None, *, threshold=.5, metal_threshold=.3,
                  event_problem=None, event_slack=_UNSET):
    """Decode independent fragments with fixed ordinary boundary events.

    Unique fragments share preparation; a sweep boundary retains its coupled
    correction. Selection and reported bins use this approximate event model.
    Model bins need not have a complete joint mapping witness. Original mapping
    constraints remain attached and govern explicit realizations without rematching.
    event_problem may supply different scoring weights for the same indexed
    atoms; the catalogue's matching constraints remain unchanged.
    The default window is the model minimum plus six events. Explicit None
    retains all model bins without a window.
    """
    if config is not None and event_slack is not _UNSET:
        raise ValueError('set event_slack in config or as a keyword, not both')
    if config is None:
        config = EventDecodeConfig(threshold, metal_threshold,
                                  6 if event_slack is _UNSET else event_slack)
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
        from .window import select_event_window
        window = select_event_window(catalogue, self.index, self.config.event_slack,
                                     local=local)
        retained = frozenset(window.retained_family_ids)
        families = []
        for fid, family in enumerate(catalogue.families):
            if fid not in retained:
                families.append(())
                continue
            responses = []
            for original, action in zip(family.actions, local.programs.program(family.actions)):
                if action not in self.actions:
                    # Self generators do not describe all D alternatives and
                    # therefore cannot certify the entire relation invariant.
                    if action[0] == 'relation':
                        responses_for_action = None
                    elif action[0] == 'group':
                        responses_for_action = self.index.target_action_responses(action[1])
                    else:
                        raise ValueError('decoding requires fragment relations or placement groups')
                    self.actions[action] = SymmetryEventResponse(original, responses_for_action)
                    if len(self.actions) > 1024:
                        self.actions.popitem(last=False)
                self.actions.move_to_end(action)
                responses.append(self.actions[action])
            families.append(tuple(responses))
        return EventInterpretation(catalogue, self.index, tuple(families),
                                   inventory, local, self.event_problem_sha256, window)
