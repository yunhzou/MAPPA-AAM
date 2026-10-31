"""Event interpretation on the original compressed atom mappings."""
from collections import OrderedDict
from .records import EventDecodeConfig, EventInterpretation, SymmetryEventResponse
from .events import SignedEventIndex
from ..final_branches import FinalBranchCatalogue


def decode_events(aam, config=None, *, threshold=.5, metal_threshold=.3):
    """Attach tolerance-dependent responses to original compressed families.

    Unique saved actions share checks. Pool generating swaps and group
    generators are checked directly; no group elements, fragment rematches,
    complete mappings or Cartesian combinations are enumerated. Invariant
    actions remain in place because later actions may depend on their motion.
    """
    config = config or EventDecodeConfig(threshold, metal_threshold)
    catalogue = (aam if isinstance(aam, FinalBranchCatalogue)
                 else FinalBranchCatalogue(aam.problem).add_aam(aam))
    return EventDecodeSession(catalogue.problem, config).decode(catalogue)


class EventDecodeSession:
    """Reuse endpoint responses across catalogues of one reaction and tolerance."""
    def __init__(self, problem, config=None):
        config = config or EventDecodeConfig()
        self.config = config
        self.problem_sha256 = FinalBranchCatalogue(problem).problem_sha256
        self.index = SignedEventIndex(problem, threshold=config.threshold,
                                     metal_threshold=config.metal_threshold)
        self.actions = OrderedDict()

    def decode(self, catalogue):
        if catalogue.problem_sha256 != self.problem_sha256:
            raise ValueError('Catalogue belongs to different endpoints')
        families = []
        for family in catalogue.families:
            responses = []
            for action in family.actions:
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
                        action, self.index.target_action_responses(generators))
                    if len(self.actions) > 1024:
                        self.actions.popitem(last=False)
                self.actions.move_to_end(action)
                responses.append(self.actions[action])
            families.append(tuple(responses))
        return EventInterpretation(catalogue, self.index, tuple(families))
