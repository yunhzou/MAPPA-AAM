"""Event interpretation on the original compressed atom mappings."""
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
    threshold, metal_threshold = config.threshold, config.metal_threshold
    catalogue = (aam if isinstance(aam, FinalBranchCatalogue)
                 else FinalBranchCatalogue(aam.problem).add_aam(aam))
    index = SignedEventIndex(catalogue.problem, threshold=threshold,
                             metal_threshold=metal_threshold)
    cache, families = {}, []
    for family in catalogue.families:
        responses = []
        for action in family.actions:
            if action not in cache:
                kind, data = action
                if kind == 'group':
                    generators = data
                elif kind == 'pool':
                    generators = []
                    for a, b in zip(data, data[1:]):
                        images = list(range(index.m))
                        images[a], images[b] = images[b], images[a]
                        generators.append(tuple(images))
                else:
                    raise ValueError(f'unknown saved action: {kind}')
                cache[action] = SymmetryEventResponse(action, index.target_action_responses(generators))
            responses.append(cache[action])
        families.append(tuple(responses))
    return EventInterpretation(catalogue, index, tuple(families))
