"""Distinct net bond changes among witnesses already retained by decoding."""
from dataclasses import dataclass, replace

from mappa.mappa_core.post_processing.decoder.records import EventInterpretation
from .api import DirectionalResults


@dataclass(frozen=True)
class ReactionCore:
    """One signed event pattern, a representative and its decoded contexts.

    Decreased bonds use base indices; increased bonds use target indices.
    They include weakening and strengthening, respectively. Model counts refer
    to the supporting contexts; physical_event_count refers to the mapping.
    """

    id: str
    mapping: tuple[tuple[int, int], ...]
    decreased_bonds: tuple[tuple[int, int], ...]
    increased_bonds: tuple[tuple[int, int], ...]
    contexts: tuple[tuple[int, int, int], ...] = ()
    model_event_counts: tuple[int, ...] = ()

    @property
    def physical_event_count(self):
        return len(self.decreased_bonds) + len(self.increased_bonds)

    @property
    def base_atoms(self):
        inverse = {target: base for base, target in self.mapping}
        return tuple(sorted({a for pair in self.decreased_bonds for a in pair}
                            | {inverse[a] for pair in self.increased_bonds
                               for a in pair if a in inverse}))

    @property
    def target_atoms(self):
        mapping = dict(self.mapping)
        return tuple(sorted({a for pair in self.increased_bonds for a in pair}
                            | {mapping[a] for pair in self.decreased_bonds
                               for a in pair if a in mapping}))


def extract_reaction_cores(decoded: EventInterpretation | DirectionalResults, *, event_count=None):
    """Group saved witnesses inside the decoded window by signed bond changes.

    Pass event_count to restrict the lookup to one model count, or 'minimum'
    for each direction's detected minimum. Otherwise use all recorded witnesses
    inside the window. Directional inputs return DirectionalResults, preserving
    each direction's indexing and event identity. Symmetry-equivalent event
    patterns share a core even when fragment compositions differ. Each core
    retains its (family, partition, cut) references and a representative from
    the lowest recorded model count in that group.

    This explicit post-processing operation performs no matching, minimum
    search or enumeration of internal shuffles. It reports observed cores, not
    every core attainable by the compressed responses or an elementary pathway.
    """
    if isinstance(decoded, DirectionalResults):
        return decoded._map(lambda value, _: extract_reaction_cores(value, event_count=event_count))
    if event_count == 'minimum':
        event_count = decoded.minimum_events
    if event_count is not None and (isinstance(event_count, bool)
                                   or not isinstance(event_count, int) or event_count < 0):
        raise ValueError("event_count must be a nonnegative integer or 'minimum'")
    window = decoded.window
    if set(window.branch_minima) != set(window.contexts):
        raise ValueError('this archive lacks saved branch counts')
    cores, contexts, counts = {}, {}, {}
    for context in sorted(window.contexts, key=lambda key: (window.branch_minima[key][0], key)):
        count = window.branch_minima[context][0]
        if not window.lower_bound <= count <= window.maximum:
            continue
        if event_count is not None and count != event_count:
            continue
        family, partition, cut = context
        mapping = decoded.saved_minimum(family, partition=partition, cut=cut).witness
        # Event identity already has a bounded cache. Retain a full mapping only
        # for each output core, not for every branch or equivalent shuffle.
        report = decoded.index.describe(decoded.index.vector(dict(mapping)))
        core_id = report['id']
        if core_id not in cores:
            cores[core_id] = ReactionCore(core_id, mapping,
                report['endpoint_events']['broken'], report['endpoint_events']['formed'])
            contexts[core_id], counts[core_id] = [], set()
        contexts[core_id].append(context)
        counts[core_id].add(count)
    return tuple(replace(core, contexts=tuple(contexts[key]), model_event_counts=tuple(sorted(counts[key])))
                 for key, core in cores.items())
