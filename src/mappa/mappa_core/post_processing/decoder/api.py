"""Anchored fragment conversion and direct event counting."""
from .records import EventDecodeConfig, EventCandidate, DecodedEvents
from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex
from mappa.mappa_core.post_processing.decoder.events import FragmentEventCounter
from mappa.mappa_core.post_processing.decoder.conversion import convert_fragments

def decode_events(aam, config=None):
    """Convert saved fragment branches and count their bond differences.

    Accept an AAMResult or a FinalBranchCatalogue with saved fragment atom sets.
    Distinct pair/anchor contexts are rematched once at the event thresholds.
    Counts and optional minimum+n filtering apply to the new retained branches;
    no hidden-permutation optimization or solver enumeration is performed.
    """
    config = config or EventDecodeConfig()
    conversion = convert_fragments(aam, threshold=config.threshold,
        metal_threshold=config.metal_threshold, branch_limit=config.branch_limit,
        random_seed=config.random_seed)
    catalogue = conversion.catalogue
    index = SignedEventIndex(catalogue.problem, threshold=config.threshold,
                             metal_threshold=config.metal_threshold)
    candidates, reports = {}, []
    counter = FragmentEventCounter(index)
    for branch in conversion.branches:
        family = catalogue.families[branch.family_id]
        total, internal = counter.count(family.mapping, branch.islands)
        complete = internal == 0 and all(conversion.rematches[i].complete
                                        for i, _ in branch.fragment_paths)
        reports.append(dict(family_id=branch.family_id, source_family_ids=branch.source_family_ids,
            complete=complete, internal_events=internal, total=total, solver_queries=0,
            reason='converted_branch' if complete else 'incomplete_fragment_conversion'))
    best = min((r['total'] for r in reports), default=None)
    bound = (best + config.event_slack if best is not None and config.event_slack is not None
             else config.max_events)
    descriptions = {}
    for branch, report in zip(conversion.branches, reports):
        if bound is not None and report['total'] > bound:
            continue
        family = catalogue.families[report['family_id']]
        if family.mapping not in descriptions:
            broken, formed = counter.events(family.mapping, branch.islands)
            descriptions[family.mapping] = index._describe_events(
                index.vector(dict(family.mapping)), broken, formed)
        pattern = descriptions[family.mapping]
        candidate = candidates.setdefault(pattern['id'], EventCandidate(pattern['id'],
            dict(family.mapping), pattern['events'], pattern['total'],
            endpoint_events=pattern['endpoint_events']))
        if report['family_id'] not in candidate.family_ids:
            candidate.family_ids.append(report['family_id'])
    selected = list(candidates.values())
    return DecodedEvents(sorted(selected, key=lambda c: (c.total, c.id)),
                         catalogue, index, reports, config, conversion, counter.statistics)
