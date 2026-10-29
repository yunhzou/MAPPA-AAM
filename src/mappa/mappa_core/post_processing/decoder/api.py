"""Anchored fragment conversion and direct event counting."""
from .records import EventDecodeConfig, EventCandidate, DecodedEvents
from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex
from mappa.frag import classify_bonds
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
    candidates, reports, counted = {}, [], {}
    problem = catalogue.problem
    for branch in conversion.branches:
        family = catalogue.families[branch.family_id]
        if family.mapping not in counted:
            broken, formed, _, _ = classify_bonds(dict(family.mapping),
                problem.reactant.wbo, problem.product.wbo,
                dwbo_threshold=config.threshold, metal_dwbo_threshold=config.metal_threshold,
                elements_R=problem.reactant.elements, elements_P=problem.product.elements)
            counted[family.mapping] = dict(broken=tuple((a,b) for a,b,*_ in broken),
                                          formed=tuple((a,b) for a,b,*_ in formed))
        events = counted[family.mapping]
        total = sum(map(len, events.values()))
        owners = dict(branch.islands)
        target_owners = {p: owners[r] for r, p in family.mapping}
        internal = sum(a in own and b in own and own[a] == own[b]
            for kind, own in (('broken', owners), ('formed', target_owners))
            for a, b in events[kind])
        complete = internal == 0 and all(conversion.rematches[i].complete
                                        for i, _ in branch.fragment_paths)
        reports.append(dict(family_id=branch.family_id, source_family_ids=branch.source_family_ids,
            complete=complete, internal_events=internal, total=total, solver_queries=0,
            reason='converted_branch' if complete else 'incomplete_fragment_conversion'))
    best = min((r['total'] for r in reports), default=None)
    bound = (best + config.event_slack if best is not None and config.event_slack is not None
             else config.max_events)
    descriptions = {}
    for report in reports:
        if bound is not None and report['total'] > bound:
            continue
        family = catalogue.families[report['family_id']]
        if family.mapping not in descriptions:
            descriptions[family.mapping] = index.describe(index.vector(dict(family.mapping)))
        pattern = descriptions[family.mapping]
        candidate = candidates.setdefault(pattern['id'], EventCandidate(pattern['id'],
            dict(family.mapping), pattern['events'], pattern['total'],
            endpoint_events=pattern['endpoint_events']))
        candidate.family_ids.append(report['family_id'])
    selected = list(candidates.values())
    return DecodedEvents(sorted(selected, key=lambda c: (c.total, c.id)),
                         catalogue, index, reports, config, conversion)
