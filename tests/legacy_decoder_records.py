"""Small legacy-format fixtures for historical proof-audit regression tests."""
from mappa.event_patterns import SignedEventIndex
from reference_event_decoder import extract_path_events


def record(catalogue, *, max_events=None, event_slack=None, seconds=30, query_log=None):
    index = SignedEventIndex(catalogue.problem)
    reports, candidates = [], {}
    for i, family in enumerate(catalogue.families):
        result = extract_path_events(family.as_path(catalogue.problem), catalogue.problem, index,
            max_events=max_events, seconds=seconds, max_patterns=None, query_log=query_log)
        reports.append(dict(family_id=i, **{k:v for k,v in result.items() if k != 'patterns'}))
        for pattern in result['patterns']:
            candidate = candidates.setdefault(pattern['id'], dict(pattern, mapping={str(r):p for r,p in index.mapping(pattern['mapping']).items()}, family_ids=[]))
            candidate['family_ids'].append(i)
    minimum = min((c['total'] for c in candidates.values()), default=None)
    bound = minimum+event_slack if minimum is not None and event_slack is not None else max_events
    complete = all(r['complete'] for r in reports)
    return dict(schema='mappa.decoded_events/v1', complete=complete,
        minimum_proven=complete and minimum is not None, observed_minimum=minimum,
        event_window=bound, incomplete_paths=0, event_policy=index.policy,
        config=dict(threshold=.5,metal_threshold=.3,max_events=max_events,event_slack=event_slack,
                    seconds_per_family=seconds,max_patterns_per_family=None),
        family_reports=reports, candidates=sorted((c for c in candidates.values() if bound is None or c['total']<=bound),key=lambda c:(c['total'],c['id'])))
