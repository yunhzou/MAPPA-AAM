"""Stage adapters: cross-cut reduction and catalogue decoding, without rematching."""
from collections import Counter

from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
from mappa.mappa_core.post_processing.decoder.api import EventDecodeSession
from mappa.mappa_core.utils.io.records import decoded_events_record
from .results import MatchingResult, CatalogueResult, DecodingResult


def reduce_matching(payload, workspace):
    matched = workspace.load(payload['matching'])
    if not isinstance(matched, MatchingResult):
        raise TypeError('Deduplication requires a MatchingResult')
    source = workspace.load(matched.input)
    catalogue = FinalBranchCatalogue(source.problem)
    stops = Counter()
    if len(matched.scope) != len(matched.cuts):
        raise ValueError('Missing matching cuts')
    for scope, reference in zip(matched.scope, matched.cuts):
        cut = workspace.load(reference)
        if cut.cut != scope.number or cut.seeds != len(scope.seeds):
            raise ValueError('Cut or seed coverage differs from the matching plan')
        if len(cut.graph.contexts) != len(scope.seeds):
            raise ValueError('Missing seed contexts')
        if any(tuple(c.cuts) != scope.bonds for c in cut.graph.contexts):
            raise ValueError('Graph belongs to another cut')
        catalogue.add_graph(cut.graph, archive=f'{matched.case}/cut/{scope.number:05d}')
        stops.update(stop.reason for stop in cut.graph.stops)
        del cut
    return CatalogueResult(matched.case, catalogue, matched.scope, tuple(sorted(stops.items())),
                           source.config, source.execution)


def decode_catalogue(payload, workspace):
    saved = workspace.load(payload['catalogue'])
    if not isinstance(saved, CatalogueResult):
        raise TypeError('Decoding requires a deduplicated CatalogueResult')
    config = payload['config']
    session = workspace.prepared(('decode', saved.catalogue.problem_sha256, config),
        lambda: EventDecodeSession(saved.catalogue.problem, config))
    interpreted = session.decode(saved.catalogue)
    return DecodingResult(saved.case, config, saved.scope, saved.stops,
                          decoded_events_record(interpreted), saved.search_config, saved.execution)
