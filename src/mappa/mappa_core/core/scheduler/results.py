"""Typed stage handoffs. These records contain no scheduling behavior."""
from dataclasses import dataclass

from mappa.mappa_core.core.search.graph import AAMSearchGraph
from mappa.mappa_core.domain import AAMProblem, AAMSearchConfig
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig
from .storage import ArtifactRef
from mappa.mappa_core.utils.io.records import endpoint_record, problem_from_record


def _problem_record(problem):
    return dict(reactant=endpoint_record(problem.reactant),
                product=endpoint_record(problem.product), name=problem.name)


def _search_input(problem, config, execution):
    return SearchInput(problem_from_record(problem), config, execution)


def _catalogue_result(case, problem, catalogue, scope, stops, search_config, execution):
    return CatalogueResult(case, FinalBranchCatalogue.from_record(problem_from_record(problem), catalogue),
                           scope, stops, search_config, execution)


@dataclass(frozen=True)
class SearchInput:
    problem: AAMProblem
    config: AAMSearchConfig
    execution: str

    def __reduce__(self):
        return _search_input, (_problem_record(self.problem), self.config, self.execution)


@dataclass(frozen=True)
class SeedBatchResult:
    graph: AAMSearchGraph
    indices: tuple[int, ...]
    counts: dict


@dataclass(frozen=True)
class CutResult:
    graph: AAMSearchGraph
    cut: int
    seeds: int
    search_counts: tuple[dict, ...]
    finalization_counts: dict


@dataclass(frozen=True)
class CutScope:
    number: int
    bonds: tuple
    seeds: tuple[int, ...]


@dataclass(frozen=True)
class MatchingResult:
    """All finalized cuts of one declared reaction/direction; graphs stay in store."""
    case: str
    input: ArtifactRef
    scope: tuple[CutScope, ...]
    cuts: tuple[ArtifactRef, ...]


@dataclass(frozen=True)
class CatalogueResult:
    case: str
    catalogue: FinalBranchCatalogue
    scope: tuple[CutScope, ...]
    stops: tuple[tuple[str, int], ...]
    search_config: AAMSearchConfig
    execution: str

    def __reduce__(self):
        return _catalogue_result, (self.case, _problem_record(self.catalogue.problem),
            self.catalogue.to_record(), self.scope, self.stops, self.search_config, self.execution)


@dataclass(frozen=True)
class DecodingResult:
    """Saved original-space event record, with its complete search scope."""
    case: str
    config: EventDecodeConfig
    scope: tuple[CutScope, ...]
    stops: tuple[tuple[str, int], ...]
    record: dict
    search_config: AAMSearchConfig
    execution: str
