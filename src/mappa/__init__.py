"""Composable atom mapping and signed-event decoding.

The main workflow is search_aam -> deduplicate -> decode_events.
TS analysis consumes an explicitly supplied ReactionContext.
"""
from mappa.mappa_core.core.search.api import search_aam, search_aam_checkpoints, AAMCheckpointResult
from mappa.mappa_core.core.search.planning import AAMSearchPlan, plan_aam_search
from mappa.mappa_core.post_processing.hierarchical_mapping import convert_fragments, FragmentConversion
from mappa.mappa_core.post_processing.decoder.api import decode_events
from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig, EventInterpretation
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue, deduplicate_final_branches as deduplicate
from mappa.mappa_core.core.search.queries.family import query_path
from mappa.mappa_core.core.matching.policy import BondTolerance
from mappa.mappa_core.core.matching.api import match_fragment
from mappa.mappa_core.core.matching.representation.records import FragmentMatchConfig, FragmentMatchContext, FragmentMatchResult, Fragment
from mappa.mappa_core.core.search.graph import AAMSearchGraph, SearchContext, SearchState, SearchStop, FragmentTransition, SearchPath, SearchBranch, PathRealization
from mappa.mappa_core.core.search.queries.partial import search_core_assignments
from mappa.mappa_core.domain import AAMProblem, AAMResult, AAMSearchConfig, AAMSearchMetrics, AtomAssignment, CoreAAMBranch, CoreAAMResult, MolecularEndpoint, ReactionContext, ResolvedMechanism, VibrationalModes
from mappa.mappa_workflow.ts_analytics.records import TSMechanismResult, TSResult, TSScore, TSScoringConfig, TransitionStateTarget
from mappa.mappa_workflow.ts_analytics.api import analyze_transition_state
from mappa.mappa_workflow.mapping_analytics.reaction_cores import ReactionCore, extract_reaction_cores
from mappa.mappa_core.utils.io.records import aam_record, aam_from_record, reaction_from_record, reaction_record
from mappa.mappa_core.utils.io.records import decoded_events_record, decoded_events_from_record
from mappa.presentation.viewers import write_aam_bundle
from mappa.mappa_workflow.ts_analytics.serialization import ts_record, write_ts_record
from .mappa_core.domain import (
    AAMHierarchy, AtomBijection, AtomPermutation,
    PermutationGroup,
)

# Stable graph and chemistry primitives are reusable public components.
from mappa.mappa_core.core.matching.graphs import WeightedGraph, WeightedNode, build_graph, build_weighted_graph
from mappa.mappa_core.post_processing.decoder.events import classify_bonds
from mappa.mappa_workflow.ts_analytics.scoring import POLICY as TS_SCORING_POLICY, subset_score
from mappa.mappa_core.core.search.queries.subgraph import SubgraphMatch, SubgraphSearchResult, match_weighted_subgraph

__all__ = [
    "deduplicate", "FinalBranchCatalogue", "EventDecodeConfig", "EventInterpretation",
    "query_path", "convert_fragments", "FragmentConversion", "decode_events",
    "AAMSearchPlan", "plan_aam_search",
    "aam_record", "aam_from_record", "write_aam_bundle",
    "decoded_events_record", "decoded_events_from_record",
    "BondTolerance", "match_fragment", "FragmentMatchConfig", "FragmentMatchContext", "FragmentMatchResult", "Fragment",
    "AAMSearchGraph", "SearchState", "FragmentTransition", "SearchPath", "SearchBranch",
    "SearchContext", "SearchStop",
    "PathRealization", "SubgraphSearchResult",
    "AAMProblem", "AAMResult", "AAMSearchConfig",
    "AAMSearchMetrics",
    "AtomAssignment", "CoreAAMBranch",
    "CoreAAMResult", "MolecularEndpoint",
    "ReactionContext", "ResolvedMechanism",
    "TSMechanismResult", "TSResult", "TSScore", "TSScoringConfig",
    "TransitionStateTarget", "VibrationalModes",
    "AAMHierarchy", "AtomBijection", "AtomPermutation",
    "PermutationGroup",
    "search_aam_checkpoints", "AAMCheckpointResult",
    "search_aam", "search_core_assignments",
    "analyze_transition_state",
    "ReactionCore", "extract_reaction_cores",
    "ts_record", "write_ts_record",
    "reaction_record", "reaction_from_record",
    "WeightedGraph", "WeightedNode", "build_graph", "build_weighted_graph",
    "classify_bonds", "TS_SCORING_POLICY", "subset_score",
    "SubgraphMatch", "match_weighted_subgraph",
]

__version__ = "0.2.0"

# Raw search helpers; signed-event decoding remains in mappa.mappa_core.post_processing.decoder.api.
from mappa.mappa_core.core.search.api import DirectedAAMResult, search_aam_directions
__all__ += ["DirectedAAMResult", "search_aam_directions"]
