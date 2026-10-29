"""Composable atom mapping and signed-event decoding.

The main workflow is search_aam -> mappa.mappa_core.post_processing.decoder.api.decode_events.
TS analysis consumes an explicitly supplied ReactionContext.
"""
from mappa.mappa_core.core.search.api import search_aam, search_aam_checkpoints, AAMCheckpointResult
from mappa.mappa_core.core.search.planning import AAMSearchPlan, plan_aam_search
from mappa.mappa_core.post_processing.decoder.conversion import convert_fragments, FragmentConversion
from mappa.mappa_workflow.mapping_analytics.family import query_path
from mappa.mappa_core.core.matching.policy import BondTolerance
from mappa.mappa_core.core.matching.api import match_fragment
from mappa.mappa_core.core.matching.representation.records import FragmentMatchConfig, FragmentMatchContext, FragmentMatchResult, FragmentPlacement
from mappa.mappa_core.core.search.graph import AAMSearchGraph, SearchContext, SearchState, SearchStop, FragmentTransition, SearchPath, SearchBranch, PathRealization
from mappa.mappa_core.core.search.branching.symmetry import finalize_graph_symmetry
from mappa.mappa_core.core.search.queries.partial import search_core_assignments
from mappa.mappa_core.domain import AAMProblem, AAMResult, AAMSearchConfig, AAMSearchMetrics, AtomAssignment, CoreAAMBranch, CoreAAMResult, MolecularEndpoint, ReactionContext, ResolvedMechanism, VibrationalModes
from mappa.mappa_workflow.ts_analytics.records import TSMechanismResult, TSResult, TSScore, TSScoringConfig, TransitionStateTarget
from mappa.mappa_workflow.ts_analytics.api import analyze_transition_state
from mappa.mappa_core.utils.io.records import aam_record, aam_from_record, reaction_from_record, reaction_record
from mappa.presentation.viewers import write_aam_bundle
from mappa.mappa_workflow.ts_analytics.serialization import ts_record, write_ts_record
from .mappa_core.domain import (
    AAMHierarchy, AtomBijection, AtomPermutation,
    FragmentMatch, PermutationGroup, SymmetryDomain,
)

# Stable graph and chemistry primitives are reusable public components.
from mappa.mappa_core.core.matching.graphs import WeightedGraph, WeightedNode, build_graph, build_weighted_graph
from mappa.mappa_core.post_processing.decoder.events import classify_bonds
from mappa.mappa_workflow.ts_analytics.scoring import POLICY as TS_SCORING_POLICY, subset_score
from mappa.mappa_core.core.search.queries.subgraph import SubgraphMatch, SubgraphSearchResult, match_weighted_subgraph

__all__ = [
    "query_path", "convert_fragments", "FragmentConversion",
    "AAMSearchPlan", "plan_aam_search",
    "aam_record", "aam_from_record", "write_aam_bundle",
    "BondTolerance", "match_fragment", "FragmentMatchConfig", "FragmentMatchContext", "FragmentMatchResult", "FragmentPlacement",
    "AAMSearchGraph", "SearchState", "FragmentTransition", "SearchPath", "SearchBranch",
    "SearchContext", "SearchStop",
    "PathRealization", "finalize_graph_symmetry", "SubgraphSearchResult",
    "AAMProblem", "AAMResult", "AAMSearchConfig",
    "AAMSearchMetrics",
    "AtomAssignment", "CoreAAMBranch",
    "CoreAAMResult", "MolecularEndpoint",
    "ReactionContext", "ResolvedMechanism",
    "TSMechanismResult", "TSResult", "TSScore", "TSScoringConfig",
    "TransitionStateTarget", "VibrationalModes",
    "AAMHierarchy", "AtomBijection", "AtomPermutation",
    "FragmentMatch", "PermutationGroup", "SymmetryDomain",
    "search_aam_checkpoints", "AAMCheckpointResult",
    "search_aam", "search_core_assignments",
    "analyze_transition_state",
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
