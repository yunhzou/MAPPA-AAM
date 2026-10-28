"""Composable atom mapping and signed-event decoding.

The main workflow is search_aam -> mappa.postprocessing.decode_events.
TS analysis consumes an explicitly supplied ReactionContext.
"""
from .aam import search_aam, search_aam_checkpoints, AAMCheckpointResult
from .search_orientation import AAMSearchPlan, plan_aam_search
from .fragment_conversion import convert_fragments, FragmentConversion
from .family_query import query_path
from .analytical import compile_mapping_families
from .matcher.bond_policy import BondTolerance
from .fragment import match_fragment, FragmentMatchConfig, FragmentMatchContext, FragmentMatchResult, FragmentPlacement
from .search_graph import AAMSearchGraph, SearchContext, SearchState, SearchStop, FragmentTransition, SearchPath, SearchBranch, PathRealization
from .search_symmetry import finalize_graph_symmetry
from .core_aam import search_core_assignments
from .domain import (
    AAMProblem,
    AAMResult,
    MappingFamilyResult,
    AAMSearchConfig,
    AAMSearchMetrics,
    AnalyticalBranch,
    AtomAssignment,
    CoreAAMBranch,
    CoreAAMResult,
    MolecularEndpoint,
    ReactionContext,
    ResolvedMechanism,
    TSMechanismResult,
    TSResult,
    TSScore,
    TSScoringConfig,
    TransitionStateTarget,
    VibrationalModes,
)
from .ts import analyze_transition_state
from .artifacts import (
    aam_record, aam_from_record, write_aam_bundle,
    reaction_from_record, reaction_record, ts_record,
    write_ts_record,
)
from .alignment.post_aam import (
    AAMBranch, AAMHierarchy, AtomBijection, AtomPermutation,
    FragmentMatch, PermutationGroup, SymmetryDomain,
)

# Stable graph and chemistry primitives are reusable public components.
from .frag import (
    WeightedGraph,
    WeightedNode,
    build_graph,
    build_weighted_graph,
    classify_bonds,
)
from .scoring import POLICY as TS_SCORING_POLICY, subset_score
from .subgraph import SubgraphMatch, SubgraphSearchResult, match_weighted_subgraph

__all__ = [
    "query_path", "convert_fragments", "FragmentConversion",
    "AAMSearchPlan", "plan_aam_search",
    "aam_record", "aam_from_record", "write_aam_bundle",
    "MappingFamilyResult",
    "BondTolerance", "match_fragment", "FragmentMatchConfig", "FragmentMatchContext", "FragmentMatchResult", "FragmentPlacement",
    "AAMSearchGraph", "SearchState", "FragmentTransition", "SearchPath", "SearchBranch",
    "SearchContext", "SearchStop",
    "PathRealization", "finalize_graph_symmetry", "SubgraphSearchResult",
    "AAMProblem", "AAMResult", "AAMSearchConfig",
    "AAMSearchMetrics", "AnalyticalBranch",
    "AtomAssignment", "CoreAAMBranch",
    "CoreAAMResult", "MolecularEndpoint",
    "ReactionContext", "ResolvedMechanism",
    "TSMechanismResult", "TSResult", "TSScore", "TSScoringConfig",
    "TransitionStateTarget", "VibrationalModes",
    "AAMBranch", "AAMHierarchy", "AtomBijection", "AtomPermutation",
    "FragmentMatch", "PermutationGroup", "SymmetryDomain",
    "search_aam_checkpoints", "AAMCheckpointResult",
    "search_aam", "compile_mapping_families", "search_core_assignments",
    "analyze_transition_state",
    "ts_record", "write_ts_record",
    "reaction_record", "reaction_from_record",
    "WeightedGraph", "WeightedNode", "build_graph", "build_weighted_graph",
    "classify_bonds", "TS_SCORING_POLICY", "subset_score",
    "SubgraphMatch", "match_weighted_subgraph",
]

__version__ = "0.2.0"

# Raw search helpers; signed-event decoding remains in mappa.postprocessing.
from .search_orientation import DirectedAAMResult, search_aam_directions
__all__ += ["DirectedAAMResult", "search_aam_directions"]
