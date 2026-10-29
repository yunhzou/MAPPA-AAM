"""Fragment matching and growth with compressed symmetry."""
from __future__ import annotations

from mappa.mappa_core.core.matching.api import match_fragment
from mappa.mappa_core.core.matching.canonical import _PartialMappingCanonicalizer
from mappa.mappa_core.core.matching.dedupe import (
    _boundary_signature,
    _dedup_sym_cands,
    _p_relation_signature,
)
from mappa.mappa_core.core.matching.extend import _extend_sym_cands
from mappa.mappa_core.core.matching.frontier import (
    _frontier_boundary_edges,
    _push_edges_from,
    _set_unique,
)
from mappa.mappa_core.core.matching.growth import grow_island
from mappa.mappa_core.core.matching.orbits import (
    _atom_tuple_orbit,
    _cand_canon_signature,
    _dedup_cands_by_orbit,
    _group_nodes_by_signature,
    _nauty_atom_generators,
    _nauty_orbits,
    _orbit_wbo_bucket,
)
from mappa.mappa_core.core.matching.policy import (
    AttributeNodeMatchPolicy,
    BondTolerance,
    CallableNodeMatchPolicy,
    DEFAULT_NODE_POLICY,
    ElementNodeMatchPolicy,
    SYM_SUPPORT_MAX_STATES,
    _edge_wbo,
    _growth_edge_supported,
    _orbit_id,
    _wbo_bucket,
    as_node_match_policy,
)
from mappa.mappa_core.core.matching.records import (
    FragmentMatchConfig,
    FragmentMatchContext,
    FragmentMatchResult,
    FragmentPlacement,
    IslandBranchLimitExceeded,
    _IsoResult,
)
from mappa.mappa_core.core.matching.state import (
    _SymBlock,
    _SymCand,
    _cand_has_open_choice,
    _cand_map,
    _cand_possible_p_atoms,
    _sym_block_assignment_expr,
    _sym_block_indexes,
    _symmetry_state,
)
from mappa.mappa_core.core.matching.support import (
    _r_compatible_with_block,
    _refine_sym_assignments,
    _support_witness_for_value,
)
