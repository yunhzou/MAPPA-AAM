"""Shared branch scheduling and typed hierarchy structures.

Public reaction search belongs to mappa.aam; conversion and event counting to
mappa.postprocessing. This package does not expose another end-to-end pipeline.
"""
from __future__ import annotations

from .sweep import cut_sweep_items
from .branch import (
    SYM_REPAIR_MAX_EVALS,
    _Branch,
    _alignment_state_signature,
    _chemistry_orbit_signature,
    _chirality_violations,
    _generate_seed_orders,
    _mapping_change_score,
    _orbit_pair,
    find_islands,
    symmetry_repair_mapping,
)
from .mapping_families import (
    AnalyticalMappingFamily,
    MappingFamilyConflict,
    MappingFamilyError,
    compile_analytical_mapping_family,
    mapping_event_signature,
)

from .post_aam import (
    AAMBranch,
    AAMHierarchy,
    AtomBijection,
    AtomPermutation,
    FragmentMatch,
    PermutationGroup,
    SymmetryDomain,
)
