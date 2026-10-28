"""Shared branch scheduling and typed hierarchy structures.

Public reaction search belongs to mappa.aam; conversion and event counting to
mappa.postprocessing. This package does not expose another end-to-end pipeline.
"""
from __future__ import annotations

from .sweep import cut_sweep_items
from .branch import (
    _Branch,
    _alignment_state_signature,
    _chemistry_orbit_signature,
    _generate_seed_orders,
    _orbit_pair,
    find_islands,
)

from .post_aam import (
    AAMHierarchy,
    AtomBijection,
    AtomPermutation,
    FragmentMatch,
    PermutationGroup,
    SymmetryDomain,
)
