"""Fragment matching and growth with compressed symmetry."""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace
from typing import Any
from mappa.mappa_core.core.matching.growth.compression import _dedup_sym_cands, _p_relation_signature
from mappa.mappa_core.core.matching.policy import DEFAULT_NODE_POLICY, as_node_match_policy
from mappa.mappa_core.core.matching.policy import _edge_wbo
from mappa.mappa_core.core.matching.representation.state import _SymCand, _sym_block_indexes
from mappa.mappa_core.core.matching.growth.support import _refine_sym_assignments, _r_compatible_with_block, _support_assignments_for_value


Node = int


Wbo = float


Support = dict[Node, Node]


SymCandidate = _SymCand | Mapping[Node, Node]


OrbitMap = Mapping[Node, int] | None


EdgeKey = tuple[Node, Node]


TargetEntry = tuple[Node, Support, bool]


def _extend_sym_cands(
    cands: Iterable[SymCandidate],
    fragment_old: set[Node],
    n: Node,
    g_R,
    g_P,
    mapping: Mapping[Node, Node],
    iso_tol: float,
    islands_R: Mapping[Node, int] | None,
    p_orbits: OrbitMap = None,
    r_orbits: OrbitMap = None,
    deferred_edges: Iterable[EdgeKey] = (),
    anchor_u: Node | None = None,
    anchor_wbo: Wbo | None = None,
    dedupe_edges: Iterable[EdgeKey] | None = None,
    node_policy=None,
    bond_tolerance=None,
) -> list[_SymCand]:
    """Symmetry-compressed incremental extension.

    This is the elementary one-atom fragment-matching step:

        old fragment F  ->  new fragment F union {n}

    Parameters
    ----------
    cands
        Parallel canonical states for the old fragment.  Each `_SymCand`
        stores one witness plus compact automorphism data; the outer list stores
        multiple canonical-distinct local possibilities for this island.
    fragment_old
        R atom set for the old fragment F.
    n
        New R atom being added to F.
    g_R, g_P
        Reactant/target WBO graphs.  They carry both node labels and a
        full WBO matrix at `graph["wbo_matrix"]`.
    mapping
        Already locked global R->P mapping from other islands.  Free extension
        derives the inverse internally only to prevent reusing locked P atoms;
        merge extension checks the pre-locked image of `n`.
    iso_tol
        WBO tolerance for the selected growth edge `anchor_u -> n`.
        A proposed image `n -> v` is valid only if some supported witness satisfies
        `WBO_P[v,image(u)] >= graph_floor` and
        `abs(WBO_R[n,u] - WBO_P[v,image(u)]) <= iso_tol`.
    islands_R
        Locked-island labels for R atoms.  If `n` is already mapped and belongs
        to an island, extension merges the whole not-yet-fragment part of that
        island, not just one atom.
    p_orbits, r_orbits
        Symmetry orbit maps, usually from `_nauty_orbits(...)` with the
        configured symmetry WBO tolerance.
        They are compression keys only; active R-pair validity is still checked
        against exact WBO values by `iso_tol`.
    deferred_edges
        Previously failed frontier edges kept as one-hop boundary evidence for
        dedupe.  They are not chemistry cuts and are not removed from the WBO
        graph.
    anchor_u, anchor_wbo
        The selected growth edge `anchor_u -> n` and its WBO. Only this edge
        adds a matching constraint; other molecular connections are retained
        for event interpretation.
    dedupe_edges
        Optional boundary-edge set to use for deduping the child states.  The
        growth loop passes current deferred edges plus current one-hop frontier
        edges so symmetry dedupe does not erase distinguishable boundary states.
    node_policy
        Node-level compatibility policy.  The default is same element.  Custom
        policies can admit target nodes by electronic or other descriptors
        before the normal WBO/``iso_tol`` edge verifier is applied.
    Returns
    -------
    list[_SymCand]
        Canonical-distinct compressed matches for `fragment_old union {n}`.
        Each returned `_SymCand` may represent many concrete injective
        assignments through symmetry blocks and exact automorphism domains.
    """
    node_policy = as_node_match_policy(node_policy)
    ctx = _make_extension_context(
        fragment_old, n, g_R, g_P, mapping, iso_tol, islands_R,
        p_orbits, r_orbits, deferred_edges, anchor_u, anchor_wbo,
        dedupe_edges, node_policy, bond_tolerance)
    if ctx is None:
        return []

    children: list[_SymCand] = []
    for raw_cand in cands:
        cand = raw_cand if isinstance(raw_cand, _SymCand) else _SymCand(raw_cand)
        children.extend(_extend_one_candidate(cand, ctx))
    return _dedupe_children(children, ctx)


@dataclass(frozen=True)
class _ExtensionContext:
    """All shared state for one ``F -> F union {n}`` extension.

    Keeping this as an explicit object makes the rest of the file easier to
    audit: helper functions either transform one `_SymCand`, or group/dedupe
    the resulting children.  Nothing here changes matching semantics.
    """

    fragment_old: frozenset[Node]
    n: Node
    n_element: str
    g_R: Any
    g_P: Any
    mapping: Mapping[Node, Node]
    locked_p_atoms: frozenset[Node]
    compatible_p_atoms: tuple[Node, ...]
    iso_tol: float
    islands_R: Mapping[Node, int] | None
    p_orbits: OrbitMap
    r_orbits: OrbitMap
    deferred_edges: tuple[EdgeKey, ...]
    anchor_u: Node | None
    anchor_wbo: Wbo | None
    dedupe_edges: tuple[EdgeKey, ...] | None
    bonded_in_frag: tuple[Node, ...]
    r_wbos: tuple[tuple[Node, Wbo], ...]
    strict_r_wbos: Mapping[Node, Wbo]
    island_atoms: tuple[Node, ...]
    node_policy: Any = DEFAULT_NODE_POLICY
    bond_tolerance: Any = None

    @property
    def is_merge(self) -> bool:
        """True when the new atom already has a locked global image."""
        return self.n in self.mapping

    @property
    def sig_fragment(self) -> set[Node]:
        """Fragment used in dedupe signatures after this extension attempt."""
        return set(self.fragment_old) | {self.n, *self.island_atoms}

    @property
    def boundary_edges(self) -> Iterable[EdgeKey]:
        """Boundary evidence used when deduping child canonical states."""
        return self.deferred_edges if self.dedupe_edges is None else self.dedupe_edges


def _locked_island_atoms(
    n: Node,
    fragment_old: set[Node],
    mapping: Mapping[Node, Node],
    islands_R: Mapping[Node, int] | None,
) -> tuple[Node, ...]:
    """Atoms that must be merged when ``n`` is already locked.

    If ``n`` belongs to an existing island, growth should absorb the whole
    not-yet-fragment part of that island.  Otherwise the merge is just ``n``.
    For a free atom, the only new atom is ``n``.
    """
    if n not in mapping:
        return (n,)
    if islands_R is None or n not in islands_R:
        return (n,)
    iid = islands_R[n]
    return tuple(
        r for r, island_id in islands_R.items()
        if island_id == iid and r not in fragment_old
    )


def _strict_growth_wbos(
    fragment_old: set[Node],
    anchor_u: Node | None,
    anchor_wbo: Wbo | None,
    g_R,
    n: Node,
) -> dict[Node, Wbo]:
    """Exact WBO for the heap edge that triggered this extension.

    This selected edge is the new matching constraint. Other connections
    from the atom to the fragment do not enter the extension test.
    """
    if anchor_u is None or anchor_u not in fragment_old:
        return {}
    return {
        anchor_u: _edge_wbo(g_R, anchor_u, n)
        if anchor_wbo is None else anchor_wbo
    }


def _make_extension_context(
    fragment_old: set[Node],
    n: Node,
    g_R,
    g_P,
    mapping: Mapping[Node, Node],
    iso_tol: float,
    islands_R: Mapping[Node, int] | None,
    p_orbits: OrbitMap,
    r_orbits: OrbitMap,
    deferred_edges: Iterable[EdgeKey],
    anchor_u: Node | None,
    anchor_wbo: Wbo | None,
    dedupe_edges: Iterable[EdgeKey] | None,
    node_policy,
    bond_tolerance=None,
) -> _ExtensionContext | None:
    """Collect repeated extension inputs into one typed context."""
    if anchor_u not in fragment_old or not g_R.has_edge(anchor_u, n):
        raise ValueError("Extension requires its selected active growth edge")
    bonded_in_frag = (anchor_u,)
    locked_p_atoms = frozenset(mapping.values())
    compatible_p_atoms = tuple(
        v for v in g_P.nodes()
        if v not in locked_p_atoms
        and node_policy.compatible(g_R, n, g_P, v)
    )
    return _ExtensionContext(
        fragment_old=frozenset(fragment_old),
        n=n,
        n_element=g_R.nodes[n].get('element'),
        g_R=g_R,
        g_P=g_P,
        mapping=mapping,
        locked_p_atoms=locked_p_atoms,
        compatible_p_atoms=compatible_p_atoms,
        iso_tol=iso_tol,
        islands_R=islands_R,
        p_orbits=p_orbits,
        r_orbits=r_orbits,
        deferred_edges=tuple(tuple(sorted(e)) for e in deferred_edges),
        anchor_u=anchor_u,
        anchor_wbo=anchor_wbo,
        dedupe_edges=(None if dedupe_edges is None
                      else tuple(tuple(sorted(e)) for e in dedupe_edges)),
        bonded_in_frag=bonded_in_frag,
        r_wbos=tuple((u, _edge_wbo(g_R, u, n)) for u in bonded_in_frag),
        strict_r_wbos=_strict_growth_wbos(
            fragment_old, anchor_u, anchor_wbo, g_R, n),
        island_atoms=_locked_island_atoms(n, fragment_old, mapping, islands_R),
        node_policy=node_policy,
        bond_tolerance=bond_tolerance,
    )


def _candidate_covers_fragment(cand: _SymCand, ctx: _ExtensionContext) -> bool:
    """A candidate can extend only if it maps every old fragment atom."""
    return all(u in cand.mapping for u in ctx.bonded_in_frag)


def _supported_value(
    cand: _SymCand,
    ctx: _ExtensionContext,
    v: Node,
    join_block_idx: int | None,
    block_indexes=None,
) -> list[Support]:
    """Return support assignments proving ``ctx.n -> v`` is valid.

    This is where active R-pair WBO validity enters the extension step.  The
    helper searches inside unresolved symmetry blocks instead of trusting the
    stored witness, so an arbitrary witness cannot incorrectly reject a valid
    correlated assignment.
    """
    return _support_assignments_for_value(
        cand,
        ctx.n,
        v,
        list(ctx.bonded_in_frag),
        list(ctx.r_wbos),
        ctx.g_P,
        ctx.iso_tol,
        join_block_idx=join_block_idx,
        block_indexes=block_indexes,
        strict_r_wbos=ctx.strict_r_wbos,
        bond_tolerance=ctx.bond_tolerance,
    )


def _extend_locked_merge(cand: _SymCand, ctx: _ExtensionContext) -> list[_SymCand]:
    """Match the joining edge, then carry the already locked island unchanged."""
    children = [cand]
    fragment = set(ctx.fragment_old)
    for atom in (ctx.n, *(r for r in ctx.island_atoms if r != ctx.n)):
        bonded = ctx.bonded_in_frag if atom == ctx.n else ()
        step = replace(ctx, n=atom, fragment_old=frozenset(fragment),
                       bonded_in_frag=bonded,
                       r_wbos=tuple((r, _edge_wbo(ctx.g_R, r, atom)) for r in bonded),
                       strict_r_wbos=ctx.strict_r_wbos if atom == ctx.n else {})
        target = ctx.mapping[atom]
        extended = []
        for child in children:
            indexes = _sym_block_indexes(child)
            for support in _supported_value(child, step, target,
                                             indexes[1].get(target), indexes):
                refined = _refine_sym_assignments(child, {**support, atom: target})
                if refined is not None:
                    extended.append(refined)
        children = extended
        fragment.add(atom)
    return children


def _target_join_info(
    cand: _SymCand,
    ctx: _ExtensionContext,
    v: Node,
    p_to_block=None,
) -> tuple[int | None, bool]:
    """Return the open block containing ``v`` and whether it can grow freely.

    The compatibility predicate is a compression rule, not the validity rule.
    If it fails, exact support checking may still prove that ``ctx.n -> v`` is
    valid.  In that case the child must refine/fix the assignment instead of
    enlarging the block as a symmetric set-to-set choice.
    """
    if p_to_block is None:
        _, p_to_block = _sym_block_indexes(cand)
    join_idx = p_to_block.get(v)
    if join_idx is None:
        return None, False
    block = cand.blocks[join_idx]
    if not block.open:
        return None, False
    can_extend = _r_compatible_with_block(
        cand, join_idx, ctx.n, set(ctx.fragment_old),
        ctx.g_R, ctx.r_orbits, bond_tolerance=ctx.bond_tolerance)
    return join_idx, can_extend


def _collect_free_target_entries(
    cand: _SymCand,
    ctx: _ExtensionContext,
) -> tuple[dict[int, list[TargetEntry]], dict[Any, list[TargetEntry]]]:
    """Find every target atom that can support free ``ctx.n`` extension.

    Results are split into two groups:

    - ``block_join``: target atom lies inside an existing open symmetry block.
      The child either extends that block or refines it under a correlated
      support assignment.
    - ``by_group``: target atom is outside existing blocks.  Entries are
      grouped by element/orbit/context signature so equivalent target atoms can
      become one new `_SymBlock`.
    """
    block_join: dict[int, list[TargetEntry]] = defaultdict(list)
    by_group: dict[Any, list[TargetEntry]] = defaultdict(list)
    block_indexes = _sym_block_indexes(cand)
    _, p_to_block = block_indexes
    for v in ctx.compatible_p_atoms:
        join_idx, can_extend = _target_join_info(
            cand, ctx, v, p_to_block=p_to_block)
        supports = _supported_value(
            cand, ctx, v, join_idx, block_indexes=block_indexes)
        if join_idx is not None:
            block_join[join_idx].extend((v, support, can_extend) for support in supports)
        elif supports:
            sig = _p_relation_signature(
                cand, v, ctx.g_P, ctx.p_orbits,
                node_policy=ctx.node_policy)
            by_group[sig].extend((v, support, True) for support in supports)
    return block_join, by_group


def _children_from_block_join(
    cand: _SymCand,
    n: Node,
    join_idx: int,
    entries: list[TargetEntry],
) -> list[_SymCand]:
    """Build children where ``n`` joins an existing open symmetry block."""
    children: list[_SymCand] = []

    free_entries = [
        (v, support)
        for v, support, can_extend in entries
        if not support and can_extend
    ]
    whole_pool = {v for v, _ in free_entries} == set(cand.blocks[join_idx].p_atoms)
    if whole_pool:
        # No old-block witness needs to be fixed.  The existing block expands
        # from, for example, "one H occupies methyl-H pool" to "two H atoms
        # occupy methyl-H pool".
        witness_v = sorted(v for v, _ in free_entries)[0]
        child = cand.with_extended_block(join_idx, n)
        if child is not None:
            child = child.with_witness({n: witness_v})
        if child is not None:
            children.append(child)

    for v, support, can_extend in entries:
        if whole_pool and not support and can_extend:
            continue
        # The join is either correlated with an old block assignment, or exact
        # support exists while the source atom is not block-compatible enough
        # to enlarge the block.  Refine/freeze instead of treating it as an
        # independent symmetric choice.
        fixed = dict(support)
        fixed[n] = v
        child = _refine_sym_assignments(cand, fixed)
        if child is not None:
            children.append(child)
    return children


def _children_from_context_group(
    cand: _SymCand,
    n: Node,
    entries: list[TargetEntry],
) -> list[_SymCand]:
    """Build children for target atoms outside existing symmetry blocks."""
    children: list[_SymCand] = []

    by_support = defaultdict(set)
    for v, support, _can_extend in entries:
        by_support[tuple(sorted(support.items()))].add(v)
    for support, values in by_support.items():
        base = _refine_sym_assignments(cand, dict(support)) if support else cand
        if base is None:
            continue
        group = tuple(sorted(values))
        child = (base.with_new_block(n, group, extendable=True) if len(group) > 1
                 else base.with_fixed(n, group[0]))
        if child is not None:
            child = child.with_witness({n: group[0]})
        if child is not None:
            children.append(child)

    return children


def _extend_free_atom(cand: _SymCand, ctx: _ExtensionContext) -> list[_SymCand]:
    """Extend by an unmapped atom using all valid same-element target atoms."""
    block_join, by_group = _collect_free_target_entries(cand, ctx)
    children: list[_SymCand] = []
    for join_idx, entries in sorted(block_join.items()):
        children.extend(_children_from_block_join(cand, ctx.n, join_idx, entries))
    for _, entries in sorted(by_group.items(), key=lambda kv: str(kv[0])):
        children.extend(_children_from_context_group(cand, ctx.n, entries))
    return children


def _extend_one_candidate(cand: _SymCand, ctx: _ExtensionContext) -> list[_SymCand]:
    """Extend one canonical state; return zero or more canonical child states."""
    if not _candidate_covers_fragment(cand, ctx):
        return []
    if ctx.is_merge:
        return _extend_locked_merge(cand, ctx)
    return _extend_free_atom(cand, ctx)


def _dedupe_children(children: list[_SymCand], ctx: _ExtensionContext) -> list[_SymCand]:
    """Compress equivalent states without changing any future extension."""
    return _dedup_sym_cands(
        children,
        ctx.g_R,
        ctx.g_P,
        ctx.r_orbits,
        ctx.p_orbits,
        fragment=ctx.sig_fragment,
        deferred_edges=ctx.boundary_edges,
        locked_mapping=ctx.mapping,
        node_policy=ctx.node_policy,
    )
