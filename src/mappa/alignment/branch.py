"""Fragment-decision scheduling and exact frontier admission."""
from __future__ import annotations

import random
from collections import Counter


from ..frag import bond_event_threshold
from ..fragment import FragmentMatchConfig, FragmentMatchContext, match_fragment
from ..search_graph import SearchContext, SearchGraphBuilder, frozen_value
from ..matcher import (
    _boundary_signature,
    _edge_wbo,
    _nauty_orbits,
    _orbit_id,
    _wbo_bucket,
    as_node_match_policy,
)


class _StateKey(tuple):
    """Canonical state tuples with a hash reused at every frontier lookup."""
    def __new__(cls, values):
        key = super().__new__(cls, values)
        key._hash = hash(values)
        return key

    def __hash__(self):
        return self._hash


class BranchLimitExceeded(RuntimeError):
    """Raised when an explicitly enumerated result set exceeds its cap."""

    def __init__(self, max_branches, *, seed=None, pass_no=None,
                 branch_count=None, stage=None):
        self.max_branches = int(max_branches)
        self.seed = seed
        self.pass_no = pass_no
        self.branch_count = branch_count
        self.stage = stage
        msg = f"alignment branch cap hit: {branch_count}/{max_branches}"
        if seed is not None:
            msg += f" seed={seed}"
        if pass_no is not None:
            msg += f" pass={pass_no}"
        if stage:
            msg += f" stage={stage}"
        super().__init__(msg)


class _Branch:
    __slots__ = ('mapping', 'islands_R', 'next_iid',
                 'deferred_edges', 'graph', 'node', '_signature_cache')
    def __init__(self, graph=None, anchor_map=()):
        self.mapping = dict(anchor_map)
        self.islands_R = {}
        self.next_iid = 1
        self.deferred_edges = set()
        self._signature_cache = None
        for r, p in sorted(self.mapping.items()):
            self.islands_R[r] = self.next_iid
            self.next_iid += 1
        self.graph = graph or SearchGraphBuilder(SearchContext((), (), ()))
        self.node = self.graph.root(self)

    @classmethod
    def from_anchor_map(cls, anchor_map, graph=None):
        return cls(graph, anchor_map)

    def fork(self):
        b = object.__new__(_Branch)
        b.mapping = dict(self.mapping)
        # Commit replaces these partitions; sibling branches can read-share
        # them until their next fragment decision.
        b.islands_R = self.islands_R
        b.next_iid = self.next_iid
        b.deferred_edges = set(self.deferred_edges)
        b.graph = self.graph
        b.node = self.node
        b._signature_cache = self._signature_cache
        return b

    @property
    def islands_P(self):
        """Derived inspection view; matching only needs the source partition."""
        return {self.mapping[atom]: island for atom, island in self.islands_R.items()}

    @property
    def symmetry_paths(self):
        """Explicit path projection for consumers; never copied during growth."""
        return [path.fragments for path in self.graph.finish().paths(self.node)]

    def merge_exact_paths(self, other):
        """Union histories after exact equality of cumulative search state."""
        self.graph.join(self, other)

    def state_key(self):
        """One canonical snapshot for both frontier dedup and DAG storage.

        Product island labels are derived from the mapping and source labels;
        storing a second copy in the dedup key adds no information.
        """
        if self._signature_cache is None:
            self._signature_cache = _StateKey((tuple(sorted(self.mapping.items())),
                                              tuple(sorted(self.islands_R.items())),
                                              tuple(sorted(self.deferred_edges))))
        return self._signature_cache

    def commit(self, iso, events=None):
        """Consume a fragment placement; its symmetry record is now DAG-owned."""
        before = self.state_key()
        if events is not None:
            touched = {self.islands_R[r] for r in iso if r in self.islands_R}
            iid = min(touched) if touched else self.next_iid
            relabeled = [(int(r), int(self.islands_R[r])) for r in iso
                         if r in self.islands_R and self.islands_R[r] != iid]
            relabeled.extend((int(r), int(k)) for r, k in self.islands_R.items()
                             if r not in iso and k in touched and k != iid)
        committed_new = []
        for r, p in iso.items():
            if r not in self.mapping:
                self.mapping[r] = p
                committed_new.append((int(r), int(p)))
        self.deferred_edges.update(iso.deferred_edges)
        ordered, islands, self.next_iid = self.graph.merged_partition(before[1], frozenset(iso))
        self.islands_R = dict(ordered)
        # A commit only extends the witness. Reuse the parent's immutable atom
        # pairs instead of allocating a complete copy for every history node.
        mapping = tuple(sorted(before[0] + tuple(committed_new))) if committed_new else before[0]
        self._signature_cache = _StateKey((mapping, islands, tuple(sorted(self.deferred_edges))))
        if before != self.state_key():
            record = {
                'island_idx': self.islands_R[next(iter(iso))],
                'fragment': sorted(int(r) for r in iso.fragment),
                'deferred_edges': self.graph.bond_record(iso.deferred_edges),
                'symmetry': iso.symmetry,
            }
            self.node = self.graph.commit(self.node, self, record, iso.preserved_bonds)
        if events is not None:
            events.append({
                'type': 'island_locked',
                'island_idx': int(iid),
                'pairs': committed_new,
                'merged_with': sorted(int(t) for t in touched - {iid}),
                'relabeled': relabeled,
                'mapped_total': len(self.mapping),
            })


def _freeze_branch_value(value):
    return frozen_value(value)


def _mapping_variation_blocks(mappings, source='interbranch'):
    """Compressed R/P components where equivalent mappings differ."""
    mappings = [
        {int(r): int(p) for r, p in dict(mapping or {}).items()}
        for mapping in mappings
    ]
    mappings = [mapping for mapping in mappings if mapping]
    if len(mappings) <= 1:
        return []

    graph = {}
    for r in sorted({r for mapping in mappings for r in mapping}):
        possible = {
            mapping[r]
            for mapping in mappings
            if r in mapping
        }
        if len(possible) <= 1:
            continue
        r_node = ('r', int(r))
        graph.setdefault(r_node, set())
        for p in possible:
            p_node = ('p', int(p))
            graph.setdefault(p_node, set()).add(r_node)
            graph[r_node].add(p_node)

    blocks = []
    seen = set()
    for start in list(graph):
        if start in seen:
            continue
        stack = [start]
        comp = set()
        seen.add(start)
        while stack:
            node = stack.pop()
            comp.add(node)
            for nb in graph.get(node, ()):
                if nb not in seen:
                    seen.add(nb)
                    stack.append(nb)
        r_atoms = sorted(v for tag, v in comp if tag == 'r')
        p_atoms = sorted(v for tag, v in comp if tag == 'p')
        if not r_atoms or len(p_atoms) <= 1:
            continue
        assignments = {
            tuple(
                (r, mapping.get(r))
                for r in r_atoms
                if r in mapping
            )
            for mapping in mappings
        }
        assignments = {item for item in assignments if item}
        if len(assignments) <= 1:
            continue
        blocks.append({
            'source': source,
            'r_atoms': r_atoms,
            'p_atoms': p_atoms,
            'extendable': False,
            'open': len(r_atoms) < len(p_atoms),
            'assignments': len(assignments),
        })
    return blocks


def _orbit_pair(a, b, orbits):
    oa, ob = _orbit_id(orbits, a), _orbit_id(orbits, b)
    return (oa, ob) if oa <= ob else (ob, oa)


def _chemistry_orbit_signature(mapping, g_R, g_P, r_orbits=None, p_orbits=None,
                               dwbo_threshold=0.5,
                               metal_dwbo_threshold=None):
    """Broken/formed bond signature in joint R-orbit/P-orbit space."""
    br_pairs = []
    fm_pairs = []
    mapped = sorted(mapping)
    elements_R = [g_R.nodes[x].get('element') for x in range(len(g_R.nodes))]
    for i, u in enumerate(mapped):
        for v in mapped[i + 1:]:
            pu, pv = mapping[u], mapping[v]
            threshold = bond_event_threshold(
                elements_R, u, v,
                default_threshold=dwbo_threshold,
                metal_threshold=metal_dwbo_threshold)
            wR = _edge_wbo(g_R, u, v)
            wP = _edge_wbo(g_P, pu, pv)
            if wR - wP >= threshold:
                br_pairs.append((
                    _orbit_pair(u, v, r_orbits),
                    _orbit_pair(pu, pv, p_orbits),
                ))
            elif wP - wR >= threshold:
                fm_pairs.append((
                    _orbit_pair(u, v, r_orbits),
                    _orbit_pair(pu, pv, p_orbits),
                ))
    return tuple(sorted(br_pairs)), tuple(sorted(fm_pairs))


def _alignment_state_signature(mapping, deferred_edges, g_R, g_P,
                               r_orbits=None, p_orbits=None, core_R=(),
                               node_policy=None):
    core_R = tuple(sorted(core_R or ()))
    core_set = set(core_R)
    fixed = []
    for r, p in sorted(mapping.items()):
        if r in core_set:
            fixed.append((('exact', r), ('exact', p)))
        else:
            fixed.append((
                ('orbit', _orbit_id(r_orbits, r)),
                ('orbit', _orbit_id(p_orbits, p)),
            ))
    internal_pairs = []
    mapped = sorted(mapping)
    for i, r1 in enumerate(mapped):
        for r2 in mapped[i + 1:]:
            p1, p2 = mapping[r1], mapping[r2]
            internal_pairs.append((
                _orbit_pair(r1, r2, r_orbits),
                _orbit_pair(p1, p2, p_orbits),
                _wbo_bucket(_edge_wbo(g_R, r1, r2)),
                _wbo_bucket(_edge_wbo(g_P, p1, p2)),
            ))
    boundary = _boundary_signature(
        mapping, g_R, g_P, fragment=set(mapping),
        deferred_edges=deferred_edges, r_orbits=r_orbits,
        p_orbits=p_orbits, locked_mapping=mapping,
        node_policy=node_policy)
    return tuple(fixed), tuple(sorted(internal_pairs)), boundary


def _normalize_anchor_map(anchor_map, g_R, g_P):
    if not anchor_map:
        return {}
    out = {}
    used_p = set()
    for raw_r, raw_p in dict(anchor_map).items():
        r, p = int(raw_r), int(raw_p)
        if r not in g_R:
            raise ValueError(f"anchor R atom is not in query graph: {r}")
        if p not in g_P:
            raise ValueError(f"anchor P atom is not in target graph: {p}")
        if r in out and out[r] != p:
            raise ValueError(f"conflicting anchors for R atom {r}")
        if p in used_p:
            raise ValueError(f"anchor target atom is used more than once: {p}")
        out[r] = p
        used_p.add(p)
    return out


def _anchor_nodes_compatible(anchor_map, g_R, g_P, node_policy):
    return all(
        node_policy.compatible(g_R, r, g_P, p)
        for r, p in anchor_map.items()
    )


def find_islands(g_R, g_P, seed_order,
                 graph_floor=0.2, iso_tol=1.0,
                 dwbo_threshold=0.5, symmetry_wbo_tol=0.2,
                 metal_dwbo_threshold=None,
                 max_branches=1_000_000, events=None,
                 orbit_dedup=True, core_R=None,
                 stop_when_core_mapped=False,
                 p_orbits=None, r_orbits=None,
                 node_policy=None,
                 anchor_map=None,
                 profile=None, cuts=(), growth_replay=None, random_seed=42,
                 bond_tolerance=None):
    """Run growth over a single seed ordering, branching on
    non-set-unique locks. Returns an AAMSearchGraph of fragment decisions.

    Optional `events` only records the FIRST (best-mapped) branch's
    trajectory — multi-branch traces would be confusing on a slider.

    orbit_dedup: when True (default), uses exact automorphism orbits on a
    tolerance-bucketed WBO-colored graph via pynauty as the compression key
    for orbit-equivalence. Callers that run many seed orders against the same
    graph can pass precomputed `p_orbits` / `r_orbits`; missing maps are
    computed here. Chemistry signatures and active R-pair extension checks
    remain the verifier.

    core_R: optional R atoms that define the scoring-relevant alignment.
    With stop_when_core_mapped=True the search
    returns once all live branches have mapped the core; spectators remain
    represented only by the compressed witness state already discovered.

    anchor_map: optional exact R->P constraints.  Anchors are preloaded as
    locked single-atom islands, but their R atoms may still be used as growth
    seeds so the surrounding subgraph can be discovered from the constraint.

    bond_tolerance: optional conditional bond rule forwarded unchanged to the
    shared fragment matcher, including during anchored internal rematching.
    """
    # Orbit equivalence and edge compatibility must describe the same
    # tolerance model.  Keep symmetry_wbo_tol in the signature for callers
    # written against older releases, but do not permit it to diverge.
    symmetry_wbo_tol = float(iso_tol)
    node_policy = as_node_match_policy(node_policy)
    anchor_map = _normalize_anchor_map(anchor_map, g_R, g_P)
    seed_order = list(dict.fromkeys(seed_order))
    graph = SearchGraphBuilder(SearchContext(
        tuple(sorted(g_R)), tuple(sorted(g_P)), tuple(seed_order),
        tuple(cuts), tuple(sorted(core_R or ())), tuple(sorted(anchor_map.items())),
        float(graph_floor), float(iso_tol), int(max_branches),
        objective='requested_core' if stop_when_core_mapped and core_R else 'full_source',
        random_seed=random_seed))
    if not _anchor_nodes_compatible(anchor_map, g_R, g_P, node_policy):
        branch = _Branch.from_anchor_map({}, graph)
        graph.stop(branch, 'incompatible_anchors')
        return graph.finish()
    if orbit_dedup:
        if p_orbits is None:
            p_orbits = _nauty_orbits(
                g_P, wbo_tol=symmetry_wbo_tol,
                node_policy=node_policy)
        if r_orbits is None:
            r_orbits = _nauty_orbits(
                g_R, wbo_tol=symmetry_wbo_tol,
                node_policy=node_policy)
    else:
        p_orbits = None
        r_orbits = None
    # Completed hierarchy/coset equivalence is richer than the live
    # mapping/island state.  Do not quotient distinct concrete live states by
    # endpoint automorphism here: two such states can carry different exact
    # fragment-assignment domains and therefore lead to different mechanisms.
    # Exact semantic family dedupe is performed after AAM has the complete
    # hierarchy available.
    core_R = tuple(sorted(set(core_R or ())))
    branches = [_Branch.from_anchor_map(anchor_map, graph)]
    progressed = True
    pass_no = 0

    def _core_complete(mapping):
        return core_R and all(r in mapping for r in core_R)

    def _mapped_anchor_can_seed(branch, seed):
        if seed not in anchor_map:
            return False
        seed_iid = branch.islands_R.get(seed)
        for nb in g_R.neighbors(seed):
            if g_R[seed][nb].get('wbo', 0.0) < graph_floor:
                continue
            if branch.islands_R.get(nb) != seed_iid:
                return True
        return False

    def _branch_signature(branch):
        return branch.state_key()

    def _merge_equivalent_paths(kept, other):
        if kept.mapping == other.mapping:
            kept.merge_exact_paths(other)
            return
        raise RuntimeError("nonidentical branches shared an exact state key")

    while progressed:
        progressed = False
        pass_no += 1
        if events is not None and branches:
            events.append({'type': 'pass_start', 'pass': pass_no,
                           'mapped': len(branches[0].mapping)})
        for seed_position, seed in enumerate(seed_order):
            graph.step = (pass_no, seed_position)
            graph.seed = int(seed)
            new_branches = []
            pending_seen = {}

            def _admit_subtree(subtree, *, made_progress=False,
                               source_branch=None):
                """Admit one parent branch's result atomically.

                A subtree that would take the live post-dedupe leaf count over
                ``max_branches`` is discarded by itself.  Sibling branches and
                other seed orders remain valid; exactly ``max_branches`` live
                leaves is allowed.
                """
                nonlocal progressed
                additions = []
                local_seen = {}
                for candidate in subtree:
                    sig = _branch_signature(candidate)
                    if sig in pending_seen:
                        _merge_equivalent_paths(
                            pending_seen[sig], candidate)
                        continue
                    if sig in local_seen:
                        _merge_equivalent_paths(local_seen[sig], candidate)
                        continue
                    local_seen[sig] = candidate
                    additions.append((sig, candidate))
                if len(new_branches) + len(additions) > max_branches:
                    for _sig, candidate in additions:
                        graph.stop(candidate, 'capped', stage='combined_live_leaves',
                                   count=len(new_branches) + len(additions), limit=max_branches)
                    if profile is not None:
                        profile.append({
                            'seed': int(seed),
                            'pass': int(pass_no),
                            'branch_index': (
                                None if source_branch is None
                                else int(source_branch)),
                            'mapped_before': (
                                0 if source_branch is None
                                else len(branches[source_branch].mapping)),
                            'result': 'subtree_branch_cap',
                            'stage': 'combined_live_leaves',
                            'subtree_branches': len(additions),
                            'live_branches_before': len(new_branches),
                            'max_branches': int(max_branches),
                            'elapsed_sec': 0.0,
                            'extend_elapsed_sec': 0.0,
                        })
                    return False
                for sig, candidate in additions:
                    pending_seen[sig] = candidate
                    new_branches.append(candidate)
                if additions and made_progress:
                    progressed = True
                return bool(additions)

            for bi, b in enumerate(branches):
                if stop_when_core_mapped and _core_complete(b.mapping):
                    _admit_subtree([b], source_branch=bi)
                    continue
                seed_is_mapped_anchor = seed in b.mapping and seed in anchor_map
                if seed in b.mapping and not seed_is_mapped_anchor:
                    _admit_subtree([b], source_branch=bi)
                    continue
                if seed_is_mapped_anchor and not _mapped_anchor_can_seed(b, seed):
                    _admit_subtree([b], source_branch=bi)
                    continue
                # Only record events for branch 0 to keep trace linear
                ev_arg = events if (events is not None and bi == 0) else None
                result = match_fragment(
                        g_R, g_P, seed=seed,
                        context=FragmentMatchContext(
                            b.mapping, b.islands_R, tuple(b.deferred_edges), r_orbits, p_orbits,
                            growth_replay),
                        config=FragmentMatchConfig(
                            graph_floor=graph_floor, iso_tolerance=iso_tol,
                            branch_limit=max_branches, node_policy=node_policy,
                            orbit_dedup=orbit_dedup, random_seed=random_seed,
                            allow_mapped_seed=seed_is_mapped_anchor,
                            bond_tolerance=bond_tolerance),
                        events=ev_arg,
                        profile=profile,
                        profile_context={
                            'pass': int(pass_no),
                            'branch_index': int(bi),
                            'mapped_before': len(b.mapping),
                        })
                if result.capped:
                    graph.stop(b, 'capped', stage='fragment_growth',
                               count=result.branch_count, limit=result.branch_limit)
                    # Only this parent branch's descendant subtree is
                    # pathological.  Remove it and continue its siblings.
                    continue
                isos = result.matches
                if not isos:
                    _admit_subtree([b], source_branch=bi)
                    continue
                # Fragment candidates were already quotiented by an exact
                # pynauty transporter inside ``grow_island``.  Distinct
                # surviving candidates are distinct symbolic branch choices;
                # do not collapse them by a mechanism/event summary here.
                subtree = []
                subtree_progressed = False
                for ii, iso in enumerate(isos):
                    before_state = _branch_signature(b)
                    b2 = b.fork()
                    b2.commit(iso,
                              events=events if (bi == 0 and ii == 0) else None)
                    after_state = _branch_signature(b2)
                    if after_state == before_state:
                        subtree.append(b)
                    else:
                        subtree.append(b2)
                        subtree_progressed = True
                _admit_subtree(
                    subtree, made_progress=subtree_progressed,
                    source_branch=bi)
            new_branches.sort(key=lambda b: -len(b.mapping))
            # Admission already enforces unique live state keys.
            branches = new_branches
            # Soft warning: pathological symmetry can blow this up. Default
            # cap is 1e6 so we should never hit it on real molecules; if we
            # do, surface it so we know.
            if len(branches) >= 10_000:
                import sys
                print(f"  [warn] alignment branch count = {len(branches)}  "
                      f"max_branches={max_branches}  "
                      f"new_branches_in={len(new_branches)}",
                      file=sys.stderr, flush=True)
            if (stop_when_core_mapped and core_R and branches and
                    all(_core_complete(b.mapping) for b in branches)):
                if events is not None:
                    events.append({'type': 'done',
                                   'mapped': len(branches[0].mapping),
                                   'stop_reason': 'core_mapped',
                                   'core_size': len(core_R)})
                for branch in branches:
                    graph.stop(branch, 'objective_met')
                return graph.finish()
    if events is not None and branches:
        events.append({'type': 'done',
                       'mapped': len(branches[0].mapping)})
    for branch in branches:
        graph.stop(branch, 'objective_met' if len(branch.mapping) == len(g_R) else 'stalled')
    return graph.finish()


def _ordered_seed_nodes(g_R, rng, common_element_threshold=3):
    nodes = list(g_R.nodes())
    element_counts = Counter(g_R.nodes[n].get('element') for n in nodes)
    threshold = max(1, int(common_element_threshold))

    def ambiguous_isolated(n):
        element = g_R.nodes[n].get('element')
        return int(g_R.degree[n]) == 0 and element_counts[element] > threshold

    contextual = [n for n in nodes if not ambiguous_isolated(n)]
    ambiguous = [n for n in nodes if ambiguous_isolated(n)]
    rng.shuffle(contextual)
    rng.shuffle(ambiguous)
    return contextual + ambiguous


def _generate_seed_orders(g_R, n_trials, rng_seed=42,
                          common_element_threshold=1, seed_selection='random'):
    """Generate at most `n_trials` deterministic seed orderings.

    Seeds are deterministic random orderings, not chemistry/core-prioritized
    rankings.  The only special case is an isolated atom whose element occurs
    multiple times: that atom has no graph context, so it is retained but kept
    behind contextual atoms and is not used as an anchor unless no contextual
    anchor exists.

    ``distance`` changes only the trial anchors: sample without replacement
    with weight 1 + distance to the nearest previously selected anchor. An
    unreachable atom has distance len(graph). Selection is adaptive to
    anchors, never to search results, and adds O(n_trials * (vertices +
    edges)) preprocessing.

    All trials draw from ONE shuffled pool. Trial ``idx`` uses the pool
    rotated to start after position ``idx`` (minus its anchor), so the seed
    that follows the first locked island differs across trials instead of
    being an independent re-shuffle that can repeat the same atom. With
    ``n_trials`` <= atoms, every atom appears exactly once at each position.
    Orders depend only on ``idx`` and ``rng_seed``, so prefixes are stable
    when ``n_trials`` grows.
    """
    if seed_selection not in ('random', 'distance'):
        raise ValueError('unknown seed selection policy')
    nodes = _ordered_seed_nodes(
        g_R,
        random.Random(rng_seed),
        common_element_threshold=common_element_threshold,
    )
    element_counts = Counter(g_R.nodes[n].get('element') for n in g_R.nodes())
    threshold = max(1, int(common_element_threshold))

    def ambiguous_isolated(n):
        element = g_R.nodes[n].get('element')
        return int(g_R.degree[n]) == 0 and element_counts[element] > threshold

    anchors = [n for n in nodes if not ambiguous_isolated(n)] or nodes
    orders = []
    n_trials = max(0, int(n_trials))
    if not anchors:
        return orders
    remaining = []
    nearest = {}
    selection_rng = random.Random(rng_seed)
    for idx in range(n_trials):
        if seed_selection == 'random':
            anchor = anchors[idx % len(anchors)]
        else:
            if not remaining:
                remaining = list(anchors)
                nearest = dict.fromkeys(anchors, len(nodes))
            anchor = (anchors[0] if idx == 0 else selection_rng.choices(
                remaining, weights=[1 + nearest[n] for n in remaining], k=1)[0])
            remaining.remove(anchor)
            # Unweighted graph distance; disconnected regions retain the largest
            # weight. No matching, automorphism expansion or geometry is needed.
            distances = {anchor: 0}
            frontier = [anchor]
            for node in frontier:
                for neighbor in g_R.neighbors(node):
                    if neighbor not in distances:
                        distances[neighbor] = distances[node] + 1
                        frontier.append(neighbor)
            for node in remaining:
                nearest[node] = min(nearest[node], distances.get(node, len(nodes)))
        # Shared pool, rotated per trial: no independent re-shuffle, so the
        # next seed after the anchor is drawn without replacement across trials.
        shift = (idx + 1) % len(nodes)
        rotated = [n for n in nodes[shift:] + nodes[:shift] if n != anchor]
        # Isolated ambiguous atoms stay behind contextual atoms.
        rest = ([n for n in rotated if not ambiguous_isolated(n)]
                + [n for n in rotated if ambiguous_isolated(n)])
        orders.append([anchor] + rest)
    return orders
