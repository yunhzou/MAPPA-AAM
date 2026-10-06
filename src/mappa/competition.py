"""Experimental, opt-in fragment competition; disabled in the default pipeline.

An explicit call may add this stage after sweep search and before decoding.
It is not part of the published MAPPA algorithm.

TODO: establish a controlled coverage/cost benefit and improve proposal
scheduling and deduplication before considering default integration.

The local-collision policy is extracted from the validated 140-case campaign.
No dataset, reference, filesystem, environment-variable, or CLI dependencies.
"""

from dataclasses import dataclass, asdict
from .domain import AAMResult
import math
import time
import collections


@dataclass(frozen=True)
class CompetitionConfig:
    operation_budget: int = 128
    seconds: float = 270.0
    parent_limit: int = 8
    depth_limit: int = 2
    queue_limit: int = 512
    dependent_component_limit: int = 8

    def __post_init__(self):
        for name in (
            "operation_budget",
            "parent_limit",
            "depth_limit",
            "queue_limit",
            "dependent_component_limit",
        ):
            if not isinstance(getattr(self, name), int) or getattr(self, name) < 1:
                raise ValueError(f"{name} must be a positive integer")
        if not math.isfinite(self.seconds) or self.seconds < 0:
            raise ValueError("seconds must be finite and nonnegative")


@dataclass
class CompetitionResult:
    baseline: AAMResult
    repairs: tuple[AAMResult, ...]
    counts: dict
    offers: tuple
    pending: int
    elapsed_seconds: float

    def final_catalogue(self):
        from .final_branches import FinalBranchCatalogue

        catalogue = FinalBranchCatalogue(self.baseline.problem).add_aam(
            self.baseline, "baseline"
        )
        for i, repair in enumerate(self.repairs):
            catalogue.add_aam(repair, f"takeover{i}")
        return catalogue


def takeover_plan(old, owners, owner, new, *, allow_b_relocation=False):
    """Exact pairs that survive B's priority; reject malformed proposals."""
    old = dict(old)
    new = dict(new)
    owners = dict(owners)
    if len(set(new.values())) != len(new):
        raise ValueError("noninjective takeover")
    base = {r: p for (r, p) in old.items() if owners[r] == owner}
    if not allow_b_relocation and (
        not all((new.get(r, p) == p for (r, p) in base.items()))
    ):
        raise ValueError("B prefix was changed")
    bmap = {} if allow_b_relocation else dict(base)
    bmap.update(new)
    used = set(bmap.values())
    retained = {r: p for (r, p) in old.items() if r not in bmap and p not in used}
    anchors = dict(retained)
    anchors.update(bmap)
    holes = sorted(set(old) - anchors.keys())
    eaten = sorted((r for r in bmap if owners[r] != owner))
    displaced = sorted((r for r in old if r not in bmap and old[r] in used))
    changed = sorted((r for r in bmap if bmap[r] != old[r]))
    return dict(
        anchors=anchors,
        retained=retained,
        bmap=bmap,
        holes=holes,
        eaten=eaten,
        displaced=displaced,
        changed=changed,
        touched_owners=sorted({owners[r] for r in eaten + displaced}),
    )


def compete_fragments(a, config=None):
    """Retain baseline and return accepted B-priority repair families.

    Requires balanced endpoints and the optional native growth engine. Budgets
    limit calls/depth, not global optimality; the wall budget is soft. Hard user
    anchors remain fixed. No event decoding is performed here. Parent ordering
    uses observed event counts with the AAM config's event thresholds.
    """
    from .domain import AAMSearchConfig, AAMResult, AAMSearchMetrics
    from .frag import build_graph
    from .fragment import match_fragment, FragmentMatchConfig, FragmentMatchContext
    from .matcher import _nauty_orbits
    from .native_search import find_islands_native
    from .search_symmetry import finalize_graph_symmetry
    from .event_patterns import SignedEventIndex
    from .growth.native import available

    if not available():
        raise RuntimeError(
            "Fragment competition requires the optional native engine; see native/README.md"
        )
    config = config or CompetitionConfig()
    problem, cfg = (a.problem, a.config)
    idx = SignedEventIndex(
        problem,
        threshold=cfg.event_threshold,
        metal_threshold=cfg.metal_event_threshold,
    )
    start = time.perf_counter()
    deadline = start + config.seconds
    BUDGET = config.operation_budget
    repairs = []
    terminals = []
    batch = []

    def flush_parents():
        if not batch:
            return
        scores = idx.counts([v for (t, v, part) in batch]).sum(axis=1)
        terminals.extend(
            ((int(score), t, part) for ((t, v, part), score) in zip(batch, scores))
        )
        batch.clear()

    for t in a.graph.terminals:
        st = a.graph.states[t]
        if len(st.mapping) != idx.n:
            continue
        groups = collections.defaultdict(list)
        for atom, g in st.islands:
            groups[g].append(atom)
        if len(groups) < 2:
            continue
        partition = tuple(sorted((tuple(v) for v in groups.values())))
        m = dict(st.mapping)
        batch.append((t, [m[r] for r in range(idx.n)], partition))
        if len(batch) >= 64:
            flush_parents()
    flush_parents()
    terminals.sort()
    selected = []
    partitions = set()
    for score, t, part in terminals:
        if part in partitions:
            continue
        partitions.add(part)
        selected.append(t)
        if len(selected) >= config.parent_limit:
            break
    target = build_graph(
        problem.product.elements, problem.product.wbo, bond_cut=cfg.graph_floor
    )
    po = _nauty_orbits(target, wbo_tol=cfg.iso_tolerance)
    from .conditioned_symmetry import ConditionedSymmetryWorkspace
    from .search_validation import RepresentativeValidationWorkspace

    symmetry_workspace = ConditionedSymmetryWorkspace(target, cfg.iso_tolerance)
    validation_workspace = RepresentativeValidationWorkspace(problem)
    sources = {}
    offers_cache = {}
    completion_cache = {}
    queue = collections.deque()
    seen_states = set()
    offers = []
    counts = collections.Counter()
    release_cache = {}
    for t in selected:
        state = a.graph.states[t]
        context = a.graph.contexts[state.context]
        queue.append(
            dict(
                mapping=dict(state.mapping),
                owners=dict(state.islands),
                cuts=tuple(context.cuts),
                order=tuple(context.seed_order),
                parent=t,
                depth=0,
            )
        )

    def prepared(cuts):
        if cuts not in sources:
            g = build_graph(
                problem.reactant.elements,
                problem.reactant.wbo,
                bond_cut=cfg.graph_floor,
            )
            g.remove_edges_from(cuts)
            sources[cuts] = g
        return sources[cuts]

    while (
        queue
        and counts["growth_calls"] < BUDGET
        and (counts["completion_calls"] < BUDGET)
        and (time.perf_counter() < deadline)
    ):
        state = queue.popleft()
        old = state["mapping"]
        owners = state["owners"]
        cuts = state["cuts"]
        order = state["order"]
        source = prepared(cuts)
        offer_source = prepared(())
        statekey = (tuple(sorted(old.items())), tuple(sorted(owners.items())), cuts)
        if statekey in seen_states:
            continue
        seen_states.add(statekey)
        groups = collections.defaultdict(list)
        for r, g in owners.items():
            groups[g].append(r)
        # Smaller established fragments challenge adjacent owners first.
        tasks = []
        for owner, atoms in sorted(
            groups.items(), key=lambda item: (len(item[1]), min(item[1]))
        ):
            if len(atoms) == idx.n:
                continue
            boundary = [
                r
                for r in atoms
                if any((owners[s] != owner for s in offer_source.neighbors(r)))
            ]
            if not boundary:
                continue
            for contact in sorted(
                {
                    n
                    for r in boundary
                    for n in offer_source.neighbors(r)
                    if owners[n] != owner
                }
            ):
                contact_seed = min(
                    (r for r in boundary if offer_source.has_edge(r, contact))
                )
                tasks.append((owner, atoms, contact_seed, contact))
        for owner, atoms, seed, contact in tasks:
            if (
                counts["growth_calls"] >= BUDGET
                or counts["completion_calls"] >= BUDGET
                or time.perf_counter() >= deadline
            ):
                break
            base = {r: old[r] for r in atoms}
            key = (cuts, tuple(sorted(base.items())), seed, contact)
            if key in offers_cache:
                result = offers_cache[key]
                counts["growth_cache_hits"] += 1
            else:
                counts["growth_calls"] += 1
                proposal_source = offer_source.subgraph(set(atoms) | {contact}).copy()
                result = match_fragment(
                    proposal_source,
                    target,
                    seed=seed,
                    context=FragmentMatchContext({}, {}, (), target_orbits=po),
                    config=FragmentMatchConfig(
                        graph_floor=cfg.graph_floor,
                        iso_tolerance=cfg.iso_tolerance,
                        branch_limit=cfg.branch_limit,
                        allow_mapped_seed=False,
                    ),
                )
                offers_cache[key] = result
            if result.capped:
                counts["growth_capped"] += 1
                continue
            for placement in result.matches:
                plan = takeover_plan(
                    old, owners, owner, dict(placement), allow_b_relocation=True
                )
                counts["offers"] += 1
                if contact not in placement:
                    counts["no_contact_reached"] += 1
                    continue
                if not plan["touched_owners"] or not plan["changed"]:
                    counts["no_changed_takeover"] += 1
                    continue
                if (
                    counts["completion_calls"] >= BUDGET
                    or time.perf_counter() >= deadline
                ):
                    break
                import networkx as nx

                # Release small dependents around changed/displaced atoms.
                affected = set(plan["changed"]) | set(plan["displaced"])
                release = set()
                for side, (graph0, centers, inverse) in enumerate(
                    (
                        (prepared(()), affected, {r: r for r in old}),
                        (
                            target,
                            {old[r] for r in affected}
                            | {plan["bmap"][r] for r in affected if r in plan["bmap"]},
                            {p: r for (r, p) in old.items()},
                        ),
                    )
                ):
                    component_key = (side, frozenset(centers))
                    if component_key not in release_cache:
                        released_nodes = set()
                        for component in nx.connected_components(
                            nx.subgraph_view(
                                graph0, filter_node=lambda n: n not in centers
                            )
                        ):
                            if len(
                                component
                            ) <= config.dependent_component_limit and any(
                                (
                                    v in centers
                                    for n in component
                                    for v in graph0.neighbors(n)
                                )
                            ):
                                released_nodes.update(component)
                        release_cache[component_key] = released_nodes
                    else:
                        counts["dependent_cache_hits"] += 1
                    release.update(
                        (
                            inverse[n]
                            for n in release_cache[component_key]
                            if n in inverse
                        )
                    )
                release -= plan["bmap"].keys()
                plan["released_dependents"] = sorted(release & plan["retained"].keys())
                for r in release:
                    plan["anchors"].pop(r, None)
                    plan["retained"].pop(r, None)
                plan["holes"] = sorted(set(old) - plan["anchors"].keys())
                if any((plan["anchors"].get(r) != p for (r, p) in cfg.anchors)):
                    counts["anchor_conflicts"] += 1
                    continue
                plan["anchors"].update(cfg.anchors)
                anchors = plan["anchors"]
                completion_key = (cuts, tuple(sorted(anchors.items())))
                if completion_key in completion_cache:
                    counts["completion_cache_hits"] += 1
                    continue
                counts["completion_calls"] += 1
                released_edges = ()
                effective_cuts = tuple(sorted(set(cuts) | set(released_edges)))
                completion_source = prepared(effective_cuts)
                completion_order = order
                # Rebuild a fresh family: displaced ancestor constraints must not return.
                graph = find_islands_native(
                    completion_source,
                    target,
                    completion_order,
                    graph_floor=cfg.graph_floor,
                    iso_tol=cfg.iso_tolerance,
                    max_branches=cfg.branch_limit,
                    p_orbits=po,
                    cuts=effective_cuts,
                    anchor_map=anchors,
                )
                graph, _ = finalize_graph_symmetry(
                    graph,
                    target,
                    iso_tolerance=cfg.iso_tolerance,
                    workspace=symmetry_workspace,
                )
                # A shared terminal witness can satisfy all ancestor constraints
                # at once; invalid unions fall back to checking each history.
                validated, invalid = validation_workspace.validate_graph(graph)
                counts["validated_paths"] += validated
                if invalid:
                    counts["rejected_repair_graphs"] += 1
                    counts["invalid_paths"] += len(invalid)
                    completion_cache[completion_key] = True
                    continue
                completion_cache[completion_key] = True
                repair_label = f"takeover{len(repairs)}"
                repaircfg = AAMSearchConfig(
                    **dict(asdict(cfg), anchors=tuple(sorted(anchors.items())))
                )
                out_aam = AAMResult(
                    problem, repaircfg, graph, AAMSearchMetrics.from_record({}, 0)
                )
                repairs.append(out_aam)
                proof = dict(
                    parent=state["parent"],
                    depth=state["depth"] + 1,
                    seed=seed,
                    contact=contact,
                    B_atoms=atoms,
                    B_prefix=sorted(base.items()),
                    B_placement=sorted(dict(placement).items()),
                    cuts=cuts,
                    released_preservation_edges=released_edges,
                    effective_cuts=effective_cuts,
                    eaten=plan["eaten"],
                    displaced=plan["displaced"],
                    holes=plan["holes"],
                    released_dependents=plan.get("released_dependents", []),
                    changed=plan["changed"],
                    retained_A_pairs=sorted(plan["retained"].items()),
                    repair_id=repair_label,
                    capped=graph.capped,
                )
                counts["completion_capped"] += int(graph.capped)
                full = 0
                for t in graph.terminals:
                    m = dict(graph.states[t].mapping)
                    if len(m) != idx.n:
                        continue
                    assert all((m[r] == p for (r, p) in anchors.items()))
                    assert len(set(m.values())) == idx.n
                    assert all(
                        (
                            problem.reactant.elements[r] == problem.product.elements[p]
                            for (r, p) in m.items()
                        )
                    )
                    full += 1
                    counts["full_witnesses"] += 1
                    if (
                        state["depth"] + 1 < config.depth_limit
                        and len(queue) < config.queue_limit
                    ):
                        newowners = {r: owners[r] for r in plan["retained"]}
                        fresh = max(owners.values()) + 1
                        newowners.update({r: fresh for r in plan["bmap"]})
                        for j, r in enumerate(plan["holes"], fresh + 1):
                            newowners[r] = j
                        queue.append(
                            dict(
                                mapping=m,
                                owners=newowners,
                                cuts=effective_cuts,
                                order=order,
                                parent=repair_label + f":{t}",
                                depth=state["depth"] + 1,
                            )
                        )
                offers.append(dict(proof, full_witnesses=full))
    for name, value in symmetry_workspace.stats().items():
        counts[f"symmetry_{name}"] = value
    for name, value in validation_workspace.stats.items():
        counts[f"validation_{name}"] = value
    return CompetitionResult(
        a,
        tuple(repairs),
        dict(counts),
        tuple(offers),
        len(queue),
        time.perf_counter() - start,
    )
