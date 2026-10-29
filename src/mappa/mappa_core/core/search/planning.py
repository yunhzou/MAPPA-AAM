"""Direction, seed and cut planning with reproducible checkpoint identity."""

from __future__ import annotations
from dataclasses import dataclass, replace, asdict
from mappa.domain import AAMProblem, AAMSearchConfig, AAMResult
import json, hashlib, random
from pathlib import Path
from collections import Counter


@dataclass(frozen=True)
class AAMSearchPlan:
    """Keep input R/P identities separate from search source/target identities.

    The raw AAM result and its symmetry actions retain search orientation.
    Convert only a concrete mapping after selecting/realizing its search path.
    Do not invert generators or reinterpret a compressed graph in place.
    """
    input_problem: AAMProblem
    problem: AAMProblem
    config: AAMSearchConfig
    reversed: bool

    def to_input_mapping(self, mapping):
        pairs = dict(mapping)
        return {r:p for p,r in pairs.items()} if self.reversed else pairs

    def to_search_mapping(self, mapping):
        return self.to_input_mapping(mapping)

    @property
    def direction(self):
        return 'P_to_R' if self.reversed else 'R_to_P'


def plan_aam_search(problem, config=None):
    """Search smaller explicit-atom endpoint first; retain R→P on ties."""
    config = config or AAMSearchConfig()
    reverse = problem.source_atom_count > problem.target_atom_count
    search_problem = AAMProblem(problem.product, problem.reactant, problem.name) if reverse else problem
    search_config = replace(config, anchors=tuple((p,r) for r,p in config.anchors)) if reverse else config
    return AAMSearchPlan(problem, search_problem, search_config, reverse)






def cut_seed(cut, root_seed=42):
    """Distinct reproducible streams, independent of scheduling and cut order."""
    edges = tuple(sorted(tuple(sorted(edge)) for edge in cut))
    if not edges:
        return root_seed
    payload = json.dumps([root_seed, edges], separators=(',', ':')).encode('ascii')
    return int.from_bytes(hashlib.blake2b(payload, digest_size=16,
                                         person=b'AAM-cut-seeds-v1').digest(), 'big')


def checkpoint_manifest(problem, config):
    """Identity required before reusing cut checkpoints (not reference labels)."""
    def endpoint(value):
        return dict(elements=list(value.elements), wbo=value.wbo.tolist(),
                    coordinates=value.coordinates.tolist())
    settings = asdict(config)
    if config.random_seed == 42:
        settings.pop('random_seed')
    if config.sweep_cuts:
        settings.pop('sweep_cuts')
    if config.seed_selection == 'random':
        settings.pop('seed_selection')  # unchanged policy can resume older cuts
    return json.loads(json.dumps(dict(schema='rxn_core.aam_checkpoints/v5',
        matching_policy='exact_fragment_target_action_v3',
        growth_policy='uniform_boundary_atom_splitmix64_v1',
        branch_cap_policy='completed_fragment_frontier_v1',
        seed_policy='independent_per_cut_blake2b_v2_shared_pool',
        reactant=endpoint(problem.reactant),product=endpoint(problem.product),
        config=settings)))


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


def _strong_edges(wboR, cut_floor):
    n = int(wboR.shape[0])
    return [(i, j) for i in range(n) for j in range(i + 1, n)
            if float(wboR[i, j]) >= cut_floor]


def cut_sweep_items(wboR, cut_floor=0.2, *, elements=None,
                    heavy_only=False):
    """Return independent no-cut plus one-edge cut work items.

    ``heavy_only`` keeps explicit hydrogens in AAM mapping and final scoring;
    it only avoids using X-H bonds as artificial graph-search cuts.
    """
    edges = _strong_edges(wboR, float(cut_floor))
    if heavy_only:
        if elements is None:
            raise ValueError("elements are required for heavy-only cut sweep")
        edges = [
            (i, j) for i, j in edges
            if str(elements[i]) != 'H' and str(elements[j]) != 'H'
        ]
    return [()] + [((int(i), int(j)),) for i, j in edges]
