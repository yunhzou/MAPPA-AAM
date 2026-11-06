"""Shared seed-policy search and its condition/result records."""

from __future__ import annotations
from collections import Counter, deque
from dataclasses import replace, dataclass
import time
from mappa.mappa_core.core.search.planning import cut_seed, growth_seed, _generate_seed_orders
from mappa.mappa_core.core.search.branching.construction import _Branch
from mappa.mappa_core.domain import AAMResult, AAMSearchMetrics
from mappa.mappa_core.core.matching.api import match_fragment
from mappa.mappa_core.core.matching.representation.records import FragmentMatchConfig, FragmentMatchContext
from mappa.mappa_core.core.search.graph import (
    SearchContext,
    SearchGraphBuilder,
    FragmentTransition,
    SearchStop,
    prune_frontier,
)
from mappa.mappa_core.core.matching.graphs import build_graph


class _SharedGraphBuilder(SearchGraphBuilder):
    """Hash-cons exact states and equal fragment edges, preserving correlations."""
    def __init__(self, context):
        super().__init__(context)
        self.nodes = {}
        self.edges = {}
        self.stop_keys = set()
        self.last_transition = None

    def state(self, branch):
        key = branch.state_key()
        if key not in self.nodes:
            self.nodes[key] = super().state(branch)
        return self.nodes[key]

    def commit(self, parent, branch, match, island, placement_generators=()):
        node = self.state(branch)
        # Index by endpoints, then compare complete records exactly. Recursively
        # freezing every record duplicated the largest payloads just for hashing.
        # The DAG already owns them; no approximate digest or dropped field is needed.
        siblings = self.edges.setdefault((parent, node), [])
        edge = next((i for i in siblings
                     if self.transitions[i].match.key == match.key
                     and self.transitions[i].island == island
                     and self.transitions[i].placement_generators == placement_generators), None)
        if edge is None:
            edge = len(self.transitions)
            siblings.append(edge)
            self.transitions.append(FragmentTransition(len(self.transitions), parent,
                node, self.seed, self.step, match, island, placement_generators))
        self.last_transition = edge
        return node

    def stop(self, branch, reason, *, stage='', count=0, limit=0):
        key = (branch.node, reason, stage, count, limit)
        if key not in self.stop_keys:
            self.stop_keys.add(key)
            super().stop(branch, reason, stage=stage, count=count, limit=limit)


class SharedSeedSearch:
    """One cut's original seed policies sharing an exact decision DAG.

The future of a seed decision depends on mapping, island partition, deferred
edges, the cut condition and the policy's random stream. Distinct policies
share retained graph structure and fragment knowledge, not growth decisions.
Each policy retains its own decision edges and prunes its hierarchy independently.
Identical growth output is never expanded into bijections.
"""
    def __init__(self, problem, config, *, condition, profile=None):
        if config.anchors:
            raise ValueError('Shared policy search currently accepts unanchored AAM only')
        self.problem, self.config, self.condition = problem, config, condition
        self.profile = profile
        self.source, self.target = condition.source, condition.target
        self.source_elements = {atom:data['element'] for atom,data in self.source.nodes(data=True)}
        self.target_elements = {atom:data['element'] for atom,data in self.target.nodes(data=True)}
        self.target_counts = Counter(self.target_elements.values())
        self.available_elements = {}
        self.orders = _generate_seed_orders(self.source, config.seed_count,
            rng_seed=cut_seed(condition.cuts, config.random_seed), seed_selection=config.seed_selection)
        # Multiple policies use this context; edge.step records policy/pass/position.
        self.builder = _SharedGraphBuilder(SearchContext(tuple(sorted(self.source)),
            tuple(sorted(self.target)), (), cuts=condition.cuts, graph_floor=config.graph_floor,
            iso_tolerance=config.iso_tolerance, branch_limit=config.branch_limit, random_seed=config.random_seed))
        self.root = _Branch(self.builder)
        self.branches = {self.root.state_key():self.root}
        self.decisions = {}
        self.work = self.growth_calls = self.reused_states = 0
        self.max_live_branches = 0
        self.elapsed = 0.
        self.frontiers = {}
        self.policy_edges = {}
        self.agenda = deque()
        for index, order in enumerate(self.orders):
            policy = self._policy(index, order)
            request = next(policy, None)
            if request is not None:
                self.agenda.append((policy, request))
        if not self.agenda:
            self._finish_policies()

    def _finish_policies(self):
        # Memoized decisions may serve another active policy. Once all policies
        # finish, only returned-path ancestry and compact stop evidence remain.
        retained = set().union(*self.policy_edges.values())
        self.builder.transitions = [replace(edge, id=index) for index, edge in enumerate(
            edge for edge in self.builder.transitions if edge.id in retained)]
        self.builder.retain((self.root,))
        self.branches.clear()
        self.decisions.clear()
        self.available_elements.clear()
        self.builder.nodes.clear()
        self.builder.edges.clear()
        self.builder.stop_keys.clear()
        self.policy_edges.clear()

    def _expand(self, branch, seed, policy):
        key = (branch.node, seed, policy)
        if key in self.decisions:
            self.reused_states += 1
            return self.decisions[key]
        self.growth_calls += 1
        result = match_fragment(self.source, self.target, seed=seed,
            context=FragmentMatchContext(branch.mapping, branch.islands_R,
                tuple(branch.deferred_edges),
                fragment_cache=self.condition.fragment_cache),
            config=FragmentMatchConfig(graph_floor=self.config.graph_floor,
                iso_tolerance=self.config.iso_tolerance, branch_limit=self.config.branch_limit,
                random_seed=growth_seed(self.condition.cuts, policy, self.config.random_seed)),
            profile=self.profile)
        edges = []
        if result.capped:
            self.builder.stop(branch, 'capped', stage='fragment_growth',
                              count=result.branch_count, limit=result.branch_limit)
            children, changed = (), False
        elif not result.matches:
            children, changed = (branch,), False
        else:
            children, changed = [], False
            for match in result.matches:
                child = branch.fork()
                self.builder.last_transition = None
                child.commit(match, deferred_edges=result.deferred_edges,
                             placement_generators=result.placement_generators)
                if self.builder.last_transition is not None:
                    edges.append(self.builder.last_transition)
                changed |= child.state_key() != branch.state_key()
                children.append(self.branches.setdefault(child.state_key(), child))
            children = tuple(children)
        self.decisions[key] = (children, changed, tuple(edges))
        return self.decisions[key]

    def _prune_policy(self, index, branches):
        # Local edge order keeps tie decisions independent of which other
        # policy first computed a shared transition.
        edges = tuple(self.policy_edges[index])
        transitions = tuple(replace(self.builder.transitions[edge], id=position)
                            for position, edge in enumerate(edges))
        kept, retained, cuts = prune_frontier(self.builder.states, transitions,
            (self.root.node,), [branch.node for branch in branches], self.config.branch_limit)
        self.policy_edges[index] = dict.fromkeys(edges[edge] for edge in sorted(retained))
        for _node, _edge, count, _size in cuts:
            self.builder._record_stop(SearchStop(None, 'capped', self.builder.seed,
                'combined_live_leaves', count, self.config.branch_limit,
                self.builder.step, context=0))
        return [branch for branch in branches if branch.node in kept]

    def _can_seed(self, branch, seed):
        if seed in branch.mapping:
            return False
        if branch.node not in self.available_elements:
            remaining = self.target_counts.copy()
            remaining.subtract(self.target_elements[p] for p in branch.mapping.values())
            self.available_elements[branch.node] = frozenset(e for e,n in remaining.items() if n>0)
        return self.source_elements[seed] in self.available_elements[branch.node]

    def _policy(self, index, order):
        """Complete each seed step, then prune its policy-local hierarchy."""
        branches, progressed, pass_no = [self.root], True, 0
        cap = self.config.branch_limit
        self.policy_edges[index] = {}
        while progressed:
            progressed = False
            pass_no += 1
            for position, seed in enumerate(order):
                # Equal elements and injectivity are mandatory. If this seed has
                # no remaining target capacity, original growth can only carry
                # the branch unchanged. Keep its place in frontier admission.
                if not any(self._can_seed(branch,seed) for branch in branches) and len(branches) <= cap:
                    continue
                admitted, seen, advanced = [], set(), set()
                for branch in branches:
                    self.frontiers[index] = tuple(branches) + tuple(admitted)
                    if not self._can_seed(branch,seed):
                        subtree, changed, edges = (branch,), False, ()
                    else:
                        subtree, changed, edges = yield (index, pass_no, position, seed, branch)
                    self.policy_edges[index].update(dict.fromkeys(edges))
                    for child in subtree:
                        if child.node not in seen:
                            admitted.append(child)
                            seen.add(child.node)
                        if changed and child.node != branch.node:
                            advanced.add(child.node)
                if len(admitted) > cap:
                    admitted = self._prune_policy(index, admitted)
                progressed |= any(branch.node in advanced for branch in admitted)
                branches = sorted(admitted, key=lambda branch:-len(branch.mapping))
                self.max_live_branches = max(self.max_live_branches, len(branches))
        branches = self._prune_policy(index, branches)
        for branch in branches:
            self.builder.stop(branch, 'objective_met' if len(branch.mapping)==len(self.source) else 'stalled')
        self.frontiers.pop(index, None)

    def advance(self):
        if not self.agenda:
            return False
        started = time.perf_counter()
        policy, (index, pass_no, position, seed, branch) = self.agenda.popleft()
        self.builder.seed, self.builder.step = seed, (index, pass_no, position)
        self.work += 1
        output = self._expand(branch, seed, index)
        try:
            request = policy.send(output)
        except StopIteration:
            pass
        else:
            self.agenda.append((policy, request))
        if not self.agenda:
            self._finish_policies()
        self.elapsed += time.perf_counter()-started
        return True

    def snapshot(self):
        pending = tuple(dict(kind='policy_frontier', policy=index, state=state)
            for index, frontier in sorted(self.frontiers.items())
            for state in sorted({branch.node for branch in frontier}))
        graph = self.builder.finish()
        if self.agenda:
            retained = set().union(*self.policy_edges.values())
            graph = replace(graph, transitions=tuple(replace(edge, id=index)
                for index, edge in enumerate(edge for edge in graph.transitions if edge.id in retained)))
        if pending:
            graph = replace(graph, stops=graph.stops + tuple(SearchStop(state, 'budget', stage='shared_policy')
                for state in sorted({p['state'] for p in pending})))
        metrics = AAMSearchMetrics.from_record(dict(cuts=1, raw_result_count=len(graph.terminals),
            retained_branch_count=len(graph.terminals),
            subtree_branch_cap_count=sum(s.occurrences for s in graph.stops
                                         if s.reason=='capped')), self.elapsed)
        return AdaptiveResult(AAMResult(self.problem, self.config, graph, metrics), self.work,
            self.growth_calls, 0, self.reused_states, pending, not self.agenda)


@dataclass(frozen=True)
class AdaptiveResult:
    aam: AAMResult
    work: int
    growth_calls: int
    closure_calls: int
    reused_states: int
    pending: tuple
    exhausted: bool


@dataclass(frozen=True)
class AdaptiveSearchCondition:
    """Prepared directed graphs and one explicit topology/seed condition."""
    source: object
    target: object
    order: tuple
    cuts: tuple = ()
    fragment_cache: object = None

    @classmethod
    def uncut(cls, problem, config):
        source, target = [build_graph(endpoint.elements, endpoint.wbo,
            bond_cut=config.graph_floor) for endpoint in (problem.reactant, problem.product)]
        order = tuple(_generate_seed_orders(source, 1, rng_seed=config.random_seed,
                                             seed_selection=config.seed_selection)[0])
        return cls(source, target, order)
