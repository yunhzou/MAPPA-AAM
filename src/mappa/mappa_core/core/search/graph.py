"""Persistent fragment-decision search graph."""

from __future__ import annotations
import gc, json
from collections import defaultdict
from dataclasses import dataclass, replace
from functools import cached_property
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from mappa.mappa_core.core.matching.representation.records import Fragment
from .pairs import SharedPairs


def frozen_value(value):
    if isinstance(value, dict):
        return tuple(sorted((str(key), frozen_value(item)) for key, item in value.items()))
    if isinstance(value, (list, tuple)):
        return tuple(frozen_value(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return tuple(sorted((frozen_value(item) for item in value), key=repr))
    return value


@dataclass(frozen=True)
class SearchContext:
    source_atoms: tuple[int, ...]
    target_atoms: tuple[int, ...]
    seed_order: tuple[int, ...]
    cuts: tuple = ()
    core_atoms: tuple = ()
    anchors: tuple = ()
    graph_floor: float = 0.2
    iso_tolerance: float = 0.5
    branch_limit: int = 100
    objective: str = 'full_source'
    random_seed: int = 42
    growth_policy: str = 'growth_edge_splitmix64_v2'



@dataclass(frozen=True)
class SearchState:
    id: int
    context: int
    mapping: tuple[tuple[int, int], ...] | SharedPairs
    islands: tuple[tuple[int, int], ...] | SharedPairs
    deferred_edges: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class FragmentTransition:
    id: int
    source: int
    target: int
    seed: int | None
    step: tuple
    match: Fragment | None
    island: int = 0
    placement_generators: tuple = ()


@dataclass(frozen=True)
class SearchStop:
    """A saved state stop, or counted cap evidence with an explicit context.

    Attached stops take their context from the state. Only detached stops use
    the context field; occurrences counts stops, not attempted branch choices.
    """
    state: int | None
    reason: str
    seed: int | None = None
    stage: str = ""
    count: int = 0
    limit: int = 0
    step: tuple[int, ...] = ()
    context: int | None = None
    occurrences: int = 1


@dataclass(frozen=True)
class SearchPath:
    graph: AAMSearchGraph
    terminal: int
    transitions: tuple[int, ...]

    @property
    def mapping(self):
        return dict(self.graph.states[self.terminal].mapping)

    @property
    def deferred_edges(self):
        return self.graph.states[self.terminal].deferred_edges

    @property
    def fragments(self):
        return tuple(self.graph.transitions[index].match for index in self.transitions
                     if self.graph.transitions[index].match is not None)

    @property
    def symmetry_paths(self):
        return (self.fragments,)

    @cached_property
    def hierarchy(self):
        from mappa.mappa_core.domain import AAMHierarchy
        return AAMHierarchy(self.fragments)

    @property
    def context(self):
        return self.graph.contexts[self.graph.states[self.terminal].context]

    def to_reference(self, graph_id):
        return {'graph': graph_id, 'terminal': self.terminal,
                'transitions': list(self.transitions)}

    @classmethod
    def from_reference(cls, record, graphs):
        return cls(graphs[record['graph']], record['terminal'], tuple(record['transitions']))

    def realize(self, fragment_choices=None, placement_choices=None):
        """Apply validated local relation choices in chronological target frames.

        Each transition choice supplies discovery and either signed self words
        or checked local permutations. The resulting target permutation moves
        later decisions together, preserving their original anchor correlations.
        """
        from mappa.mappa_core.core.symmetry.actions import apply_action
        choices = {} if fragment_choices is None else dict(fragment_choices)
        placements = {} if placement_choices is None else dict(placement_choices)
        if (set(choices) | set(placements)) - set(self.transitions):
            raise ValueError('fragment choice is not on this path')
        action = {atom: atom for atom in self.context.target_atoms}
        degree = max(self.context.target_atoms, default=-1) + 1
        for edge_id in self.transitions:
            if edge_id not in choices and edge_id not in placements:
                continue
            transition = self.graph.transitions[edge_id]
            delta = list(range(degree))
            if edge_id in choices:
                if transition.match is None:
                    raise ValueError('transition has no fragment relation')
                delta = apply_action(delta, ('relation', transition.match), choices[edge_id], degree)
            if edge_id in placements:
                delta = apply_action(delta, ('group', transition.placement_generators),
                                     placements[edge_id], degree)
            action = {atom: action[int(delta[atom])] for atom in action}
        return PathRealization(self,
            tuple(sorted((atom, action[image]) for atom, image in self.mapping.items())),
            tuple(sorted(action.items())), tuple(sorted(choices.items())),
            tuple(sorted(placements.items())))


@dataclass(frozen=True)
class PathRealization:
    path: SearchPath
    mapping: tuple
    target_action: tuple
    fragment_choices: tuple
    placement_choices: tuple = ()

    @property
    def hierarchy(self):
        return self.path.hierarchy.relabel_target(self.target_action)


@dataclass(frozen=True)
class SearchBranch:
    """One unordered final fragment-pair combination and its saved alternatives.

    representative/hierarchy describe the first witness only; paths contains
    the complete union. Consumers needing one literal relation must explicitly
    use AAMSearchGraph.literal_branches().
    """

    paths: tuple[SearchPath, ...]

    @property
    def fragment_pairs(self):
        """Atom-set identity; mapping choices and matching settings are contents."""
        from mappa.mappa_core.post_processing.final_branches import final_fragment_pairs
        path = self.paths[0]
        return final_fragment_pairs(path.graph.states[path.terminal])

    @property
    def representative(self):
        from mappa.mappa_core.domain import AtomAssignment
        return AtomAssignment.from_mapping(self.paths[0].mapping)

    @property
    def hierarchy(self):
        return self.paths[0].hierarchy


class SearchHistories:
    """Shared path references containing only transition IDs and parent links.

    Reduction needs final placement programs and path counts, not expanded
    growth paths. Explicit replay can still recover every original path.
    """
    def __init__(self):
        self.parents = [()]
        self.counts = [1]

    def add(self, parents):
        node = len(self.parents)
        parents = tuple(parents)
        self.parents.append(parents)
        self.counts.append(self.counts[parents[0][1]] if len(parents) == 1
                           else sum(self.counts[parent] for _, parent in parents))
        return node

    def paths(self, node):
        stack = [(node, ())]
        while stack:
            current, suffix = stack.pop()
            if not current:
                yield suffix
            else:
                stack.extend((parent, suffix if edge is None else (edge,) + suffix)
                             for edge, parent in reversed(self.parents[current]))

    def to_record(self):
        return self.parents

    @classmethod
    def from_record(cls, record):
        result = cls()
        for parents in record[1:]:
            result.add(tuple(pair) for pair in parents)
        return result

    @classmethod
    def collect(cls, graph, paths=None):
        """Group histories by their ordered placement program once per state."""
        histories = cls()
        placements, placement_ids, objects = [], {}, {}
        programs, program_ids = [None], {}

        def append(program, edge):
            generators = edge.placement_generators
            if not generators:
                return program
            identity = id(generators)
            if identity not in objects:
                group = placement_ids.get(generators)
                if group is None:
                    group = placement_ids[generators] = len(placements)
                    placements.append(generators)
                objects[identity] = group
            key = objects[identity], program
            if key not in program_ids:
                program_ids[key] = len(programs)
                programs.append(key)
            return program_ids[key]

        if paths is not None:
            # An explicitly selected path scope uses the same compact records.
            endings, links = {}, {}
            for path in paths:
                program = node = 0
                for eid in path.transitions:
                    program = append(program, graph.transitions[eid])
                    key = eid, node
                    if key not in links:
                        links[key] = histories.add((key,))
                    node = links[key]
                endings.setdefault(path.terminal, {}).setdefault(program, []).append(node)
            endings = {end: {program: histories.add(tuple((None, node) for node in nodes))
                              for program, nodes in choices.items()}
                       for end, choices in endings.items()}
        else:
            incoming = [[] for _ in graph.states]
            for edge in graph.transitions:
                incoming[edge.target].append(edge)
            roots, ready, order = set(graph.roots), set(), []
            for end in graph.terminals:
                stack = [(end, False)]
                while stack:
                    state, complete = stack.pop()
                    if state in ready:
                        continue
                    if complete or state in roots:
                        ready.add(state)
                        order.append(state)
                    else:
                        stack.append((state, True))
                        stack.extend((edge.source, False) for edge in reversed(incoming[state]))
            choices = [None] * len(graph.states)
            for state in order:
                if state in roots:
                    choices[state] = {0: 0}
                    continue
                parents = {}
                for edge in incoming[state]:
                    for program, node in choices[edge.source].items():
                        parents.setdefault(append(program, edge), []).append((edge.id, node))
                choices[state] = {program: histories.add(parents)
                                  for program, parents in parents.items()}
            endings = {end: choices[end] for end in graph.terminals}

        expanded = {}
        for choices in endings.values():
            for program in choices:
                if program not in expanded:
                    sequence, current = [], program
                    while current:
                        group, current = programs[current]
                        sequence.append(placements[group])
                    expanded[program] = tuple(sequence)
        return histories, tuple((end, expanded[program], node)
                                for end, choices in endings.items()
                                for program, node in choices.items())


@dataclass(frozen=True)
class AAMSearchGraph:
    contexts: tuple[SearchContext, ...]
    roots: tuple[int, ...]
    states: tuple[SearchState, ...]
    transitions: tuple[FragmentTransition, ...]
    stops: tuple[SearchStop, ...]

    def stop_context(self, stop):
        """Attached historical stops derive context from their saved state."""
        if stop.state is not None:
            return self.states[stop.state].context
        if stop.context is None:
            raise ValueError('Detached stop requires its search context')
        return stop.context

    @property
    def terminals(self):
        return tuple(dict.fromkeys(stop.state for stop in self.stops
                                   if stop.reason in {"objective_met", "stalled"}))

    @property
    def capped(self):
        return any(stop.reason == "capped" for stop in self.stops)

    def ancestor_transitions(self, states):
        """All incoming decisions for selected states, without unfolding paths."""
        incoming = [[] for _state in self.states]
        for edge in self.transitions:
            incoming[edge.target].append(edge)
        pending, seen, edges = list(states), set(), set()
        while pending:
            state = pending.pop()
            if state in seen:
                continue
            seen.add(state)
            for edge in incoming[state]:
                edges.add(edge.id)
                pending.append(edge.source)
        return frozenset(edges)

    def paths(self, terminal=None):
        """Lazily unfold recorded decision paths, never atom permutations."""
        incoming = [[] for _state in self.states]
        for edge in self.transitions:
            incoming[edge.target].append(edge)
        roots = set(self.roots)
        for end in self.terminals if terminal is None else (terminal,):
            stack = [(end, ())]
            while stack:
                state, suffix = stack.pop()
                if state in roots:
                    yield SearchPath(self, end, suffix)
                else:
                    for edge in reversed(incoming[state]):
                        stack.append((edge.source, (edge.id,) + suffix))

    def branches(self):
        """Group final fragment pairs independently of seed/growth history.

        All original paths remain available: grouping never selects only one
        mapping or assumes that an intrinsic group exhausts a saved family.
        """
        from mappa.mappa_core.post_processing.final_branches import final_fragment_pairs
        groups, keys = {}, {}
        for path in self.paths():
            if path.terminal not in keys:
                keys[path.terminal] = final_fragment_pairs(self.states[path.terminal])
            groups.setdefault(keys[path.terminal], []).append(path)
        return tuple(SearchBranch(tuple(paths)) for paths in groups.values())

    def literal_branches(self):
        """Group exactly equal ordered relations without merging different histories."""
        relations = {}
        for path in self.paths():
            fragments = tuple((None if edge.match is None else edge.match.key,
                               edge.placement_generators,
                               self.states[edge.target].deferred_edges)
                              for index in path.transitions
                              if (edge := self.transitions[index]).match is not None
                              or edge.placement_generators)
            key = (tuple(sorted(path.mapping.items())), fragments)
            relations.setdefault(key, []).append(path)
        return tuple(SearchBranch(tuple(paths)) for paths in relations.values())

    def to_record(self, *, copy=True):
        """Return detached data by default; borrow nested data for immediate I/O.

        The read-only serialization view avoids recursively copying shared
        fragment/group records. It must not be retained or mutated by callers.
        """
        from copy import deepcopy
        from dataclasses import fields
        def rows(items):
            return tuple({field.name: getattr(item, field.name)
                          for field in fields(item)} for item in items)
        record = dict(schema='mappa.aam_search_graph/3',
                      contexts=rows(self.contexts), roots=self.roots,
                      states=rows(self.states), transitions=rows(self.transitions),
                      stops=rows(self.stops))
        # Expansion is explicit at the JSON/export boundary, never during search
        # or binary checkpoint publication. Binary storage preserves block sharing.
        for edge in record['transitions']:
            if edge['match'] is not None:
                edge['match'] = edge['match'].to_record()
        for state in record['states']:
            state['mapping'] = tuple(state['mapping'])
            state['islands'] = tuple(state['islands'])
        return deepcopy(record) if copy else record

    def to_json(self):
        """Encode an acyclic archive without collecting its temporary views."""
        enabled=gc.isenabled()
        try:
            gc.disable()
            return json.dumps(self.to_record(copy=False))+'\n'
        finally:
            if enabled:gc.enable()

    @classmethod
    def initial_fragment_search(cls, source, target, seed, result, config):
        """Record an unanchored one-seed call, without inventing continuation."""
        context = SearchContext(tuple(sorted(source)), tuple(sorted(target)),
            (int(seed),), graph_floor=config.graph_floor,
            iso_tolerance=config.iso_tolerance, branch_limit=config.branch_limit,

            objective='seeded_fragment')
        states = [SearchState(0, 0, (), (), ())]
        edges, stops = [], []
        for placement in result.matches:
            node = len(states)
            states.append(SearchState(node, 0, placement.witness,
                tuple((atom, 1) for atom in placement.source_atoms),
                tuple(sorted(result.deferred_edges))))
            edges.append(FragmentTransition(len(edges), 0, node, seed, (1, 0),
                                            placement, 1, result.placement_generators))
            stops.append(SearchStop(node, 'objective_met', seed, 'seeded_fragment'))
        if result.capped:
            stops.append(SearchStop(None, 'capped', seed, 'fragment_growth',
                                    result.branch_count, result.branch_limit, context=0))
        elif not result.matches:
            stops.append(SearchStop(0, 'stalled', seed, 'seeded_fragment'))
        return cls((context,), (0,), tuple(states), tuple(edges), tuple(stops))

    @classmethod
    def from_record(cls, record, *, copy=True, tuple_pool=None):
        from mappa.mappa_core.core.matching.representation.records import Fragment
        if record["schema"] != "mappa.aam_search_graph/3":
            raise ValueError("unsupported AAM search graph schema")
        tuple_pool = {} if tuple_pool is None else tuple_pool
        def pairs(items):
            values = tuple(map(tuple, items))
            shared = tuple_pool.get(values)
            if shared is not None:
                return shared
            # Atom-index pairs and unchanged mapping/island sequences recur
            # throughout the DAG. Intern immutable values, never alternatives.
            values = tuple(tuple_pool.setdefault(pair,pair) for pair in values)
            return tuple_pool.setdefault(values,values)
        contexts = tuple(SearchContext(
            **{**item,
               "source_atoms": tuple(item["source_atoms"]),
               "target_atoms": tuple(item["target_atoms"]),
               "seed_order": tuple(item["seed_order"]), "cuts": pairs(item["cuts"]),
               "anchors": pairs(item["anchors"]), "core_atoms": tuple(item["core_atoms"])})
            for item in record["contexts"])
        states = tuple(SearchState(**{**item, "mapping": pairs(item["mapping"]),
            "islands": pairs(item["islands"]), "deferred_edges": pairs(item["deferred_edges"])})
            for item in record["states"])
        self_pool = {}
        transitions = tuple(FragmentTransition(**{**item,
            'match': None if item['match'] is None else Fragment.from_record(item['match'], self_pool=self_pool),
            'placement_generators': pairs(item['placement_generators']),
            'step': tuple(item['step'])}) for item in record['transitions'])
        return cls(contexts, tuple(record["roots"]), states, tuple(transitions),
                   tuple(SearchStop(**{**item, 'step': tuple(item.get('step', ()))})
                         for item in record["stops"]))

    @classmethod
    def combine(cls, graphs):
        contexts, roots, states, transitions, stops = [], [], [], [], []
        for graph in graphs:
            c, s, t = len(contexts), len(states), len(transitions)
            contexts.extend(graph.contexts)
            roots.extend(root + s for root in graph.roots)
            states.extend(SearchState(state.id+s, state.context+c, state.mapping,
                state.islands, state.deferred_edges) for state in graph.states)
            transitions.extend(FragmentTransition(edge.id+t, edge.source+s, edge.target+s,
                edge.seed, edge.step, edge.match, edge.island, edge.placement_generators)
                for edge in graph.transitions)
            stops.extend(replace(stop,
                state=None if stop.state is None else stop.state+s,
                context=graph.stop_context(stop)+c
                if stop.state is None or stop.context is not None else None)
                for stop in graph.stops)
        return cls(tuple(contexts), tuple(roots), tuple(states), tuple(transitions), tuple(stops))


def prune_frontier(states, transitions, roots, frontier_ids, cap):
    """Limit descendant continuations bottom-up without unfolding DAG paths.

    A rejected decision edge does not remove a shared child reached through
    another parent. The final root walk retains only surviving histories.
    """
    frontier = set(frontier_ids)
    outgoing = [[] for _ in states]
    incoming_count = [0] * len(states)
    ordered = True
    for edge in transitions:
        outgoing[edge.source].append(edge)
        incoming_count[edge.target] += 1
        ordered &= edge.source < edge.target
    order = list(range(len(states))) if ordered else [
        node for node, count in enumerate(incoming_count) if not count]
    if not ordered:
        # Shared policies can reuse a state created before its new parent.
        for node in order:
            for edge in outgoing[node]:
                incoming_count[edge.target] -= 1
                if not incoming_count[edge.target]:
                    order.append(edge.target)

    descendants = [frozenset()] * len(states)
    rejected, cuts = set(), []
    for node in reversed(order):
        children = [edge for edge in outgoing[node] if descendants[edge.target]]
        own = frozenset((node,)) if node in frontier else frozenset()
        if len(children) == 1 and not own:
            descendants[node] = descendants[children[0].target]
            continue
        live = own.union(*(descendants[edge.target] for edge in children))
        if len(live) > cap:
            counts = {leaf: 1 for leaf in own}
            for edge in children:
                for leaf in descendants[edge.target]:
                    counts[leaf] = counts.get(leaf, 0) + 1
            children.sort(key=lambda edge: (len(descendants[edge.target]), edge.id), reverse=True)
            for edge in children:
                if len(counts) <= cap:
                    break
                removed = descendants[edge.target]
                cuts.append((node, edge.id, len(counts), len(removed)))
                rejected.add(edge.id)
                for leaf in removed:
                    counts[leaf] -= 1
                    if not counts[leaf]:
                        del counts[leaf]
            live = frozenset(counts)
        descendants[node] = live

    reachable, kept_edges = set(roots), set()
    pending = list(roots)
    while pending:
        for edge in outgoing[pending.pop()]:
            if edge.id in rejected or not descendants[edge.target]:
                continue
            kept_edges.add(edge.id)
            if edge.target not in reachable:
                reachable.add(edge.target)
                pending.append(edge.target)
    return frontier & reachable, kept_edges, cuts


class SearchGraphBuilder:
    def __init__(self, context):
        self.context = context
        self.states = []
        self.transitions = []
        self.stops = []
        self.roots = []
        self.step = (0, 0)
        self.seed = None
        self._partitions = {}
        self._cap_stops = {}

    def merged_partition(self, previous, atoms):
        """Merge touched source islands once for all sibling target placements."""
        key = (previous, atoms)
        partition = self._partitions.get(key)
        if partition is None:
            labels = dict(previous)
            touched = {labels[atom] for atom in atoms if atom in labels}
            groups = defaultdict(list)
            merged = set(atoms)
            for atom, island in previous:
                if island in touched:
                    merged.add(atom)
                else:
                    groups[island].append(atom)
            components = list(groups.values())
            if merged:
                components.append(sorted(merged))
            ordered = tuple((atom, index) for index, group in enumerate(sorted(components), 1)
                            for atom in group)
            partition = ordered, tuple(sorted(ordered)), len(components) + 1
            self._partitions[key] = partition
        return partition

    def state(self, branch):
        node = SearchState(len(self.states), 0, *branch.state_key())
        self.states.append(node)
        return node.id

    def root(self, branch):
        node = self.state(branch)
        self.roots.append(node)
        return node

    def commit(self, parent, branch, match, island, placement_generators=()):
        """Reference the immutable completed fragment directly."""
        node = self.state(branch)
        self.transitions.append(FragmentTransition(len(self.transitions), parent,
            node, self.seed, self.step, match, island, placement_generators))
        return node

    def join(self, kept, other):
        if kept.node == other.node:
            return
        node = self.state(kept)
        for parent in (kept.node, other.node):
            self.transitions.append(FragmentTransition(len(self.transitions), parent,
                                                       node, self.seed, self.step, None))
        kept.node = node

    def stop(self, branch, reason, *, stage="", count=0, limit=0):
        self._record_stop(SearchStop(None if reason == 'capped' else branch.node,
            reason, self.seed, stage, count, limit, self.step,
            context=0 if reason == 'capped' else None))

    def _record_stop(self, stop):
        if stop.reason == 'capped':
            key = (stop.context, stop.reason, stop.seed, stop.step,
                   stop.stage, stop.count, stop.limit)
            index = self._cap_stops.get(key)
            if index is not None:
                previous = self.stops[index]
                self.stops[index] = replace(previous,
                    occurrences=previous.occurrences+stop.occurrences)
                return
            self._cap_stops[key] = len(self.stops)
        self.stops.append(stop)

    def cap_frontier(self, branches, cap, profile=None):
        kept, edges, cuts = prune_frontier(self.states, self.transitions, self.roots,
                                          (branch.node for branch in branches), cap)
        for node, _edge, count, removed in cuts:
            self._record_stop(SearchStop(None, 'capped', self.seed,
                'combined_live_leaves', count, cap, self.step, context=0))
            if profile is not None:
                profile.append(dict(seed=self.seed, **{'pass': self.step[0]},
                    result='subtree_branch_cap', stage='combined_live_leaves',
                    mapped_before=len(self.states[node].mapping),
                    subtree_branches=removed, live_branches_before=count,
                    max_branches=cap, elapsed_sec=0.0, extend_elapsed_sec=0.0))
        self.transitions = [edge for edge in self.transitions if edge.id in edges]
        return [branch for branch in branches if branch.node in kept]

    def retain(self, branches):
        """Release discarded history while preserving every live incoming path.

        Matching supplies all branch objects whose node IDs it will reuse.
        Loading and inspecting saved graphs never prunes them.
        """
        branches = tuple({id(branch): branch for branch in branches}.values())
        required = set(self.roots)
        required.update(branch.node for branch in branches)
        required.update(stop.state for stop in self.stops
                        if stop.reason != 'capped' and stop.state is not None)
        incoming = [[] for _state in self.states]
        for edge in self.transitions:
            incoming[edge.target].append(edge)
        pending, retained = list(required), set()
        while pending:
            node = pending.pop()
            for edge in incoming[node]:
                retained.add(edge.id)
                if edge.source not in required:
                    required.add(edge.source)
                    pending.append(edge.source)
        indices = {node: index for index, node in enumerate(sorted(required))}
        states = [self.states[node] if node == index else replace(self.states[node], id=index)
                  for node, index in indices.items()]
        transitions = []
        for edge in self.transitions:
            if edge.id not in retained:
                continue
            index, source, target = len(transitions), indices[edge.source], indices[edge.target]
            transitions.append(edge if (edge.id, edge.source, edge.target) == (index, source, target)
                               else replace(edge, id=index, source=source, target=target))
        stops, self.stops = self.stops, []
        self._cap_stops.clear()
        for stop in stops:
            if stop.reason == 'capped':
                context = self.states[stop.state].context if stop.state is not None else stop.context
                self._record_stop(replace(stop, state=None, context=context))
            else:
                self._record_stop(replace(stop, state=indices[stop.state]))
        self.roots = [indices[node] for node in self.roots]
        self.states, self.transitions = states, transitions
        for branch in branches:
            branch.node = indices[branch.node]
        self._partitions.clear()

    def finish(self):
        return AAMSearchGraph((self.context,), tuple(self.roots), tuple(self.states),
                              tuple(self.transitions), tuple(self.stops))
