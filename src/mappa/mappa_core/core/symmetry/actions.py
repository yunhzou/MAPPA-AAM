"""Shared ordered structural actions and joint endpoint-pair projection."""
from collections import defaultdict, OrderedDict
from collections.abc import Mapping
from functools import cached_property
import heapq
from numbers import Integral

import numpy as np


def apply_action(images, action, choice, degree):
    """Apply one checked joint action to atom images.

    The full permutation acts jointly on all supplied images. None selects
    the saved witness; structural group words use signed generator numbers.
    A fragment choice specifies discovery plus its two self rearrangements.
    """
    images = np.asarray(images, dtype=int)
    if choice is None:
        return images
    kind, data = action
    if kind == 'relation':
        if not isinstance(choice, Mapping):
            raise ValueError('fragment choice must name discovery and self rearrangements')
        words = {'source_word', 'target_word'}.intersection(choice)
        permutations = {'source_permutation', 'target_permutation'}.intersection(choice)
        if words and permutations:
            raise ValueError('use either self words or self permutations for one fragment choice')
        realize = data.relation.realize_permutations if permutations else data.relation.realize
        mapping = realize(**choice)
        permutation = np.arange(degree)
        for source, target in data.witness:
            permutation[target] = mapping[source]
        return permutation[images]
    if kind == 'pool':
        selected = tuple(choice)
        if (any(not isinstance(p, Integral) for p in selected)
                or sorted(selected) != sorted(data)):
            raise ValueError('pool choice must permute exactly its saved atoms')
        permutation = np.arange(degree)
        permutation[list(data)] = selected
        return permutation[images]
    if kind == 'group':
        if isinstance(choice, Mapping):
            from .groups import permutation_group
            if set(choice) != {'permutation'}:
                raise ValueError('placement choice must supply one permutation or a generator word')
            selected = tuple(choice['permutation'])
            if (len(selected) != degree or any(not isinstance(p, Integral) for p in selected)
                    or sorted(selected) != list(range(degree))):
                raise ValueError('placement permutation is outside the saved group')
            # The unchanged placement belongs to every group. Checking it must
            # not reconstruct a stabilizer chain for each composed branch.
            if selected == tuple(range(degree)):
                return images
            if not permutation_group(data, degree).contains(selected):
                raise ValueError('placement permutation is outside the saved group')
            return np.asarray(selected)[images]
        for number in choice:
            if not isinstance(number, Integral) or not 1 <= abs(number) <= len(data):
                raise ValueError('group word must use signed 1-based generator numbers')
            permutation = np.asarray(data[abs(number) - 1], dtype=int)
            images = (permutation if number > 0 else np.argsort(permutation))[images]
        return images
    raise ValueError(f'unknown saved action: {kind}')


class _PreparedAction(tuple):
    """An unchanged saved action whose potentially large hash is computed once."""
    def __new__(cls, action):
        value = super().__new__(cls, action)
        value._hash = tuple.__hash__(value)
        return value

    def __hash__(self):
        return self._hash

    @cached_property
    def support(self):
        kind, value = self
        if kind == 'relation':
            return frozenset(value.target_atoms) - {p for _, p in value.relation.anchors}
        if kind == 'pool':
            return frozenset(value)
        if kind == 'group':
            return frozenset(a for generator in value for a, b in enumerate(generator) if a != b)
        raise ValueError(f'unsupported action: {kind}')

    @cached_property
    def moves(self):
        """Changed endpoint pairs of a structural placement, shared across uses."""
        kind, generators = self
        if kind != 'group':
            raise ValueError('only structural placements have generator moves')
        return tuple({(a, b) for generator in generators for a, b in enumerate(generator) if a != b})

    def __reduce__(self):
        """Publish scientific tuples; owners rebuild derived preparation."""
        return tuple, (tuple(self),)


def _support(action):
    return (action if isinstance(action, _PreparedAction) else _PreparedAction(action)).support


def decision_blocks(actions):
    """Keep a local relation and its following placement action together."""
    blocks = []
    for action in actions:
        if action[0] == 'group' and blocks and blocks[-1][-1][0] == 'relation':
            blocks[-1] += (action,)
        else:
            blocks.append((action,))
    return tuple(blocks)


def _block_key(block):
    return tuple((kind, value.key if kind == 'relation' else value) for kind, value in block)


def canonical_actions(actions):
    """Normalize independent decisions without separating their local/placement order."""
    blocks = decision_blocks(actions)
    support = tuple(frozenset().union(*(_support(a) for a in block)) for block in blocks)
    outgoing = [[] for _ in blocks]
    incoming = [0] * len(blocks)
    for j in range(len(blocks)):
        for i in range(j):
            if support[i] & support[j]:
                outgoing[i].append(j)
                incoming[j] += 1
    ready = [(_block_key(blocks[i]), i) for i, degree in enumerate(incoming) if not degree]
    heapq.heapify(ready)
    result = []
    while ready:
        _, i = heapq.heappop(ready)
        result.extend(blocks[i])
        for j in outgoing[i]:
            incoming[j] -= 1
            if not incoming[j]:
                heapq.heappush(ready, (_block_key(blocks[j]), j))
    return tuple(result)


class _ActionSequence(tuple):
    """Shared ordered actions; branch bindings remain separate from preparation."""
    def __new__(cls, actions, owner, number):
        value = super().__new__(cls, actions)
        value._hash = tuple.__hash__(value)
        value._owner, value._number = owner, number
        value.native = tuple(action.native for action in value)
        return value

    def __hash__(self):
        return self._hash


class PreparedActionPrograms:
    """Reaction-owned exact programs; input references make identity lookup safe."""
    def __init__(self, index):
        self.index = index
        self.actions, self.programs = {}, {}
        self.action_references, self.program_references = {}, {}
        self.sequences = {}
        self._owner = object()
        from mappa._group_ops import DomainPrograms
        self.native = DomainPrograms(index.m)

    def _action(self, action):
        reference = self.action_references.get(id(action))
        if reference is None:
            prepared = _PreparedAction(action)
            if prepared not in self.actions:
                prepared.domains = {image: tuple(sorted(domain | {image})) for image, domain in
                    enumerate(self.index.action_domains(prepared)) if len(domain) > 1}
                from mappa._group_ops import DomainAction
                prepared.native = DomainAction(prepared.domains)
                self.actions[prepared] = prepared
            reference = self.action_references[id(action)] = action, self.actions[prepared]
        return reference[1]

    def program(self, actions):
        if isinstance(actions, _ActionSequence) and actions._owner is self._owner:
            return actions
        reference = self.program_references.get(id(actions))
        if reference is None:
            prepared = tuple(self._action(action) for action in actions)
            if prepared not in self.programs:
                self.programs[prepared] = _ActionSequence(prepared, self._owner, len(self.programs))
            reference = self.program_references[id(actions)] = actions, self.programs[prepared]
        return reference[1]

    def select(self, program, positions):
        positions = tuple(positions)
        if positions == tuple(range(len(program))):
            return program
        key = program._number, positions
        if key not in self.sequences:
            prepared = tuple(program[i] for i in positions)
            if prepared not in self.programs:
                self.programs[prepared] = _ActionSequence(prepared, self._owner, len(self.programs))
            self.sequences[key] = self.programs[prepared]
        return self.sequences[key]

    def bind(self, actions, mapping):
        """Bind one ordered branch using shared native domains and bounded reuse."""
        return self.native.bind(self.program(actions).native, mapping)


def _action_orbits(action, degree):
    """Atom orbits from generator edges; never construct permutation elements."""
    kind, data = action
    parent = list(range(degree))

    def root(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    edges = zip(data, data[1:]) if kind == 'pool' else (
        (a, b) for generator in data for a, b in enumerate(generator) if a != b)
    for a, b in edges:
        parent[root(a)] = root(b)
    cells = defaultdict(list)
    for a in range(degree):
        cells[root(a)].append(a)
    result = [None] * degree
    for cell in cells.values():
        orbit = frozenset(cell)
        for a in cell:
            result[a] = orbit
    return tuple(result)


def _relation_domains(fragment, degree):
    """Conservative endpoint bounds; they do not permit independent choices.

    A source atom can move under alpha inverse, then D, then beta. Structural
    self orbits bound those motions without closing the relaxed relationship.
    All joint predicates remain on the saved relation and its compiler.
    """
    relation = fragment.relation
    source, target = relation.source, relation.target
    source_orbits = _action_orbits(('group', source.generators), len(source.atoms))
    target_orbits = _action_orbits(('group', target.generators), len(target.atoms))
    target_positions = {atom: i for i, atom in enumerate(target.atoms)}
    discovery = dict(relation.domains)
    discovery.update((r, (p,)) for r, p in relation.anchors)
    bounds = {}
    for orbit in set(source_orbits):
        images = {p for i in orbit for p in discovery[source.atoms[i]]}
        bounds[orbit] = frozenset(target.atoms[j] for p in images
                                 for j in target_orbits[target_positions[p]])
    result = [frozenset((i,)) for i in range(degree)]
    source_positions = {atom: i for i, atom in enumerate(source.atoms)}
    for r, p in relation.witness:
        result[p] = bounds[source_orbits[source_positions[r]]]
    return tuple(result)


class ActionIndex:
    """Prepared actions and bounded projections in one atom-index frame."""
    def __init__(self, degree):
        self.m = degree
        self._domain_tables = OrderedDict()
        self._endpoint_pairs = OrderedDict()
        self._compiled_pair_actions = OrderedDict()
        self._pair_programs = OrderedDict()

    def action_domains(self, action):
        """Endpoint bounds: exact structural orbits, conservative relation domains."""
        if action not in self._domain_tables:
            self._domain_tables[action] = (_relation_domains(action[1], self.m)
                if action[0] == 'relation' else _action_orbits(action, self.m))
            if len(self._domain_tables) > 1024:
                self._domain_tables.popitem(last=False)
        return self._domain_tables[action]

    def endpoint_pairs(self, pair, actions):
        """Inspect endpoint bounds without expanding complete mappings.

        Structural actions have exact pair projections. A program containing a
        relaxed relation yields conservative pairs lazily, because separate
        endpoint bounds cannot express its joint constraints. Membership always
        requires a checked complete choice or the finite relation compiler.
        """
        if not actions:
            return (pair,)
        if any(kind == 'relation' for kind, _ in actions):
            domains = [frozenset((p,)) for p in pair]
            for action in actions:
                bounds = self.action_domains(action)
                domains = [frozenset(q for p in domain
                    for q in (bounds[p] if p >= 0 else (-1,))) for domain in domains]
            return ((x, y) for x in sorted(domains[0]) for y in sorted(domains[1])
                    if x != y or x == -1)
        key = pair, actions
        if key not in self._endpoint_pairs:
            from mappa._group_ops import PairAction, EndpointPairProgram
            if actions not in self._pair_programs:
                compiled = []
                for action in actions:
                    if action not in self._compiled_pair_actions:
                        kind, data = action
                        factory = {'pool': PairAction.pool, 'group': PairAction.group}[kind]
                        self._compiled_pair_actions[action] = factory(self.m, data)
                        if len(self._compiled_pair_actions) > 1024:
                            self._compiled_pair_actions.popitem(last=False)
                    self._compiled_pair_actions.move_to_end(action)
                    compiled.append(self._compiled_pair_actions[action])
                self._pair_programs[actions] = EndpointPairProgram(self.m, compiled)
                if len(self._pair_programs) > 1024:
                    self._pair_programs.popitem(last=False)
            self._pair_programs.move_to_end(actions)
            self._endpoint_pairs[key] = self._pair_programs[actions].project(pair)
            if len(self._endpoint_pairs) > 1024:
                self._endpoint_pairs.popitem(last=False)
        self._endpoint_pairs.move_to_end(key)
        return self._endpoint_pairs[key]
