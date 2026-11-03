"""Shared ordered structural actions and joint endpoint-pair projection."""
from collections import defaultdict, OrderedDict
import sys
from numbers import Integral

import numpy as np


def apply_action(images, action, choice, degree):
    """Apply one pool permutation or signed generator word to atom images.

    The full permutation acts jointly on all supplied images. None selects
    identity; a group word uses signed 1-based generator numbers.
    """
    images = np.asarray(images, dtype=int)
    if choice is None:
        return images
    kind, data = action
    if kind == 'pool':
        selected = tuple(choice)
        if (any(not isinstance(p, Integral) for p in selected)
                or sorted(selected) != sorted(data)):
            raise ValueError('pool choice must permute exactly its saved atoms')
        permutation = np.arange(degree)
        permutation[list(data)] = selected
        return permutation[images]
    if kind == 'group':
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


class _AtomRows:
    """Bounded final atom rows; numeric keys do not retain their programs."""
    def __init__(self):
        self.rows = OrderedDict()
        self.bytes = 0
        self.next_program = 0

    def atom(self, program, image):
        key = program._number, image
        entry = self.rows.get(key)
        if entry is not None:
            self.rows.move_to_end(key)
            return entry[0]
        row = program._prepare_atom(image)
        size = sys.getsizeof(key) + sys.getsizeof(row) + sum(map(sys.getsizeof, row))
        self.rows[key] = row, size
        self.bytes += size
        while len(self.rows) > 4096 or self.bytes > 16 * 2**20:
            _, (_, removed) = self.rows.popitem(last=False)
            self.bytes -= removed
        return row


class _ActionSequence(tuple):
    """Exact ordered action values with bounded final atom preparation."""
    def __new__(cls, actions, rows):
        value = super().__new__(cls, actions)
        value._hash = tuple.__hash__(value)
        value._rows = rows
        value._number = rows.next_program
        rows.next_program += 1
        return value

    def __hash__(self):
        return self._hash

    def atom(self, image):
        return self._rows.atom(self, image)

    def _prepare_atom(self, image):
        domain, dependencies = frozenset((image,)), []
        for position, action in enumerate(self):
            moved = domain.intersection(action.orbits)
            if moved:
                dependencies.append(position)
                reached = {action.orbits[p] for p in moved}
                domain = domain.union(*reached)
        return tuple(sorted(domain)), tuple(dependencies)


class PreparedActionPrograms:
    """Reaction-owned exact programs; input references make identity lookup safe."""
    def __init__(self, index):
        self.index = index
        self.actions, self.programs = {}, {}
        self.action_references, self.program_references = {}, {}
        self.sequences = {}

    def _action(self, action):
        reference = self.action_references.get(id(action))
        if reference is None:
            prepared = _PreparedAction(action)
            if prepared not in self.actions:
                prepared.orbits = {image: orbit for image, orbit in
                    enumerate(self.index.action_orbits(prepared)) if len(orbit) > 1}
                self.actions[prepared] = prepared
            reference = self.action_references[id(action)] = action, self.actions[prepared]
        return reference[1]

    def program(self, actions):
        if isinstance(actions, _ActionSequence) and actions._rows is self.index._prepared_atoms:
            return actions
        reference = self.program_references.get(id(actions))
        if reference is None:
            prepared = tuple(self._action(action) for action in actions)
            if prepared not in self.programs:
                self.programs[prepared] = _ActionSequence(prepared, self.index._prepared_atoms)
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
                self.programs[prepared] = _ActionSequence(prepared, self.index._prepared_atoms)
            self.sequences[key] = self.programs[prepared]
        return self.sequences[key]


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


class ActionIndex:
    """Prepared actions and bounded projections in one atom-index frame."""
    def __init__(self, degree):
        self.m = degree
        self._orbit_tables = OrderedDict()
        self._prepared_atoms = _AtomRows()
        self._endpoint_pairs = OrderedDict()
        self._compiled_pair_actions = OrderedDict()
        self._pair_programs = OrderedDict()

    def action_orbits(self, action):
        if action not in self._orbit_tables:
            self._orbit_tables[action] = _action_orbits(action, self.m)
            if len(self._orbit_tables) > 1024:
                self._orbit_tables.popitem(last=False)
        return self._orbit_tables[action]

    def endpoint_pairs(self, pair, actions):
        """Project the ordered action program onto two endpoints, never mappings.

        Each group closes only its pair orbit before the next action starts.
        A pool contributes at most k*(k-1) pairs, not k! permutations. These
        projections retain correlation between endpoints; constraints involving
        other atoms are still enforced when a complete choice is realized.
        """
        if not actions:
            return (pair,)
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
