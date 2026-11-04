"""Unordered final fragment pairs with a lossless union of mapping families.

Branch identity excludes growth history. Families retain complete local fragment
relationships and global placement actions in their composition order. Only
disjoint decisions may be reordered; a constrained local relationship is never
absorbed as a placement group.
"""
from collections import defaultdict
from dataclasses import dataclass, field
import heapq
import hashlib
import json

from mappa.mappa_core.core.symmetry.actions import _PreparedAction


def _problem_fingerprint(problem):
    import numpy as np
    digest = hashlib.sha256()
    for endpoint in (problem.reactant, problem.product):
        digest.update(json.dumps(endpoint.elements, separators=(',', ':')).encode())
        digest.update(np.asarray(endpoint.wbo, dtype='<f8').tobytes())
    return digest.hexdigest()


def final_fragment_pairs(state):
    mapping = dict(state.mapping)
    owners = dict(state.islands)
    groups = defaultdict(list)
    for r in sorted(mapping):
        # Unlabelled anchors are singleton fragments.
        groups[('island', owners[r]) if r in owners else ('single', r)].append(r)
    return tuple(sorted((tuple(atoms), tuple(sorted(mapping[r] for r in atoms)))
                        for atoms in groups.values()))


def _support(action):
    return (action if isinstance(action, _PreparedAction) else _PreparedAction(action)).support


def _decision_blocks(actions):
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
    blocks = _decision_blocks(actions)
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


@dataclass
class FinalFamily:
    mapping: tuple
    required_edges: tuple
    actions: tuple
    policy: tuple
    provenance: list = field(default_factory=list)

    @property
    def fixed_assignments(self):
        """Read direct anchors from the saved, already conditioned actions."""
        moving = set()
        for action in self.actions:
            moving.update(_support(action))
        return tuple((r, p) for r, p in self.mapping if p not in moving)

    def validate_representative(self, problem):
        mapping = dict(self.mapping)
        assert set(mapping) <= set(range(problem.source_atom_count))
        assert set(mapping.values()) <= set(range(problem.target_atom_count))
        assert len(set(mapping.values())) == len(mapping)
        assert all(problem.reactant.elements[r] == problem.product.elements[p] for r, p in mapping.items())
        floor, iso = self.policy
        assert all(problem.product.wbo[mapping[a], mapping[b]] >= floor and
                   abs(problem.reactant.wbo[a, b] - problem.product.wbo[mapping[a], mapping[b]]) <= iso + 1e-9
                   for a, b in self.required_edges)

    @property
    def fragments(self):
        return tuple(fragment for kind, fragment in reversed(self.actions) if kind == 'relation')

    @property
    def decisions(self):
        """Chronological local/placement blocks, including explicit placement-only actions."""
        return tuple((next((value for kind, value in block if kind == 'relation'), None),
                      next((value for kind, value in block if kind == 'group'), ()))
                     for block in reversed(_decision_blocks(self.actions)))

    def as_path(self, problem, source_edges=None):
        """Expose the saved local/placement program through the hierarchy API."""
        from mappa.mappa_core.core.search.graph import AAMSearchGraph, SearchContext, SearchState, FragmentTransition, SearchStop
        floor, iso = self.policy
        cuts = tuple(sorted({tuple(pair) for origin in self.provenance
                             for pair in origin.get('cuts', ())}))
        context = SearchContext(tuple(range(problem.source_atom_count)),
            tuple(range(problem.target_atom_count)), (), cuts=cuts,
            graph_floor=floor, iso_tolerance=iso, objective='saved_family')
        decisions = self.decisions
        states = tuple(SearchState(i, 0, self.mapping, (), ()) for i in range(len(decisions) + 1))
        edges = tuple(FragmentTransition(i, i, i + 1, None, (i, 0), fragment, i,
                                        placement_generators=generators)
                      for i, (fragment, generators) in enumerate(decisions))
        graph = AAMSearchGraph((context,), (0,), states, edges,
                              (SearchStop(len(decisions), 'objective_met'),))
        return next(graph.paths())

    @classmethod
    def from_record(cls, record, inventory, placements):
        inventories = {'relation': inventory, 'group': placements}
        actions = tuple((kind, inventories[kind][index])
                        for kind, index in record['actions'])
        return cls(tuple(tuple(pair) for pair in record['mapping']),
                   tuple(sorted({edge for kind, fragment in actions if kind == 'relation'
                                 for edge in fragment.accepted_edges})), actions,
                   tuple(record['policy']), list(record.get('provenance', ())))

    def to_record(self, fragment_ids, placement_ids):
        inventories = {'relation': fragment_ids, 'group': placement_ids}
        return dict(mapping=self.mapping,
                    actions=tuple((kind, inventories[kind][value])
                                  for kind, value in self.actions),
                    policy=self.policy, provenance=self.provenance)


class FinalBranchCatalogue:
    """One problem, many saved searches; all histories remain as references."""
    schema = 'mappa.final_fragment_branches/v3'

    def __init__(self, problem):
        self.problem = problem
        self.problem_sha256 = _problem_fingerprint(problem)
        self.branches = {}
        self.coupled_branches = {}
        self.families = []
        self._family_ids = {}
        self._required = {}
        self._programs = {}
        self._fragments = {}
        self._placements = {}
        self._actions = {}
        self.path_count = 0
        self.incomplete_paths = 0
        self.mode = 'lossless_saved_family_union'
        self.reconstruction = {}

    def _intern(self, cache, value):
        return cache.setdefault(value, value)

    def intern_actions(self, actions):
        """Share immutable local relations and placement generators across families."""
        prepared = []
        for kind, value in actions:
            if kind == 'relation':
                value = self._fragments.setdefault(value.key, value)
            elif kind == 'group':
                # Subgraph searches retain original labels in a shorter frame.
                if value and len(value[0]) < self.problem.target_atom_count:
                    value = tuple(tuple(g) + tuple(range(len(g), self.problem.target_atom_count))
                                  for g in value)
                value = self._intern(self._placements, value)
            else:
                raise ValueError(f'unsupported family action: {kind}')
            prepared.append(self._intern(self._actions, _PreparedAction((kind, value))))
        return tuple(prepared)

    def __reduce__(self):
        """Binary transport retains shared scientific objects, not JSON lists."""
        from mappa.mappa_core.utils.io.records import endpoint_record
        problem = dict(name=self.problem.name,
            reactant=endpoint_record(self.problem.reactant),
            product=endpoint_record(self.problem.product))
        state = {name: getattr(self, name) for name in (
            'branches', 'coupled_branches', 'families', 'path_count',
            'incomplete_paths', 'mode', 'reconstruction')}
        return type(self)._from_binary, (self.schema, problem, state)

    @classmethod
    def _from_binary(cls, schema, problem, state):
        from mappa.mappa_core.utils.io.records import problem_from_record
        if schema != cls.schema:
            raise ValueError('Unsupported final-branch schema')
        result = cls(problem_from_record(problem))
        for name, value in state.items():
            setattr(result, name, value)
        result._index_families()
        return result

    def _index_families(self):
        """Rebuild derived lookup indexes while preserving saved joint choices."""
        mappings, pairs = {}, {}
        for index, family in enumerate(self.families):
            family.mapping = self._intern(mappings, tuple(
                self._intern(pairs, pair) for pair in family.mapping))
            family.required_edges = self._intern(self._required, tuple(
                self._intern(pairs, pair) for pair in family.required_edges))
            family.actions = self._intern(self._programs, self.intern_actions(family.actions))
            self._family_ids[family.policy, family.mapping, family.actions] = index

    def add_aam(self, aam, archive='in_memory'):
        if _problem_fingerprint(aam.problem) != self.problem_sha256:
            raise ValueError('Cannot combine different endpoint problems')
        return self.add_graph(aam.graph, archive)

    def add_graph(self, graph, archive='in_memory', *, paths=None):
        pairs = {}
        for path in graph.paths() if paths is None else paths:
            state = graph.states[path.terminal]
            self.path_count += 1
            context = path.context
            policy = context.graph_floor, context.iso_tolerance
            if path.terminal not in pairs:
                pairs[path.terminal] = final_fragment_pairs(state)
            branch = self.branches.setdefault((policy, pairs[path.terminal]), set())
            program = []
            for eid in reversed(path.transitions):
                transition = graph.transitions[eid]
                fragment = transition.match
                if fragment is not None:
                    fragment = self._fragments.setdefault(fragment.key, fragment)
                    program.append(('relation', fragment))
                if transition.placement_generators:
                    program.append(('group', transition.placement_generators))
            program = self._intern(self._programs, canonical_actions(self.intern_actions(program)))
            required = self._intern(self._required, tuple(sorted({edge
                for kind, fragment in program if kind == 'relation' for edge in fragment.accepted_edges})))
            mapping = tuple(sorted(state.mapping))
            key = policy, mapping, program
            if key not in self._family_ids:
                self._family_ids[key] = len(self.families)
                self.families.append(FinalFamily(mapping, required, program, policy))
            family_id = self._family_ids[key]
            branch.add(family_id)
            self.families[family_id].provenance.append(dict(
                archive=str(archive), context=state.context, terminal=path.terminal,
                transitions=path.transitions, cuts=context.cuts,
                source_partition=tuple(atoms for atoms, _ in pairs[path.terminal])))
        return self

    @property
    def branch_count(self):
        return len(self.branches) + len(self.coupled_branches)

    @classmethod
    def from_record(cls, problem, record):
        if record.get('schema') != cls.schema:
            raise ValueError('Unsupported final-branch schema')
        result = cls(problem)
        if record.get('problem_sha256') != result.problem_sha256:
            raise ValueError('Final branches belong to a different endpoint problem')
        result.mode = record.get('mode', 'lossless_saved_family_union')
        result.path_count = int(record['input_paths'])
        result.incomplete_paths = int(record.get('incomplete_paths', 0))
        result.reconstruction = dict(record.get('reconstruction', {}))
        from mappa.mappa_core.core.matching.representation.records import Fragment
        self_pool = {}
        inventory = tuple(Fragment.from_record(fragment, self_pool=self_pool)
                          for fragment in record['fragments'])
        placements = tuple(tuple(tuple(generator) for generator in group)
                           for group in record['placements'])
        result.families = [FinalFamily.from_record(family_record, inventory, placements)
                           for family_record in record['families']]
        result._index_families()
        atom_sets, fragment_pairs = {}, {}
        def restore_pairs(pairs):
            restored = []
            for rs, ps in pairs:
                pair = (result._intern(atom_sets, tuple(rs)), result._intern(atom_sets, tuple(ps)))
                restored.append(result._intern(fragment_pairs, pair))
            return tuple(restored)
        for branch in record['branches']:
            pairs = restore_pairs(branch['fragment_pairs'])
            result.branches[tuple(branch['policy']), pairs] = set(branch['families'])
        for branch in record.get('coupled_branches', ()):
            result.coupled_branches[branch['family']] = {
                restore_pairs(pairs) for pairs in branch['source_fragment_pairs']}
        for ids in (*result.branches.values(), result.coupled_branches):
            if any(i < 0 or i >= len(result.families) for i in ids):
                raise ValueError('Invalid family reference')
        return result

    def to_record(self):
        fragment_ids, placement_ids = {}, {}
        for family in self.families:
            for fragment in family.fragments:
                fragment_ids.setdefault(fragment, len(fragment_ids))
            for kind, generators in family.actions:
                if kind == 'group':
                    placement_ids.setdefault(generators, len(placement_ids))
        self_pool = {}
        return dict(schema=self.schema, problem_sha256=self.problem_sha256, mode=self.mode, reconstruction=self.reconstruction, input_paths=self.path_count,
                    incomplete_paths=self.incomplete_paths,
                    branches=[dict(policy=policy, fragment_pairs=pairs, families=sorted(families))
                              for (policy, pairs), families in self.branches.items()],
                    coupled_branches=[dict(family=k, source_fragment_pairs=list(v)) for k, v in self.coupled_branches.items()],
                    fragments=[fragment.to_record(self_pool=self_pool) for fragment in fragment_ids],
                    placements=list(placement_ids),
                    families=[f.to_record(fragment_ids, placement_ids) for f in self.families])


def deduplicate_final_branches(aam):
    return FinalBranchCatalogue(aam.problem).add_aam(aam)

