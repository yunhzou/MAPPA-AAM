"""Unique final atom-set fragments and their branch compositions.

Growth proofs remain in the search archive. A completed fragment owns the
internal mapping search; ordered global placements remain branch context.
"""
from collections.abc import Sequence
from dataclasses import dataclass, field
from itertools import chain
import hashlib
import json

from mappa.mappa_core.core.symmetry.actions import _PreparedAction, decision_blocks


def _problem_fingerprint(problem):
    import numpy as np
    digest = hashlib.sha256()
    for endpoint in (problem.reactant, problem.product):
        digest.update(json.dumps(endpoint.elements, separators=(',', ':')).encode())
        digest.update(np.asarray(endpoint.wbo, dtype='<f8').tobytes())
    return digest.hexdigest()


def final_fragment_pairs(state):
    from mappa._group_ops import fragment_compositions
    records, compositions = fragment_compositions((state,), (0,), (0,) * (state.context + 1))
    return tuple(records[index][1:3] for index in compositions[0][1])


class FamilyOrigins(Sequence):
    """Compact archive references; expanding growth paths is an explicit view."""
    def __init__(self, origins=()):
        self.references = []
        self.count = 0
        for origin in origins:
            self.append(origin)

    def append(self, origin):
        self.references.append((None, (), origin))
        self.count += 1

    def extend(self, origins):
        for origin in origins:
            self.append(origin)

    def add(self, histories, node, origin):
        terminal = origin.pop('terminal')
        if (self.references and self.references[-1][0] is histories
                and self.references[-1][2] == origin):
            self.references[-1][1].append((terminal, node))
        else:
            self.references.append((histories, [(terminal, node)], origin))
        self.count += histories.counts[node]

    def __len__(self):
        return self.count

    def __iter__(self):
        from itertools import groupby
        from heapq import merge
        for histories, nodes, origin in self.references:
            if histories is None:
                yield origin
            else:
                for terminal, entries in groupby(nodes, key=lambda entry: entry[0]):
                    paths = [histories.paths(node) for _, node in entries]
                    ordered = paths[0] if len(paths) == 1 else merge(
                        *paths, key=lambda path: tuple(reversed(path)))
                    for transitions in ordered:
                        yield dict(origin, terminal=terminal, transitions=transitions)

    def __getitem__(self, index):
        from itertools import islice
        if isinstance(index, slice):
            return list(self)[index]
        index = index + len(self) if index < 0 else index
        if not 0 <= index < len(self):
            raise IndexError(index)
        return next(islice(self, index, None))

    def __eq__(self, other):
        return len(self) == len(other) and all(a == b for a, b in zip(self, other))

    @property
    def metadata(self):
        return (origin for _, _, origin in self.references)

    def to_record(self, history_ids):
        return [dict(history=None if histories is None else history_ids[histories],
                     nodes=nodes, origin=origin)
                for histories, nodes, origin in self.references]

    @classmethod
    def from_record(cls, record, histories):
        result = cls()
        for entry in record:
            history = entry['history']
            if history is None:
                result.append(entry['origin'])
            else:
                history = histories[history]
                nodes = [tuple(pair) for pair in entry['nodes']]
                result.references.append((history, nodes, entry['origin']))
                result.count += sum(history.counts[node] for _, node in nodes)
        return result


@dataclass
class FinalFamily:
    mapping: tuple
    required_edges: tuple
    actions: tuple
    policy: tuple
    provenance: FamilyOrigins = field(default_factory=FamilyOrigins)

    def __post_init__(self):
        if not isinstance(self.provenance, FamilyOrigins):
            self.provenance = FamilyOrigins(self.provenance)

    @property
    def fixed_assignments(self):
        """Read direct anchors from the saved, already conditioned actions."""
        moving = set()
        for action in self.actions:
            moving.update(_PreparedAction(action).support)
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
                     for block in reversed(decision_blocks(self.actions)))

    def as_path(self, problem, source_edges=None):
        """Expose the saved local/placement program through the hierarchy API."""
        from mappa.mappa_core.core.search.graph import AAMSearchGraph, SearchContext, SearchState, FragmentTransition, SearchStop
        floor, iso = self.policy
        cuts = tuple(sorted({tuple(pair) for origin in self.provenance.metadata
                             for pair in origin.get('cuts', ())}))
        anchors = tuple(sorted(pair for fragment in self.fragments
                               for pair in fragment.fixed_assignments))
        context = SearchContext(tuple(range(problem.source_atom_count)),
            tuple(range(problem.target_atom_count)), (), cuts=cuts, anchors=anchors,
            graph_floor=floor, iso_tolerance=iso, objective='saved_family')
        decisions = self.decisions
        islands = tuple(sorted((atom, i) for i, fragment in enumerate(self.fragments)
                               for atom in fragment.source_atoms))
        states = tuple(SearchState(i, 0, self.mapping, islands, ()) for i in range(len(decisions) + 1))
        edges = tuple(FragmentTransition(i, i, i + 1, None, (i, 0), fragment, i,
                                        placement_generators=generators)
                      for i, (fragment, generators) in enumerate(decisions))
        graph = AAMSearchGraph((context,), (0,), states, edges,
                              (SearchStop(len(decisions), 'objective_met'),))
        return next(graph.paths())

    @classmethod
    def from_record(cls, record, inventory, placements, histories):
        inventories = {'relation': inventory, 'group': placements}
        actions = tuple((kind, inventories[kind][index])
                        for kind, index in record['actions'])
        return cls(tuple(tuple(pair) for pair in record['mapping']),
                   tuple(sorted({edge for kind, fragment in actions if kind == 'relation'
                                 for edge in fragment.accepted_edges})), actions,
                   tuple(record['policy']), FamilyOrigins.from_record(record['provenance'], histories))

    def to_record(self, fragment_ids, placement_ids, history_ids):
        inventories = {'relation': fragment_ids, 'group': placement_ids}
        return dict(mapping=self.mapping,
                    actions=tuple((kind, inventories[kind][value])
                                  for kind, value in self.actions),
                    policy=self.policy, provenance=self.provenance.to_record(history_ids))


class FinalBranchCatalogue:
    """One problem, many saved searches; all histories remain as references."""
    schema = 'mappa.final_fragment_branches/v5'

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
        self._self_pool = {}
        self._completed_fragments = {}
        self._fragment_ids = {}
        self._fragment_records = []
        self._fragment_actions = []
        self._compositions = {}
        self._composition_owners = {}
        self.path_count = 0
        self.incomplete_paths = 0
        self.mode = 'final_atom_set_compositions'
        self.reconstruction = {}

    def _intern(self, cache, value):
        return cache.setdefault(value, value)

    def intern_actions(self, actions):
        """Share immutable local relations and placement generators across families."""
        prepared = []
        for kind, value in actions:
            if kind == 'relation':
                value = self._fragments.setdefault(value, value)
            elif kind == 'group':
                # Subgraph searches retain original labels in a shorter frame.
                if value and len(value[0]) < self.problem.target_atom_count:
                    value = tuple(tuple(g) + tuple(range(len(g), self.problem.target_atom_count))
                                  for g in value)
                value = self._intern(self._placements, value)
            else:
                raise ValueError(f'unsupported family action: {kind}')
            key = kind, value
            if key not in self._actions:
                self._actions[key] = _PreparedAction(key)
            prepared.append(self._actions[key])
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
            groups = tuple(action for action in family.actions if action[0] == 'group')
            anchors = tuple(sorted(pair for fragment in family.fragments
                                   for pair in fragment.fixed_assignments))
            ids = tuple(self._fragment_id(fragment)
                        for fragment in sorted(family.fragments, key=lambda f: f.atom_sets))
            self._family_ids[(family.policy, anchors, ids), groups] = index
            for fragment in family.fragments:
                self._completed_fragments[family.policy, anchors, fragment.atom_sets] = fragment

    def _fragment_id(self, fragment):
        if fragment not in self._fragment_ids:
            self._fragment_ids[fragment] = len(self._fragment_records)
            self._fragment_records.append(fragment)
            self._fragment_actions.append(self.intern_actions((('relation', fragment),))[0])
        return self._fragment_ids[fragment]

    def _composition(self, policy, ids, anchors):
        key = policy, anchors, ids
        if key not in self._compositions:
            fragments = tuple(self._fragment_records[index] for index in ids)
            pairs = tuple(fragment.atom_sets for fragment in fragments)
            actions = tuple(self._fragment_actions[index] for index in ids)
            witness = tuple(sorted(chain.from_iterable(fragment.witness for fragment in fragments)))
            branch = self.branches.setdefault((policy, pairs), set())
            self._compositions[key] = actions, witness, pairs, tuple(rs for rs, _ in pairs), branch
        return key, self._compositions[key]

    def add_aam(self, aam, archive='in_memory'):
        if _problem_fingerprint(aam.problem) != self.problem_sha256:
            raise ValueError('Cannot combine different endpoint problems')
        return self.add_graph(aam.graph, archive)

    def add_graph(self, graph, archive='in_memory', *, paths=None):
        from mappa._group_ops import fragment_compositions
        from mappa.mappa_core.core.matching.representation.records import Fragment
        from mappa.mappa_core.core.search.graph import SearchHistories
        histories, endings = SearchHistories.collect(graph, paths)
        conditions, context_ids = {}, []
        for context in graph.contexts:
            condition = (context.graph_floor, context.iso_tolerance), tuple(sorted(context.anchors))
            context_ids.append(conditions.setdefault(condition, len(conditions)))
        conditions = tuple(conditions)
        records, raw_compositions = fragment_compositions(graph.states,
            tuple(dict.fromkeys(end for end, _, _ in endings)), context_ids)
        local_ids = []
        for condition, source, target, witness in records:
            policy, anchors = conditions[condition]
            key = policy, anchors, (source, target)
            if key not in self._completed_fragments:
                self._completed_fragments[key] = Fragment.from_atom_sets(
                    self.problem, source, target, witness, anchors=anchors,
                    graph_floor=policy[0], self_pool=self._self_pool)
            local_ids.append(self._fragment_id(self._completed_fragments[key]))
        ids_by_terminal = {end: tuple(local_ids[index] for index in ids)
                           for end, ids in raw_compositions}
        compositions, placements = {}, {}
        for terminal, sequence, node in endings:
            state = graph.states[terminal]
            self.path_count += histories.counts[node]
            context = graph.contexts[state.context]
            policy = context.graph_floor, context.iso_tolerance
            anchors = tuple(sorted(context.anchors))
            if terminal not in compositions:
                compositions[terminal] = self._composition(policy, ids_by_terminal[terminal], anchors)
            composition, (actions, witness, pairs, partition, branch) = compositions[terminal]
            identity = id(sequence)
            if identity not in placements:
                placements[identity] = self.intern_actions(('group', group) for group in sequence)
            groups = placements[identity]
            # Before a placement moves a whole atom set, operations confined to
            # its pieces already belong to their target binary symmetries.
            first = 0
            if groups:
                if composition not in self._composition_owners:
                    owners = [-1] * self.problem.target_atom_count
                    for owner, (_, target) in enumerate(pairs):
                        for atom in target:
                            owners[atom] = owner
                    for _, atom in anchors:
                        owners[atom] = len(pairs) + atom
                    self._composition_owners[composition] = owners
                owners = self._composition_owners[composition]
                while first < len(groups) and all(owners[a] == owners[b] for a, b in groups[first].moves):
                    first += 1
            groups = groups[first:]
            key = composition, groups
            if key not in self._family_ids:
                program = self._intern(self._programs, actions + groups)
                self._family_ids[key] = len(self.families)
                self.families.append(FinalFamily(witness, (), program, policy))
            family_id = self._family_ids[key]
            branch.add(family_id)
            self.families[family_id].provenance.add(histories, node, dict(
                archive=str(archive), context=state.context, terminal=terminal,
                cuts=context.cuts, fragment_partition=partition))
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
        result.mode = record['mode']
        result.path_count = int(record['input_paths'])
        result.incomplete_paths = int(record.get('incomplete_paths', 0))
        result.reconstruction = dict(record.get('reconstruction', {}))
        from mappa.mappa_core.core.matching.representation.records import Fragment
        self_pool = {}
        inventory = tuple(Fragment.from_record(fragment, self_pool=self_pool)
                          for fragment in record['fragments'])
        placements = tuple(tuple(tuple(generator) for generator in group)
                           for group in record['placements'])
        from mappa.mappa_core.core.search.graph import SearchHistories
        histories = tuple(SearchHistories.from_record(history) for history in record['histories'])
        result.families = [FinalFamily.from_record(family_record, inventory, placements, histories)
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
        fragment_ids, placement_ids, history_ids = {}, {}, {}
        for family in self.families:
            for histories, _, _ in family.provenance.references:
                if histories is not None:
                    history_ids.setdefault(histories, len(history_ids))
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
                    histories=[history.to_record() for history in history_ids],
                    families=[f.to_record(fragment_ids, placement_ids, history_ids) for f in self.families])


def deduplicate_final_branches(aam):
    return FinalBranchCatalogue(aam.problem).add_aam(aam)

