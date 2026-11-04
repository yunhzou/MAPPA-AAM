"""Collect unique fragments, then share their conditional local event results."""
from dataclasses import dataclass, replace

from mappa.mappa_core.utils.bonds import bond_event_threshold
from mappa.mappa_core.core.symmetry.actions import PreparedActionPrograms
from .records import BondEventTerm


@dataclass(frozen=True, slots=True)
class PartitionPiece:
    """Disjoint bond ownership in one final branch, not a matching record."""
    source_atoms: tuple[int, ...]
    target_atoms: tuple[int, ...]


class FragmentInventory:
    """Collect actual completed fragments before preparing event responses.

    Historical fragments can overlap. Final partition pieces give each bond
    one owner; the original fragment relations retain all joint permissions.
    """
    def __init__(self, catalogue):
        fragments = {}
        for family in catalogue.families:
            for fragment in family.fragments:
                fragments.setdefault(fragment.key, fragment)
        self.fragments = tuple(fragments.values())
        pieces, partitions = {}, [set() for _ in catalogue.families]

        def intern(key):
            if key not in pieces:
                pieces[key] = PartitionPiece(*key)
            return pieces[key]

        def collect(pairs, family_ids):
            source = tuple(sorted(tuple(sorted(rs)) for rs, _ in pairs))
            for rs, ps in pairs:
                key = tuple(sorted(rs)), tuple(sorted(ps))
                intern(key)
            for fid in family_ids:
                partitions[fid].add(source)

        for (_, pairs), family_ids in catalogue.branches.items():
            collect(pairs, family_ids)
        for fid, alternatives in catalogue.coupled_branches.items():
            for pairs in alternatives:
                collect(pairs, (fid,))
        self.partitions = tuple(tuple(sorted(parts)) for parts in partitions)
        self.uses = []
        for family, alternatives in zip(catalogue.families, self.partitions):
            mapping = dict(family.mapping)
            uses = []
            for partition in alternatives:
                refs = []
                for atoms in partition:
                    key = atoms, tuple(sorted(mapping[a] for a in atoms))
                    refs.append(intern(key))
                uses.append(tuple(refs))
            self.uses.append(tuple(uses))
        self.uses = tuple(self.uses)
        self.pieces = tuple(pieces[key] for key in sorted(pieces))


@dataclass(frozen=True, slots=True)
class FragmentResponse:
    piece: PartitionPiece
    atom_domains: dict
    atom_actions: dict
    terms: tuple[BondEventTerm, ...]
    projected_pair_states: int
    neighboring_pair_count: int


class FragmentDecoder:
    """Owned decoded results shared by fragment uses and boundary references."""
    def __init__(self, index, inventory):
        self.index, self.inventory = index, inventory
        self.programs = PreparedActionPrograms(index)
        self.fragments = {}
        self.bonds = {}
        self.terms = {}
        self.bindings = {}

    def bond(self, atoms, pair, actions, domains, action_ids, manual_cut):
        """Decode one exact local relationship, retaining its ordered actions."""
        xs, ys = domains
        key = atoms, pair, actions, xs, ys
        if key not in self.bonds:
            a, b = atoms
            tau = bond_event_threshold(self.index.elements, a, b,
                default_threshold=self.index.threshold, metal_threshold=self.index.metal_threshold)
            weight = self.index.r[a, b]
            response = self.index.domain_response(weight, tau, xs, ys)
            response, projected = self.index.refine_response(weight, tau, pair, actions, response)
            self.bonds[key] = response, projected
        response, projected = self.bonds[key]
        if response.constant == 0 and not manual_cut:
            return None, projected
        term_key = key, action_ids, manual_cut
        if term_key not in self.terms:
            moving = int(len(xs) > 1) + int(len(ys) > 1)
            self.terms[term_key] = BondEventTerm(atoms, action_ids, response, pair,
                actions, self.index, ('fixed_fixed', 'fixed_shuffle', 'shuffle_shuffle')[moving], manual_cut)
        return self.terms[term_key], projected

    def fragment(self, piece, mapping, actions, domains, dependencies, cut_bonds):
        """Reuse a partition's local mapping and conditional event responses.

        Local action positions are bound back to each family's original positions;
        a shared action affecting two fragments therefore remains one choice.
        """
        atoms = piece.source_atoms
        binding = tuple(sorted({i for a in atoms for i in dependencies[a]}))
        local_actions = self.programs.select(self.programs.program(actions), binding)
        pairs = tuple((a, mapping[a]) for a in atoms)
        cuts = tuple(sorted((a, b) for a, b in cut_bonds if a in atoms and b in atoms))
        key = piece, pairs, local_actions, cuts
        if key not in self.fragments:
            positions = {aid: i for i, aid in enumerate(binding)}
            local_domains = {a: domains[a] for a in atoms}
            local_dependencies = {a: tuple(positions[i] for i in dependencies[a]) for a in atoms}
            neighbors = self.index.neighboring_pairs(local_domains)
            neighbors = {(a, b) for a, b in neighbors if a in local_domains and b in local_domains}
            neighbors.update(cuts)
            terms, projected = [], 0
            sequences = {}
            for a, b in sorted(neighbors):
                ids = tuple(sorted(set(local_dependencies[a]) | set(local_dependencies[b])))
                if ids not in sequences:
                    sequences[ids] = self.programs.select(local_actions, ids)
                term, count = self.bond((a, b), (mapping[a], mapping[b]), sequences[ids],
                    (local_domains[a], local_domains[b]), ids, (a, b) in cuts)
                projected += count
                if term is not None:
                    terms.append(term)
            self.fragments[key] = FragmentResponse(piece, local_domains,
                local_dependencies, tuple(terms), projected, len(neighbors))
        return self.fragments[key], binding

    def bind(self, terms, binding):
        if binding == tuple(range(len(binding))):
            return terms
        bound = []
        for term in terms:
            ids = tuple(binding[i] for i in term.action_ids)
            if ids == term.action_ids:
                bound.append(term)
                continue
            key = term, ids
            if key not in self.bindings:
                self.bindings[key] = replace(term, action_ids=ids)
            bound.append(self.bindings[key])
        return tuple(bound)
