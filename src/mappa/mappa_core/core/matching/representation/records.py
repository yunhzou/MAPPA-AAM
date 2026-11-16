"""Fragment matching inputs and compressed growth results."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
from functools import cached_property
from mappa.mappa_core.core.matching.policy import BondTolerance
from mappa.mappa_core.core.symmetry.relations import RelaxedIsomorphism
from mappa.mappa_core.core.symmetry.relations import DiscoveryEdge
from mappa.mappa_core.core.symmetry.structure import BinarySymmetry


@dataclass
class FragmentGrowth:
    """One live accepted subgraph and its current existence witness.

    Deferred edges are search continuation, not constraints on the completed
    mapping. The frontier contains proposals; it does not duplicate the graph.
    """
    source_atoms: set
    witness: dict
    accepted_edges: set = field(default_factory=set)
    deferred_edges: set = field(default_factory=set)
    frontier: list = field(default_factory=list)
    used_edges: set = field(default_factory=set, repr=False)


class IslandBranchLimitExceeded(RuntimeError):
    """A completed fragment exceeded its distinct-placement branch limit."""

    def __init__(self, count, limit, *, seed=None):
        self.count = int(count)
        self.limit = int(limit)
        self.seed = seed
        message = f"fragment subtree branch cap hit: {count}>{limit}"
        if seed is not None:
            message += f" seed={seed}"
        super().__init__(message)


@dataclass(frozen=True, eq=False)
class Fragment:
    """One completed atom-set placement, shared by search and decoding.

    The relation owns the atom lists, discovery evidence and both self
    relationships. The properties expose that data without a second payload.
    """
    relation: RelaxedIsomorphism

    @cached_property
    def atom_sets(self):
        return tuple(sorted(self.source_atoms)), tuple(sorted(self.target_atoms))

    def __reduce__(self):
        return type(self), (self.relation,)

    @property
    def source_atoms(self):
        return self.relation.source_atoms

    @property
    def target_atoms(self):
        return self.relation.target_atoms

    @property
    def accepted_edges(self):
        return tuple(edge.atoms for edge in self.relation.accepted_edges)

    @property
    def source_symmetry(self):
        return self.relation.source

    @property
    def target_symmetry(self):
        return self.relation.target

    @property
    def witness(self):
        return self.relation.witness

    @property
    def fixed_assignments(self):
        return self.relation.anchors

    @cached_property
    def key(self):
        relation = self.relation
        return (relation.source.key, relation.target.key, relation.witness,
                relation.anchors, relation.domains,
                tuple((edge.atoms, edge.weight, edge.tolerance, edge.inclusive)
                      for edge in relation.accepted_edges),
                relation.target_edges, relation.graph_floor)

    def __eq__(self, other):
        return isinstance(other, Fragment) and self.key == other.key

    @cached_property
    def _hash(self):
        return hash(self.key)

    def __hash__(self):
        return self._hash

    def relabel_target(self, action):
        """Transport coordinates and constraints together, without rematching."""
        action = dict(action)
        image = lambda atom: action.get(atom, atom)
        relation = self.relation
        target = replace(relation.target,
            atoms=tuple(image(p) for p in relation.target_atoms),
            anchors=tuple(sorted(image(p) for p in relation.target.anchors)))
        return Fragment(replace(relation, target=target,
            witness=tuple((r, image(p)) for r, p in relation.witness),
            anchors=tuple((r, image(p)) for r, p in relation.anchors),
            domains=tuple((r, tuple(image(p) for p in values)) for r, values in relation.domains),
            target_edges=tuple((image(a), image(b)) for a, b in relation.target_edges)))

    def to_record(self, *, self_pool=None, array_lists=True):
        relation = self.relation
        return dict(schema='mappa.fragment/2', source=relation.source.to_record(pool=self_pool, array_lists=array_lists),
                    target=relation.target.to_record(pool=self_pool, array_lists=array_lists), witness=relation.witness,
                    anchors=relation.anchors, domains=relation.domains,
                    accepted_edges=tuple(asdict(edge) for edge in relation.accepted_edges),
                    target_edges=relation.target_edges, graph_floor=relation.graph_floor)

    @classmethod
    def from_record(cls, record, *, self_pool=None):
        if record.get('schema') != 'mappa.fragment/2':
            raise ValueError('unsupported fragment record; historical choices cannot be reinterpreted')
        return cls(RelaxedIsomorphism(
            BinarySymmetry.from_record(record['source'], pool=self_pool),
            BinarySymmetry.from_record(record['target'], pool=self_pool),
            tuple(tuple(pair) for pair in record['witness']),
            tuple(tuple(pair) for pair in record['anchors']),
            tuple((r, tuple(values)) for r, values in record['domains']),
            tuple(DiscoveryEdge(**edge) for edge in record['accepted_edges']),
            tuple(tuple(pair) for pair in record['target_edges']), record['graph_floor']))


def fragment_from_native(record, source_nodes, target_nodes, self_pool):
    """Translate native fragment frames once, preserving shared preparation."""
    def prepared(data, nodes, side):
        key = side, int(data['self_id'])
        if key not in self_pool:
            self_pool[key] = BinarySymmetry.from_record({**data,
                'atoms': tuple(nodes[int(a)] for a in data['atoms']),
                'anchors': tuple(nodes[int(a)] for a in data['anchors'])})
        return self_pool[key]

    source = prepared(record['source'], source_nodes, 'source')
    target = prepared(record['target'], target_nodes, 'target')
    source_atoms, target_atoms = source.atoms, target.atoms
    edges = tuple(DiscoveryEdge((source_atoms[a], source_atoms[b]),
                              float(source.weights[a, b]), tolerance,
                              record['inclusive'])
                  for (a, b), tolerance in zip(record['accepted_edges'], record['edge_tolerances']))
    return Fragment(RelaxedIsomorphism(source, target,
        tuple((r, target_atoms[p]) for r, p in zip(source_atoms, record['witness'])),
        tuple((source_atoms[r], target_atoms[p]) for r, p in record['anchors']),
        tuple((r, tuple(target_atoms[p] for p in values))
              for r, values in zip(source_atoms, record['domains'])), edges,
        tuple((target_atoms[a], target_atoms[b]) for a, b in record['target_edges']),
        record['target_bond_floor']))

@dataclass(frozen=True)
class FragmentMatchContext:
    locked_mapping: dict = field(default_factory=dict)
    islands: dict | None = None
    deferred_edges: tuple = ()
    fragment_cache: object = None
    placement_fixed_targets: tuple = ()

    @classmethod
    def from_saved_fragment(cls, fragment, *, deferred_edges=()):
        """Use recorded fixed assignments as anchors in the fragment's frame.

        Within the supplied R/P fragment graphs, the matcher rediscovers free
        assignments and symmetry under its new bond tolerance.
        """
        anchors = dict(fragment.fixed_assignments)
        return cls(locked_mapping=anchors, islands={r: r for r in anchors},
                   deferred_edges=tuple(deferred_edges))


@dataclass(frozen=True)
class FragmentMatchConfig:
    graph_floor: float = 0.2
    iso_tolerance: float = 0.5
    minimum_size: int = 1
    branch_limit: int = 100  # Representative placement continuations after completion.
    node_policy: object = None
    allow_mapped_seed: bool = False
    random_seed: int = 42
    bond_tolerance: BondTolerance | None = None


@dataclass(frozen=True)
class FragmentMatchResult:
    matches: tuple
    capped: bool
    branch_count: int
    branch_limit: int
    elapsed_seconds: float
    deferred_edges: frozenset = frozenset()
    placement_generators: tuple = ()
