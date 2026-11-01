"""Construct each compressed fragment placement from its final constraints."""
from collections import defaultdict
from time import perf_counter

from mappa.mappa_core.core.matching.growth.existence import FragmentQuery
from mappa.mappa_core.core.matching.policy import ElementNodeMatchPolicy, _orbit_id, as_node_match_policy
from mappa.mappa_core.core.matching.representation.state import _SymBlock, _SymCand
from mappa.mappa_core.core.matching.symmetry.canonical import _CandidateAutomorphismCanonicalizer


def _twin_classes(atoms, colors, rows, width):
    """Partition vertices whose transpositions preserve every edge predicate."""
    slot = (1 << width) - 1
    positions = {atom: slot << (index * width) for index, atom in enumerate(atoms)}
    by_color = defaultdict(list)
    classes = []
    for atom in atoms:
        for index in by_color[colors[atom]]:
            representative = classes[index][0]
            outside = ~(positions[atom] | positions[representative])
            if not (rows[atom] ^ rows[representative]) & outside:
                classes[index].append(atom)
                break
        else:
            by_color[colors[atom]].append(len(classes))
            classes.append([atom])
    return tuple(tuple(group) for group in classes)


class FragmentCompletion:
    """Allocate source atoms to interchangeable target classes, up to symmetry.

    A target class is interchangeable for every accepted matching predicate.
    A class allocation therefore represents all injective assignments inside
    its pools. The exact target group prevents equivalent class allocations
    before any completed candidate is constructed.
    """

    def __init__(self, source, target, fragment, bonds, locked, target_orbits,
                 tolerance, node_policy, bond_tolerance=None):
        node_policy = as_node_match_policy(node_policy)
        self.query = FragmentQuery(source, target, fragment, bonds, locked,
                                   node_policy, tolerance, bond_tolerance)
        self.canonical = _CandidateAutomorphismCanonicalizer(
            target, target_orbits, locked, node_policy)
        self.target_orbits = target_orbits
        self.custom_policy = type(node_policy) is not ElementNodeMatchPolicy
        self.compatibility = {}
        if self.custom_policy:
            self.compatibility = {p: tuple(r for r in sorted(source)
                if node_policy.compatible(source, r, target, p)) for p in target}
        self.target_generators = None
        self.locked = locked
        self.target_atoms = tuple(sorted(target))
        self.nodes = self.rejected = 0
        self.results = []
        self.classes = ()
        self.class_of = {}

    def prepare_classes(self, domains):
        allowed = {r: frozenset(values) for r, values in domains.items()}
        eligible = {p: tuple(r for r in self.query.atoms if p in allowed[r])
                    for p in self.target_atoms}
        active = tuple(p for p in self.target_atoms if eligible[p])
        predicates = {identity: bit for bit, identity in enumerate(
            dict.fromkeys(id(relation) for relation in self.query.supports.values()))}
        width = max(1, len(predicates))
        positions = {p: index * width for index, p in enumerate(active)}
        rows = dict.fromkeys(active, 0)
        constraints = {(id(relation), allowed[a], allowed[b]): relation
                       for (a, b), relation in self.query.supports.items()}
        for (identity, sources, destinations), relation in constraints.items():
            bit = predicates[identity]
            for p in sources:
                for q in relation[p] & destinations:
                    rows[p] |= 1 << (positions[q] + bit)
        # Earlier fragment actions also move later assignments. Preserve the
        # custom predicate for every source atom, including future fragments.
        self.context_roles = {p: (('eligible', eligible[p]),
            ('policy', self.compatibility.get(p, ()))) for p in self.target_atoms}
        locked_roles = {p: r for r, p in self.locked.items()}
        colors = {p: (self.context_roles[p], _orbit_id(self.target_orbits, p), locked_roles.get(p),
                      self.canonical.node_policy.key(self.query.target, p))
                  for p in self.target_atoms}
        self.classes = _twin_classes(active, colors, rows, width)
        self.class_of = {p: index for index, group in enumerate(self.classes) for p in group}

    def generators(self, assigned):
        sources = defaultdict(list)
        for r, index in sorted(assigned.items()):
            sources[index].append(r)
        roles = {p: (*self.context_roles[p], ('pool', tuple(sources[self.class_of[p]])))
                 if self.class_of.get(p) in sources else self.context_roles[p]
                 for p in self.target_atoms}
        coloring = self.canonical._colored_vertices(None, candidate_roles=roles)
        actions = self.canonical.atom_generators(None, colored_vertices=coloring)
        if self.custom_policy and not assigned:
            self.target_generators = actions
        result = []
        identity = tuple(range(len(self.classes)))
        for action in actions:
            projected = tuple(self.class_of[action[group[0]]] for group in self.classes)
            if any(self.class_of[action[p]] != projected[index]
                   for index, group in enumerate(self.classes) for p in group):
                raise ValueError('supplied target orbit partition is incompatible '
                                 'with exact target symmetry')
            if projected != identity:
                result.append(projected)
        return result

    def representatives(self, values, generators):
        remaining = set(values)
        for first in values:
            if first not in remaining:
                continue
            remaining.remove(first)
            orbit = [first]
            for index in orbit:
                for action in generators:
                    image = action[index]
                    if image in remaining:
                        remaining.remove(image)
                        orbit.append(image)
            yield first

    def emit(self, assigned):
        sources = defaultdict(list)
        for r, index in sorted(assigned.items()):
            sources[index].append(r)
        image, pools = {}, []
        for index, atoms in sorted(sources.items()):
            values = self.classes[index]
            image.update(zip(atoms, values))
            if len(values) > 1:
                pools.append(_SymBlock(tuple(atoms), values, extendable=False))
        self.results.append(_SymCand(image, pools,
            exact_fixed=tuple(r for r in self.query.atoms if r in self.locked),
            target_generators=self.target_generators))

    def search(self, domains, assigned, generators):
        self.nodes += 1
        if not self.query._reduce(domains):
            self.rejected += 1
            return
        choices = {r: tuple(dict.fromkeys(self.class_of[p] for p in domains[r]))
                   for r in self.query.atoms if r not in assigned}
        forced = {r: values[0] for r, values in choices.items() if len(values) == 1}
        if forced:
            assigned = {**assigned, **forced}
            choices = {r: values for r, values in choices.items() if r not in forced}
        if not choices:
            self.emit(assigned)
            return
        if generators is None or any(g[index] != index
                for index in set(assigned.values()) for g in generators):
            generators = self.generators(assigned)
        atom = min(choices, key=lambda r: (len(choices[r]), -len(self.query.neighbors[r]), r))
        for index in self.representatives(choices[atom], generators):
            placement = {**assigned, atom: index}
            self.search({**domains, atom: self.classes[index]}, placement, generators)

    def run(self):
        domains = dict(self.query.domains)
        if not self.query._reduce(domains):
            self.nodes = self.rejected = 1
            return []
        self.prepare_classes(domains)
        generators = self.generators({}) if self.custom_policy else None
        self.search(domains, {}, generators)
        return self.results


def complete_fragment(fragment, bonds, source, target, locked, target_orbits,
                      tolerance, node_policy, bond_tolerance=None, profile=None):
    started = perf_counter()
    query = FragmentCompletion(source, target, fragment, bonds, locked,
                               target_orbits, tolerance, node_policy, bond_tolerance)
    result = query.run()
    if not result:
        raise RuntimeError('completed fragment lost its proven matching witness')
    if profile is not None:
        profile.update(completion_seconds=perf_counter() - started,
                       completion_queries=query.nodes, completion_pruned=query.rejected,
                       completion_generated=len(result), completion_peak_candidates=len(result))
    return result
