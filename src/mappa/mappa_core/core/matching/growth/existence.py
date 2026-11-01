"""First-solution matching of accepted growth constraints."""
from mappa.mappa_core.core.matching.policy import _edge_wbo, _growth_edge_supported


def extend_witness(source, target, witness, anchor, atom, enlarged, bonds, locked,
                   node_policy, tolerance, bond_tolerance=None):
    """Reuse a proven witness; search the whole fragment only to reassign it.

    Earlier accepted edges are already satisfied. Only the proposed edge and
    newly joined locked atoms need checking when the old assignment is retained.
    Returns the witness, visited search nodes, and whether reassignment was tried.
    """
    limit = tolerance if bond_tolerance is None else bond_tolerance.limit(anchor, atom)
    inclusive = bond_tolerance is None or bond_tolerance.inclusive
    weight = _edge_wbo(source, anchor, atom)
    floor = target.graph.get('bond_cut', .2)
    image = witness[anchor]

    def supports(value):
        return target.has_edge(image, value) and _growth_edge_supported(
            weight, _edge_wbo(target, image, value), limit, floor, inclusive=inclusive)

    used = set(witness.values())
    if atom in locked:
        added = {r: locked[r] for r in enlarged if r not in witness}
        if (supports(locked[atom]) and not used.intersection(added.values())
                and len(set(added.values())) == len(added)
                and all(node_policy.compatible(source, r, target, p) for r, p in added.items())):
            return {**witness, **added}, 1, False
    else:
        occupied = used | set(locked.values())
        for value in sorted(target.neighbors(image)):
            if (value not in occupied and node_policy.compatible(source, atom, target, value)
                    and supports(value)):
                return {**witness, atom: value}, 1, False

    query = FragmentQuery(source, target, enlarged, [*bonds, (anchor, atom)], locked,
                          node_policy, tolerance, bond_tolerance)
    result = query.find()
    return result, query.nodes, True


class FragmentQuery:
    """A fixed fragment query; a previous witness is a hint, never a constraint."""

    def __init__(self, source, target, atoms, edges, locked, node_policy,
                 tolerance, bond_tolerance=None):
        self.atoms = tuple(sorted(atoms))
        self.target = target
        self.calls = self.nodes = self.rejected = 0
        occupied = set(locked.values())
        self.domains = {
            r: tuple(p for p in target if
                     (p == locked[r] if r in locked else p not in occupied)
                     and node_policy.compatible(source, r, target, p))
            for r in self.atoms
        }
        self.neighbors = {r: [] for r in self.atoms}
        self.supports = {}
        relations = {}
        floor = target.graph.get('bond_cut', .2)
        for a, b in edges:
            weight = _edge_wbo(source, a, b)
            limit = tolerance if bond_tolerance is None else bond_tolerance.limit(a, b)
            inclusive = bond_tolerance is None or bond_tolerance.inclusive
            self.neighbors[a].append(b)
            self.neighbors[b].append(a)
            key = (weight, limit, inclusive)
            if key not in relations:
                relations[key] = {
                    p: frozenset(q for q in target.neighbors(p) if
                        _growth_edge_supported(weight, _edge_wbo(target, p, q),
                                               limit, floor, inclusive=inclusive))
                    for p in target}
            self.supports[a, b] = self.supports[b, a] = relations[key]
        self.domains = {r: tuple(p for p in values if target.degree[p] >= len(self.neighbors[r]))
                        for r, values in self.domains.items()}

    def _distinct_destinations(self, domains):
        """A bipartite matching rejects competing pools without enumerating maps."""
        owner = {values[0]: r for r, values in domains.items() if len(values) == 1}

        def augment(r, seen):
            for p in domains[r]:
                if p in seen:
                    continue
                seen.add(p)
                if p not in owner or augment(owner[p], seen):
                    owner[p] = r
                    return True
            return False

        return all(augment(r, set()) for r in sorted(self.atoms, key=lambda r: len(domains[r]))
                   if len(domains[r]) > 1)

    def _preferred_extension(self, preferred, domains):
        mapping = {r: preferred[r] for r in self.atoms if r in preferred}
        used = set(mapping.values())
        if len(used) != len(mapping) or any(p not in domains[r] for r, p in mapping.items()):
            return None
        if any(a in mapping and b in mapping and mapping[b] not in allowed[mapping[a]]
               for (a, b), allowed in self.supports.items()):
            return None
        for r in self.atoms:
            if r in mapping:
                continue
            p = next((p for p in domains[r] if p not in used and all(
                mapping[b] in self.supports[r, b][p]
                for b in self.neighbors[r] if b in mapping)), None)
            if p is None:
                return None
            mapping[r] = p
            used.add(p)
        return mapping

    def _reduce(self, domains):
        """Propagate bonds and occupied singleton domains to a fixed point."""
        changed = True
        while changed:
            changed = False
            singles = [values[0] for values in domains.values() if len(values) == 1]
            if len(set(singles)) != len(singles):
                return False
            occupied = set(singles)
            for r, values in domains.items():
                if len(values) > 1:
                    remaining = tuple(p for p in values if p not in occupied)
                    changed |= remaining != values
                    domains[r] = remaining
                if not domains[r]:
                    return False
            for (a, b), allowed in self.supports.items():
                values = tuple(p for p in domains[a]
                               if not allowed[p].isdisjoint(domains[b]))
                changed |= values != domains[a]
                domains[a] = values
                if not values:
                    return False
        return (len(set().union(*domains.values())) >= len(self.atoms)
                and self._distinct_destinations(domains))

    def _variable(self, domains):
        return min((r for r in self.atoms if len(domains[r]) > 1),
                   key=lambda r: (len(domains[r]), -len(self.neighbors[r]), r), default=None)

    def find(self, *, preferred=None):
        """Return one witness; only an exhaustive failure means no match."""
        self.calls += 1
        domains = dict(self.domains)
        if preferred is not None:
            result = self._preferred_extension(preferred, domains)
            if result is not None:
                self.nodes += 1
                return result

        def search(domains):
            self.nodes += 1
            if not self._reduce(domains):
                return None
            r = self._variable(domains)
            if r is None:
                return {r: values[0] for r, values in domains.items()}
            for p in domains[r]:
                found = search({**domains, r: (p,)})
                if found is not None:
                    return found
            return None

        result = search(domains)
        self.rejected += result is None
        return result
