"""First-solution matching of accepted growth constraints, without symmetry expansion."""
from mappa.mappa_core.core.matching.policy import _edge_wbo, _growth_edge_supported


class FragmentQuery:
    """A fixed fragment query; witnesses are hints, never additional constraints."""

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
            allowed = relations[key]
            self.supports[a, b] = allowed
            self.supports[b, a] = allowed

    def find(self, *, candidate=None, preferred=None):
        """Return one full witness, or None only after proving absence.

        A candidate restricts this query to its witness/pool relation for final
        compressed matching. Its exact target-group images have the same
        feasibility because the graph, constraints and earlier locks are fixed.
        """
        self.calls += 1
        result = self._find(candidate, preferred)
        self.rejected += result is None
        return result

    def _find(self, candidate, preferred):
        domains = dict(self.domains)
        if candidate is not None:
            restrictions = {r: (p,) for r, p in candidate.mapping.items()}
            for block in candidate.blocks:
                for r in block.r_atoms:
                    restrictions[r] = block.p_atoms
            for r, allowed in restrictions.items():
                domains[r] = tuple(p for p in domains[r] if p in allowed)
        # All final constraints participate in pruning, including edges whose
        # atoms have not yet been reached by compressed candidate construction.
        changed = True
        while changed:
            changed = False
            for (a, b), allowed in self.supports.items():
                values = tuple(p for p in domains[a]
                               if not allowed[p].isdisjoint(domains[b]))
                changed |= len(values) != len(domains[a])
                domains[a] = values
                if not values:
                    return None
        if any(not values for values in domains.values()):
            return None
        if len(set().union(*(set(d) for d in domains.values()))) < len(self.atoms):
            return None
        mapping, used = {}, set()
        preferred = preferred or (candidate.mapping if candidate is not None else {})

        def search():
            self.nodes += 1
            if len(mapping) == len(self.atoms):
                return dict(mapping)
            selected, values, rank = None, (), None
            for r in self.atoms:
                if r in mapping:
                    continue
                available = [p for p in domains[r] if p not in used and all(
                    mapping[b] in self.supports[r, b][p]
                    for b in self.neighbors[r] if b in mapping)]
                if not available:
                    return None
                score = (len(available), -len(self.neighbors[r]), r)
                if rank is None or score < rank:
                    selected, values, rank = r, available, score
            hint = preferred.get(selected)
            if hint in values:
                values.remove(hint)
                values.insert(0, hint)
            for p in values:
                mapping[selected] = p
                used.add(p)
                found = search()
                if found is not None:
                    return found
                used.remove(p)
                del mapping[selected]
            return None

        return search()
