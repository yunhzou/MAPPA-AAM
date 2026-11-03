"""First-solution matching of accepted growth constraints."""
from functools import cached_property

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
        self.global_generators = self.global_images = self.global_pruned = 0
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

    @cached_property
    def _image_components(self):
        components, forbidden = {}, set()
        # A connected query fits within one eligible target component. A
        # component with exactly enough atoms has only one possible image set.
        connected = False
        if self.atoms:
            reached = {self.atoms[0]}
            pending = [self.atoms[0]]
            for r in pending:
                for neighbor in self.neighbors[r]:
                    if neighbor not in reached:
                        reached.add(neighbor)
                        pending.append(neighbor)
            connected = len(reached) == len(self.atoms)
        if connected:
            eligible = set().union(*self.domains.values())
            for first in sorted(eligible):
                if first in components:
                    continue
                reached, pending = {first}, [first]
                for p in pending:
                    for q in self.target.neighbors(p):
                        if q in eligible and q not in reached:
                            reached.add(q)
                            pending.append(q)
                component = frozenset(reached)
                components.update((p, component) for p in component)
                if len(component) < len(self.atoms):
                    forbidden.update(component)

        return connected, components, frozenset(forbidden)

    def needs_global_symmetry(self):
        """Skip preparation when ordinary image exclusions already close the search."""
        if len(set().union(*self.domains.values())) <= len(self.atoms):
            return False
        connected, components, _ = self._image_components
        return not connected or any(len(component) > len(self.atoms)
                                    for component in components.values())

    def _global_actions(self, prepared):
        """Certify a subgroup preserving discovery and original self tolerance.

        Preserving true support pairs suffices because each generator bijects
        both domains. Self tolerance needs a whole pair-orbit proof: checking
        generators separately would incorrectly assume transitivity.
        """
        domains = {r: frozenset(values) for r, values in self.domains.items()}
        actions, generators = [], []
        for generator in prepared.generators:
            action = {atom: prepared.atoms[image]
                      for atom, image in zip(prepared.atoms, generator)}
            if any({action[p] for p in values} != values for values in domains.values()):
                continue
            if any(action[q] not in allowed[action[p]]
                   for (a, b), allowed in self.supports.items() if a < b
                   for p in domains[a] for q in allowed[p] if q in domains[b]):
                continue
            generators.append(generator)
            actions.append(action)
        if not actions or not prepared.constraints:
            return tuple(actions)
        weights, tolerance = prepared.weights, prepared.tolerance
        active = {(a, b) for a in range(len(prepared.atoms))
                  for b in range(a + 1, len(prepared.atoms)) if weights[a, b] > tolerance}
        for first in sorted(active):
            if first not in active:
                continue
            orbit, pending = {first}, [first]
            low = high = weights[first]
            for a, b in pending:
                for generator in generators:
                    pair = tuple(sorted((generator[a], generator[b])))
                    if pair in orbit:
                        continue
                    orbit.add(pair)
                    pending.append(pair)
                    weight = weights[pair]
                    low, high = min(low, weight), max(high, weight)
                    if high - low > tolerance:
                        return ()
            active.difference_update(orbit)
        return tuple(actions)

    def images(self, limit, *, preferred=None, global_self=None):
        """Return one witness per distinct target atom set, up to ``limit``.

        Emitted images become search exclusions. A certified global subgroup
        shares discovery across their complete set orbit, so exclusions remain
        invariant while internal assignment choices are pruned.
        """
        if limit < 1:
            raise ValueError('image witness limit must be positive')
        self.calls += 1
        results, excluded = [], []
        connected, components, initially_forbidden = self._image_components
        forbidden = set(initially_forbidden)
        actions = () if global_self is None else self._global_actions(global_self)
        self.global_generators = len(actions)

        def emit_one(witness):
            if connected and len(components[witness[self.atoms[0]]]) == len(self.atoms):
                forbidden.update(components[witness[self.atoms[0]]])
            else:
                excluded.append(frozenset(witness.values()))
            results.append(witness)
            return len(results) < limit

        def emit(witness):
            first = len(results)
            if not emit_one(witness):
                return False
            if not actions:
                return True
            reached = {frozenset(witness.values())}
            index = first
            while index < len(results):
                current = results[index]
                image = frozenset(current.values())
                for action in actions:
                    moved = frozenset(action[p] for p in image)
                    if moved in reached:
                        continue
                    reached.add(moved)
                    self.global_images += 1
                    if not emit_one({r: action[p] for r, p in current.items()}):
                        return False
                index += 1
            return True

        def choices(domains, atom):
            if not actions:
                return domains[atom]
            generators = tuple(action for action in actions if all(
                {action[p] for p in values} == set(values) for values in domains.values()))
            if not generators:
                return domains[atom]
            reached, representatives = set(), []
            for first in domains[atom]:
                if first in reached:
                    continue
                representatives.append(first)
                reached.add(first)
                pending = [first]
                for p in pending:
                    for action in generators:
                        image = action[p]
                        if image not in reached:
                            reached.add(image)
                            pending.append(image)
            self.global_pruned += len(domains[atom]) - len(representatives)
            return representatives

        def constrain(domains):
            if not excluded and not forbidden:
                return self._reduce(domains)
            while True:
                before = sum(map(len, domains.values()))
                if forbidden:
                    for r, values in domains.items():
                        domains[r] = tuple(p for p in values if p not in forbidden)
                        if not domains[r]:
                            return False
                for image in excluded:
                    outside = [r for r in self.atoms if not image.issuperset(domains[r])]
                    if not outside:
                        return False
                    if len(outside) == 1:
                        r = outside[0]
                        domains[r] = tuple(p for p in domains[r] if p not in image)
                if not self._reduce(domains):
                    return False
                if before == sum(map(len, domains.values())):
                    return True

        if preferred is not None:
            witness = self._preferred_extension(preferred, self.domains)
            if witness is not None:
                self.nodes += 1
                if not emit(witness):
                    return results

        def search(domains):
            self.nodes += 1
            if not constrain(domains):
                return True
            r = self._variable(domains)
            if r is None:
                return emit({r: values[0] for r, values in domains.items()})
            for p in choices(domains, r):
                if p not in domains[r]:
                    continue
                before = len(results)
                if not search({**domains, r: (p,)}):
                    return False
                # A descendant's image exclusion constrains the remaining
                # siblings without restarting the search from its root.
                if len(results) != before and not constrain(domains):
                    return True
            return True

        search(dict(self.domains))
        self.rejected += not results
        return results
