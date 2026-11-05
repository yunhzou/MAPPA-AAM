"""Compact permutation groups, independent factors and event invariance."""
import numpy as np
from sympy.combinatorics import Permutation, PermutationGroup

from .structure import structural_generators


def permutation_group(generators, size):
    operations = [Permutation(list(g), size=size) for g in generators]
    return PermutationGroup(operations or [Permutation(list(range(size)))])


def group_key(group):
    return tuple(sorted(tuple(g.array_form) for g in group.generators if not g.is_Identity))


def support(group):
    return tuple(sorted({i for g in group.generators
                         for i, image in enumerate(g.array_form) if image != i}))


def product_group(left, right):
    """Independent operations on two disjoint copies of the atom positions."""
    n = left.degree
    return permutation_group(
        [g.array_form + list(range(n, 2*n)) for g in left.generators]
        + [list(range(n)) + [n+i for i in g.array_form] for g in right.generators], 2*n)


def pair_orbits(group, size, pairs=None):
    """Only unordered atom pairs are expanded, never complete shuffles."""
    remaining = set(zip(*np.triu_indices(size, 1))) if pairs is None else set(pairs)
    operations = [g.array_form for g in group.generators]
    while remaining:
        first = min(remaining)
        orbit, pending = {first}, [first]
        while pending:
            i, j = pending.pop()
            for operation in operations:
                pair = tuple(sorted((operation[i], operation[j])))
                if pair not in orbit:
                    orbit.add(pair)
                    pending.append(pair)
        remaining.difference_update(orbit)
        yield tuple(sorted(orbit))


def constant_response(left, right, threshold):
    """Classify all weight comparisons without building their Cartesian table."""
    if max(abs(left[0]-right[-1]), abs(left[-1]-right[0])) <= threshold:
        return False
    if left[0]-right[-1] > threshold or right[0]-left[-1] > threshold:
        return True
    position = np.searchsorted(right, left)
    nearest = np.minimum(np.abs(left-right[np.minimum(position, len(right)-1)]),
                         np.abs(left-right[np.maximum(position-1, 0)]))
    return True if np.all(nearest > threshold) else None


def binary_group(size, edges, colors):
    return permutation_group(structural_generators(size, edges, colors), size)


def invariant_event_count(source, target, source_group, target_group, threshold):
    """Certify a constant internal count over a fixed witness's two-sided orbit.

    Matrices and groups use the witness-aligned frame. Every generator must
    preserve each pair's event response against every weight on the other side.
    Compositions then preserve those same responses. None means inconclusive;
    it is neither a minimum nor evidence that no allowed mapping exists.
    """
    for weights, other, group in ((source, target, source_group),
                                  (target, source, target_group)):
        responses = np.abs(weights[:, :, None] - np.unique(other)) > threshold
        for operation in group.generators:
            image = operation.array_form
            if not np.array_equal(responses, responses[np.ix_(image, image)]):
                return None
    return int(np.count_nonzero(np.triu(np.abs(source - target) > threshold, 1)))


def _merge_blocks(blocks, joins):
    blocks = [set(block) for block in blocks]
    for joined in joins:
        touched = [i for i, block in enumerate(blocks) if block.intersection(joined)]
        if len(touched) > 1:
            combined = set().union(*(blocks[i] for i in touched))
            blocks = [block for i, block in enumerate(blocks) if i not in touched] + [combined]
    return sorted((tuple(sorted(block)) for block in blocks), key=lambda block: block[0])


def independent_factors(group, adjacency=None):
    """Extract certified factors; keep any unresolved correlations together."""
    moved = support(group)
    if not moved:
        return ()
    generator_supports = [{i for i, image in enumerate(g.array_form) if image != i}
                          for g in group.generators]
    blocks = _merge_blocks([(i,) for i in moved], generator_supports)
    if len(blocks) > 1:
        factors = []
        for block in blocks:
            generators = [g.array_form for g, affected in zip(group.generators, generator_supports)
                          if affected and affected.issubset(block)]
            factors.extend(independent_factors(permutation_group(generators, group.degree), adjacency))
        return tuple(factors)
    orbits = sorted((sorted(orbit) for orbit in group.orbits() if len(orbit) > 1), key=len)
    if len(orbits) > 1:
        for orbit in orbits:
            if adjacency is not None:
                outside_atoms = sorted(set(range(group.degree)).difference(orbit))
                crossing = adjacency[np.ix_(outside_atoms, orbit)]
                # An independent transitive orbit must look identical to every
                # fixed outside atom. Otherwise its movement necessarily carries
                # other atoms with it, and the expensive group-order test is moot.
                if np.any(crossing != crossing[:, :1]):
                    continue
            positions = {atom: i for i, atom in enumerate(orbit)}
            outside = sorted(set(moved).difference(orbit))
            isolated = group.pointwise_stabilizer(outside, incremental=False)
            if all(g.is_Identity for g in isolated.generators):
                continue
            projected = permutation_group(
                [[positions[g.array_form[atom]] for atom in orbit] for g in group.generators], len(orbit))
            if isolated.order() == projected.order():
                remaining = group.pointwise_stabilizer(orbit, incremental=False)
                return independent_factors(isolated, adjacency) + independent_factors(remaining, adjacency)
    return (group,)


def event_components(factors, pairs):
    """A variable event connecting two factors prevents independent minimization."""
    if len(factors) <= 1:
        return factors
    owner = {atom: i for i, factor in enumerate(factors) for atom in support(factor)}
    joins = [set(owner[atom] for atom in pair if atom in owner) for pair in pairs]
    blocks = _merge_blocks([(i,) for i in range(len(factors))], joins)
    return tuple(factors[block[0]] if len(block) == 1 else permutation_group(
        [g.array_form for i in block for g in factors[i].generators], factors[0].degree)
        for block in blocks)
