"""Compact permutation groups, independent factors and event invariance."""
from mappa._group_ops import PermutationGroup

from .structure import structural_generators


def permutation_group(generators, size):
    return PermutationGroup(size, list(generators))


def group_key(group):
    return tuple(sorted(group.generators))


def support(group):
    return tuple(group.support())


def product_group(left, right):
    """Independent operations on two disjoint copies of the atom positions."""
    n = left.degree
    return permutation_group(
        [list(g) + list(range(n, 2*n)) for g in left.generators]
        + [list(range(n)) + [n+i for i in g] for g in right.generators], 2*n)


def binary_group(size, edges, colors):
    return permutation_group(structural_generators(size, edges, colors), size)


def independent_factors(group, adjacency=None):
    """Extract certified native factors, retaining unresolved correlations."""
    from mappa._group_ops import independent_factors as split
    return tuple(split(group, adjacency))


def event_components(factors, pairs):
    """Join native factors whenever an event or endpoint condition couples them."""
    from mappa._group_ops import event_components as join
    return tuple(join(factors, [tuple(pair) for pair in pairs]))
