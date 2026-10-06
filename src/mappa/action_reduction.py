"""Exact absorption of adjacent nested action subgroups, without enumeration.

For H <= G, both GH and HG equal G. Final-family preservation constraints
remain imposed on the resulting mapping; no intermediate condition is erased.
"""
from functools import lru_cache

@lru_cache(maxsize=1024)
def _support(action):
    kind, data = action
    if kind == 'pool':
        return frozenset(data) if len(data) > 1 else frozenset()
    return frozenset(i for g in data for i, j in enumerate(g) if i != j)

@lru_cache(maxsize=128)
def _generators(action, degree):
    kind, data = action
    if kind == 'group':
        return data
    result = []
    for a, b in zip(data, data[1:]):
        g = list(range(degree))
        g[a], g[b] = g[b], g[a]
        result.append(tuple(g))
    return tuple(result)

@lru_cache(maxsize=128)
def _group(action, degree):
    from sympy.combinatorics import Permutation, PermutationGroup
    generators = _generators(action, degree)
    return PermutationGroup([Permutation(list(g), size=degree) for g in generators]
                            or [Permutation(list(range(degree)))])

@lru_cache(maxsize=2048)
def _contains(container, member, degree):
    if container == member:
        return True
    if not _support(member) <= _support(container):
        return False
    if container[0] == 'pool':
        # A pool is the full symmetric group on its free atoms.
        return True
    from sympy.combinatorics import Permutation
    group = _group(container, degree)
    # Schreier-Sims membership checks generators, not group elements.
    return all(group.contains(Permutation(list(g), size=degree))
               for g in _generators(member, degree))

@lru_cache(maxsize=2048)
def absorb_subgroups(actions, degree):
    """Reduce a group product by exact adjacent containments only."""
    stack = []
    for action in actions:
        while stack:
            if _contains(stack[-1], action, degree):
                break
            if _contains(action, stack[-1], degree):
                stack.pop()
            else:
                stack.append(action)
                break
        else:
            stack.append(action)
    return tuple(stack)
