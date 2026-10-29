"""Fragment matching and growth with compressed symmetry."""
from __future__ import annotations

from mappa.mappa_core.core.matching.state import _cand_has_open_choice, _cand_map


MASK64 = (1 << 64) - 1


class GrowthRandom:
    """Portable SplitMix64 stream shared with the native growth engine."""
    def __init__(self, seed):
        self.state = int(seed) & MASK64

    def below(self, n):
        # Rejection removes modulo bias. Even singleton choices consume a draw.
        threshold = (1 << 64) % n
        while True:
            self.state = (self.state + 0x9E3779B97F4A7C15) & MASK64
            z = self.state
            z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
            z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK64
            z ^= z >> 31
            if z >= threshold:
                return z % n


def _pop_random_atom(pool, used_edges, fragment, rng):
    """Uniform boundary atom, then uniform incident edge; WBO never ranks."""
    pool[:] = [e for e in pool if e[2] not in fragment
               and frozenset(e[1:]) not in used_edges]
    if not pool:
        return None
    atoms = sorted({e[2] for e in pool})
    atom = atoms[rng.below(len(atoms))]
    edges = sorted((e for e in pool if e[2] == atom), key=lambda e: e[1])
    edge = edges[rng.below(len(edges))]
    pool.remove(edge)
    return edge


def _set_unique(cands):
    """True iff there is at most one fully resolved distinct mapping.

    Earlier this used `frozenset(c.values())` as the equivalence key,
    which collapses every same-value-set witness to one, silently dropping
    symmetric variants like (R87->IG87, R88->IG88) vs
    (R87->IG88, R88->IG87).  The right equivalence for a resolved candidate is
    the witness pair set itself.  Compressed symmetry candidates with an open
    target pool are not lockable yet; they must either absorb enough context or
    reach saturation first."""
    if not cands:
        return False
    if any(_cand_has_open_choice(c) for c in cands):
        return False
    if len(cands) == 1:
        return True
    sig0 = tuple(sorted(_cand_map(cands[0]).items()))
    return all(tuple(sorted(_cand_map(c).items())) == sig0 for c in cands[1:])


def _push_edges_from(heap, used_edges, g_R, atom, fragment, graph_floor):
    """Add not-yet-seen outgoing edges to the unranked boundary pool."""
    for nb in g_R.neighbors(atom):
        if nb in fragment:
            continue
        edge = frozenset({atom, nb})
        if edge in used_edges:
            continue
        w = g_R[atom][nb]['wbo']
        if w >= graph_floor:
            heap.append((-w, atom, nb))


def _frontier_boundary_edges(g_R, fragment, graph_floor):
    """One-hop live frontier edges that may still distinguish symmetry."""
    fragment = set(fragment)
    edges = set()
    for atom in fragment:
        for nb in g_R.neighbors(atom):
            if nb in fragment:
                continue
            if g_R[atom][nb]['wbo'] >= graph_floor:
                edges.add(tuple(sorted((atom, nb))))
    return edges
