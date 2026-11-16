"""Fragment matching and growth with compressed symmetry."""
from __future__ import annotations



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
