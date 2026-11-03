"""Fragment matching and growth with compressed symmetry."""
from __future__ import annotations

from collections import defaultdict


def cands_sample(witnesses, k=10):
    return [dict(witness) for witness in witnesses[:k]]


def cand_possible_values(witness, atom):
    return {witness[atom]} if atom in witness else set()


def represented_assignment_expr(witnesses, max_terms=6):
    return str(len(witnesses))


def cands_pattern_sample(witnesses, k=5):
    return [dict(witness=dict(witness)) for witness in witnesses[:k]]


def heap_snapshot(heap, used_edges, fragment, mapping, k=None):
    """Pending random-pool entries in atom-index display order, not priority order."""
    peek = list(heap)
    peek.sort(key=lambda e: (e[1], e[2]))
    live = [(w, uu, nn) for (w, uu, nn) in peek
            if frozenset({uu, nn}) not in used_edges and nn not in fragment]
    if k is not None:
        live = live[:k]
    return [{'frag_atom': int(uu), 'ext_atom': int(nn),
             'wbo': round(-w, 3),
             'ext_status': ('mapped' if nn in mapping else 'free')}
            for w, uu, nn in live]


def pool_by_frag_atom(heap, used_edges, fragment, mapping, g_R):
    """Live propagation pool grouped by fragment atom."""
    peek = list(heap)
    peek.sort(key=lambda e: (e[1], e[2]))
    by_u = defaultdict(list)
    for w, uu, nn in peek:
        if frozenset({uu, nn}) in used_edges:
            continue
        if nn in fragment:
            continue
        by_u[int(uu)].append({
            'ext_atom': int(nn),
            'wbo': round(-w, 3),
            'ext_status': ('mapped' if nn in mapping else 'free'),
            'ext_element': g_R.nodes[nn]['element'],
        })
    return [{'frag_atom': int(u),
             'frag_element': g_R.nodes[u]['element'],
             'edges': sorted(by_u[u], key=lambda x: x['ext_atom'])}
            for u in sorted(by_u.keys())]
