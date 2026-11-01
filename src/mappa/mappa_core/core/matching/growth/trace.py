"""Fragment matching and growth with compressed symmetry."""
from __future__ import annotations

from collections import Counter, defaultdict
from mappa.mappa_core.core.matching.representation.state import _SymCand, _cand_map, _sym_block_assignment_expr


def cands_sample(cands, k=10):
    return [{int(a): int(b) for a, b in _cand_map(c).items()}
            for c in cands[:k]]


def cand_possible_values(cand, r):
    if isinstance(cand, _SymCand):
        for block in cand.blocks:
            if r in block.r_atoms:
                return set(block.p_atoms)
    cm = _cand_map(cand)
    return {cm[r]} if r in cm else set()


def cand_assignment_expr(cand):
    if not isinstance(cand, _SymCand):
        return '1'
    factors = []
    if cand.multiplicity != 1:
        factors.append(str(cand.multiplicity))
    for block in cand.blocks:
        expr = _sym_block_assignment_expr(block)
        if expr != '1':
            factors.append(expr)
    return ' * '.join(factors) if factors else '1'


def represented_assignment_expr(cands, max_terms=6):
    counts = Counter(cand_assignment_expr(c) for c in cands)
    if not counts:
        return '0'
    terms = []
    for expr, count in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:max_terms]:
        if expr == '1':
            terms.append(str(count))
        elif count == 1:
            terms.append(expr)
        else:
            terms.append(f'{count}*({expr})')
    remaining = max(0, len(counts) - max_terms)
    if remaining:
        terms.append(f'... + {remaining} expression groups')
    return ' + '.join(terms)


def cands_pattern_sample(cands, k=5):
    out = []
    for cand in cands[:k]:
        item = {
            'witness': {int(a): int(b) for a, b in _cand_map(cand).items()},
            'blocks': [],
        }
        if isinstance(cand, _SymCand):
            item['compression_version'] = cand.compression_version
            item['multiplicity'] = int(cand.multiplicity)
            item['automorph_domains'] = len(cand.automorph_blocks)
            for block in cand.blocks:
                item['blocks'].append({
                    'r_atoms': [int(x) for x in block.r_atoms],
                    'p_atoms': [int(x) for x in block.p_atoms],
                    'extendable': bool(block.extendable),
                    'assignments': _sym_block_assignment_expr(block),
                })
        out.append(item)
    return out


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
