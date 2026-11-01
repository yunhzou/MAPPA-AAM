"""Fragment matching and growth with compressed symmetry."""
from __future__ import annotations

from dataclasses import dataclass


def candidate_from_record(state):
    """Reconstruct the compressed native/Python candidate, not its bijections."""
    # Unversioned archives describe the historical group, not a new search.
    witness = {int(r): int(p) for r, p in state['witness'].items()}
    automorph_blocks = tuple(_SymBlock(
        tuple(block['r_atoms']), tuple(block['p_atoms']), extendable=False)
        for block in state.get('automorph_blocks', ()))
    blocks = tuple(_SymBlock(tuple(block['r_atoms']), tuple(block['p_atoms']),
                             extendable=bool(block.get('extendable', False)))
        for block in state.get('blocks', ())
        if block.get('source') != 'exact_automorph_group')
    generators = state.get('automorph_generators')
    return _SymCand(witness, blocks, exact_fixed=tuple(state.get('exact_fixed', ())),
                    multiplicity=int(state.get('multiplicity', 1)),
                    automorph_blocks=automorph_blocks,
                    compression_version=state.get('compression_version', 1),
                    target_generators=None if generators is None else
                    tuple(tuple(generator) for generator in generators))


@dataclass(frozen=True)
class _SymBlock:
    """Source atoms occupy distinct members of an interchangeable target pool."""
    r_atoms: tuple
    p_atoms: tuple
    extendable: bool = True

    def __post_init__(self):
        object.__setattr__(self, 'r_atoms', tuple(sorted(set(self.r_atoms))))
        object.__setattr__(self, 'p_atoms', tuple(sorted(set(self.p_atoms))))

    @property
    def open(self):
        return len(self.r_atoms) < len(self.p_atoms)

    @property
    def complete(self):
        return len(self.r_atoms) >= len(self.p_atoms)


def pool_witness(mapping, blocks):
    """Keep valid supplied assignments and fill unassigned pool members."""
    block_r = {r for b in blocks for r in b.r_atoms}
    raw = dict(mapping or {})
    m = {r: p for r, p in raw.items() if r not in block_r}
    used = set(m.values())
    for b in blocks:
        if len(b.r_atoms) > len(b.p_atoms):
            raise ValueError("symmetry block has more R atoms than P atoms")
        for r in b.r_atoms:
            if r not in raw:
                continue
            p = raw[r]
            if p not in b.p_atoms or p in used:
                raise ValueError("symmetry block witness conflicts with fixed atoms")
            m[r] = p
            used.add(p)
        available = [p for p in b.p_atoms if p not in used]
        missing = [r for r in b.r_atoms if r not in m]
        if len(available) < len(missing):
            raise ValueError("symmetry block witness conflicts with fixed atoms")
        for r, p in zip(missing, available):
            m[r] = p
            used.add(p)
    return m


class _SymCand:
    """A completed placement: one witness, injective pools and exact fixed roles.

    Legacy group domains remain readable for historical saved relations.
    New results use the conditioned global target action directly.
    """
    __slots__ = (
        'mapping', 'blocks', 'exact_fixed', 'multiplicity',
        'automorph_blocks', 'compression_version', 'target_generators',
    )

    def __init__(self, mapping=None, blocks=(), exact_fixed=(), multiplicity=1,
                 automorph_blocks=(), compression_version=3, target_generators=None):
        self.compression_version = compression_version
        # None defers group construction; an empty tuple is the exact identity
        # group. Completion shares this immutable group across its placements.
        self.target_generators = target_generators
        blocks = tuple(blocks)
        m = pool_witness(mapping, blocks)
        self.mapping = m
        self.blocks = blocks
        self.exact_fixed = frozenset(exact_fixed)
        self.multiplicity = int(multiplicity)
        expanded_automorph_blocks = []
        for block in automorph_blocks:
            p_atoms = set(block.p_atoms)
            r_atoms = set(block.r_atoms)
            r_atoms.update(r for r, p in m.items() if p in p_atoms)
            expanded_automorph_blocks.append(_SymBlock(
                tuple(r_atoms), tuple(p_atoms), extendable=False))
        self.automorph_blocks = tuple(expanded_automorph_blocks)

    def __contains__(self, r):
        return r in self.mapping

    def __getitem__(self, r):
        return self.mapping[r]

    def get(self, r, default=None):
        return self.mapping.get(r, default)

    def items(self):
        return self.mapping.items()

    def values(self):
        return self.mapping.values()

    def materialize(self):
        return dict(self.mapping)

    def has_open_choice(self):
        return any(b.open for b in self.blocks)


def _cand_map(cand):
    return cand.materialize() if isinstance(cand, _SymCand) else dict(cand)


def _cand_possible_p_atoms(cand):
    atoms = set(_cand_map(cand).values())
    if isinstance(cand, _SymCand):
        for block in cand.blocks:
            atoms.update(block.p_atoms)
    return atoms


def _cand_has_open_choice(cand):
    return isinstance(cand, _SymCand) and cand.has_open_choice()


def _sym_block_assignment_expr(block):
    n = len(block.p_atoms)
    k = len(block.r_atoms)
    if k <= 0 or n <= 1:
        return '1'
    if k == 1:
        return str(n)
    if k == n:
        return f'{n}!'
    return f'P({n},{k})'


def _symmetry_state(cand):
    """Serialize only symmetry carried by the candidate itself."""
    item = {
        'witness': {int(a): int(b) for a, b in _cand_map(cand).items()},
        'blocks': [],
    }
    if isinstance(cand, _SymCand):
        item['compression_version'] = cand.compression_version
        item['exact_fixed'] = [int(x) for x in sorted(cand.exact_fixed)]
        if cand.target_generators is not None:
            item['automorph_generators'] = cand.target_generators
            item['automorph_group_source'] = 'conditioned_search_transition'
        if cand.compression_version < 3:
            item['multiplicity'] = int(cand.multiplicity)
            item['automorph_blocks'] = [{
                'r_atoms': [int(x) for x in block.r_atoms],
                'p_atoms': [int(x) for x in block.p_atoms],
                'extendable': False,
                'open': False,
                'assignments': 'exact_group',
                'source': 'exact_automorph_group',
            } for block in cand.automorph_blocks]
        for block in cand.blocks:
            item['blocks'].append({
                'r_atoms': [int(x) for x in block.r_atoms],
                'p_atoms': [int(x) for x in block.p_atoms],
                'extendable': bool(block.extendable),
                'open': bool(block.open),
                'assignments': _sym_block_assignment_expr(block),
            })
        if cand.compression_version < 3:
            item['blocks'].extend(item['automorph_blocks'])
    return item
