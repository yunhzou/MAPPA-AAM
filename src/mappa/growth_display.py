"""Display domains derived from recorded candidates, never global atom orbits."""
from copy import deepcopy
import json

from .frag import build_graph
from .matcher.state import candidate_from_record
from .search_symmetry import SymmetryWorkspace


def with_represented_targets(document):
    """Annotate per-atom marginal domains without enumerating full mappings.

    Local pools act first, then the conditioned candidate group, as in the
    family decoder. Marginal domains must never be shown as independently
    selectable assignments. Raw automorph-block unions are not target domains.
    """
    import numpy as np
    doc = deepcopy(document)
    ep = doc['input']['product']
    cache = {}
    for run in doc['runs']:
        cfg = run['config']
        target = build_graph(ep['elements'], np.asarray(ep['wbo']), bond_cut=cfg['graph_floor'])
        workspace = SymmetryWorkspace(target, cfg['iso_tolerance'])
        for path in run['paths']:
            for frame in path['frames']:
                for record in frame['candidates']:
                    if ((record.get('automorph_domains', 0) and 'automorph_blocks' not in record)
                            or (frame['kind'] == 'seed_end' and not record.get('symmetry_complete'))):
                        record['target_options'] = None
                        continue
                    key = json.dumps([cfg['graph_floor'], cfg['iso_tolerance'], frame['locked'], record], sort_keys=True)
                    if key not in cache:
                        cand = candidate_from_record(record)
                        locked = tuple(sorted((int(r), int(p)) for r, p in frame['locked'].items()))
                        generators = workspace.conditioned_generators(cand, locked)
                        block_r = {r for b in cand.blocks for r in b.r_atoms}
                        fixed = {p for _, p in locked} | {
                            p for r, p in cand.mapping.items() if r not in block_r or r in cand.exact_fixed}
                        options = {}
                        for r, p in cand.mapping.items():
                            domain = {p}
                            for block in cand.blocks:
                                free = set(block.p_atoms) - fixed
                                if domain & free:
                                    domain |= free
                            queue = list(domain)
                            while queue:
                                atom = queue.pop()
                                for gen in generators:
                                    image = gen[atom]
                                    if image not in domain:
                                        domain.add(image)
                                        queue.append(image)
                            options[str(r)] = sorted(domain)
                        cache[key] = options
                    record['target_options'] = cache[key]
    return doc
