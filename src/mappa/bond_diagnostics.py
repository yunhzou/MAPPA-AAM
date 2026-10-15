"""Descriptive endpoint-to-target WBO changes; never a scoring gate."""
from itertools import combinations
import numpy as np
from .frag import bond_event_threshold


def endpoint_bond_changes(reaction, mechanism, target, assignment):
    r, p = reaction.problem.reactant, reaction.problem.product
    t = target.molecule
    core_mapping = assignment.as_dict()
    declared = target.reactant_atom_mapping
    if declared is not None:
        images = tuple(declared)
        if (len(images) != r.atom_count or t.atom_count != r.atom_count
                or set(images) != set(range(t.atom_count))):
            raise ValueError('diagnostic reactant_atom_mapping must be a full bijection')
        mapping = dict(enumerate(images))
        if any(r.elements[i] != t.elements[j] for i, j in mapping.items()):
            raise ValueError('diagnostic atom mapping changes element identity')
        agreement = all(mapping[i] == j for i, j in core_mapping.items())
    else:
        mapping = core_mapping
        agreement = True
    rp = mechanism.mapping.as_dict()
    intended = {tuple(sorted(pair)) for pair in (*mechanism.broken_bonds, *mechanism.formed_bonds)}
    core = set(mechanism.core_atoms)
    reports = {}
    for name, endpoint, to_endpoint in [('reactant_to_ts', r, {i:i for i in mapping}),
                                       ('product_to_ts', p, rp)]:
        changes = []
        for i, j in combinations(sorted(mapping), 2):
            a, b = to_endpoint[i], to_endpoint[j]
            x, y = mapping[i], mapping[j]
            before, after = float(endpoint.wbo[a,b]), float(t.wbo[x,y])
            threshold = bond_event_threshold(r.elements, i, j,
                default_threshold=reaction.config.event_threshold,
                metal_threshold=reaction.config.metal_event_threshold)
            delta = after-before
            if abs(delta) < threshold:
                continue
            changes.append(dict(kind='forming' if delta > 0 else 'breaking',
                reactant_atoms=[i,j], endpoint_atoms=[a,b], ts_atoms=[x,y],
                elements=[r.elements[i],r.elements[j]], wbo_endpoint=before, wbo_ts=after,
                delta_wbo=delta, threshold=threshold,
                distance_endpoint_angstrom=float(np.linalg.norm(endpoint.coordinates[a]-endpoint.coordinates[b])),
                distance_ts_angstrom=float(np.linalg.norm(t.coordinates[x]-t.coordinates[y])),
                intended_reaction_pair=(i,j) in intended,
                touches_outside_reaction_core=not {i,j}.issubset(core)))
        reports[name] = changes
    return dict(coverage='all_atom_pairs' if declared is not None else 'mapped_core_pairs_only',
                mapping_source='explicit_atom_identity' if declared is not None else 'core_aam',
                core_assignment_agrees_with_atom_identity=agreement,
                outside_core_assessed=declared is not None, atom_indices='zero_based',
                interpretation='Directional WBO increases/decreases above existing reaction thresholds; not proof of a formed intermediate or an IRC match.',
                **reports)
