"""Viewer must distinguish a completed relation from its growth witness."""
from copy import deepcopy

import numpy as np

from mappa import AAMProblem, MolecularEndpoint
from mappa.presentation.trajectory import with_represented_targets
from growth_cases import saved_fragment


def test_target_options_follow_local_relationship_and_fixed_anchors():
    # An external chlorine distinguishes the target hydrogens globally. Inside
    # the selected H-C-H fragment they interchange until an anchor fixes one.
    source = np.zeros((3, 3))
    target = np.zeros((4, 4))
    for a, b in ((0, 1), (1, 2)):
        source[a, b] = source[b, a] = target[a, b] = target[b, a] = 1.
    target[0, 3] = target[3, 0] = 1.
    problem = AAMProblem(
        MolecularEndpoint(('H', 'C', 'H'), np.zeros((3, 3)), source),
        MolecularEndpoint(('H', 'C', 'H', 'Cl'), np.zeros((4, 3)), target))
    free = saved_fragment(problem, accepted_edges=((0, 1), (1, 2)))
    anchored = saved_fragment(problem, anchors=((0, 0),),
                             accepted_edges=((0, 1), (1, 2)))
    candidates = [dict(witness={1: 1}),
                  dict(witness=dict(free.witness), relation=free.to_record()),
                  dict(witness=dict(anchored.witness), relation=anchored.to_record())]
    document = dict(runs=[dict(paths=[dict(frames=[dict(candidates=[candidate])
                                                    for candidate in candidates])])])
    original = deepcopy(document)
    frames = with_represented_targets(document)['runs'][0]['paths'][0]['frames']
    assert frames[0]['candidates'][0]['target_options'] == {'1': [1]}
    assert frames[1]['candidates'][0]['target_options'] == {
        '0': [0, 2], '1': [1], '2': [0, 2]}
    assert frames[2]['candidates'][0]['target_options'] == {
        '0': [0], '1': [1], '2': [2]}
    assert document == original
