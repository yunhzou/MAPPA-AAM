"""Viewer must distinguish a compressed relation from its witness."""
from copy import deepcopy

from mappa.growth_display import with_represented_targets


def test_target_options_follow_saved_groups_not_global_symmetry():
    # Two identical O-H-C components. A fixed newly grown carbon restricts
    # the recorded group; the viewer must display that actual restriction.
    wbo = [[0.0]*6 for _ in range(6)]
    for a, b in [(0, 1), (0, 2), (3, 4), (3, 5)]:
        wbo[a][b] = wbo[b][a] = 1.0
    pair = dict(witness={0: 3, 1: 4}, blocks=[],
                automorph_blocks=[dict(r_atoms=[0, 1], p_atoms=[0, 1, 3, 4])])
    triple = deepcopy(pair); triple['witness'][2] = 5
    seed = dict(witness={0: 0}, blocks=[dict(r_atoms=[0], p_atoms=[0, 3])])
    doc = dict(input=dict(product=dict(elements=['O', 'H', 'C']*2, wbo=wbo)),
               runs=[dict(config=dict(graph_floor=.2, iso_tolerance=1), paths=[dict(frames=[
                   dict(kind='seed_start', locked={}, candidates=[seed]),
                   dict(kind='commit', locked={}, candidates=[pair]),
                   dict(kind='commit', locked={}, candidates=[triple]),
                   dict(kind='seed_end', locked={}, candidates=[dict(witness={0: 3}, blocks=[])])])])])
    original = deepcopy(doc)
    frames = with_represented_targets(doc)['runs'][0]['paths'][0]['frames']
    assert frames[0]['candidates'][0]['target_options']['0'] == [0, 3]
    assert frames[1]['candidates'][0]['target_options'] == {'0': [0, 3], '1': [1, 4]}
    assert frames[2]['candidates'][0]['target_options'] == {'0': [3], '1': [4], '2': [5]}
    assert frames[3]['candidates'][0]['target_options'] is None
    assert doc == original
