import numpy as np
import pytest
import mappa
from mappa.scoring import subset_score, POLICY
from mappa.domain import TSScoringConfig

def case():
    return np.array([[0.,0,0],[1.4,0,0],[5.,0,0],[6.4,0,0]]), ['C']*4, [dict(kind='broken',target_bond=(0,1)),dict(kind='broken',target_bond=(2,3))], [(0,1),(2,3)]

def test_async_one_active_event_remains_substantial():
    xyz,e,events,pairs=case();q=np.zeros((4,3));q[0,0]=-1;q[1,0]=1
    d=subset_score(xyz,q,e,events,pairs)
    assert d['score']==pytest.approx(.95)
    assert d['subset_strength']==1
    assert d['effective_event_count']==pytest.approx(1)
    assert subset_score(xyz,-4*q,e,events,pairs)['score']==pytest.approx(d['score'])

def test_competing_contact_motion_is_penalized():
    xyz,e,events,pairs=case();q=np.zeros((4,3));q[0,0]=-1;q[1,0]=1;q[2,0]=2;q[3,0]=-2
    assert subset_score(xyz,q,e,events[:1],pairs)['score'] < .5

def test_invalid_modes_fail_without_other_formula():
    xyz,e,events,pairs=case()
    for q in [np.zeros((4,3)),np.full((4,3),np.nan)]:
        with pytest.raises(ValueError):subset_score(xyz,q,e,events,pairs)
    with pytest.raises(ValueError):subset_score(xyz,np.ones((4,3)),e,[],pairs)

def test_old_scorers_and_configuration_are_removed():
    assert mappa.TS_SCORING_POLICY == POLICY == 'significant_subset_motion_v1'
    assert not hasattr(mappa,'bond_overlap_per_mode')
    assert not hasattr(mappa,'rxn_overlap_per_mode')
    with pytest.raises(TypeError):TSScoringConfig(wbo_progress_power=0)
    with pytest.raises(TypeError):TSScoringConfig(event_weight_power=1)

def test_saved_async_scores_are_preserved_exactly():
    import json
    from pathlib import Path
    cases=json.loads((Path(__file__).parent/'data/significant_subset_gold.json').read_text())
    for c in cases:
        result=subset_score(c['coordinates'],c['mode'],c['elements'],c['events'],c['pairs'])
        assert result['score']==pytest.approx(c['expected_score'],abs=1e-12)
