"""Event windows filter converted branches; they do not invoke minimization."""
from dataclasses import replace

import pytest

from mappa.mappa_core.post_processing.hierarchical_mapping import RefinementConfig
from mappa.mappa_core.post_processing.hierarchical_mapping import refine_mapping
from test_fragment_conversion import endpoint, saved


def catalogue():
    ep=endpoint('OCC',[(0,1,1.),(0,2,2.)])
    result=saved(ep,ep,{0:0,1:2,2:1},actions=(('pool',(1,2)),))
    result.families.append(replace(result.families[0],actions=()))
    next(iter(result.branches.values())).add(1)
    return result


@pytest.mark.parametrize('slack,expected', [(0,[0]),(1,[0]),(2,[0,2]),(5,[0,2])])
def test_window_filters_direct_branch_counts(slack, expected):
    result=refine_mapping(catalogue(),RefinementConfig(event_slack=slack))
    assert result.complete and result.minimum_proven
    assert result.observed_minimum==0 and result.event_window==slack
    assert [c.total for c in result.candidates]==expected
    assert sorted(r['total'] for r in result.family_reports)==[0,2]


def test_absolute_window_and_empty_catalogue():
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    cat=catalogue()
    assert [c.total for c in refine_mapping(cat,RefinementConfig(max_events=1)).candidates]==[0]
    empty=refine_mapping(FinalBranchCatalogue(cat.problem),RefinementConfig(event_slack=5))
    assert empty.observed_minimum is None and not empty.minimum_proven


def test_invalid_or_ambiguous_configuration_rejected():
    for kwargs in ({'event_slack':-1},{'event_slack':1.5},{'event_slack':5,'max_events':9},
                   {'branch_limit':0},{'threshold':0},{'random_seed':.1}):
        with pytest.raises(ValueError):RefinementConfig(**kwargs)


def test_event_ids_are_only_built_for_the_requested_window(monkeypatch):
    from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex
    original = SignedEventIndex._describe_events
    observed = []
    def describe(index, vector, broken, formed):
        count = len(broken) + len(formed)
        observed.append(count)
        return original(index, vector, broken, formed)
    monkeypatch.setattr(SignedEventIndex, '_describe_events', describe)
    result = refine_mapping(catalogue(), RefinementConfig(event_slack=0))
    assert sorted(r['total'] for r in result.family_reports) == [0,2]
    assert observed == [0]
