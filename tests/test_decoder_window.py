"""The window uses local constraints and only the selected sweep coupling."""
import pytest

from mappa import AAMProblem, decode_events
from mappa.mappa_core.post_processing.decoder.fragments import FragmentDecoder
from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
from mappa.mappa_core.utils.io.records import decoded_events_record
from growth_cases import saved_fragment
from test_event_interpretation import endpoint, saved, swap


def test_boundary_bound_rejects_composition_before_local_preparation(monkeypatch):
    ep = endpoint('CCCC', [(0, 1, 1.), (1, 2, 1.), (2, 3, 1.)])
    source = saved(ep, ep, dict(enumerate(range(4))))
    source.coupled_branches[0] = {tuple(((a,), (a,)) for a in range(4))}
    prepare = FragmentDecoder.fragment
    calls = []
    def checked(self, piece, *args):
        calls.append(piece.source_atoms)
        assert len(piece.source_atoms) == 4, 'fixed boundary already excludes singleton partition'
        return prepare(self, piece, *args)
    monkeypatch.setattr(FragmentDecoder, 'fragment', checked)
    result = decode_events(source, event_slack=0)
    record = decoded_events_record(result)
    assert result.minimum_events == result.event_window == 0
    assert len(record['tables']) == 1
    assert record['tables'][0]['fragments'] == ((0, 1, 2, 3),)
    assert len(result.local.fragments) == 1 and calls
    assert record['catalogue'] == source.to_record()


def test_sweep_resolves_only_its_endpoint_fragments_and_preserves_count_gaps():
    r = endpoint('CCC', [(0, 1, 1.)])
    p = endpoint('CCC', [(0, 2, 1.)])
    fragment = saved_fragment(AAMProblem(r, p), source_atoms=(1, 2))
    source = saved(r, p, {0: 0, 1: 1, 2: 2}, fragments=(fragment,),
                   pairs=(((0,), (0,)), ((1, 2), (1, 2))))
    source.families[0].provenance = [dict(cuts=((0, 1),))]
    result = decode_events(source, event_slack=1)
    assert result.minimum_events == 0 and result.event_window == 1
    original = result.contributions(0)
    moved = result.contributions(0, [swap(fragment)])
    assert (original['total'], moved['total']) == (2, 0)
    assert not original['within_window'] and moved['within_window']
    assert [row.sum() for row in result.counts(0, [None, [swap(fragment)]])] == [2, 0]
    record = decoded_events_record(result)
    assert record['tables'][0]['event_window'] == 1
    assert record['minimum_scope'] == 'independent_fragment_model'
    assert record['event_counts_approximate'] and not record['global_minimum_proven']


def test_separate_sweep_histories_do_not_become_a_two_cut_query():
    ep = endpoint('CCC', [(0, 1, 1.), (1, 2, 1.)])
    source = saved(ep, ep, {0: 0, 1: 1, 2: 2},
                   pairs=(((0,), (0,)), ((1,), (1,)), ((2,), (2,))))
    source.families[0].provenance = [dict(cuts=((0, 1),)), dict(cuts=((1, 2),))]
    result = decode_events(source)
    assert result.cut_contexts == (((((0, 1),), ((1, 2),)),),)
    assert result.minimum_events == 2
    assert result.tables(0, cut=0).cut_bonds == frozenset(((0, 1),))
    assert result.tables(0, cut=1).cut_bonds == frozenset(((1, 2),))
    assert len(decoded_events_record(result)['tables']) == 2


def test_ordinary_boundary_approximation_is_separate_from_physical_report():
    ep = endpoint('CC', [(0, 1, 1.)])
    source = saved(ep, ep, {0: 0, 1: 1}, pairs=(((0,), (0,)), ((1,), (1,))))
    result = decode_events(source)
    assert result.minimum_events == result.contributions(0)['total'] == 2
    assert result.evaluate(0)['total'] == 0
    assert result.realize(0) == {0: 0, 1: 1}


def test_empty_catalogue_has_no_minimum_or_window():
    ep = endpoint('C', [])
    result = decode_events(FinalBranchCatalogue(AAMProblem(ep, ep)))
    assert result.minimum_events is None and result.event_window is None
    assert decoded_events_record(result)['tables'] == []


def test_fixed_sweep_response_needs_no_joint_fragment_query(monkeypatch):
    from mappa.mappa_core.post_processing.decoder import local_counts
    r, p = endpoint('CCCC', [(0, 2, 1.)]), endpoint('CCCC', [])
    problem = AAMProblem(r, p)
    fragments = tuple(saved_fragment(problem, source_atoms=atoms) for atoms in ((0, 1), (2, 3)))
    source = saved(r, p, dict(enumerate(range(4))), fragments=fragments,
                   pairs=(((0, 1), (0, 1)), ((2, 3), (2, 3))))
    source.families[0].provenance = [dict(cuts=((0, 2),))]
    def forbidden(*args):
        pytest.fail('the sweep costs one event for every allowed endpoint pair')
    monkeypatch.setattr(local_counts, 'compile_local', forbidden)
    result = decode_events(source)
    assert result.minimum_events == 1
    for choice in ([None, None], [swap(fragment) for fragment in fragments]):
        assert result.contributions(0, choice)['total'] == 1


def test_partition_provenance_keeps_its_own_sweep_context():
    ep = endpoint('CCC', [(0, 1, 1.), (1, 2, 1.)])
    source = saved(ep, ep, {0: 0, 1: 1, 2: 2},
                   pairs=(((0,), (0,)), ((1, 2), (1, 2))))
    source.coupled_branches[0] = {(((0, 1), (0, 1)), ((2,), (2,)))}
    source.families[0].provenance = [
        dict(cuts=((0, 1),), fragment_partition=((0,), (1, 2))),
        dict(cuts=((1, 2),), fragment_partition=((0, 1), (2,)))]
    result = decode_events(source)
    assert result.cut_contexts == (((((0, 1),),), (((1, 2),),)),)
    assert len(decoded_events_record(result)['tables']) == 2


@pytest.mark.parametrize('slack', [-1, 1.5, None, True])
def test_window_has_one_nonnegative_integer_setting(slack):
    with pytest.raises(ValueError, match='event_slack'):
        EventDecodeConfig(event_slack=slack)
