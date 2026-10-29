from copy import deepcopy
import pytest
from rebenchmark.decoder_metrics import summarize, feature_variants

REF={k:'ref' for k in ['strict','without_stereo','without_charge','without_both']}

def row(key='a', minimum=2, bound=7, complete=True, reference=True, tier=(2,3)):
    candidate=dict(event_id='event'+str(minimum),total=minimum,**{k:('ref' if reference else 'other') for k in REF})
    return dict(key=key,coverage=tier,candidates=[candidate],decoded_complete=complete,
                catalogue_complete=True,bound=bound)


def test_missing_output_prevents_minimum_and_tier_claim():
    result=summarize([row(minimum=0)],['a','missing'],REF)
    assert result['observed_minimum']==0
    assert not result['minimum_proven'] and not result['coverage_tier_known']
    assert not result['decoded_window_complete'] and result['missing_rows']==['missing']


def test_bounds_gate_minimum_and_full_window_separately():
    result=summarize([row(bound=3)],['a'],REF)
    assert result['minimum_proven']
    assert not result['decoded_window_complete'] and result['insufficient_window_rows']==['a']
    result=summarize([row(bound=1)],['a'],REF)
    assert not result['minimum_proven']


def test_incomplete_catalogue_cannot_certify_zero_at_unknown_tier():
    value=row(minimum=0);value['catalogue_complete']=False
    result=summarize([value],['a'],REF)
    assert not result['minimum_proven'] and not result['coverage_tier_known']


def test_best_tier_and_prefix_union_keep_saved_empty_windows():
    first=row('first',minimum=2)
    empty=row('empty',bound=7);empty['candidates']=[]
    lower=row('lower',minimum=0,tier=(1,1))
    result=summarize([first,empty,lower],['first','empty','lower'],REF)
    assert result['minimum_proven'] and result['decoded_window_complete']
    assert result['observed_minimum']==2 and result['selected_rows']==2
    assert result['event_class_count']==1


def test_one_event_can_have_several_strict_representatives():
    first=row();second=row('b',reference=False)
    result=summarize([first,second],['a','b'],REF)
    assert result['event_class_count']==1 and result['strict_representative_classes']==2
    assert result['representative_hits']['strict']


def test_no_reference_promotion_from_relaxed_labels():
    value=row(reference=False);value['candidates'][0]['without_stereo']='ref'
    result=summarize([value],['a'],REF)
    assert result['representative_hits']['without_stereo']
    assert not result['representative_hits']['strict']
    assert result['reference_compatible_witness_offset'] is None


def test_repeated_or_unexpected_outputs_rejected():
    with pytest.raises(ValueError,match='Duplicate'):summarize([row(),row()],['a'],REF)
    with pytest.raises(ValueError,match='unexpected'):summarize([row()],['b'],REF)


def test_relax_only_named_stereo_and_charge_features():
    f=[dict(heavy=[0],colors=[('C',1,13,3,'R')],bonds=[(0,1,(1,'E'))])]
    v=feature_variants(f)
    assert v['without_stereo'][0]['colors'][0]==['C',1,13,3,'']
    assert v['without_charge'][0]['colors'][0]==['C',0,13,3,'R']
    assert v['without_charge'][0]['bonds']==f[0]['bonds']
    assert f[0]['colors'][0][-1]=='R'


def test_reverse_concrete_witness_normalizes_before_event_union():
    from mappa import AAMProblem
    from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex
    from test_partial_event_decoding import endpoint
    problem=AAMProblem(endpoint('CCC',[(0,1,1)]),endpoint('CCC',[(1,2,1)]))
    forward={0:1,1:0,2:2};reverse={b:a for a,b in forward.items()}
    index=SignedEventIndex(problem)
    expected=index.describe(index.vector(forward))
    actual=index.describe(index.vector({b:a for a,b in reverse.items()}))
    native=SignedEventIndex(AAMProblem(problem.product,problem.reactant)).describe(
        SignedEventIndex(AAMProblem(problem.product,problem.reactant)).vector(reverse))
    assert actual['id']==expected['id'] and native['total']==expected['total']
