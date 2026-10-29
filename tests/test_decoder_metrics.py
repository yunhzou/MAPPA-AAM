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


def test_membership_prefixes_preserve_unknown_and_early_witness():
    from rebenchmark.decoder_metrics import membership_summary
    rows=[dict(key='a',seed=1,cut=0,status='not_recovered'),
          dict(key='b',seed=2,cut=1,status='recovered')]
    assert membership_summary(rows,['a'])['status']=='not_recovered'
    assert membership_summary(rows,['a','missing'])['status']=='unknown'
    assert membership_summary(rows,['a','b','missing'])['status']=='recovered'


def test_saved_bidirectional_membership_needs_no_search_or_conversion_rerun(tmp_path, monkeypatch):
    import hashlib,json
    from dataclasses import asdict
    from golden_evaluation import prepare
    from mappa import AAMProblem,AAMSearchConfig,search_aam
    from mappa.mappa_core.core.scheduler import BundleWriter
    from mappa.mappa_core.core.scheduler.runtime import Workspace
    from mappa.mappa_core.utils.io.records import endpoint_record
    from rebenchmark.decoding import decode_graph
    from rebenchmark.decoder_metrics import case_membership
    problem,features,reference=prepare('[CH3:1][OH:2]>>[CH3:1][OH:2]')
    writer=BundleWriter(tmp_path/'saved',tmp_path)
    payload=dict(index=0,expected=[],outputs=[],search_outputs=[],missing_tasks=[],
        directions={'0':'R_to_P','1':'P_to_R'})
    for name,value in [('input',dict(name='test',reactant=endpoint_record(problem.reactant),product=endpoint_record(problem.product))),
                       ('reference',dict(features=features,mapping=reference))]:
        path=tmp_path/(name+'.json');path.write_text(json.dumps(value,default=lambda a:a.tolist()))
        payload[name]=dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    for direction in (0,1):
        directed=AAMProblem(problem.product,problem.reactant) if direction else problem
        graph=search_aam(directed,AAMSearchConfig(seed_count=2,iso_tolerance=1.)).graph
        parent=f'{direction:05d}/cut/00000/finalize'
        payload['search_outputs'].append(dict(key=parent,blob=asdict(writer.put(dict(graph=graph)))))
        for seed in range(2):
            value=decode_graph(directed,graph,dict(seed=seed,slack=5,branch_limit=100,random_seed=42))
            value.update(task_id=direction,cut=0)
            from rebenchmark.check_conversion_outputs import verify_conversion
            checked=verify_conversion(directed,value,{},'test')
            assert checked['branches']==len(value['decoded']['converted_branches'])
            key=f'{direction:05d}/cut/00000/seed/{seed:02d}'
            payload['expected'].append(key)
            payload['outputs'].append(dict(key=key,blob=asdict(writer.put(value))))
    def forbidden(*args,**kwargs):raise AssertionError('Saved-result membership must not rerun science')
    monkeypatch.setattr('mappa.mappa_core.post_processing.decoder.api.decode_events',forbidden)
    monkeypatch.setattr('mappa.mappa_core.core.search.api.search_aam',forbidden)
    result=case_membership(payload,Workspace(tmp_path,tmp_path/'reports'))
    assert len(result['modes'])==80
    assert all(m[k]['status']=='recovered' for m in result['modes'].values() for k in ('aam','converted'))
    assert result['reads']==dict(aam=2,converted=2)
    assert all(r['seed']==1 for rows in result['rows'].values() for r in rows)
    assert all(r['witness']['input_orientation_mapping'] for rows in result['rows'].values() for r in rows)


def test_original_timing_keeps_batches_publication_and_checking_separate():
    from rebenchmark.timings import case_timing,distribution
    resource=dict(cpu_seconds=4.,compute_seconds=5.,peak_tree_rss_bytes=100,
                  publish_cpu_seconds=2.,publish_seconds=3.)
    result=case_timing(dict(conversion_checks=[dict(key='00000/cut/00000/seed/00',complete=True,
        timing=dict(decoding=dict(cpu_seconds=1.,wall_seconds=2.)),resource=resource)],
        search_resources=[dict(key='00000/cut/00000/seeds/00000',**resource),
                          dict(key='00000/cut/00000/finalize',**resource)],
        calculation=dict(cpu_seconds=999.)))
    assert result['phases']['conversion_decoding']['cpu']==[1.]
    assert result['phases']['matching_batch']['cpu']==[4.]
    assert result['phases']['conversion_publication']['wall']==[3.]
    assert result['prefixes']['uncut/R_to_P/seeds1']['matching_batch_cpu']==4.
    assert distribution([1.,3.])['median']==2.


def test_claim_concurrency_is_not_summed_latency(tmp_path):
    import sqlite3
    from rebenchmark.timings import execution_timeline
    path=tmp_path/'timing.sqlite';db=sqlite3.connect(path)
    db.execute('CREATE TABLE attempts(started REAL,finished REAL,state TEXT)')
    db.executemany('INSERT INTO attempts VALUES (?,?,?)',[(0.,2.,'complete'),(1.,3.,'complete'),(10.,11.,'complete')])
    db.commit();db.close()
    value=execution_timeline(path)
    assert value['elapsed_seconds']==11.
    assert value['claimed_task_seconds']==5.
    assert value['seconds_with_claimed_work']==4.
    assert value['peak_claimed_tasks']==2
