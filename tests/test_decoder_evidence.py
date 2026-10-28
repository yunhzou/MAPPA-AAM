from copy import deepcopy
import hashlib
import pytest
from mappa import AAMProblem
from mappa.artifacts import decoded_events_record
from legacy_decoder_records import record as legacy_record
from rebenchmark.decoder_evidence import ProofReplay, validate_record
from test_partial_event_decoding import endpoint, catalogue


def example(symbolic=False):
    problem = AAMProblem(endpoint('CCC', [(0,1,1),(1,2,1)]),
                        endpoint('CCCC', [(0,1,1),(1,2,.6),(2,3,1.8),(0,3,.4)]))
    cat = catalogue(problem, {0:0,2:1}, range(4) if symbolic else ())
    queries = []
    def log(text, status, seconds):
        queries.append(dict(smt2=text,status=status,sha256=hashlib.sha256(text.encode()).hexdigest()))
        return len(queries)-1
    result = legacy_record(cat, max_events=10, query_log=log)
    return problem, dict(decoded=result, families=[f.to_record() for f in cat.families],selected_families=1), queries


@pytest.mark.parametrize('symbolic', [False, True])
def test_saved_candidates_and_complete_proofs_are_checked(symbolic):
    p,v,q = example(symbolic)
    result = validate_record(v,p,q.__getitem__,ProofReplay())
    assert result['verified'] and result['decoded_complete']
    assert result['candidates'] == len(v['decoded']['candidates'])
    if symbolic:
        assert q and result['proof_checks']


def test_changed_event_count_is_rejected():
    p,v,q = example()
    v['decoded']['candidates'][0]['total'] += 1
    with pytest.raises(ValueError,match='witness'):
        validate_record(v,p,q.__getitem__,ProofReplay())


def test_missing_query_or_forged_completeness_is_rejected():
    p,v,q = example(True)
    bad = deepcopy(v);bad['decoded']['family_reports'][0]['query_evidence'] = []
    with pytest.raises(ValueError,match='missing'):
        validate_record(bad,p,q.__getitem__,ProofReplay())
    bad = deepcopy(v);bad['decoded']['family_reports'][0]['reason'] = 'time_budget'
    with pytest.raises(ValueError,match='certificate'):
        validate_record(bad,p,q.__getitem__,ProofReplay())


def test_incomplete_decoding_is_not_promoted():
    p,v,q = example()
    v['decoded']['family_reports'][0].update(complete=False,reason='time_budget')
    v['decoded'].update(complete=False,minimum_proven=False)
    result = validate_record(v,p,q.__getitem__,ProofReplay())
    assert result['verified'] and not result['decoded_complete']
    assert not result['minimum_proven']


def test_replay_cache_checks_expected_outcome():
    replay = ProofReplay()
    text='(declare-const a Int) (assert (= a 1))'
    assert replay.check(text,'sat')['status']=='sat'
    with pytest.raises(ValueError,match='Conflicting'):
        replay.check(text,'unsat')


def test_inconclusive_replay_leaves_audit_unresolved():
    class Unresolved:
        def check(self, text, expected):
            return dict(status='unknown',expected=expected,reason='audit budget')
    p,v,q = example(True)
    result = validate_record(v,p,q.__getitem__,Unresolved())
    assert result['decoded_complete']
    assert not result['verified'] and result['unresolved_proof_checks'] > 0


def test_original_task_metadata_is_not_reconstructed_from_output():
    from dataclasses import replace
    from mappa.execution import Task
    from rebenchmark.decoder_evidence import verify_task
    from rebenchmark.decoder_inventory import signature
    _,v,_ = example()
    v.update(task_id=0,cut=0,seed=1,search_rerun=False,source={},graph={},raw={},finalized=None,
             window_basis=dict(slack=5))
    t = Task('00000/cut/00000/seed/00','rebenchmark.decoding:decode_seed',dict(
        task_id=0,cut=0,seed=0,source={},graph={},raw={},finalized=True,slack=5,
        family_seconds=30,cut_seconds=3000.),case='00000')
    assert verify_task(v,t,signature(t)) == t.key
    with pytest.raises(ValueError,match='identity'):
        verify_task(v,replace(t,payload=dict(t.payload,cut_seconds=2999.)),signature(t))
    v['window_basis']['slack']=6
    with pytest.raises(ValueError,match='scientific'):
        verify_task(v,t,signature(t))


def test_relative_pilot_window_has_explicit_policy_and_same_proofs():
    from rebenchmark.decoder_evidence import window_slack, verify_window
    problem = AAMProblem(endpoint('CCC', [(0,1,1),(1,2,1)]),
                        endpoint('CCCC', [(0,1,1),(1,2,.6),(2,3,1.8),(0,3,.4)]))
    cat = catalogue(problem, {0:0,2:1}, range(4))
    queries=[]
    def log(text,status,seconds):
        queries.append(dict(smt2=text,status=status,sha256=hashlib.sha256(text.encode()).hexdigest()))
        return len(queries)-1
    result=legacy_record(cat,event_slack=5,query_log=log)
    value=dict(decoded=result,families=[f.to_record() for f in cat.families],selected_families=1)
    assert window_slack(value)==5
    verify_window(value,problem,None,None)
    assert validate_record(value,problem,queries.__getitem__,ProofReplay())['verified']
    changed=deepcopy(value);changed['decoded']['event_window']+=1
    with pytest.raises(ValueError,match='Relative'):
        verify_window(changed,problem,None,None)
    changed=deepcopy(value);changed['window_basis']=dict(slack=5)
    with pytest.raises(ValueError,match='conflicting'):
        window_slack(changed)
