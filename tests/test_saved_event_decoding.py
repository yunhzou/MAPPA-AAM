"""Saved seed/cut windows must agree with direct public decoding of their union."""
import z3

from mappa.artifacts import decoded_events_record
from mappa.execution import Task
from mappa.postprocessing import decode_events, EventDecodeConfig
from mappa.final_branches import FinalBranchCatalogue
from rebenchmark.decoding import seed_tasks, merge_windows, coverage
from test_partial_event_decoding import endpoint, catalogue
from mappa import AAMProblem


def test_each_seed_window_preserves_earlier_prefixes_and_cut_unions():
    problem = AAMProblem(endpoint('CCC', [(0,1,1),(1,2,1)]),
                         endpoint('CCCCC', [(0,1,1),(1,2,1),(2,3,1),(3,4,1)]))
    groups = [catalogue(problem, m) for m in ({0:0,1:2,2:4}, {0:0,1:1,2:2}, {0:2,1:3,2:4})]
    rows = [dict(coverage=coverage(problem,c.families[0]),
                 decoded=decoded_events_record(decode_events(c,EventDecodeConfig(event_slack=1)))) for c in groups]
    assert rows[0]['decoded']['observed_minimum'] > rows[1]['decoded']['observed_minimum']
    for count in (1,2,3):
        cat=FinalBranchCatalogue(problem)
        cat.families=[f for c in groups[:count] for f in c.families]
        direct=decode_events(cat,EventDecodeConfig(event_slack=1))
        merged=merge_windows(rows[:count],slack=1)
        assert merged['complete'] and merged['minimum_proven']
        assert merged['minimum']==direct.observed_minimum
        assert {c['id'] for c in merged['candidates']}=={c.id for c in direct.candidates}


def test_original_coverage_tier_and_unknown_are_preserved():
    problem=AAMProblem(endpoint('CC',[(0,1,1)]),endpoint('CC',[]))
    full=catalogue(problem,{0:0,1:1})
    empty=catalogue(problem,{})
    rows=[dict(coverage=coverage(problem,c.families[0]),
         decoded=decoded_events_record(decode_events(c,EventDecodeConfig(event_slack=0)))) for c in (empty,full)]
    assert merge_windows(rows,slack=0)['minimum']==1
    rows[1]['decoded']['complete']=False
    assert not merge_windows(rows,slack=0)['complete']
    assert not merge_windows(rows,slack=0)['minimum_proven']
    assert not merge_windows([],slack=0)['complete']


def test_public_decoder_query_evidence_replays():
    problem=AAMProblem(endpoint('CCC',[(0,1,1),(1,2,1)]),
        endpoint('CCCC',[(0,1,1),(1,2,.6),(2,3,1.8),(0,3,.4)]))
    cat=catalogue(problem,{0:0,2:1},range(4))
    evidence=[]
    def log(text,status,seconds):
        evidence.append((text,status,seconds))
        return len(evidence)-1
    result=decode_events(cat,EventDecodeConfig(event_slack=1,seconds_per_family=30),query_log=log)
    assert result.complete and evidence
    assert result.family_reports[0]['query_evidence']==list(range(len(evidence)))
    for text,status,seconds in evidence:
        solver=z3.Solver();solver.from_string(text)
        assert str(solver.check())==status and seconds>=0


def test_tasks_share_direction_budget_and_keep_seed_cut_identity():
    cut=Task('old','unused',dict(task_id=12,cut=3))
    tasks=list(seed_tasks([cut],seeds=10))
    assert len({t.key for t in tasks})==10
    assert {t.case for t in tasks}=={'00012'}
    assert {t.payload['seed'] for t in tasks}==set(range(10))
    assert all(t.case_seconds==3600 for t in tasks)
