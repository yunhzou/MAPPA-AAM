"""Saved seed/cut windows must agree with direct public decoding of their union."""
import pytest

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
        for c in groups[:count]:
            offset=len(cat.families)
            cat.families.extend(c.families)
            for key,ids in c.branches.items():cat.branches.setdefault(key,set()).update(offset+i for i in ids)
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


def test_converted_record_does_not_emit_solver_evidence():
    problem=AAMProblem(endpoint('CO',[(0,1,1)]),endpoint('CO',[]))
    result=decode_events(catalogue(problem,{0:0,1:1}))
    assert result.complete
    assert all(r['solver_queries']==0 and 'query_evidence' not in r for r in result.family_reports)


def test_tasks_share_direction_budget_and_keep_seed_cut_identity():
    cut=Task('old','unused',dict(task_id=12,cut=3,source={},raw={},graph={},finalized=True,
                               family_seconds=10,cut_seconds=60))
    tasks=list(seed_tasks([cut],seeds=10))
    assert len({t.key for t in tasks})==10
    assert {t.case for t in tasks}=={'00012'}
    assert {t.payload['seed'] for t in tasks}==set(range(10))
    assert all(t.case_seconds==3600 for t in tasks)
    assert all('family_seconds' not in t.payload and 'cut_seconds' not in t.payload for t in tasks)
    assert cut.payload['family_seconds'] == 10


def test_saved_worker_checks_native_pairs_and_publishes_partial_mapping(monkeypatch,tmp_path):
    from types import SimpleNamespace
    from mappa.execution.runtime import Workspace
    from rebenchmark.decoding import decode_seed
    problem=AAMProblem(endpoint('CO',[(0,1,1)]),endpoint('CON',[(0,2,1)]))
    cat=catalogue(problem,{0:0,1:1})
    graph=cat.families[0].as_path(problem).graph
    monkeypatch.setattr('rebenchmark.decoding.prepare_cut',lambda payload,workspace:
        (SimpleNamespace(problem=problem),{},graph,'checked',None))
    workspace=Workspace(tmp_path,tmp_path/'bundles')
    value=decode_seed(dict(seed=0,slack=5,branch_limit=100,random_seed=42,decoder='anchored_fragment_conversion',task_id=0,cut=0,
        source={},raw={},graph={}),workspace)
    assert value['decoded']['complete'] and value['decoded']['observed_minimum']==2
    assert value['coverage']==(2,2) and value['search_rerun'] is False


def test_prefix_representatives_bound_later_seed_without_minimum_stage():
    from types import SimpleNamespace
    from rebenchmark.decoding import representative_bound
    problem=AAMProblem(endpoint('CCC',[(0,1,1),(1,2,1)]),
        endpoint('CCCCC',[(0,1,1),(1,2,1),(2,3,1),(3,4,1)]))
    maps=[{0:0,1:2,2:4},{0:0,1:1,2:2},{0:0,1:2,2:4}]
    graph=SimpleNamespace(terminals=range(3),states=[SimpleNamespace(context=i,mapping=tuple(m.items())) for i,m in enumerate(maps)])
    rows=[]
    bounds=[]
    for seed,mapping in enumerate(maps):
        cat=catalogue(problem,mapping)
        tier=coverage(problem,cat.families[0])
        bound=representative_bound(problem,graph,seed,tier,1)
        bounds.append(bound)
        result=decode_events(cat,EventDecodeConfig(max_events=bound))
        rows.append(dict(coverage=tier,decoded=decoded_events_record(result)))
        merged=merge_windows(rows,slack=1)
        direct=FinalBranchCatalogue(problem)
        for m in maps[:seed+1]:
            c=catalogue(problem,m)
            offset=len(direct.families)
            direct.families.extend(c.families)
            for key,ids in c.branches.items():direct.branches.setdefault(key,set()).update(offset+i for i in ids)
        expected=decode_events(direct,EventDecodeConfig(event_slack=1))
        assert merged['complete'] and merged['minimum']==expected.observed_minimum
        assert {r['id'] for r in merged['candidates']}=={r.id for r in expected.candidates}
    assert bounds[2]==bounds[1]<bounds[0]


def test_retired_solver_tasks_cannot_run_as_conversion(tmp_path):
    from rebenchmark.decoding import decode_seed
    with pytest.raises(ValueError, match='retired solver tasks'):
        decode_seed(dict(family_seconds=10), None)
