"""Pipeline validation: identical graphs, independent prefixes and durable evidence."""
from dataclasses import asdict,replace
import json
from pathlib import Path
import sys

import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bench'))
from rebenchmark.evidence import save,read,sha,publish,restore
from rebenchmark.search import search,select_contexts
from golden_evaluation import prepare
from mappa import AAMSearchConfig,AAMProblem,search_aam_checkpoints
from mappa.artifacts import read_raw_cut


@pytest.mark.parametrize('reverse',[False,True])
def test_parallel_prefixes_and_resume(tmp_path,reverse,monkeypatch):
    problem,features,reference=prepare('[CH3:1][OH:2]>>[CH3:1][OH:2].[Cl-]')
    if reverse:problem=AAMProblem(problem.product,problem.reactant)
    config=AAMSearchConfig(seed_count=10,branch_limit=100,iso_tolerance=1)
    parallel=search(problem,config,tmp_path/'parallel',2)
    serial=search_aam_checkpoints(problem,config,intermediate_dir=tmp_path/'serial',workers=1)
    for a,b in zip(parallel.cut_paths,serial.cut_paths):
        assert read_raw_cut(a)==read_raw_cut(b)
    rows=[json.loads(line) for p in (tmp_path/'parallel/seed-timing').glob('*.jsonl') for line in p.read_text().splitlines()]
    assert len(rows)==10*len(parallel.cut_paths)
    for n in [1,2,3]:
        prefix=search_aam_checkpoints(problem,replace(config,seed_count=n),intermediate_dir=tmp_path/f'prefix{n}')
        for a,b in zip(parallel.cut_paths,prefix.cut_paths):
            assert select_contexts(read_raw_cut(a),range(n))==read_raw_cut(b)
    import mappa.aam as core
    monkeypatch.setattr(core,'_search_cut',lambda *a:pytest.fail('Completed cut searched again'))
    resumed=search_aam_checkpoints(problem,config,intermediate_dir=tmp_path/'serial',resume=True)
    assert resumed.metrics.retained_branch_count==serial.metrics.retained_branch_count


def test_archive_integrity_and_corruption(tmp_path):
    source=tmp_path/'source';save(source/'a.json',{'a':[1,2]})
    receipt=publish(source,tmp_path/'case.tar')
    assert receipt['files']==1
    restore(tmp_path/'case.tar',tmp_path/'restored')
    assert sha(source/'a.json')==sha(tmp_path/'restored/a.json')
    with pytest.raises(AssertionError):restore(tmp_path/'case.tar',tmp_path/'restored')


def test_logged_resume_searches_only_missing_cut(tmp_path,monkeypatch):
    problem,_,_=prepare('[CH3:1][OH:2]>>[CH3:1][OH:2]')
    config=AAMSearchConfig(seed_count=2,branch_limit=100)
    full=search(problem,config,tmp_path,1)
    expected={p.name:read_raw_cut(p) for p in full.cut_paths}
    missing=full.cut_paths[-1];missing.unlink()
    from mappa import aam as core
    original=core._search_cut;calls=[]
    def observed(cut):
        calls.append(cut)
        return original(cut)
    monkeypatch.setattr(core,'_search_cut',observed)
    result=search(problem,config,tmp_path,1)
    assert len(calls)==1
    assert read(tmp_path/'search.json')['resumed_cuts']==len(full.cut_paths)-1
    assert {p.name:read_raw_cut(p) for p in result.cut_paths}==expected


def test_worker_continuation_preserves_inputs_and_prior_budget(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from rebenchmark import campaign
    run=tmp_path/'run';seed=tmp_path/'seed'
    task=dict(index=0,direction='R_to_P',config={})
    for name in ('input.json','reference.json','source.json'):
        save(run/'inputs/0'/name,{})
        save(seed/name,{})
    save(seed/'task.json',task);save(seed/'cuts/manifest.json',{})
    save(seed/'cuts/saved.json',{'unchanged':True})
    receipt=publish(seed,tmp_path/'seed.tar')
    budgets=[]
    monkeypatch.setattr(campaign,'CASE_SECONDS',10)
    monkeypatch.setattr(campaign.time,'perf_counter',lambda:100.0)
    def stage(command,folder,name,seconds,memory):
        assert read(folder/'cuts/saved.json')=={'unchanged':True}
        budgets.append(seconds)
        return dict(status='complete')
    monkeypatch.setattr(campaign,'guarded',stage)
    manifest=dict(tasks=[task],workers=1,stage_seconds=dict(search=100,analysis=100),
        memory_mib=100,resume_archives={'0':dict(path=receipt['path'],sha256=receipt['sha256'],charged_seconds=3)})
    campaign.worker(SimpleNamespace(run=run,task=0),manifest)
    assert budgets==[7,7]
    state=read(run/'progress/00000.json')
    assert state['complete'] and state['prior_charged_seconds']==3
    restore(run/'archives/00000.tar',tmp_path/'restored')
    assert read(tmp_path/'restored/continuation.json')['sha256']==receipt['sha256']


@pytest.mark.parametrize('reverse',[False,True])
def test_full_analysis_keeps_seed_evidence(tmp_path,reverse):
    from rebenchmark.analysis import analyze
    original,features,reference=prepare('[CH3:1][OH:2]>>[CH3:1][OH:2].[Cl-]')
    config=AAMSearchConfig(seed_count=2,branch_limit=100)
    def ep(e):return dict(elements=e.elements,coordinates=e.coordinates.tolist(),wbo=e.wbo.tolist(),label=e.label,metadata=dict(e.metadata))
    save(tmp_path/'input.json',dict(name='control',reactant=ep(original.reactant),product=ep(original.product)))
    save(tmp_path/'reference.json',dict(features=features,mapping=sorted(reference.items())))
    save(tmp_path/'task.json',dict(direction='P_to_R' if reverse else 'R_to_P',config=asdict(config)))
    problem=AAMProblem(original.product,original.reactant) if reverse else original
    result=search(problem,config,tmp_path,2)
    analyze(tmp_path,2)
    records=list((tmp_path/'analysis').glob('*/*/evaluation.json'))
    assert len(records)==2*len(result.cut_paths)
    assert all(read(p)['reference_recovery']=='recovered' for p in records)
    assert all(read(p)['raw_sha256'] for p in records)


def test_saved_exact_queries_replay(tmp_path):
    import z3
    from rebenchmark.analysis import exact_query
    from mappa import search_aam
    problem,features,reference=prepare('[CH3:1][OH:2]>>[CH3:1][OH:2]')
    result=search_aam(problem,AAMSearchConfig(seed_count=1,sweep_cuts=False),execution='reused_native')
    path=next(result.graph.paths())
    mapping=path.mapping
    query=exact_query(tmp_path)
    status,witness=query(path,problem,mapping,source_atoms=list(mapping),complete_reference=True)
    assert status=='recovered'
    wrong=dict(mapping);atoms=list(mapping);wrong[atoms[0]]=wrong[atoms[1]]
    status,_=query(path,problem,wrong,source_atoms=list(mapping),complete_reference=True)
    assert status=='not_recovered'
    statuses=[]
    for p in tmp_path.glob('*.smt2'):
        solver=z3.Solver();solver.from_string(p.read_text())
        result=solver.check();statuses.append(str(result))
        record=read(p.with_suffix('.json'))
        assert record['outcome']==('recovered' if result==z3.sat else 'not_recovered')
    assert sorted(statuses)==['sat','unsat']


def test_report_missing_seed_remains_unknown(tmp_path):
    from rebenchmark.report import summarize
    run=tmp_path/'run'
    _,features,reference=prepare('[CH3:1][OH:2]>>[CH3:1][OH:2]')
    save(run/'inputs/0/reference.json',dict(features=features,mapping=sorted(reference.items())))
    save(run/'manifest.json',dict(total=1,inputs=[dict(index=0,status='prepared',reactant_atoms=6,product_atoms=6)],
        tasks=[dict(index=0,direction=d) for d in ('R_to_P','P_to_R')],timing='test'))
    for task in range(2):
        folder=tmp_path/f'worker{task}'
        save(folder/'search.json',dict(metrics=dict(cut_count=1)))
        save(folder/'analysis/00000/00/evaluation.json',dict(cut=0,seed=1,reference_recovery='not_recovered'))
        save(folder/'analysis/00000/00/classes.json',[])
        publish(folder,run/'archives'/f'{task:05d}.tar')
    result=summarize(run,tmp_path/'report')
    assert result['modes']['bidirectional/seeds1']['not_recovered']==1
    assert result['modes']['bidirectional/seeds2']['unknown']==1
    assert result['modes']['bidirectional/seeds2']['ranking_complete']==0


def test_case_deadline_spans_stages_and_preserves_evidence(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from rebenchmark import campaign
    run=tmp_path/'run'
    for name in ('input.json','reference.json','source.json'):
        save(run/'inputs/0'/name,{})
    clock=[100.0];budgets=[]
    monkeypatch.setattr(campaign,'CASE_SECONDS',10)
    monkeypatch.setattr(campaign.time,'perf_counter',lambda:clock[0])
    def stage(command,folder,name,seconds,memory):
        budgets.append(seconds)
        save(folder/(name+'-saved.json'),dict(evidence='retained'))
        clock[0]+=8 if name=='search' else seconds
        return dict(status='complete' if name=='search' else 'timeout')
    monkeypatch.setattr(campaign,'guarded',stage)
    manifest=dict(tasks=[dict(index=0,direction='R_to_P',config={})],
        workers=1,stage_seconds=dict(search=100,analysis=100),memory_mib=100)
    campaign.worker(SimpleNamespace(run=run,task=0),manifest)
    assert budgets==[10,2]
    state=read(run/'progress/00000.json')
    assert not state['complete'] and state['active_stage'] is None
    assert state['stages']['analysis']['limit_scope']=='case'
    restore(run/'archives/00000.tar',tmp_path/'restored')
    assert read(tmp_path/'restored/search-saved.json')['evidence']=='retained'
    assert read(tmp_path/'restored/analysis-saved.json')['evidence']=='retained'


def test_report_uncut_complete_does_not_make_incomplete_sweep_negative(tmp_path):
    from rebenchmark.report import summarize
    run=tmp_path/'run'
    _,features,reference=prepare('[CH3:1][OH:2]>>[CH3:1][OH:2]')
    save(run/'inputs/0/reference.json',dict(features=features,mapping=sorted(reference.items())))
    save(run/'manifest.json',dict(total=1,inputs=[dict(index=0,status='prepared',reactant_atoms=6,product_atoms=6)],
        tasks=[dict(index=0,direction=d) for d in ('R_to_P','P_to_R')],timing='test'))
    for task in range(2):
        folder=tmp_path/f'worker{task}'
        save(folder/'search.json',dict(metrics=dict(cut_count=2)))
        save(folder/'analysis/00000/00/evaluation.json',dict(cut=0,seed=1,reference_recovery='not_recovered'))
        save(folder/'analysis/00000/00/classes.json',[])
        publish(folder,run/'archives'/f'{task:05d}.tar')
    serial=summarize(run,tmp_path/'serial')
    parallel=summarize(run,tmp_path/'parallel',workers=2)
    assert serial==parallel
    assert serial['modes']['uncut/bidirectional/seeds1']['not_recovered']==1
    assert serial['modes']['bidirectional/seeds1']['unknown']==1
    assert serial['modes']['bidirectional/seeds1']['ranking_complete']==0


def test_report_continuation_requires_same_task_and_preserved_cuts(tmp_path):
    import pytest
    from rebenchmark.report import resolve_archives
    from rebenchmark.evidence import sha
    run=tmp_path/'run'
    task=dict(index=17,direction='P_to_R',config={'seed_count':10})
    manifest={'tasks':[task]}
    save(run/'manifest.json',manifest)
    old=tmp_path/'old';new=tmp_path/'new'
    for folder in (old,new):
        save(folder/'task.json',task)
        save(folder/'cuts/cut_00000.raw.pkl.gz',{'saved':'unchanged'})
    publish(old,run/'archives/00000.tar')
    target=tmp_path/'continuation.tar';publish(new,target)
    spec={'source_manifest_sha256':sha(run/'manifest.json'),
          'replacements':[{'task':0,'archive':str(target),'sha256':sha(target)}]}
    assert resolve_archives(run,manifest,spec)[0]==target
    save(new/'cuts/cut_00000.raw.pkl.gz',{'saved':'changed'})
    publish(new,target);spec['replacements'][0]['sha256']=sha(target)
    with pytest.raises(AssertionError):resolve_archives(run,manifest,spec)
    save(new/'cuts/cut_00000.raw.pkl.gz',{'saved':'unchanged'})
    save(new/'task.json',{**task,'direction':'R_to_P'})
    publish(new,target);spec['replacements'][0]['sha256']=sha(target)
    with pytest.raises(AssertionError):resolve_archives(run,manifest,spec)


def test_report_overlay_rejects_replacement_of_completed_rows(tmp_path):
    from rebenchmark.report import archive_results, load_overlays
    folder = tmp_path/'case'
    save(folder/'search.json', dict(metrics=dict(cut_count=1)))
    save(folder/'input.json', {})
    save(folder/'reference.json', {})
    save(folder/'task.json', {})
    save(folder/'cuts/cut_00000.raw.pkl.gz', {})
    save(folder/'finalized/00000.pkl.gz', {})
    digest = sha(folder/'cuts/cut_00000.raw.pkl.gz')
    row = dict(cut=0,seed=1,reference_recovery='recovered',raw_sha256=digest,finalized_sha256=digest)
    overlay = dict(state='complete',evaluation=row,classes=[],raw=dict(sha256=digest),
        source={n:dict(sha256=sha(folder/n)) for n in ('input.json','reference.json','task.json')})
    publish(folder,tmp_path/'before.tar')
    assert archive_results(tmp_path/'before.tar',[overlay])['seeds'][0,1]['evaluation']==row
    save(folder/'analysis/00000/00/evaluation.json',row)
    save(folder/'analysis/00000/00/classes.json',[])
    publish(folder,tmp_path/'after.tar')
    with pytest.raises(AssertionError):
        archive_results(tmp_path/'after.tar',[overlay])
    save(tmp_path/'delta.json',{'00000/cut/00000/seed/00':overlay})
    specification=[dict(path=str(tmp_path/'delta.json'),sha256=sha(tmp_path/'delta.json'))]
    assert len(load_overlays(specification)[0])==1
    with pytest.raises(AssertionError):load_overlays(specification*2)


def test_incremental_report_matches_complete_rebuild(tmp_path):
    from collections import Counter
    from rebenchmark.report import summarize, refresh
    run = tmp_path/'run'
    _, features, reference = prepare('[OH2:1]>>[OH2:1]')
    save(run/'inputs/0/reference.json',dict(features=features,mapping=sorted(reference.items())))
    save(run/'manifest.json',dict(total=1,inputs=[dict(index=0,status='prepared',reactant_atoms=3,product_atoms=3)],
        tasks=[dict(index=0,direction=d) for d in ('R_to_P','P_to_R')],timing='test'))
    for task in range(2):
        folder = tmp_path/f'worker{task}'
        save(folder/'search.json', dict(metrics=dict(cut_count=1)))
        for name in ('input.json','reference.json','task.json','cuts/cut_00000.raw.pkl.gz','finalized/00000.pkl.gz'):
            save(folder/name,{})
        for seed in range(10):
            if task == 0 and seed == 0:continue
            save(folder/f'analysis/00000/{seed:02d}/evaluation.json',dict(cut=0,seed=seed+1,reference_recovery='not_recovered'))
            save(folder/f'analysis/00000/{seed:02d}/classes.json',[])
        publish(folder,run/'archives'/f'{task:05d}.tar')
    summarize(run,tmp_path/'old')
    publish(tmp_path/'old',tmp_path/'base-report.tar')
    folder = tmp_path/'worker0';digest=sha(folder/'input.json')
    row = dict(cut=0,seed=1,reference_recovery='recovered',raw_sha256=digest,finalized_sha256=digest)
    delta = dict(state='complete',evaluation=row,classes=[],raw=dict(sha256=digest),
        source={n:dict(sha256=digest) for n in ('input.json','reference.json','task.json')})
    save(tmp_path/'delta.json',{'00000/cut/00000/seed/00':delta})
    actual = refresh(tmp_path/'base-report.tar',run,tmp_path/'updated',
        [dict(path=str(tmp_path/'delta.json'),sha256=sha(tmp_path/'delta.json'))])
    save(folder/'analysis/00000/00/evaluation.json',row)
    save(folder/'analysis/00000/00/classes.json',[])
    publish(folder,run/'archives/00000.tar')
    expected = summarize(run,tmp_path/'complete')
    assert {k:Counter(v) for k,v in actual['modes'].items()} == {k:Counter(v) for k,v in expected['modes'].items()}
