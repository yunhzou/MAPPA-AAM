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
