import gc
from dataclasses import replace
import json
from pathlib import Path
import sys
import weakref

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bench'))
import golden_checkpoint_evaluation as streaming
from golden_evaluation import prepare, evaluate_planned
from mappa import search_aam, search_aam_checkpoints, AAMSearchConfig, AAMSearchPlan
from mappa.aam import checkpoint_manifest
from mappa.artifacts import read_raw_cut
from mappa.search_graph import AAMSearchGraph
from mappa.search_symmetry import finalize_graph_symmetry
from mappa.conditioned_symmetry import ConditionedSymmetryWorkspace
from mappa.frag import build_graph


@pytest.mark.parametrize('rows,expected', [([], 'unknown'),
    ([(0, 'not_recovered')], 'unknown'),
    ([(0, 'not_recovered'), (1, 'not_recovered')], 'not_recovered'),
    ([(1, 'not_recovered'), (0, 'not_recovered')], 'not_recovered'),
    ([(0, 'unknown'), (1, 'not_recovered')], 'unknown'),
    ([(0, 'unknown'), (1, 'recovered')], 'recovered')])
def test_complete_union_required_for_absence(rows, expected):
    assert streaming.aggregate_cut_verdicts(2, [dict(cut_index=i, reference_recovery=s) for i,s in rows]) == expected


@pytest.mark.parametrize('indices', [[0,0], [2], [-1]])
def test_duplicate_or_unexpected_verdict_rejected(indices):
    with pytest.raises(ValueError):
        streaming.aggregate_cut_verdicts(2, [dict(cut_index=i, reference_recovery='not_recovered') for i in indices])


@pytest.fixture
def reaction():
    p, f, r = prepare('[CH3:1][OH:2]>>[CH3:1][OH:2]')
    cfg = AAMSearchConfig(seed_count=2, branch_limit=10)
    return AAMSearchPlan(p, p, cfg, False), f, r


def test_disk_search_matches_full_graph_and_resume_without_matching(tmp_path, reaction, monkeypatch):
    plan, features, reference = reaction
    full = search_aam(plan.problem, plan.config, execution='reused_native')
    disk = search_aam_checkpoints(plan.problem, plan.config, intermediate_dir=tmp_path)
    target = build_graph(plan.problem.product.elements, plan.problem.product.wbo, bond_cut=plan.config.graph_floor)
    raw = AAMSearchGraph.combine([read_raw_cut(p) for p in disk.cut_paths])
    combined, _ = finalize_graph_symmetry(raw, target, iso_tolerance=plan.config.iso_tolerance,
        workspace=ConditionedSymmetryWorkspace(target, plan.config.iso_tolerance))
    assert combined == full.graph
    assert disk.metrics.retained_branch_count == full.metrics.retained_branch_count
    assert not (tmp_path/'aam.pkl.gz').exists()
    actual = streaming.evaluate_checkpoints(tmp_path, plan, features, reference)
    expected = evaluate_planned(full, plan, features, reference)
    assert actual['reference_recovery'] == expected['reference_recovery'] == 'recovered'
    import mappa.aam as aam
    monkeypatch.setattr(aam, '_search_cut', lambda *a: pytest.fail('Resume searched an existing cut'))
    monkeypatch.setattr(AAMSearchGraph, 'combine', lambda *a: pytest.fail('Resume combined histories'))
    resumed = search_aam_checkpoints(plan.problem, plan.config, intermediate_dir=tmp_path, resume=True)
    assert resumed.metrics.retained_branch_count == disk.metrics.retained_branch_count


def test_serial_producer_does_not_retain_previous_cut(tmp_path, reaction, monkeypatch):
    import mappa.aam as aam
    original = aam._search_cut
    previous = []
    def observed(cut):
        gc.collect()
        assert all(ref() is None for ref in previous)
        graph, counts = original(cut)
        previous.append(weakref.ref(graph))
        return graph, counts
    monkeypatch.setattr(aam, '_search_cut', observed)
    plan, _, _ = reaction
    search_aam_checkpoints(plan.problem, plan.config, intermediate_dir=tmp_path)
    assert len(previous) > 1
    assert all(ref() is None for ref in previous)


@pytest.mark.parametrize('sweep', [False, True])
def test_verifier_releases_cut_and_rejects_missing_or_wrong_manifest(tmp_path, reaction, monkeypatch, sweep):
    plan, features, reference = reaction
    plan = replace(plan, config=replace(plan.config, sweep_cuts=sweep))
    result = search_aam_checkpoints(plan.problem, plan.config, intermediate_dir=tmp_path)
    live = []
    reader = streaming.read_raw_cut
    def observed(path):
        gc.collect()
        assert all(ref() is None for ref in live)
        graph = reader(path)
        live.append(weakref.ref(graph))
        return graph
    monkeypatch.setattr(streaming, 'read_raw_cut', observed)
    def negative(aam, *args, **kwargs):
        live.append(weakref.ref(aam.graph))
        return {'reference_recovery':'not_recovered'}
    monkeypatch.setattr(streaming, 'evaluate_planned', negative)
    verdict = streaming.evaluate_checkpoints(tmp_path, plan, features, reference)
    assert verdict['reference_recovery'] == 'not_recovered'
    assert verdict['checkpoint_verification']['all_expected_cuts_checked']
    result.cut_paths[-1].unlink()
    assert streaming.evaluate_checkpoints(tmp_path, plan, features, reference)['reference_recovery'] == 'unknown'
    assert streaming.evaluate_checkpoints(tmp_path, plan, features, reference, seconds=0)['reference_recovery'] == 'unknown'
    manifest = checkpoint_manifest(plan.problem, plan.config)
    manifest['config']['branch_limit'] += 1
    (tmp_path/'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='configuration'):
        streaming.evaluate_checkpoints(tmp_path, plan, features, reference)


def test_recorded_failed_cases_now_have_definitive_negative_union():
    # Original cut verdicts, retained with evidence hashes by the audit.
    root = Path(__file__).resolve().parents[1]/'reports/checkpoint_memory_fix_20260913'
    audit = json.loads((root/'audit.json').read_text())
    for direction in audit['directions']:
        path = root/f"case{direction['case']}-{direction['direction']}.json"
        proof = json.loads(path.read_text())['checkpoint_verification']
        rows = [dict(cut_index=int(Path(r['cut']).name.split('_')[1].split('.')[0]),
                     reference_recovery=r['reference_recovery']) for r in proof['checked']]
        assert streaming.aggregate_cut_verdicts(direction['expected_cuts'], rows) == 'not_recovered'


def test_parallel_checkpoint_producers_preserve_cut_graphs(tmp_path, reaction):
    plan, _, _ = reaction
    serial = search_aam_checkpoints(plan.problem, plan.config, intermediate_dir=tmp_path/'serial')
    parallel = search_aam_checkpoints(plan.problem, plan.config, intermediate_dir=tmp_path/'parallel', workers=2)
    assert len(serial.cut_paths) == len(parallel.cut_paths)
    assert all(read_raw_cut(a) == read_raw_cut(b) for a,b in zip(serial.cut_paths,parallel.cut_paths))


def test_wrong_cut_context_is_rejected(tmp_path, reaction):
    from dataclasses import replace
    from mappa.artifacts import write_raw_cut
    plan, features, reference = reaction
    result = search_aam_checkpoints(plan.problem, plan.config, intermediate_dir=tmp_path)
    graph = read_raw_cut(result.cut_paths[0])
    graph = replace(graph, contexts=tuple(replace(c,cuts=((999,1000),)) for c in graph.contexts))
    write_raw_cut(graph, result.cut_paths[0])
    with pytest.raises(ValueError, match='cut context'):
        streaming.evaluate_checkpoints(tmp_path, plan, features, reference)
