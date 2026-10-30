"""Saved continuation never reruns search or overwrites completed checks."""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import tarfile

import pytest

from mappa import AAMSearchConfig, search_aam
from mappa.mappa_core.utils.io.checkpoints import write_graph_checkpoint, write_raw_cut
from mappa.mappa_core.core.scheduler.runtime import Workspace
from golden_evaluation import prepare
from rebenchmark.evidence import save, publish, read
from rebenchmark.saved_analysis import plan_missing, read_member, evaluate
from rebenchmark.analysis import analyze_seed
from rebenchmark.inputs import load_plan


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('finalized', [False, True])
def test_saved_continuation_matches_direct_evaluation(tmp_path, monkeypatch, reverse, finalized):
    def endpoint(e):
        return dict(elements=e.elements, coordinates=e.coordinates.tolist(), wbo=e.wbo.tolist(),
                    label=e.label, metadata=dict(e.metadata))
    folder = tmp_path/'source'
    problem, features, reference = prepare('[CH3:1][OH:2]>>[CH3:1][OH:2].[Cl-]')
    config = AAMSearchConfig(seed_count=2, branch_limit=100, iso_tolerance=1)
    save(folder/'input.json', dict(name='test', reactant=endpoint(problem.reactant), product=endpoint(problem.product)))
    save(folder/'reference.json', dict(features=features, mapping=sorted(reference.items())))
    save(folder/'task.json', dict(index=0, direction='P_to_R' if reverse else 'R_to_P', config=asdict(config)))
    save(folder/'search.json', dict(metrics=dict(cut_count=1)))
    plan = load_plan(folder)
    graph = search_aam(plan.problem, config).graph
    (folder/'cuts').mkdir()
    write_raw_cut(graph, folder/'cuts/cut_00000.raw.pkl.gz')
    if finalized:
        (folder/'finalized').mkdir()
        write_graph_checkpoint(graph, folder/'finalized/00000.pkl.gz')
    # One evaluation exists and must never be dispatched again.
    save(folder/'analysis/00000/00/evaluation.json', dict(immutable=True))
    receipt = publish(folder, tmp_path/'case.tar')
    tasks = plan_missing(tmp_path/'case.tar', 42)
    assert [task.key for task in tasks] == ['00042/cut/00000/seed/01']
    direct = tmp_path/'direct'
    analyze_seed(plan, read(folder/'reference.json'), graph, 1, direct, 0)
    import mappa.mappa_core.core.scheduler.search
    monkeypatch.setattr(mappa.mappa_core.core.scheduler.search, '_search_cut', lambda *args: pytest.fail('Search must not run'))
    workspace = Workspace(tmp_path, tmp_path/'bundles')
    actual = evaluate(tasks[0].payload, workspace)
    expected = read(direct/'evaluation.json')
    assert actual['evaluation']['reference_recovery'] == expected['reference_recovery']
    assert actual['classes'] == read(direct/'classes.json')
    assert actual['search_rerun'] is False
    assert len(list((tmp_path/'bundles').iterdir())) <= 1
    assert hashlib.sha256((tmp_path/'case.tar').read_bytes()).hexdigest() == receipt['sha256']


def test_archive_byte_ranges_detect_corruption(tmp_path):
    folder = tmp_path/'source'
    save(folder/'item.json', dict(value=1))
    publish(folder, tmp_path/'case.tar')
    with tarfile.open(tmp_path/'case.tar') as tar:
        entry = tar.getmember('item.json')
        data = tar.extractfile(entry).read()
        member = dict(archive=str(tmp_path/'case.tar'), name=entry.name,
                      offset=entry.offset_data, bytes=entry.size, sha256=hashlib.sha256(data).hexdigest())
    assert json.loads(read_member(member)) == {'value': 1}
    with (tmp_path/'case.tar').open('r+b') as stream:
        stream.seek(member['offset']); stream.write(b'x')
    with pytest.raises(ValueError, match='checksum'):
        read_member(member)
