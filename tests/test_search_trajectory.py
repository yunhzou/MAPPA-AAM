"""Scientific replay integrity and a manifest pipeline independent of PR7."""
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from mappa import AAMProblem,AAMSearchConfig,MolecularEndpoint,search_aam
from mappa.artifacts import write_aam_checkpoint, read_aam_checkpoint, aam_json
from mappa.search_trajectory import build_trajectory,growth,extension
from mappa.viewers import growth_trace_html


def toy_archive(folder, suffix='', translation=0):
    elements=['C','O','O','H'];coords=np.array([[0.,0.,0.],[1.2,0,0],[-1,0,0],[-1,1,0]])
    w=np.array([[0,1.8,.8,0],[1.8,0,0,0],[.8,0,0,.9],[0,0,.9,0.]])
    order=[2,0,3,1]
    r=MolecularEndpoint(elements=elements,coordinates=coords,wbo=w)
    p=MolecularEndpoint(elements=[elements[i] for i in order],coordinates=coords[order]+translation,wbo=w[np.ix_(order,order)])
    result=search_aam(AAMProblem(r,p,name='four_atom_example'),
        AAMSearchConfig(iso_tolerance=.1,graph_floor=.4,branch_limit=32,seed_count=1),workers=1)
    archive=folder/f'toy{suffix}.pkl.gz';write_aam_checkpoint(result,archive)
    return archive


def test_replay_uses_archive_config_and_preserves_inputs(tmp_path):
    archive=toy_archive(tmp_path);before=archive.read_bytes()
    original=growth._extend_sym_cands,extension._dedupe_children
    doc=build_trajectory([dict(archive=str(archive),context=0)])
    assert archive.read_bytes()==before
    assert (growth._extend_sym_cands,extension._dedupe_children)==original
    assert doc['runs'][0]['config']['iso_tolerance']==.1
    assert doc['runs'][0]['config']['graph_floor']==.4
    for path in doc['runs'][0]['paths']:
        assert path['complete'] and len(path['mapping'])==4
        assert path['events']==[]
        assert any(f['kind']=='commit' for f in path['frames'])
        assert all(check['archived_fragment_matches']=='exact' for check in doc['checks'])
    html=growth_trace_html(doc)
    assert 'PR7' not in html and 'data-viewer-layout="growth_trace"' in html


def test_manifest_cli_and_render_only(tmp_path,monkeypatch):
    archive=toy_archive(tmp_path)
    interchange=tmp_path/'toy.json';interchange.write_text(aam_json(read_aam_checkpoint(archive)))
    manifest=tmp_path/'manifest.json'
    manifest.write_text(json.dumps(dict(selections=[dict(archive=interchange.name,context=0)])))
    root=Path(__file__).resolve().parents[1];output=tmp_path/'output'
    subprocess.run([sys.executable,str(root/'tools/build_search_trajectory.py'),str(manifest),str(output)],
                   check=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},capture_output=True)
    doc=json.loads((output/'trace.json').read_text())
    def fail(*a,**k):raise AssertionError('Rendering must not replay or search')
    monkeypatch.setattr(growth,'grow_island',fail)
    assert growth_trace_html(doc)==(output/'algorithm_trajectory.html').read_text()
    assert json.loads((output/'validation.json').read_text())['status']=='passed'


def test_mismatched_endpoints_are_rejected(tmp_path):
    a,b=toy_archive(tmp_path),toy_archive(tmp_path,'_other',translation=2)
    with pytest.raises(ValueError,match='identical endpoint'):
        build_trajectory([dict(archive=str(a)),dict(archive=str(b))])


def test_in_memory_replay_matches_archive_and_keeps_result(tmp_path, monkeypatch):
    from mappa import search_trajectory
    from mappa.viewers import aam_growth_html
    archive = toy_archive(tmp_path)
    aam = read_aam_checkpoint(archive)
    before = aam_json(aam)
    archived = build_trajectory([dict(archive=str(archive), context=0)])
    def fail(*a, **k):
        raise AssertionError('In-memory viewing must not read an archive')
    monkeypatch.setattr(search_trajectory, 'read_archive', fail)
    direct = build_trajectory([dict(aam=aam, context=0)])
    assert direct['runs'][0]['paths'] == archived['runs'][0]['paths']
    assert direct['runs'][0]['traces'] == archived['runs'][0]['traces']
    assert direct['runs'][0]['provenance']['source'] == 'in_memory'
    assert aam_json(aam) == before
    page = aam_growth_html(aam, context=0)
    assert 'data-viewer-layout="growth_trace"' in page
    assert 'Sweep 0:' in page
    with pytest.raises(ValueError, match='one context'):
        aam_growth_html(aam, terminals=[aam.graph.terminals[0]])
    with pytest.raises(ValueError, match='context must index'):
        aam_growth_html(aam, context=-1)
    with pytest.raises(ValueError, match='exactly one'):
        build_trajectory([dict(aam=aam, archive=str(archive))])


def test_event_overlay_uses_decoder_thresholds_and_target_indexing():
    from mappa.search_trajectory import bond_events
    # Inclusive 0.5 boundary, metal-only 0.3 event, subthreshold edge loss,
    # and a formed bond under a nonidentity atom mapping.
    elements = ['C', 'O', 'Fe', 'H']
    r = np.zeros((4, 4)); p = np.zeros((4, 4))
    for a, b, x, y in [(0, 1, 1., .5), (0, 2, .3, 0.),
                        (0, 3, .25, 0.), (1, 3, 0., .5)]:
        r[a, b] = r[b, a] = x
        p[a, b] = p[b, a] = y
    order = [3, 2, 1, 0]
    mapping = {i: order.index(i) for i in range(4)}
    raw = dict(reactant=dict(elements=elements, wbo=r.tolist()),
               product=dict(elements=[elements[i] for i in order],
                            wbo=p[np.ix_(order, order)].tolist()))
    events = bond_events(raw, mapping)
    assert {(e['kind'], tuple(e['r'])) for e in events} == {
        ('weakened', (0, 1)), ('broken', (0, 2)), ('formed', (1, 3))}
    assert all(e['p'] == [mapping[a] for a in e['r']] for e in events)
    uniform = bond_events(raw, mapping, metal_tolerance=None)
    assert {tuple(e['r']) for e in uniform} == {(0, 1), (1, 3)}
