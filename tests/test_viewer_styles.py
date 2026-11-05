"""Keep presentation shared, and prevent display fitting from changing geometry."""
from pathlib import Path
import json
import re
import shutil
import subprocess

import numpy as np
import pytest

from mappa.presentation.viewers import ASSETS, align_product, comparison_document, stylesheet


def test_only_two_shared_viewer_styles():
    assert sorted(p.name for p in ASSETS.glob('*viewer.css')) == [
        'batch_viewer.css', 'reaction_viewer.css']
    assert '#ffd700' in stylesheet('reaction')
    assert 'background:#fafafa' in stylesheet('reaction')


def test_report_generators_do_not_define_independent_skins():
    root = Path(__file__).resolve().parents[1]
    for folder in ('bench', 'tools'):
        for path in (root / folder).glob('*.py'):
            assert not re.search(r'<style(?:>|\s)', path.read_text()), path
    for name in ('elementary_comparison_viewer.html', 'missing_patterns_viewer.html'):
        assert not (root / 'bench' / name).exists()


def test_display_fit_is_a_proper_rigid_transform():
    r = np.array([[0., 0., 0.], [2., 0., 0.], [.2, 1.3, 0.], [.1, .3, 1.7]])
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    p = r @ rotation + [9., -4., 2.]
    aligned = np.asarray(align_product(r, p))
    np.testing.assert_allclose(aligned, r, atol=1e-12)
    mirrored = p * [-1, 1, 1]
    fitted = np.asarray(align_product(r, mirrored))
    np.testing.assert_allclose(
        np.linalg.norm(fitted[:, None] - fitted[None, :], axis=-1),
        np.linalg.norm(mirrored[:, None] - mirrored[None, :], axis=-1), atol=1e-12)
    assert np.sign(np.linalg.det(fitted[1:] - fitted[0])) == np.sign(np.linalg.det(mirrored[1:] - mirrored[0]))


def test_batch_style_has_one_navigation_layout():
    css = stylesheet('batch')
    assert 'data-viewer-layout="batch"' in css
    assert 'data-viewer-layout="typed_batch"' not in css
    assert 'data-viewer-layout="catalog"' not in css
    assert 'data-viewer-layout="precursor"' not in css
    with pytest.raises(ValueError):
        stylesheet('catalog')


def test_embedded_navigation_restores_documents_and_switches_cases():
    if not shutil.which('node'):
        pytest.skip('Node is needed to execute the viewer interactions')
    from mappa.presentation.viewers import embedded_collection_html
    # Dollar replacement tokens must survive library insertion literally.
    library = '<script>const labels = "$& $` $\'";</script>'
    documents = [f'<html><head>{library}</head><body>Case {i} — C–O</body></html>'
                 for i in range(2)]
    cases = [dict(id=str(i), atoms=i+2, badge='small', meta='saved', detail=f'Case {i}',
                  document=value.encode()) for i, value in enumerate(documents)]
    page = embedded_collection_html(cases, title='Saved <cases>', summary='Two cases')
    assert '<title>Saved &lt;cases&gt;</title>' in page
    script = re.search(r'<script>\s*const CASES=.*?</script>', page, re.S).group()[8:-9]
    harness = '''
const nodes=new Map();
const document={
 getElementById(id){if(!nodes.has(id))nodes.set(id,{children:[],classList:{add(){},remove(){}},
  replaceChildren(){this.children=[]},appendChild(value){this.children.push(value)}});return nodes.get(id)},
 createElement(){return {children:[],appendChild(value){this.children.push(value)}}}
};
const location={hash:''},window={open(){}};
'''
    checks = '''
await openCase(CASES[1]);
const second=frame.srcdoc;
document.getElementById('search').oninput({target:{value:'0'}});
await list.children[0].onclick();
console.log(JSON.stringify({second,first:frame.srcdoc,count:list.children.length,
 enabled:!standalone.disabled,detail:current.textContent,hash:location.hash}));
'''
    result = subprocess.run(['node', '-'], input='(async()=>{'+harness+script+checks+'})()',
                            text=True, capture_output=True, check=True)
    value = json.loads(result.stdout)
    assert value == dict(second=documents[1], first=documents[0], count=1,
                         enabled=True, detail='Case 0', hash='0')


def test_embedded_collection_rejects_different_libraries():
    from mappa.presentation.viewers import embedded_collection_html
    cases = [dict(id=str(i), document=f'<script>library{i}</script>'.encode()) for i in range(2)]
    with pytest.raises(ValueError, match='Shared library differs for 1'):
        embedded_collection_html(cases, title='Cases', summary='')


@pytest.mark.parametrize('typed', [False, True])
def test_batch_adapters_preserve_saved_case_payloads(tmp_path, typed):
    import base64
    import gzip
    from importlib import import_module
    source = b'<html><head><script>shared library</script></head><body>saved case</body></html>'
    case = tmp_path / 'cases' / 'case1'
    case.mkdir(parents=True)
    if typed:
        (case / 'ts_scores.json').write_text(json.dumps(dict(case='case1', atom_count=2, rankings=[])))
        (case / 'view.html').write_bytes(source)
        module = 'build_typed_ts_integrated_viewer'
    else:
        (tmp_path / 'inventory.json').write_text('{"case_count":1}')
        (tmp_path / 'batch_summary.json').write_text('{"error_count":0,"total_elapsed_case_seconds":1}')
        (tmp_path / 'manifests').mkdir()
        for tier in ('small', 'medium', 'large'):
            records = [dict(source_index=0, step_id='case1', atom_count=2)] if tier == 'small' else []
            (tmp_path / 'manifests' / f'{tier}.json').write_text(json.dumps(dict(cases=records)))
        (case / 'summary.json').write_text('{"elapsed_seconds":1}')
        (case / 'rp_stage.json').write_text('{"mechanisms":[]}')
        view = tmp_path / 'views' / 'case1'
        view.mkdir(parents=True)
        (view / 'view.html').write_bytes(source)
        module = 'build_integrated_batch_viewer'
    output = tmp_path / 'view.html'
    result = import_module(f'tools.{module}').build(tmp_path, output)
    page = output.read_text()
    records = json.loads(re.search(r'const CASES=(.*);', page).group(1))
    library = re.search(r'COMMON_B64="([^"]+)"', page).group(1)
    body = gzip.decompress(base64.b64decode(records[0]['payload']))
    restored = body.replace(b'<!-- MAPPA_SHARED_LIBRARY -->', gzip.decompress(base64.b64decode(library)), 1)
    assert restored == source
    assert records[0]['id'] == 'case1' and records[0]['atoms'] == 2
    assert result['cases'] == 1


def _partial_comparison():
    return dict(name='Partial reference', index=7, endpoints=[
        dict(elements=['C', 'O', 'N'], coordinates=[[0, 0, 0], [1, 0, 0], [2, 0, 0]],
             wbo=[[0, 1, 0], [1, 0, 1], [0, 1, 0]], atom_indices=[2, 5, 8]),
        dict(elements=['O', 'C', 'F', 'N'], coordinates=[[0, 1, 0], [1, 1, 0], [2, 1, 0], [3, 1, 0]],
             wbo=[[0, 2, 0, 0], [2, 0, 1.5, 0], [0, 1.5, 0, 0], [0, 0, 0, 0]],
             atom_indices=[4, 6, 9, 12])], records=[dict(name='Reference', mapping=[[0, 1], [1, 0]],
        events=[dict(kind='broken', r=[1, 2], p=[0], wbo=[1, None]),
                dict(kind='formed', r=[0], p=[1, 2], wbo=[None, 1])])])


def test_partial_comparison_preserves_unmatched_atoms_and_boundary_events():
    case = _partial_comparison()
    document = comparison_document(case)
    mechanism = document['mechanisms'][0]
    assert mechanism['mapping_RP'] == {0: 1, 1: 0}
    assert mechanism['product_xyz_in_R'] is None
    assert mechanism['product_xyz_in_R_aligned'] is None
    assert document['product']['coords'] == case['endpoints'][1]['coordinates']
    assert document['product']['atom_indices'] == [4, 6, 9, 12]
    assert mechanism['broken_bonds_R'] == [[1, 2]]
    assert mechanism['formed_bonds_R'] == []
    assert mechanism['formed_bonds_P'] == [[1, 2]]
    assert mechanism['event_counts']['formed'] == mechanism['event_counts']['broken'] == 1


@pytest.mark.parametrize('mapping', [
    [[0, 1], [0, 1]], [[0, 1], [2, 1]], [[0, 4]], [[0, 0]],
])
def test_comparison_rejects_invalid_correspondences(mapping):
    case = _partial_comparison()
    case['records'][0]['mapping'] = mapping
    with pytest.raises(ValueError):
        comparison_document(case)


def test_complete_comparison_keeps_r_order_and_rigid_alignment():
    r = dict(elements=['C', 'O', 'N'], coordinates=[[0, 0, 0], [1, 0, 0], [0, 2, 0]],
             wbo=[[0, 1, 1], [1, 0, 0], [1, 0, 0]])
    p = dict(elements=['O', 'N', 'C'], coordinates=[[6, 3, 0], [5, 5, 0], [5, 3, 0]],
             wbo=[[0, 0, 1], [0, 0, 1], [1, 1, 0]])
    document = comparison_document(dict(name='Complete', index=1, endpoints=[r, p],
        records=[dict(name='Saved', mapping=[2, 0, 1], events=[])]))
    mechanism = document['mechanisms'][0]
    assert mechanism['product_xyz_in_R'] == [p['coordinates'][i] for i in (2, 0, 1)]
    np.testing.assert_allclose(mechanism['product_xyz_in_R_aligned'], r['coordinates'], atol=1e-12)


def test_native_partial_viewer_uses_endpoint_bonds_and_original_atom_labels():
    if not shutil.which('node'):
        pytest.skip('Node is needed to execute the viewer interactions')
    document = comparison_document(_partial_comparison())
    template = (ASSETS / 'reaction_viewer.html').read_text()
    script = re.search(r'<script>\s*const DATA =.*?</script>', template, re.S).group()[8:-9]
    script = script.replace('__DATA__', json.dumps(document))
    harness = '''
const nodes = new Map();
const document = {
 getElementById(id) { if (!nodes.has(id)) nodes.set(id, {style:{},children:[],
  appendChild(child){this.children.push(child)}, replaceChildren(){this.children=[]}}); return nodes.get(id); },
 createElement(){return {dataset:{}}}
};
const window = {addEventListener(){}};
'''
    checks = '''
const mech = findMech(currentMechId), frame = productFrame(mech);
let atoms;
addEndpointModel({addModel(){return {addAtoms(value){atoms=value}}},addCylinder(){}},
                 'vw_P', frame.elements, frame.xyz);
describeWitness(mech);
const events = document.getElementById('event-table').children.map(row=>row.textContent);
mech.event_records=[]; describeWitness(mech);
console.log(JSON.stringify({aligned:frame.aligned, bonds:atoms.map(a=>a.bonds), orders:atoms.map(a=>a.bondOrder),
 labels:productAtomLabels(mech,frame.aligned), reactantLabels:reactantAtomLabels(),
 events, mappingVisible:!document.getElementById('event-details').hidden}));
'''
    result = subprocess.run(['node', '-'], input=harness + script + checks, text=True,
                            capture_output=True, check=True)
    value = json.loads(result.stdout)
    assert value['aligned'] is False  # The checkbox starts enabled, but no full alignment exists.
    assert value['bonds'] == [[1], [0, 2], [1], []]
    assert value['orders'] == [[2], [2, 1], [1], []]
    assert value['labels'] == ['P4 ← R5', 'P6 ← R2', 'P9 (unmatched)', 'P12 (unmatched)']
    assert value['reactantLabels'] == ['R2', 'R5', 'R8']
    assert 'unmatched' in value['events'][0] and 'unmatched' in value['events'][1]
    assert value['mappingVisible'] is True
