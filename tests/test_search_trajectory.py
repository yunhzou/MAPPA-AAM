"""Scientific replay integrity and a manifest pipeline independent of PR7."""
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from mappa import AAMProblem,AAMSearchConfig,MolecularEndpoint,search_aam
from mappa.mappa_core.utils.io.checkpoints import write_aam_checkpoint, read_aam_checkpoint
from mappa.mappa_core.utils.io.records import aam_json
from mappa.presentation.trajectory import build_trajectory, growth
from mappa.presentation.viewers import growth_trace_html


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
    original=growth.extend_witness
    doc=build_trajectory([dict(archive=str(archive),context=0)])
    assert archive.read_bytes()==before
    assert growth.extend_witness==original
    assert doc['runs'][0]['config']['iso_tolerance']==.1
    assert doc['runs'][0]['config']['graph_floor']==.4
    for path in doc['runs'][0]['paths']:
        assert path['complete'] and len(path['mapping'])==4
        assert path['events']==[]
        assert any(f['kind']=='commit' for f in path['frames'])
        growth_frames=[f for f in path['frames'] if f.get('candidate_kind')=='witness']
        assert growth_frames and all(len(f['candidates'])==1 for f in growth_frames)
        assert all(f['candidate_kind']=='compressed_match' for f in path['frames'] if f['kind']=='seed_end')
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
    from mappa.presentation import trajectory as search_trajectory
    from mappa.presentation.viewers import aam_growth_html
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
    from mappa.presentation.trajectory import bond_events
    # Inclusive 0.5/0.3 matching boundaries, subthreshold edge loss,
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
    assert bond_events(raw, mapping) == []
    r[0, 1] = r[1, 0] = 1.1
    r[0, 2] = r[2, 0] = .4
    p[1, 3] = p[3, 1] = .6
    raw['reactant']['wbo'] = r.tolist()
    raw['product']['wbo'] = p[np.ix_(order, order)].tolist()
    events = bond_events(raw, mapping)
    assert {(e['kind'], tuple(e['r'])) for e in events} == {
        ('weakened', (0, 1)), ('broken', (0, 2)), ('formed', (1, 3))}
    assert all(e['p'] == [mapping[a] for a in e['r']] for e in events)
    uniform = bond_events(raw, mapping, metal_tolerance=None)
    assert {tuple(e['r']) for e in uniform} == {(0, 1), (1, 3)}


def test_event_overlay_keeps_mapped_unmapped_boundary_bonds():
    from mappa.presentation.trajectory import bond_events
    raw = dict(reactant=dict(elements=['C', 'O'], wbo=[[0., 0.], [0., 0.]]),
               product=dict(elements=['C', 'O', 'O'], wbo=[[0., 0., .8], [0., 0., 0.], [.8, 0., 0.]]))
    events = bond_events(raw, {0: 0, 1: 1})
    assert events == [dict(kind='formed', r=[0], p=[0, 2], wbo=[None, .8])]


def test_partial_saved_mapping_replays_boundary_loss_and_gain():
    source=MolecularEndpoint(['C','O'],np.zeros((2,3)),_matrix(2,[(0,1)]))
    target=MolecularEndpoint(['C','N'],np.zeros((2,3)),_matrix(2,[(0,1)]))
    result=search_aam(AAMProblem(source,target),
                      AAMSearchConfig(seed_count=1,branch_limit=10),workers=1)
    document=build_trajectory([dict(aam=result)])
    path,=document['runs'][0]['paths']
    assert not path['complete'] and path['mapping']=={0:0}
    assert path['events']==[
        dict(kind='broken',r=[0,1],p=[0],wbo=[1.,None]),
        dict(kind='formed',r=[0],p=[0,1],wbo=[None,1.])]
    assert document['checks'][0]['events']==2


def test_render_geometry_uses_trace_metadata_without_changing_raw_trace(tmp_path):
    weights=_matrix(5,[(0,i) for i in range(1,5)])
    endpoint=dict(elements=['C','H','H','H','H'],coordinates=np.zeros((5,3)).tolist(),
                  wbo=weights.tolist(),metadata=dict(scoring_features=dict(
                      heavy=[0],colors=[['C',0,0,4,'']],bonds=[])))
    source=MolecularEndpoint(**endpoint)
    result=search_aam(AAMProblem(source,source),
                      AAMSearchConfig(seed_count=1,branch_limit=10),workers=1)
    document=build_trajectory([dict(aam=result)])
    trace=tmp_path/'trace.json'
    trace.write_text(json.dumps(document))
    original=trace.read_bytes()
    root=Path(__file__).resolve().parents[1]
    command=[sys.executable,str(root/'tools/render_search_trajectory.py'),str(trace),
             str(tmp_path/'algorithm_trajectory.html'),'--geometry-cache',str(tmp_path/'geometry')]
    subprocess.run(command,check=True,capture_output=True)
    displayed=json.loads((tmp_path/'display-trace.json').read_text())
    assert trace.read_bytes()==original
    assert displayed['runs']==json.loads(original)['runs']
    assert np.linalg.matrix_rank(np.array(displayed['input']['reactant']['coordinates']))==3
    assert displayed['display_geometry']['endpoints']['reactant']['source']=='generated'
    assert len(list((tmp_path/'geometry').iterdir()))==1
    for side in ('reactant','product'):
        assert displayed['input'][side]['elements']==document['input'][side]['elements']
        assert displayed['input'][side]['wbo']==document['input'][side]['wbo']
    # Old traces can supply the same indexed chemistry separately.
    for side in ('reactant','product'):
        del document['input'][side]['metadata']
    trace.write_text(json.dumps(document))
    chemistry=tmp_path/'input.json'
    chemistry.write_text(json.dumps(dict(reactant=endpoint,product=endpoint)))
    subprocess.run(command+['--geometry-input',str(chemistry)],check=True,capture_output=True)
    legacy_display=json.loads((tmp_path/'display-trace.json').read_text())
    for side in ('reactant','product'):
        expected={key:value for key,value in displayed['input'][side].items() if key!='metadata'}
        assert legacy_display['input'][side]==expected


def test_display_geometry_rejects_a_different_indexed_endpoint():
    from mappa.presentation.trajectory import with_display_geometry
    endpoint=dict(elements=['C','O'],coordinates=[[0.,0.,0.],[1.,0.,0.]],
                  wbo=[[0.,1.],[1.,0.]])
    document=dict(input=dict(reactant=endpoint,product=endpoint),runs=[])
    changed=dict(endpoint,elements=['O','C'])
    with pytest.raises(ValueError,match='atom order and WBOs'):
        with_display_geometry(document,cache=Path('unused'),
                              endpoints=dict(reactant=changed,product=endpoint))


def test_selected_terminal_diagram_excludes_other_branches_but_keeps_candidates():
    source=MolecularEndpoint(['C','O'],np.zeros((2,3)),_matrix(2,[(0,1)]))
    target=MolecularEndpoint(['C','O','F','C','O','N'],np.zeros((6,3)),
                              _matrix(6,[(0,1),(1,2),(3,4),(4,5)]))
    result=search_aam(AAMProblem(source,target),
                      AAMSearchConfig(seed_count=1,branch_limit=16,sweep_cuts=False),workers=1)
    assert len(result.graph.terminals)==2
    terminal=result.graph.terminals[0]
    path=next(result.graph.paths(terminal))
    selected={terminal}|{node for edge in path.transitions for node in (
        result.graph.transitions[edge].source,result.graph.transitions[edge].target)}
    document=build_trajectory([dict(aam=result,terminals=[terminal])])
    run=document['runs'][0]
    assert {state['id'] for state in run['graph']['states']}==selected
    assert len(selected)<len(result.graph.states)
    assert {edge['id'] for edge in run['graph']['transitions']}==set(path.transitions)
    saturated=[frame for frame in run['paths'][0]['frames'] if frame['kind']=='seed_end']
    assert len(saturated)==1 and len(saturated[0]['candidates'])==2


def test_display_marginals_use_local_relation_despite_external_asymmetry():
    import z3
    from mappa import match_fragment
    from mappa.mappa_core.core.matching.graphs import build_graph
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints, compile_relaxed_relation
    from mappa.presentation.trajectory import with_represented_targets

    source = build_graph(['O', 'C', 'O'], _matrix(3, [(0, 1), (1, 2)]))
    target = build_graph(['H', 'O', 'C', 'O'], _matrix(4, [(0, 1), (1, 2), (2, 3)]))
    fragment, = match_fragment(source, target, seed=1).matches
    assert fragment.target_atoms == (1, 2, 3)
    candidate = dict(witness=dict(fragment.witness), relation=fragment.to_record())
    document = dict(runs=[dict(paths=[dict(frames=[dict(candidates=[candidate])])])])
    projected = with_represented_targets(document)
    actual = projected['runs'][0]['paths'][0]['frames'][0]['candidates'][0]
    assert actual['target_options'] == {'0': [1, 3], '1': [2], '2': [1, 3]}
    assert 'target_options' not in candidate
    # The direct Boolean query must retain the previous exact joint relation.
    system = FiniteConstraints()
    compiled = compile_relaxed_relation(system, fragment.relation, 'comparison')
    solver = z3.Solver()
    system.emit(solver)
    previous = {str(atom): [target for target in fragment.target_atoms
                if solver.check(system.expression(value) == target) == z3.sat]
                for atom, (value, _) in compiled.images.items()}
    assert actual['target_options'] == previous
    assert solver.check(system.expression(compiled.images[0][0]) == 1,
                        system.expression(compiled.images[2][0]) == 1) == z3.unsat


def test_saturation_shows_other_components_and_rejects_lost_placement_actions(monkeypatch):
    from dataclasses import replace
    from mappa.presentation.trajectory import with_represented_targets

    source = MolecularEndpoint(elements=['C', 'O'], coordinates=np.zeros((2, 3)),
                               wbo=_matrix(2, [(0, 1)]))
    target = MolecularEndpoint(elements=['C', 'O', 'C', 'O'], coordinates=np.zeros((4, 3)),
                               wbo=_matrix(4, [(0, 1), (2, 3)]))
    result = search_aam(AAMProblem(source, target),
                       AAMSearchConfig(seed_count=1, branch_limit=32), workers=1)
    document = build_trajectory([dict(aam=result)])
    projected = with_represented_targets(document)
    saturation = [frame for frame in projected['runs'][0]['paths'][0]['frames']
                  if frame['kind'] == 'seed_end']
    assert len(saturation) == 1 and len(saturation[0]['candidates']) == 1
    assert saturation[0]['placement_generators']
    assert saturation[0]['candidates'][0]['target_options'] == {'0': [0, 2], '1': [1, 3]}

    # The same local fragment has narrower marginals without its shared action.
    canonical_frame = dict(saturation[0], placement_generators=())
    document['runs'][0]['paths'][0]['frames'].append(canonical_frame)
    projected = with_represented_targets(document)
    canonical = projected['runs'][0]['paths'][0]['frames'][-1]['candidates'][0]
    assert canonical['target_options'] == {str(a): [b] for a, b in canonical['witness'].items()}

    original = growth.grow_island_python
    monkeypatch.setattr(growth, 'grow_island_python',
                        lambda *a, **k: replace(original(*a, **k), placement_generators=()))
    with pytest.raises(ValueError, match='Replay differs'):
        build_trajectory([dict(aam=result)])


def test_replay_compares_placement_groups_instead_of_generator_lists():
    from mappa.presentation.trajectory import _same_placement_group

    rotation = (1, 2, 0)
    inverse = (2, 0, 1)
    assert _same_placement_group((rotation,), (inverse,), 3)
    assert _same_placement_group((), ((0, 1, 2),), 3)
    assert not _same_placement_group((rotation,), (), 3)


def test_saved_cut_cli_preserves_history_and_checks_source_ranges(tmp_path):
    from dataclasses import asdict
    from mappa.mappa_core.core.scheduler.bundles import BundleWriter
    from mappa.mappa_core.core.scheduler.results import CutResult, SearchInput

    archive=toy_archive(tmp_path)
    result=read_aam_checkpoint(archive)
    writer=BundleWriter(tmp_path,None)
    source=writer.put(SearchInput(result.problem,result.config,'native'))
    cut=writer.put(CutResult(result.graph,0,1,()))
    terminal=result.graph.terminals[0]
    transitions=list(next(result.graph.paths(terminal)).transitions)
    selection=dict(saved_cut=dict(input=asdict(source),cut=asdict(cut)),context=0,
                   terminals=[terminal],expected_transitions={str(terminal):transitions})
    manifest=tmp_path/'manifest.json'
    manifest.write_text(json.dumps(dict(selections=[selection])))
    root=Path(__file__).resolve().parents[1]
    command=[sys.executable,str(root/'tools/build_search_trajectory.py'),str(manifest),str(tmp_path/'view')]
    subprocess.run(command,check=True,capture_output=True)
    document=json.loads((tmp_path/'view/trace.json').read_text())
    run=document['runs'][0]
    assert run['paths'][0]['mapping']=={str(a):b for a,b in result.graph.states[terminal].mapping}
    assert run['provenance']['saved_cut']['cut']['sha256']==cut.sha256
    assert run['provenance']['metrics']=='Unmeasured archive adapter'
    selection['expected_transitions'][str(terminal)]=transitions+[999]
    manifest.write_text(json.dumps(dict(selections=[selection])))
    failed=subprocess.run(command,capture_output=True,text=True)
    assert failed.returncode and 'differs from saved provenance' in failed.stderr
    selection['expected_transitions'][str(terminal)]=transitions
    selection['saved_cut']['cut']['sha256']='0'*64
    manifest.write_text(json.dumps(dict(selections=[selection])))
    failed=subprocess.run(command,capture_output=True,text=True)
    assert failed.returncode and 'checksum differs' in failed.stderr


def _matrix(size, bonds):
    weights = np.zeros((size, size))
    for a, b in bonds:
        weights[a, b] = weights[b, a] = 1.
    return weights
