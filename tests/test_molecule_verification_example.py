"""A whole-fragment hit still requires checking extra candidate connections."""
from pathlib import Path
import runpy
import numpy as np
import pytest

pytest.importorskip("rdkit")
from mappa import AAMProblem, AAMSearchConfig, MolecularEndpoint, search_aam

EXAMPLE = Path(__file__).resolve().parents[1] / "examples/molecule_verification"
API = runpy.run_path(str(EXAMPLE / "verify.py"))


def problem(extra=False):
    r=np.zeros((4,4))
    for a,b in [(0,1),(1,2),(2,3)]:r[a,b]=r[b,a]=1
    p=r.copy()
    if extra:p[0,3]=p[3,0]=1
    xyz=np.array([[i,0,0] for i in range(4)],float)
    return AAMProblem(MolecularEndpoint(tuple("CCCC"),xyz,r),
                      MolecularEndpoint(tuple("CCCC"),xyz,p))


def test_single_fragment_does_not_certify_extra_candidate_edges():
    value=problem(extra=True)
    aam=search_aam(value,AAMSearchConfig(sweep_cuts=False,iso_tolerance=.1),workers=1)
    paths=list(aam.graph.paths())
    assert paths
    terminal=paths[0].terminal
    placements=[aam.graph.transitions[i] for i in paths[0].transitions
                if aam.graph.transitions[i].match is not None]
    assert len(placements)==1
    check=API["connection_check"](value,dict(aam.graph.states[terminal].mapping))
    assert check==dict(complete=True,missing=0,extra=1,same_connectivity=False)


def test_missing_connection_is_rejected():
    value=problem(extra=True)
    value=AAMProblem(value.product,value.reactant)
    check=API["connection_check"](value,{i:i for i in range(4)})
    assert check==dict(complete=True,missing=1,extra=0,same_connectivity=False)


def test_partial_or_noninjective_assignment_is_not_verification():
    value=problem()
    for mapping in ({0:0,1:1},{0:0,1:1,2:2,3:2}):
        assert not API["connection_check"](value,mapping)["same_connectivity"]


def test_bundled_135_atom_xyz_candidate_passes(tmp_path):
    report=API["verify"](EXAMPLE/"target.xyz",EXAMPLE/"candidate.xyz",tmp_path)
    assert report['status']=='verified_connectivity'
    assert report['atoms']==135 and report['fragments']==1
    assert report['target_edges']==148 and report['missing']==report['extra']==0


def test_broken_xyz_maps_all_atoms_but_fails_connectivity(tmp_path):
    report=API["verify"](EXAMPLE/"target.xyz",EXAMPLE/"candidate-broken.xyz",tmp_path,capture=True)
    assert report['status']=='different_connectivity'
    assert report['complete'] and len(report['mapping'])==135
    assert report['target_components']==1 and report['candidate_components']==2
    assert report['target_edges']==148 and report['candidate_edges']==147
    assert report['fragments']==2 and report['missing']==1 and report['extra']==0
    assert (tmp_path/'trajectory.json.gz').is_file()


def test_broken_preparation_removes_only_the_documented_bridge():
    import json
    original=API["endpoint_from_xyz"](EXAMPLE/"candidate.xyz")
    broken=API["endpoint_from_xyz"](EXAMPLE/"candidate-broken.xyz")
    prep=json.loads((EXAMPLE/"break-preparation.json").read_text())
    expected=original.wbo.copy()
    a,b=prep['candidate_edge'];expected[a,b]=expected[b,a]=0
    assert original.elements==broken.elements and np.array_equal(broken.wbo,expected)
