"""Generated coordinates preserve indexed chemistry and declared stereochemistry."""

import copy
import json

import numpy as np
import pytest
from rdkit import Chem

from mappa.presentation.geometry import display_endpoint


def endpoint(smiles):
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    heavy = [a.GetIdx() for a in molecule.GetAtoms() if a.GetSymbol() != "H"]
    features = dict(heavy=heavy, colors=[
        [a.GetSymbol(), a.GetFormalCharge(), a.GetIsotope(),
         a.GetTotalNumHs(includeNeighbors=True),
         a.GetProp("_CIPCode") if a.HasProp("_CIPCode") else ""]
        for a in molecule.GetAtoms() if a.GetSymbol() != "H"], bonds=[
            [b.GetBeginAtomIdx(), b.GetEndAtomIdx(),
             [b.GetBondTypeAsDouble(), str(b.GetStereo())]]
            for b in molecule.GetBonds() if b.GetBeginAtomIdx() in heavy
            and b.GetEndAtomIdx() in heavy])
    return molecule, dict(elements=[a.GetSymbol() for a in molecule.GetAtoms()],
                          coordinates=np.zeros((molecule.GetNumAtoms(), 3)).tolist(),
                          wbo=Chem.GetAdjacencyMatrix(molecule, useBO=True).tolist(),
                          metadata=dict(scoring_features=features))


def test_native_geometry_is_unchanged_without_generation(tmp_path):
    source = dict(elements=["C", "O"], coordinates=[[1, 2, 3], [2, 3, 4]],
                  wbo=[[0, 2], [2, 0]], atom_indices=[5, 8])
    saved = copy.deepcopy(source)
    result = display_endpoint(source, cache=tmp_path / "cache")
    assert source == saved
    assert {key: result[key] for key in source} == source
    assert result["geometry"] == dict(source="native", display_only=False)
    assert not (tmp_path / "cache").exists()


@pytest.mark.parametrize("smiles", ["[13CH3][C@H](O)F", "C/C(c1ccccc1)=N/N",
                                    "C/C(c1ccccc1)=N\\N", "[C-]#[N+]C(C)(C)C"])
def test_generated_geometry_preserves_stereo_charge_isotope_and_indices(smiles, tmp_path):
    original, source = endpoint(smiles)
    saved = copy.deepcopy(source)
    result = display_endpoint(source, cache=tmp_path)
    assert source == saved
    assert result["elements"] == source["elements"]
    assert result["wbo"] == source["wbo"]
    assert result["atom_indices"] == list(range(original.GetNumAtoms()))
    assert result["geometry"]["source"] == "generated"
    assert np.isfinite(result["coordinates"]).all()
    assert np.linalg.matrix_rank(np.array(result["coordinates"]) -
                                 np.mean(result["coordinates"], axis=0)) == 3
    observed = Chem.Mol(original)
    conformer = Chem.Conformer(original.GetNumAtoms())
    for atom, position in enumerate(result["coordinates"]):
        conformer.SetAtomPosition(atom, position)
    observed.AddConformer(conformer)
    Chem.AssignStereochemistryFrom3D(observed, replaceExistingTags=True)
    Chem.AssignStereochemistry(observed, cleanIt=True, force=True)
    for atom in original.GetAtoms():
        if atom.HasProp("_CIPCode"):
            assert observed.GetAtomWithIdx(atom.GetIdx()).GetProp("_CIPCode") == atom.GetProp("_CIPCode")
    for bond in original.GetBonds():
        if bond.GetStereo() != Chem.BondStereo.STEREONONE:
            assert observed.GetBondWithIdx(bond.GetIdx()).GetStereo() == bond.GetStereo()


def test_cache_reuses_geometry_and_distinguishes_stereochemistry(tmp_path, monkeypatch):
    _, source = endpoint("C/C=C/C")
    result = display_endpoint(source, cache=tmp_path)
    _, other = endpoint("C/C=C\\C")
    different = display_endpoint(other, cache=tmp_path)
    assert result["geometry"]["cache_key"] != different["geometry"]["cache_key"]
    monkeypatch.setattr("mappa.presentation.geometry._generate",
                        lambda _: pytest.fail("cached geometry was regenerated"))
    assert display_endpoint(source, cache=tmp_path) == result
    assert len(list(tmp_path.iterdir())) == 2


def test_invalid_cache_and_inconsistent_indexed_chemistry_fail(tmp_path):
    _, source = endpoint("CO")
    result = display_endpoint(source, cache=tmp_path)
    path = tmp_path / f"{result['geometry']['cache_key']}.json"
    cached = json.loads(path.read_text())
    cached["coordinates"] = [[0, 0, 0]]
    path.write_text(json.dumps(cached))
    with pytest.raises(ValueError, match="cached display coordinates"):
        display_endpoint(source, cache=tmp_path)
    source["elements"][0] = "N"
    with pytest.raises(ValueError, match="indexed endpoint"):
        display_endpoint(source, cache=tmp_path)
