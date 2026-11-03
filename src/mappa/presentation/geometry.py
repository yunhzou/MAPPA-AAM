"""Deterministic display conformers, separate from scientific endpoint data."""

import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np


_CONFIG = dict(method="ETKDGv3 + UFF", random_seed=42,
               conformers_per_component=3, embedding_timeout_seconds=30,
               embedding_max_iterations=1000, uff_max_iterations=2000,
               num_threads=1, component_gap=3.0)


def _indexed_molecule(endpoint):
    from rdkit import Chem

    features = endpoint["metadata"]["scoring_features"]
    heavy = features["heavy"]
    if list(heavy) != list(range(len(heavy))):
        raise ValueError("Generated geometry requires heavy atoms before hydrogens")
    molecule = Chem.RWMol()
    for symbol, charge, isotope, hydrogens, stereo in features["colors"]:
        atom = Chem.Atom(symbol)
        atom.SetFormalCharge(int(charge))
        atom.SetIsotope(int(isotope))
        atom.SetNumExplicitHs(int(hydrogens))
        atom.SetNoImplicit(True)
        if stereo:
            if stereo not in {"R", "S"}:
                raise ValueError(f"Unsupported atom stereochemistry: {stereo}")
            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
        molecule.AddAtom(atom)
    bond_types = {1: Chem.BondType.SINGLE, 1.5: Chem.BondType.AROMATIC,
                  2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
    for left, right, (order, _) in features["bonds"]:
        molecule.AddBond(int(left), int(right), bond_types[float(order)])
    molecule = molecule.GetMol()
    Chem.SanitizeMol(molecule)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True,
                              flagPossibleStereoCenters=True)
    for atom, color in zip(molecule.GetAtoms(), features["colors"]):
        stereo = color[4]
        if stereo and atom.GetProp("_CIPCode") != stereo:
            atom.InvertChirality()
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True,
                              flagPossibleStereoCenters=True)
    for left, right, (_, stereo) in features["bonds"]:
        if stereo == "STEREONONE":
            continue
        if stereo not in {"STEREOE", "STEREOZ"}:
            raise ValueError(f"Unsupported bond stereochemistry: {stereo}")
        bond = molecule.GetBondBetweenAtoms(int(left), int(right))
        neighbors = [max((a for a in molecule.GetAtomWithIdx(int(center)).GetNeighbors()
                          if a.GetIdx() != other), key=lambda a: a.GetIntProp("_CIPRank"))
                     for center, other in ((left, right), (right, left))]
        bond.SetStereoAtoms(*(a.GetIdx() for a in neighbors))
        bond.SetStereo(getattr(Chem.BondStereo, stereo))
    molecule = Chem.AddHs(molecule)
    if ([a.GetSymbol() for a in molecule.GetAtoms()] != list(endpoint["elements"])
            or not np.array_equal(Chem.GetAdjacencyMatrix(molecule, useBO=True),
                                  np.asarray(endpoint["wbo"]))):
        raise ValueError("Declared chemistry does not reproduce the indexed endpoint")
    return molecule


def _check_stereo(molecule, conformer):
    from rdkit import Chem

    expected_atoms = {a.GetIdx(): a.GetProp("_CIPCode") for a in molecule.GetAtoms()
                      if a.HasProp("_CIPCode")}
    expected_bonds = {b.GetIdx(): b.GetStereo() for b in molecule.GetBonds()
                      if b.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ)}
    observed = Chem.Mol(molecule)
    Chem.AssignStereochemistryFrom3D(observed, confId=conformer, replaceExistingTags=True)
    Chem.AssignStereochemistry(observed, cleanIt=True, force=True)
    if (any(not observed.GetAtomWithIdx(i).HasProp("_CIPCode")
            or observed.GetAtomWithIdx(i).GetProp("_CIPCode") != value
            for i, value in expected_atoms.items())
            or any(observed.GetBondWithIdx(i).GetStereo() != value
                   for i, value in expected_bonds.items())):
        raise ValueError("Generated conformer does not preserve declared stereochemistry")


def _generate(endpoint):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    molecule = _indexed_molecule(endpoint)
    indices = []
    components = Chem.GetMolFrags(molecule, asMols=True, fragsMolAtomMapping=indices)
    coordinates = np.empty((molecule.GetNumAtoms(), 3))
    records = []
    offset = 0.0
    for component, atoms in zip(components, indices):
        if not AllChem.UFFHasAllMoleculeParams(component):
            raise ValueError("UFF does not support this endpoint chemistry")
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = _CONFIG["random_seed"]
        parameters.numThreads = _CONFIG["num_threads"]
        parameters.maxIterations = _CONFIG["embedding_max_iterations"]
        parameters.timeout = _CONFIG["embedding_timeout_seconds"]
        conformers = list(AllChem.EmbedMultipleConfs(
            component, numConfs=_CONFIG["conformers_per_component"], params=parameters))
        if not conformers:
            raise ValueError("ETKDG did not generate a conformer")
        energies = AllChem.UFFOptimizeMoleculeConfs(
            component, numThreads=1, maxIters=_CONFIG["uff_max_iterations"])
        converged = [(energy, conformer) for conformer, (status, energy)
                     in zip(conformers, energies) if status == 0]
        if not converged:
            raise ValueError("UFF did not converge for this endpoint")
        energy, conformer = min(converged)
        _check_stereo(component, conformer)
        positions = np.array(component.GetConformer(conformer).GetPositions())
        positions -= positions.mean(axis=0)
        positions[:, 0] += offset - positions[:, 0].min()
        offset = float(positions[:, 0].max()) + _CONFIG["component_gap"]
        coordinates[list(atoms)] = positions
        records.append(dict(atoms=list(atoms), converged=True, energy=energy))
    coordinates -= coordinates.mean(axis=0)
    return dict(coordinates=coordinates.tolist(), components=records)


def display_endpoint(endpoint, *, cache: Path) -> dict:
    """Keep supplied geometry or cache a separately labelled 3D display conformer.

    Graph-only endpoints require indexed formal chemistry in ``scoring_features``.
    Components are packed for visibility, not modeled as an intermolecular complex.
    """
    coordinates = np.asarray(endpoint["coordinates"], dtype=float)
    if coordinates.shape != (len(endpoint["elements"]), 3) or not np.isfinite(coordinates).all():
        raise ValueError("Endpoint coordinates must be finite and have shape (atoms, 3)")
    result = dict(elements=list(endpoint["elements"]), coordinates=coordinates.tolist(),
                  wbo=np.asarray(endpoint["wbo"]).tolist(),
                  atom_indices=list(endpoint.get("atom_indices", range(len(coordinates)))))
    if np.any(coordinates):
        result["geometry"] = dict(source="native", display_only=False)
        return result

    from rdkit import rdBase

    chemistry = dict(elements=result["elements"], wbo=result["wbo"],
                     features=endpoint["metadata"]["scoring_features"])
    identity = dict(schema=1, chemistry=chemistry, config=_CONFIG, rdkit=rdBase.rdkitVersion)
    key = hashlib.sha256(json.dumps(identity, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()
    cache = Path(cache)
    path = cache / f"{key}.json"
    if path.exists():
        saved = json.loads(path.read_text())
        if saved["identity"] != identity:
            raise ValueError("Display geometry cache identity differs")
    else:
        saved = dict(identity=identity, **_generate(endpoint))
        cache.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode="w", dir=cache, delete=False) as stream:
            temporary = Path(stream.name)
            json.dump(saved, stream, separators=(",", ":"), allow_nan=False)
        os.replace(temporary, path)
    positions = np.asarray(saved["coordinates"], dtype=float)
    if positions.shape != coordinates.shape or not np.isfinite(positions).all():
        raise ValueError("Invalid cached display coordinates")
    result["coordinates"] = positions.tolist()
    result["geometry"] = dict(source="generated", display_only=True, cache_key=key,
                              method=_CONFIG["method"], config=dict(_CONFIG),
                              rdkit=rdBase.rdkitVersion, components=saved["components"],
                              note="Generated 3D conformer; rigid component packing is for display only.")
    return result
