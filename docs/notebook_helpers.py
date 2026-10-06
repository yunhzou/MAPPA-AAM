"""Display-only helpers for AAM_SIMPLE.ipynb (formal bond orders, explicit H)."""

from rdkit import Chem
from rdkit.Chem import AllChem
from mappa import MolecularEndpoint


def endpoint(smiles, label):
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(molecule, randomSeed=7) != 0:
        raise ValueError("Example coordinate embedding failed")
    # RDKit embeds disconnected components near the origin; separate them for display.
    import numpy as np

    coordinates = molecule.GetConformer().GetPositions()
    cursor = 0.0
    for fragment in Chem.GetMolFrags(molecule):
        atoms = list(fragment)
        coordinates[atoms, 0] += cursor - coordinates[atoms, 0].min()
        cursor = coordinates[atoms, 0].max() + 2.5
    coordinates -= np.mean(coordinates, axis=0)
    for atom, xyz in enumerate(coordinates):
        molecule.GetConformer().SetAtomPosition(atom, xyz)
    return MolecularEndpoint(
        tuple(atom.GetSymbol() for atom in molecule.GetAtoms()),
        molecule.GetConformer().GetPositions(),
        Chem.GetAdjacencyMatrix(molecule, useBO=True),
        label=label,
        metadata={
            "molblock": Chem.MolToMolBlock(molecule),
            "fragment_count": len(Chem.GetMolFrags(molecule)),
        },
    )


def show_mapping(problem, mapping):
    """Original atom indices; matched atoms have the same colors in both panes."""
    import py3Dmol

    palette = ["#4e79a7", "#e15759", "#59a14f", "#b07aa1", "#f28e2b", "#76b7b2"]
    inverse = {p: r for r, p in mapping.items()}
    view = py3Dmol.view(width=800, height=360, viewergrid=(1, 2))
    for side, ep in enumerate((problem.reactant, problem.product)):
        pane = (0, side)
        view.addModel(ep.metadata["molblock"], "mol", viewer=pane)
        view.setBackgroundColor("white", viewer=pane)
        for atom, (symbol, position) in enumerate(zip(ep.elements, ep.coordinates)):
            source = atom if side == 0 else inverse.get(atom)
            color = palette[source % len(palette)] if source is not None else "#aaaaaa"
            view.setStyle(
                {"index": atom},
                {"stick": {"color": color}, "sphere": {"scale": 0.25, "color": color}},
                viewer=pane,
            )
            view.addLabel(
                f"{symbol}{atom}",
                {
                    "position": dict(zip("xyz", map(float, position))),
                    "fontColor": "#222222",
                    "backgroundColor": "white",
                    "backgroundOpacity": 0.7,
                    "fontSize": 12,
                },
                viewer=pane,
            )
        view.addLabel(
            ep.label,
            {
                "useScreen": True,
                "position": {"x": 10, "y": 10},
                "fontColor": "black",
                "backgroundColor": "white",
            },
            viewer=pane,
        )
        view.zoomTo(viewer=pane)
        view.zoom(1.5 if ep.metadata["fragment_count"] == 1 else 1.0, viewer=pane)
    return view
