# Verify a candidate XYZ against a target

This self-contained example checks whether two XYZ files have the same **inferred connectivity**, using an uncut MAPPA search. It is intended as a starting point for checking generated coordinates.

[Watch the 28-second 3D film](../../manuscript/animations/molecule_verification/molecule-verification.mp4) · [Animated preview and evidence](../../manuscript/animations/molecule_verification/README.md)

## Run

Install the repository with its `smiles` extra (for RDKit), then:

```bash
python examples/molecule_verification/verify.py --output verification-output
```

For your own structures:

```bash
python examples/molecule_verification/verify.py \
  --target target.xyz --candidate model_output.xyz --output verification-output
```

The importable example function is `verify(target, candidate, output, capture=False)`. It uses `MolecularEndpoint`, `AAMProblem`, `AAMSearchConfig` and `search_aam` from `mappa`. Add `--capture` to save the actual atom-by-atom trajectory for a successful match or a retained complete diagnostic witness.

A result of `verified_connectivity` requires all of the following for the same retained witness:

- Every target and candidate atom participates in an element-preserving bijection.
- The uncut search places the entire connected target in one fragment.
- Every inferred connection is preserved, and the candidate has no extra connections.

A single complete fragment **alone is not sufficient**, because the growth matcher need not constrain source nonedges. The final adjacency check covers both edges and nonedges. A result of `different_connectivity` means element counts, edge counts or component counts differ, which independently rules out graph isomorphism. A result of `not_verified` means this search did not establish the criterion; it is not a proof that no isomorphism exists. This example is for a connected target, rather than a multi-component mixture.

## What the demonstration contains

`target.xyz` is the 135-atom reactant of coordinate holdout case 68, `pr17.carbene.ins_ts6a`, including 85 heavy atoms and two Rh centers. `candidate.xyz` is a **controlled demonstration, not an actual generative-model output**. It contains four accepted torsional changes, small coordinate perturbations, a rigid rotation and a shuffled atom order. The original atom permutation is retained in the preparation record only for provenance; the matcher never reads it or uses anchors.

The result is one complete 135-atom fragment with 148 inferred connections, zero missing connections and zero extra connections. Its replay contains one seed placement and 134 single-atom growth events. Symmetry-related assignments remain compressed; the film follows one actual recorded representative.

Regenerate the candidate deterministically from the bundled target:

```bash
python examples/molecule_verification/prepare_candidate.py --output /tmp/mappa-candidate
```

The preparation script uses torsion indices specific to this example. It is not a general conformer generator. [Preparation provenance](preparation.json) records the changes, seed, target hash and RDKit version.

## Scope

XYZ files contain coordinates and elements, but no bond orders. Here, RDKit infers a binary connection matrix **independently from each XYZ**, using covalent radii with `useVdw=True` and factor 1.25. That matrix is supplied to `MolecularEndpoint.wbo`; it contains binary connectivity, not computed Wiberg bond orders. Metal coordination and close contacts can be sensitive to the inference rule. Use a chemically appropriate connectivity or bond-order model for your application.

The search uses `iso_tolerance=0.1`, `graph_floor=0.5`, one seed, branch cap 100, no cut sweep and one worker. Cuts are disabled because this task asks whether the whole connected structure is preserved; it is not a reaction-alternative search. The final adjacency comparison is exact on the inferred binary graphs.

Passing establishes connectivity under this inference rule. It does not certify bond orders, charge, stereochemistry, conformation quality or energetic stability. The reported search time excludes graph inference, trajectory replay and video rendering.

## Negative control: break one connection

[Watch the broken-molecule film](../../manuscript/animations/molecule_verification_broken/README.md).

`candidate-broken.xyz` keeps the same 135 atoms but moves a 21-atom group away, removing exactly one inferred C–C connection and adding none. Connectivity is inferred again from the changed coordinates. The target has one component and 148 connections; the broken candidate has two components and 147 connections.

```bash
python examples/molecule_verification/verify.py \
  --candidate examples/molecule_verification/candidate-broken.xyz \
  --output verification-broken-output --capture
```

The recorded diagnostic witness maps all 135 atoms in **two fragments (114 + 21)** and has **one missing connection**. The result is `different_connectivity`: complete atom coverage alone does not establish the same molecule. The component and edge counts establish the mismatch independently of the selected mapping or any search cap.

Regenerate the negative input with `python examples/molecule_verification/prepare_broken.py --output /tmp/mappa-broken`. The script checks that only the selected bridge is lost. Its supplied bridge indices are specific to this demonstration. [Preparation evidence](break-preparation.json) records the translation and moved atoms. Equivalent ligand attachments can be shuffled by symmetry; the film's final red marker follows the returned MAPPA witness, whose target atom indices can differ from the construction indices.
