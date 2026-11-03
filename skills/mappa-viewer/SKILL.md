---
name: mappa-viewer
description: Open saved reactions and atom mappings in the established MAPPA 3D viewer, including Golden benchmark cases. Use for requests such as “show this case on MAPPA.” For fragment-growth replay use mappa-growth-viewer; for Reactivo pathway canvases use reactivo-test-verify-view.
---

# MAPPA reaction viewer

Resolve this skill's symlink, then go up two directories for the active MAPPA
checkout. Use its shared molecular renderer; do not build another UI or write
inline geometry-generation scripts. Viewing saved evidence does not run matching,
decoding, reference verification or a chemistry benchmark.

## Fast path

If the requested view already exists, open its HTML directly. For a collection,
`#case=1654` selects that case. Use `reactivo view DIRECTORY` for a directory
containing `view.html`, or open the HTML in the browser.

For selected Golden reference reactions, run the maintained command from the
checkout:

```sh
PYTHONPATH=src:bench python tools/render_mapping_comparison.py \
  private/local-golden-inputs/index.json \
  private/golden-viewer/view.html --golden-cases 1654 --open
```

Use the project Python environment. On this Mac it is
`/Users/yunhengz/Documents/Codex/2026-09-11/aa/work/aam-event-env/bin/python`.
Supply only the requested case IDs, in the desired display order. The command
loads their checksummed input/reference byte ranges, generates missing 3D display
geometry once, caches it by chemistry and generation settings, renders the shared
viewer and opens it. Repeating the command reuses the geometry. The source index
and scientific results remain unchanged. Bundle relocation, if required, uses
the existing `MAPPA_BUNDLE_ROOTS` mechanism with the saved relocation map.

For saved MAPPA alternatives already exported as comparison JSON:

```sh
PYTHONPATH=src python tools/render_mapping_comparison.py cases.json view.html --open
```

This displays the supplied records. Do not describe one witness or a ground-truth
reference as the complete decoded mapping relation. Keep reference and MAPPA
alternatives explicitly labelled. Missing records require an adapter to the saved
format, not a new mapping run.

Show bond changes by default: red breaking/weakening on R, green
forming/strengthening on P, with the event list visible. Derive overlays from the
displayed mapping and original bond weights, retaining mapped/unmatched boundaries.

## Geometry and delivery

Preserve native geometry, atom identities, mappings and unmatched atoms. Graph-only
Golden inputs use deterministic ETKDGv3 + UFF display conformers with hydrogens;
the reference view displays heavy atoms because hydrogen identities are not
annotated. Generated conformers are not experimental structures. Preserve specified
stereo and label unspecified stereochemistry as illustrative. Errors remain
explicit; never replace failed 3D generation with a flat drawing.

Deliver the viewer link immediately after generation and opening. Report the
measured render time when investigating workflow speed. Ordinary viewing needs no
repo-wide search, benchmark audit, exhaustive validation or new server. After a
renderer/adapter change, check the selected reaction, bond highlights and rotation
in a browser once. The shared API and comparison format are documented in
`docs/VIEWERS.md`; geometry belongs to `mappa.presentation.geometry`, and Golden
input selection belongs to `bench/golden.py`.
