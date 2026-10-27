---
name: mappa-growth-viewer
description: Open or build the established MAPPA fragment-growth trajectory viewer from saved search checkpoints or traces. Use for growth replay, seed decisions, candidate branching, fragment locking, and inspecting a saved branch. For ordinary final-mapping or Reactivo pathway viewing, use the normal MAPPA/Reactivo viewer workflow instead.
---

# MAPPA growth viewer

Use the existing interactive growth viewer, with two rotatable molecular endpoints,
playback, frame stepping, candidate and target inspection, branch navigation,
persistent locked-fragment highlighting, and final bond-event overlays.

The canonical implementation is in the MAPPA checkout containing this skill:
resolve this file's symlink, then go up two directories. The current checkout is
`/Users/yunhengz/Desktop/AAM Writing`. Read that checkout's
`docs/SEARCH_TRAJECTORIES.md` for the manifest schema and replay contract.

## Select the saved evidence

- Identify the requested case, cut/no-cut context, seed ordering, terminal branch,
  and event count from saved records. Distinguish a retained-family minimum from
  a global optimum. Do not select an unrelated convenient example.
- If a matching verified `trace.json` exists, reuse it and re-render only. Viewing
  does not authorize a new search, cut sweep, chemistry calculation, or monitor.
- Otherwise select the original AAM checkpoint/JSON and its producing engine.
  For archived native searches, use the matching Python replay and require exact
  agreement with saved placements, fragments, and deferred boundaries. A differing
  replay must be diagnosed, not presented as the original trajectory.
- Raw search graphs without endpoints/configuration need an explicit AAM archive
  adapter using the original saved inputs and settings. Do not infer a trajectory
  from only the final mapping. Preserve originals and record adaptation provenance.

## Build or refresh

Run from the canonical checkout using its Python environment:

```sh
python tools/build_search_trajectory.py manifest.json output_directory
```

For an existing verified trace, skip capture:

```sh
PYTHONPATH=src python tools/render_search_trajectory.py output_directory/trace.json output_directory/algorithm_trajectory.html
```

The renderer is `mappa.viewers.growth_trace_html`; its shared template is
`src/mappa/static/growth_trace.html`. For an in-memory AAM result use
`mappa.viewers.aam_growth_html` with an explicit context. Extend these components
when necessary; do not build a replacement molecular UI or duplicate the renderer
inside the skill. Keep dataset-specific selection data in the manifest.

Keep all saved compressed candidates and symmetry blocks; do not cap the display
to match the search branch cap. Temporary embedding counts and completed-frontier
branch counts describe different stages.

Default to all represented target options, including alternative product
components. Never present a single stored witness as the full matching relation.
Show per-atom options as correlated domains, not independently selectable pairs.
Use conditioned candidate groups and local pools; raw automorph-domain unions
and global endpoint orbits do not prove a current assignment is represented.
Preserve complete symmetry records at saturation. Label missing records explicitly.
Verify the visible atom highlights and mapping labels, not only candidate counts.

## Geometry and persistent locks

- Preserve endpoint atom order, graph, mappings, and fixed coordinates throughout
  playback. Reuse genuine saved 3D geometry where available.
- For zero/flat SMILES input coordinates, use a separately validated 3D display
  conformer in a rendering copy; retain the raw trace. Record geometry source/hash,
  verify atom identity, stereochemistry, bonds, and clashes. Display geometry does
  not participate in matching and must not alter the scientific input.
- Give each fragment a distinct color and show a numbered color legend. Use the
  same color on R and P during growth and after locking. Previously locked atoms
  retain their original fragment color if revisited by a later growth stage. The
  **Locked fragments** toggle defaults on. Derive locks from each frame's saved
  mapping, not accumulated playback state or a speculative candidate preview.
- Keep final bond losses red on R and gains green on P, independently of locking.

## Verify and show

Run `PYTHONPATH=src python tools/check_search_trajectory.py output_directory`
with Playwright/Chromium available; `MANUSCRIPT_BROWSER` can select installed
Chrome. If the local runtime uses Node Playwright instead, execute the same
behavioral checks in that runtime rather than claiming the Python checker ran.

Check locks at commit and the next fragment, backward seek/restart, candidate
selection, playback, final correspondence/event count, fixed geometry, rotation,
offline loading, and no browser errors. For multi-branch bundles also switch
branches and verify their own saved locks. Inspect a screenshot after a lock
while the next fragment grows; a screenshot of the initial frame cannot prove
persistence. Preserve the validation receipt beside the output.

Open `algorithm_trajectory.html` and give its absolute clickable path with a
short description of the selected saved branch and controls. No server is needed.
