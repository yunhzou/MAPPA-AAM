---
name: mappa-growth-viewer
description: Open or build the existing MAPPA 3D fragment-growth replay from saved search archives, work-pool cuts or traces. Use for seed and extension decisions, saturation candidates, fragment locking and saved-branch inspection. For final mappings use mappa-viewer.
---

# MAPPA growth viewer

Resolve this skill's symlink, then go up two directories for the active checkout.
Use its established UI and `docs/SEARCH_TRAJECTORIES.md`; do not create another
viewer or reconstruct growth from a final mapping. On this Mac the project Python
is `/Users/yunhengz/Documents/Codex/2026-09-11/aa/work/aam-event-env/bin/python`.

## Fast reusable workflow

1. Open the requested existing `algorithm_trajectory.html` directly when available.
2. With a saved verified trace, render only. Use `display-trace.json` for generated
   conformers, or `trace.json` when it contains native 3D coordinates:
   `PYTHONPATH=src python tools/render_search_trajectory.py display-trace.json algorithm_trajectory.html`.
   Do not render an all-zero raw trace without `--geometry-cache DIRECTORY`.
3. Otherwise select the original saved branch in a manifest and run:
   `python tools/build_search_trajectory.py manifest.json output_directory --open`.

The build command replays only the selected fragment calls and checks their saved
fragments, deferred edges and shared placement relationship. It does not run a
seed search or cut sweep. A mismatch is an error, not permission to substitute a
new trajectory. Select explicit terminals for large archives rather than replaying
all branches. Keep case-specific data in the private manifest.

A manifest accepts AAM checkpoint/JSON `archive` selections or the original typed
work-pool `saved_cut` pair: `input` and `cut` Blob records with path, offset, size
and SHA-256. The maintained command performs this adaptation; do not write a
one-off conversion script. For a catalogue branch, use its family provenance to
find cut/context/terminal/transitions and the completed `MatchingResult` to find
the input and cut ranges. Set `expected_transitions` when selecting that history.
The guide documents the exact manifest. Only copy needed ranges from remote results,
retain checksums, and keep infrastructure details in ignored private storage.

## Current representation and display

The capture owner is `mappa.presentation.trajectory`; rendering is
`mappa.presentation.viewers.growth_trace_html` and
`src/mappa/presentation/static/growth_trace.html`. Extend those components rather
than implementing a parallel UI. In-memory results use `aam_growth_html` with an
explicit context and terminals.

Growth shows one witness. Saturation shows complete compressed `Fragment`
relations, including relaxed source/target self-relationships and the shared global
placement action. Target highlights are exact correlated per-atom domains, not
independent choices. Preserve all candidates from the replayed operation; do not
truncate alternatives or replace local relations with global endpoint orbits.
The displayed committed mapping is one saved representative, not the whole family.

Keep the two rotatable 3D endpoints, playback/stepping, candidate selection,
numbered fragment colors and persistent locks. **Growth fragment** jumps directly
to a stage; **Final mapping** shows its terminal and bond overlays. Stage labels
count newly assigned atoms; joined stages need not equal final fragment partitions.
Explain hydrogen-only fragments separately when interpreting large counts.
Red indicates losses on R, green gains on P, including mapped/unmatched boundaries.
A sweep cut is a search condition, not automatically a bond event.

## Geometry and validation

Reuse native 3D coordinates. For graph-only input set manifest `geometry_cache`
(or build `--geometry-cache`) to reuse the ordinary viewer's conformer cache.
The shared geometry module preserves indexed atoms, bonds and specified stereo.
Keep `trace.json` unchanged and render `display-trace.json`; label generated
coordinates as display conformers. Older traces can supply original endpoint
chemistry using render `--geometry-input`. Never replace failed 3D generation with
flat coordinates or change scientific inputs to make replay agree.

After adapting capture or UI, run the focused trajectory tests and:
`PYTHONPATH=src python tools/check_search_trajectory.py output_directory`.
Add `--trace display-trace.json` when that rendering copy is used.
`MANUSCRIPT_BROWSER` can select installed Chrome; on this Mac Playwright is
available via `PYTHONPATH=src:/private/tmp/mappa-browser-deps` and Chrome is
`/Applications/Google Chrome.app/Contents/MacOS/Google Chrome`.

Check candidate choice, persistent locks through the next stage and reverse seeks,
fragment navigation, playback, final mapping/bond colors, fixed geometry, rotation,
offline loading and browser errors. Inspect the saved screenshot after a lock.
Routine reopening needs no repeated audit. Deliver the absolute viewer link and
briefly identify the saved case/branch. When speed is at issue, report measured
capture/render times separately; opening existing HTML needs neither step.
