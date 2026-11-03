# Viewer styles

There are **two shared viewer styles**. New reports provide data or a layout
adapter; they do not introduce another stylesheet or benchmark-specific skin.

| Style | Purpose | Authoritative files |
| --- | --- | --- |
| `reaction` | R/P comparison, TS structures/modes, mapping diagnostics and recorded growth | `src/mappa/presentation/static/reaction_viewer.css`; original R/P/TS layout in `reaction_viewer.html` |
| `batch` | AAM and TS result navigation | `src/mappa/presentation/static/batch_viewer.css` |

The molecular presentation restores the original white panels, yellow mechanism
buttons, thin molecular sticks, R-order/spatial-alignment controls and TS mode
animation from `stable/legacy-aam:src/rxn_core/pipeline.py`. The separate dark
`elementary_comparison_viewer.html` skin and its derived missing-pattern viewer
were rejected and removed. Its reports now use the shared original renderer.

`mappa.presentation.viewers.reaction_html` renders one reaction; `collection_html` adds a
case selector with one shared copy of the renderer and libraries. Existing
diagnostic layouts use `style_document` / `viewer_style`. Saved native geometries
remain unchanged. When generated 3D conformers are requested for graph-only
inputs, label their generation method, preserve atom identities and chemical
stereochemistry, and retain them separately from the scientific inputs. Rigid
placement of separate molecules is a display choice, not a predicted complex.
Batch layouts retain their existing result-navigation behavior.
Different layouts and chemical highlight colors are not independent skins.

The reusable **algorithm search trajectory** pipeline is documented in
[SEARCH_TRAJECTORIES.md](SEARCH_TRAJECTORIES.md). Use `from mappa.presentation.viewers import aam_growth_html` for direct in-memory
results; archive capture and `growth_trace_html` remain supported. The executed
[AAM notebook](AAM_SIMPLE.ipynb) shows a separate animation for each sweep.
The PR7 example is `reports/pr7_search_trajectory_20260911/algorithm_trajectory.html`.
Its endpoint coordinates remain fixed; playback follows growth and branching.

Bond-event overlays use **red for breaking or weakening** and **green for forming
or strengthening**, always in the R → P direction. Do not add separate colors for
bond-order changes. In the R/P comparison, **R shows only red changes** and
**P shows only green changes**; each event is highlighted once across the two
panels. Event names and numerical WBO changes remain in the details.

## Saved comparisons

Open selected Golden references in one command, without loading catalogues:

```sh
PYTHONPATH=src:bench python tools/render_mapping_comparison.py \
  private/local-golden-inputs/index.json private/golden-viewer/view.html \
  --golden-cases 1654 824 --open
```

The input index supplies checksummed original endpoints and reference mappings.
Native coordinates are preserved. Missing coordinates use deterministic 3D
display conformers through `mappa.presentation.geometry`; its chemistry-keyed
cache is reused on subsequent views. The reference adapter displays heavy atoms
and retains unmatched atoms. The HTML and adjacent JSON identify the reference
and geometry provenance. No matching, decoding or verification is invoked.
Use `--geometry-cache PATH` to share display geometry across report directories.
The [MAPPA viewer skill](../skills/mappa-viewer/SKILL.md) provides the fast workflow.

`tools/render_mapping_comparison.py INPUT.json OUTPUT.html` renders saved full
mapping records without starting a benchmark search. Historical comparison
exports remain with their saved reports; their retired generators are available
at the reports’ recorded source revisions.

The active missing-pattern comparison is
`reports/holdout_missing_pattern_seeds_20260910/viewer.html`.
Product fitting is a rigid display transform of the **recorded mapping**;
switching to native P restores its original coordinates. Explicit WBO connectivity
is used when provided. No search result is replaced by display fitting.

## Maintenance and validation

Edit the shared renderer/styles and generate new output through the public viewer
functions. Preserve source mappings, metrics and provenance. Static manuscript
figures and molecular highlight colors do not define additional viewer skins.
Standalone pages embed the shared CSS for offline viewing.

`tools/check_viewer_consolidation.py` checks saved witnesses, geometry,
batch/growth layouts and offline operation. Its output goes under
`reports/viewer_validation/`. Historical migrations are recorded in Git; no
migration script is part of the current rendering workflow.
