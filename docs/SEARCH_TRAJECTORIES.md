# Reusable search-trajectory pipeline

This pipeline animates the **algorithm's search history**: seed selection,
one-witness growth attempts, accepted and deferred edges, final compressed
placements, and committed fragments. R and P keep
their original coordinates throughout. It uses the shared white 3D reaction
viewer style.

The stages are independent:

1. Select saved archives and contexts in a JSON manifest.
2. Replay the selected fragment calls with their producing matcher and verify
   the placements, fragments, and deferred boundaries against the archive.
3. Store the recorded events and candidates in `trace.json`.
4. Render `trace.json` into a self-contained `algorithm_trajectory.html`.

No cut sweep is run. One recorded history is followed per selected terminal.
Growth frames show one feasible witness. Saturation frames include every
retained compressed placement and its symmetry blocks, without enumerating atom
permutations. A growth witness can change when an extension needs reassignment. The selected saved branch determines
which placement is committed before the next fragment starts.

## Build any saved case

From the repository root, using its Python environment:

```bash
python tools/build_search_trajectory.py manifest.json output_directory
```

A minimal manifest is:

```json
{
  "selections": [
    {"archive": "saved/aam.pkl.gz", "context": 0}
  ]
}
```

Archives can be AAM checkpoints, AAM JSON records, or gzipped AAM JSON records.
Paths resolve relative to the manifest. Endpoint coordinates must be present.
Multiple selections must refer to exactly the same R/P endpoints and atom order.
The matching tolerance, graph threshold, branch cap, seed order, cuts, and prior
committed state are read from the selected archive context.

Optional manifest fields:

| Field | Purpose |
| --- | --- |
| `title` | Page title |
| `engine_source` | Producing engine's `src` directory; defaults to this checkout |
| `engine_commit` | Human-supplied reference revision; runtime source hashes are recorded separately |
| `key_atoms` | Reactant atom indices labeled by default |
| `watch_targets` | Product atom indices whose occupancy is shown below the structures |
| `event_tolerance` | Bond-order event display threshold, default 0.5 |
| `selections[].terminals` | Specific terminals; defaults to all terminals in that context |
| `selections[].label` | Run button label |
| `selections[].path_labels` | Object mapping terminal IDs to display labels |
| `selections[].default_terminal` | Branch selected when opening this run |
| `selections[].focus` | Optional jump point: `{"source_edge": [a, b], "target": p}` |

Capture runs in a separate process and reads the matcher's own trace events.
The actual matching functions decide compatibility; a replay
that differs from the saved placement fails. The bundle records archive hashes,
runtime source hashes, and validation results. Native archives use the matching
Python trace implementation and must pass the same placement verification.

## Re-render an existing trace

```bash
PYTHONPATH=src python tools/render_search_trajectory.py output_directory/trace.json output_directory/algorithm_trajectory.html
```

This stage does not replay matching or load the original archives. The generated
HTML embeds its molecular renderer, shared stylesheet, and trace, so it works
offline. Playback supports frame stepping, scrubbing, candidate selection,
witness inspection, branch selection, focus/reset, and final bond overlays.
Final events use red losses on R and green gains on P. Black/gray dashed bonds
mark the tested edge, sweep cuts, or deferred search boundaries.

```bash
PYTHONPATH=src python tools/check_search_trajectory.py output_directory
```

The browser check needs Playwright and Chromium. `MANUSCRIPT_BROWSER` can select
an installed Chromium executable. It checks playback, candidates, fixed endpoint
coordinates, final event colors, offline operation, and mobile width.

`tests/test_search_trajectory.py` also exercises a separate four-atom example,
non-default matching settings, JSON/checkpoint loading, input preservation,
manifest orchestration, render-only behavior, and rejection of mixed endpoints.
