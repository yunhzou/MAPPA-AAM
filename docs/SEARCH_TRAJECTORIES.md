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
canonical placement, its local relation and the shared global placement action,
without enumerating atom permutations. Displayed target options include the
other components represented by that action; these are correlated alternatives,
not independent atom choices. A growth witness can change when an extension needs reassignment. The selected saved branch determines
which placement is committed before the next fragment starts.

## Build any saved case

From the repository root, using its Python environment:

```bash
python tools/build_search_trajectory.py manifest.json output_directory --open
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
| `geometry_cache` | Optional display-conformer cache, relative to the manifest |
| `selections[].terminals` | Specific terminals; defaults to all terminals in that context |
| `selections[].label` | Run button label |
| `selections[].path_labels` | Object mapping terminal IDs to display labels |
| `selections[].default_terminal` | Branch selected when opening this run |
| `selections[].expected_transitions` | Optional terminal-ID to transition-list map; rejects a different recorded history |
| `selections[].focus` | Optional jump point: `{"source_edge": [a, b], "target": p}` |

Capture runs in a separate process and reads the matcher's own trace events.
The actual matching functions decide compatibility; a replay
that differs from the saved fragment or shared placement relationship fails. The bundle records archive hashes,
runtime source hashes, and validation results. Native archives use the matching
Python trace implementation and must pass the same placement verification.

## Saved work-pool cuts

Instead of `archive`, a selection can reference the original typed `SearchInput`
and `CutResult` byte ranges. Paths are relative to the manifest; the loader checks
each range's SHA-256 before reading it. No new search or catalogue reconstruction
is needed:

```json
{
  "geometry_cache": "../geometry",
  "selections": [{
    "saved_cut": {
      "input": {"path": "saved.bundle", "offset": 0, "size": 123, "sha256": "INPUT_SHA256"},
      "cut": {"path": "saved.bundle", "offset": 123, "size": 456, "sha256": "CUT_SHA256"}
    },
    "context": 0,
    "terminals": [42],
    "expected_transitions": {"42": [0, 3, 8]}
  }]
}
```

The values above illustrate the schema; use the saved references and provenance.
For benchmark catalogues, a family's provenance identifies its cut, context,
terminal and transitions. The completed `MatchingResult` supplies that cut's
reference and its original input reference. Retrieve only those ranges when the
results live remotely, preserving checksums and keeping operational paths private.
Catalogue families alone cannot reconstruct a growth trajectory.

The adapter preserves the original graph, configuration and endpoints. Its
placeholder metrics are explicitly unmeasured. Capture records the source ranges
and validates the selected history, fragments, deferred edges and placement groups.
The branch diagram includes the selected histories; saturation retains all
compressed candidates of each replayed operation, including local relations and
shared placement actions.

## Graph-only inputs and 3D display

Use `geometry_cache` in the manifest, or `--geometry-cache DIRECTORY` on the build
command. The existing `mappa.presentation.geometry` owner generates and caches
display conformers from indexed chemistry metadata. It preserves atom order,
bond weights and specified stereochemistry. Coordinates remain fixed in playback.
`trace.json` retains scientific coordinates; `display-trace.json` records the
rendering copy and geometry provenance. These conformers are illustrative, not
experimental structures, and do not participate in matching.

For an older trace without indexed chemistry metadata, render with
`--geometry-input original-endpoints.json --geometry-cache DIRECTORY`. The
endpoint file must contain matching `reactant` and `product` records; mismatched
atom order or weights are rejected.

## Re-render an existing trace

```bash
PYTHONPATH=src python tools/render_search_trajectory.py output_directory/trace.json output_directory/algorithm_trajectory.html
```

This stage does not replay matching or load the original archives. The generated
HTML embeds its molecular renderer, shared stylesheet, and trace, so it works
offline. Playback supports frame stepping, scrubbing, candidate selection,
witness inspection, branch selection, focus/reset, and final bond overlays.
Use **Growth fragment** to jump to a locking stage and **Final mapping** to inspect
the recorded terminal. Fragment labels count newly assigned atoms, including
hydrogens; joined growth stages need not equal the final fragment partition.
Final events use red losses on R and green gains on P, including mapped/unmatched
boundary bonds in partial mappings. Black/gray dashed bonds
mark the tested edge, sweep cuts, or deferred search boundaries.

```bash
PYTHONPATH=src python tools/check_search_trajectory.py output_directory
```

When rendering generated geometry, add `--trace display-trace.json` so the browser
check compares against the actual display coordinates. Opening an existing HTML
requires neither capture nor rendering. `validation.json` records capture and
render times; browser checks produce a separate receipt.

The browser check needs Playwright and Chromium. `MANUSCRIPT_BROWSER` can select
an installed Chromium executable. It checks playback, candidates, fixed endpoint
coordinates, final event colors, offline operation, and mobile width.

`tests/test_search_trajectory.py` also exercises a separate four-atom example,
non-default matching settings, JSON/checkpoint loading, input preservation,
manifest orchestration, render-only behavior, and rejection of mixed endpoints.
