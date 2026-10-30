# Importable AAM API

Start with [AAM_SIMPLE.ipynb](AAM_SIMPLE.ipynb): a self-contained six-atom reaction with embedded coordinates, bond orders, and display helpers. It inspects the raw AAM, unique event candidates, certified shuffles, anchors, both search directions, conditional matching, serialization, and one growth animation per sweep. It uses one worker and can run outside the checkout after installation. Saved tables and static molecule drawings render on GitHub; interactive views require a trusted Jupyter notebook. The longer [TUTORIAL.ipynb](TUTORIAL.ipynb) covers R/P geometry and TS analysis.

Install `python -m pip install -e ".[notebook]"` for the notebook, or `.[postprocessing]` for optional symbolic mapping queries without notebook packages. Event conversion itself needs no solver. A C++17 compiler builds the bookkeeping extensions. The optional fast growth engine requires `python src/mappa/mappa_core/native/build_engine.py`; see [native installation](../src/mappa/mappa_core/native/README.md). Ordinary `search_aam` can use the Python reference growth implementation. `compete_fragments` and `execution="reused_native"` / `"shared_policies"` require that optional engine and fail explicitly if unavailable. No xTB is needed when coordinates and bond-order matrices are supplied.

## Default pipeline

```python
from mappa import AAMProblem, MolecularEndpoint, AAMSearchConfig, search_aam, decode_events
from mappa.mappa_core.utils.io.records import decoded_events_record

problem = AAMProblem(MolecularEndpoint(elements_R, xyz_R, wbo_R),
                     MolecularEndpoint(elements_P, xyz_P, wbo_P))
aam = search_aam(problem, AAMSearchConfig(iso_tolerance=1.0))
decoded = decode_events(aam, threshold=0.5, metal_threshold=0.3)
witness = decoded.realize(0)
events = decoded.evaluate(0)
contributions = decoded.contributions(0)
record = decoded_events_record(decoded)
```

Decoding keeps the original compressed mappings, their correlated actions and
required growth edges. Only event interpretation changes with the thresholds.
Sparse local tables retain internal and boundary responses, including manual cuts.
Requested choices are applied in saved order. No rematching, solver search,
Cartesian branch combination or full symmetry enumeration is performed.

Records preserve the original catalogue, shared tables and representative events.
Representatives are not all symmetry alternatives; their minimum is not a global
minimum. Partial mappings count bonds touching mapped atoms and ignore wholly
unmatched bonds. Both ordinary and metal thresholds treat equality as unchanged.

## Experimental hierarchical refinement

```python
from mappa.mappa_core.post_processing.hierarchical_mapping import refine_mapping, RefinementConfig
refined = refine_mapping(aam, RefinementConfig(event_slack=5))
```

This separate operation rematches fragment pairs at tighter tolerances. It can
change the retained mapping space. It remains available for hierarchical mapping
experiments and is never called by the normal decoder. The historical tutorial
notebooks currently demonstrate this explicitly named experimental operation.

## Inspect growth directly from a result

```python
from html import escape
from IPython.display import HTML, display
from mappa.presentation.viewers import aam_growth_html

page = aam_growth_html(aam, context=None, event_threshold=0.5,
                       metal_event_threshold=0.3)
display(HTML('<iframe style="width:100%;height:1020px" srcdoc="'
             + escape(page, quote=True) + '"></iframe>'))
```

No archive is required. `context=None` includes all saved cut/seed contexts with terminals; an integer selects one. With a selected context, `terminals=[...]` selects its terminal paths. The HTML embeds its JavaScript and can also be saved as a standalone file. Each tab provides playback, candidate inspection, compressed symmetry, and the search graph.

This is verified local replay, not another sweep search or full decoding. It replays one recorded history per terminal and checks fragment results against the saved decisions. It does not enumerate symmetry permutations or impose another branch cap. Coordinates remain fixed. Capture temporarily observes matcher calls, so run it serially within a process. Event overlays use the same inclusive WBO-change thresholds and metal rule as decoding; the graph floor only distinguishes a broken/formed edge from an order decrease/increase in the drawing.

## Directions and anchors

```python
from mappa import search_aam_directions
config = AAMSearchConfig(anchors=((r0, p0),), seed_count=1)
runs = search_aam_directions(problem, config, direction="both", workers=1)
for run in runs:
    directed_decoding = decode_events(run.aam)
    for fid in range(len(directed_decoding.catalogue.families)):
        input_mapping = run.to_input_mapping(directed_decoding.realize(fid))
```

`direction` accepts `forward`, `reverse`, `both`, `smaller_first`, and `larger_first`. Size counts explicit atoms, including H. On ties smaller-first uses input R→P; larger-first uses the reverse. The example converts each direction separately; partial mappings use the same bond-difference rule. Both directions run sequentially under the same worker ceiling. Each direction has its own checkpoint directory when `intermediate_dir` is supplied. Anchors are always supplied in input R→P indexing and reversed automatically. Anchors constrain search, not just the displayed witness. The `reference` and `reused_native` backends accept anchors; the optional `shared_policies` backend currently rejects anchored searches.

Keep compressed results in their search direction. Invert only concrete witnesses. Event IDs are direction-local: merging already reduced forward/reverse event classes is **not** an exact common-frame bidirectional decoder. That joint event quotient is not implemented; use the two directed results for audit. `plan_aam_search` remains the lightweight smaller-first planning API.

## Conditional fragment and fixed-query isomorphism

`match_fragment(source_graph, target_graph, seed=..., context=FragmentMatchContext(locked_mapping=...), config=FragmentMatchConfig(...))` grows a seeded saturated fragment subject to existing assignments. Each extension constrains only the selected growth edge. Returned placements expose those accepted edges as `preserved_bonds`; other bonds between included atoms remain available in the original molecular graphs for event interpretation. Its result includes compressed placements and an explicit cap flag. It does not promise a globally largest common subgraph.

`match_weighted_subgraph(query, target, anchor_map=..., iso_tol=..., symmetry_wbo_tol=..., graph_floor=..., max_branches=..., node_policy=...)` matches a fixed query, with complete validated query placements and unsuccessful/capped search evidence. These two operations serve different purposes. Both accept `WeightedGraph` or prepared graphs from `build_graph`.

## Adjustable parameters

| API/configuration | Controls |
|---|---|
| `AAMSearchConfig` | `iso_tolerance=1.0`; `graph_floor=0.2`; `cut_floor=0.2`; `sweep_cuts=True`; `seed_count=1`; `random_seed=42`; `seed_selection="random"` or `"distance"`; `branch_limit=100`; `anchors=()`; `task_chunksize=1` |
| `search_aam` / directional helper | `workers`; `execution="reference"`, `"reused_native"`, or `"shared_policies"`; `intermediate_dir`; `resume`; `archive_format="json"` or `"checkpoint"` |
| `RefinementConfig` | `threshold=0.5`; `metal_threshold=0.3` (`None` disables the special metal rule); `max_events=None`; `event_slack=None`; `branch_limit=100`; `random_seed=42` |
| `query_assignments` | `assignments`; `same_events=True`; `seconds=10` |
| `FragmentMatchConfig` | `graph_floor`; `iso_tolerance`; `minimum_size`; `branch_limit`; `node_policy`; `allow_mapped_seed`; `orbit_dedup` |
| `FragmentMatchContext` | `locked_mapping`; `islands`; `deferred_edges`; optional reusable `source_orbits`, `target_orbits`, `growth_replay` |
| `CompetitionConfig` | `operation_budget=128`; `seconds=270`; `parent_limit=8`; `depth_limit=2`; `queue_limit=512`; `dependent_component_limit=8` |

The library and this notebook use the branch-limit default of 100. A cap bounds a growth call's live alternatives, not the number of saved histories or all permutations. The default pipeline uses one seed and cut sweep; no-sweep is an explicit ablation (`sweep_cuts=False`). The random seed controls reproducible, independent per-cut streams and is separate from the number of seed orderings.

Raw search does not classify or repair events. `decode_events` interprets original mappings; `refine_mapping` is separate experimental hierarchical matching. `AAMSearchConfig.event_threshold` and `metal_event_threshold` remain inputs to competition ordering; they do not configure decoding. TS analysis requires a `ReactionContext` supplied by the caller. The old mechanism-grouping API, its compiler and symmetry-repair options are removed. Full-result checkpoints containing those retired options must be read with their producing revision; current readers reject them rather than discard settings. Saved raw/finalized graph-only checkpoints remain readable.

## Experimental fragment competition (optional; default off)

```python
from mappa.mappa_core.core.search.branching.competition import CompetitionConfig, compete_fragments
augmented = compete_fragments(aam, CompetitionConfig(operation_budget=128))
decoded = refine_mapping(augmented.final_catalogue())
```

This exposes the bounded local B-priority takeover policy used in the coordinate campaign. It returns the untouched baseline plus independently validated repair families, counters, and proposal provenance. It performs no decoding or reference comparison. Growth/completion use the supplied AAM branch limit and tolerances; hard user anchors remain fixed. Parent selection uses the AAM configuration's observed event counts. Competition requires an explicit call and the native growth engine. It is excluded from the published method and its reported results.

**TODO (experimental):** establish a controlled coverage/cost benefit and improve proposal scheduling and deduplication before considering this extension for the default pipeline. The historical experiments are retained under `reports/`.

## Chirality and audit scope

**TODO: a general stereochemistry-preserving AAM search/decoder contract.** Element/WBO inputs do not encode CIP labels, E/Z geometry, formal charge, or isotope labels as first-class constraints. Raw AAM and signed-event decoding do not establish a general stereochemical guarantee. Fixed-query matching can accept richer node policies, but the full AAM endpoint contract still needs that extension.

The September 2026 Python audit added importable direction orchestration, signed-event decoding with witness-linked symmetry queries, explicit random/cut controls, optional dependency groups, and the competition adapter. It leaves the computational stages separate. Focused tests verify anchors in both orientations, threshold independence, partial/timeout status, same-class versus broader-family queries, checkpoint identity, and baseline retention under competition. This is software validation using small examples, not a new benchmark campaign.

## Removed orientation-selection API

`mappa.rp`, `align_reaction`, `select_rp_mappings`, and
`mappa.alignment.index_chirality` have been removed, including the legacy
orientation-graph selector. There are no deprecated wrappers, compatibility
aliases, automatic substitutions, or fallback options. Callers must use
`search_aam` and `refine_mapping` explicitly. Coordinate orientation is not
selected or certified by the default mapping workflow.

The CLI accepts `--stage aam`, `--stage events` (default), and `--stage ts`.
The former `rp`/`full` modes and `--post-workers` argument are rejected.
TS calls require a `ReactionContext`; no hidden R/P selection is performed.

The cached TS utility (`tools/run_typed_ts_case.py`) likewise requires
`--reaction-root`, containing `cases/<case_id>/reaction.json` records produced
with `reaction_record`. Its batch launcher requires `TS_REACTION_ROOT`.
Neither entry point computes missing R/P selections.

## Event decoding

```python
from mappa.mappa_core.post_processing.hierarchical_mapping import RefinementConfig, refine_mapping
from mappa.mappa_core.utils.io.records import refined_mapping_record

result = refine_mapping(aam, RefinementConfig(event_slack=0))
print(result.observed_minimum, result.complete)
for candidate in result.candidates:
    print(candidate.total, candidate.mapping, candidate.endpoint_events)
record = refined_mapping_record(result)
```

`aam` can be an `AAMResult` or a `FinalBranchCatalogue` with recorded fragment pairs. A flattened family list alone does not contain the required fragment atom sets and is rejected. `convert_fragments(aam)` exposes the same reusable conversion without event grouping. The experiment and shared execution adapter call this public implementation; there is no cluster-specific decoder.

`result.conversion.rematches` retains the compressed rematched graphs. `result.conversion.branches` links each recombined branch to its saved family and fragment paths. `result.catalogue` contains the converted mapping families. Event records use schema `mappa.decoded_events/v2` and preserve this representation and the explicit scope.

The event count is the number of bond differences touching mapped atoms, including mapped/unmapped boundaries. Bonds wholly between unmatched atoms are ignored. No atoms are invented. Identical literal mappings share classification work; branch partitions and provenance remain distinct. Event-class grouping does not enumerate symmetry bijections.

`event_slack=0` returns the minimum-count event classes. Every converted branch and its direct count remains in `conversion` and `family_reports`. Counts are computed before symmetry-class IDs, so high-event branches outside the requested window require no canonicalization. Requesting all event classes can cost substantially more than conversion and counting.

`max_events` and `event_slack` filter the converted branch results; they are not solver search bounds. The minimum is the lowest count among these branches. `complete` requires completed uncapped rematches, retained mapped atoms, and no counted events internal to a new fragment. `minimum_proven` uses this same converted-branch scope. Neither flag certifies exhaustive original-family decoding or global AAM search.

`result.symmetry(candidate)` returns compressed choices from the converted catalogue. `query_assignments(result, ...)` is a separate, explicitly requested constraint query over those choices. Only that query requires the optional solver dependency. Conversion does not call it or use it as a fallback.

Historical benchmark proofs apply only to their recorded source revisions; they do not certify converted-branch outputs.

### Matching implementation ownership

`mappa.match_fragment` and its public configuration/result types remain available.
Their implementation is now under `mappa.mappa_core.core.matching`, together with
fragment growth and compressed symmetry. Old internal `fragment`, `matcher` and
`growth` import paths are removed. The checkpoint reader explicitly relocates
saved class names to their new modules; it does not retain a second implementation.

### Search implementation ownership

Public `mappa.search_aam`, `search_aam_checkpoints`, direction planning, partial
mapping and subgraph APIs retain their signatures. Implementations now live under
`mappa.mappa_core.core.search`: `api`, `planning`, `session` and `graph`, with branch construction under
`branching/` and partial/subgraph objectives under `queries/`.
Old internal search module paths are removed. Local workers use
`mappa.mappa_core.core.scheduler.search`; the shared task queue retains its current API.

Trusted checkpoint and bundle readers explicitly relocate saved search class
names. Rebuild the native engine after updating: the unused experimental
`FragmentChoices` export has been removed. Shared seed search and optional
fragment competition remain supported.

## Interpret events without rematching

```python
from mappa import decode_events

events = decode_events(aam, threshold=0.5, metal_threshold=0.3)
pattern = events.evaluate(0)  # Original representative of saved family 0.
actions = events.family_actions[0]
```

The returned object retains the original catalogue and all its correlated
choices. Each action has a `response_invariant` certificate and per-generator
changed-pair counts. An uncertified action remains compressed; it is not rejected
or enumerated. This API does not list all patterns or assert a minimum.

To request a realization, provide one choice per original action to
`events.evaluate(family_id, choices)`. `None` selects the identity. For a pool,
list the images of its saved atoms; for a group, provide a word of signed 1-based
generator numbers (negative numbers apply the inverse). The original matching
constraints still apply. `events.counts(family_id, choices_iterable)` streams
counts in bounded vector batches for explicitly supplied choices only.

Internal and boundary contributions are available without expanding combinations:

```python
tables = events.tables(0)  # Lazy, shared response tables for one saved partition.
tables.internal          # Fragment ID -> original source-pair terms.
tables.boundary          # Fragment pair -> source-pair terms; -1 is unmatched.
tables.product_boundary_pairs  # Product bonds with dynamic unmatched occupancy.
result = events.contributions(0, choices=None, include_events=True)
result["counts"]         # breaking, forming, strengthening, weakening
result["internal"]       # Contributions by fragment.
result["boundary"]       # Contributions by pair of fragments / unmatched material.
result["events"]         # Original R/P endpoint pairs and weights.
rows = events.counts(0, requested_choices, categories=True)
```

Bond presence uses the same pair threshold as event counting: WBO > threshold.
An event additionally needs absolute change > threshold. Equality produces no
event. The four categories are exclusive; their sum equals the historical total.
Without `categories=True`, `counts` retains decrease/increase counts.

Terms retain the IDs of original actions affecting their endpoints; these are
not new independent choices. Atom domains and table entries can overapproximate
jointly reachable pairs. Only requested choices validated against the saved
family may be composed. Product-to-unmatched contributions are conditional on
occupied target atoms and are added once, not once per possible target image.
If a family has multiple saved partitions, select one with `partition=` in
`tables` / `contributions`; counts are independent of that bookkeeping choice.
Tables never select independent per-bond minima or claim a global optimum.

`convert_fragments` and `refine_mapping` retain their existing hierarchical
rematching behavior. They do not replace this lossless interpretation.

## Local event contributions and sweep cuts

`decode_events(aam).tables(family_id)` exposes internal and boundary bond
terms. Each term retains `target_pairs`, original `action_ids`, `endpoint_case`
(`fixed_fixed`, `fixed_shuffle`, or `shuffle_shuffle`) and `manual_cut` provenance.
Pair projection preserves correlated moves and saved action order without
expanding full permutations. A projected pair alone is not permission to combine
it with another term; use `contributions` or `counts` with a joint saved choice.

`tables.components` contains `(fragment_ids, action_ids)` contribution groups.
Shared choices, required-edge constraints and partial target occupancy keep
related fragments together. Constant responses are counted once. Original
sweep cuts are recorded even when no bond event occurs; fragment separation
alone never creates an event. All four event categories use the original R/P
bond weights and the existing ordinary/metal thresholds.

Preparation follows source/product bond neighborhoods, including product bonds
reachable through allowed choices. Manual cuts add their own endpoint pair;
unrelated source nonbonds are not visited. `tables.statistics['neighboring_pairs']`
reports the number of pairs examined. This preserves internal variations and
possible forming bonds without repeating an all-atom-pair scan per family.
