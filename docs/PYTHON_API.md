# Importable AAM API

Start with [AAM_SIMPLE.ipynb](AAM_SIMPLE.ipynb): a self-contained six-atom reaction with embedded coordinates, bond orders, and display helpers. It inspects the raw AAM, unique event candidates, certified shuffles, anchors, both search directions, conditional matching, serialization, and one growth animation per sweep. It uses one worker and can run outside the checkout after installation. Saved tables and static molecule drawings render on GitHub; interactive views require a trusted Jupyter notebook. The longer [TUTORIAL.ipynb](TUTORIAL.ipynb) covers R/P geometry and TS analysis.

Install `python -m pip install -e ".[notebook]"` for the notebook, or `.[postprocessing]` for optional symbolic mapping queries without notebook packages. Event conversion itself needs no solver. A C++17 compiler builds the bookkeeping extensions. The optional fast growth engine requires `python src/mappa/mappa_core/native/build_engine.py`; see [native installation](../src/mappa/mappa_core/native/README.md). Ordinary `search_aam` can use the Python reference growth implementation. `compete_fragments` and `execution="reused_native"` / `"shared_policies"` require that optional engine and fail explicitly if unavailable. No xTB is needed when coordinates and bond-order matrices are supplied.

## Default pipeline; separate post-processing

```python
from mappa import AAMProblem, MolecularEndpoint, AAMSearchConfig, search_aam
from mappa.mappa_core.post_processing.decoder import EventDecodeConfig, decode_events

problem = AAMProblem(
    MolecularEndpoint(elements_R, xyz_R, wbo_R),
    MolecularEndpoint(elements_P, xyz_P, wbo_P))
aam = search_aam(problem, AAMSearchConfig(iso_tolerance=1.0), workers=1)
decoded = decode_events(aam,
                        EventDecodeConfig(threshold=0.5, metal_threshold=0.3, event_slack=0))
candidates = decoded.candidates
witness = candidates[0].mapping if candidates else None
```

Default MAPPA uses one-seed cut sweep followed by separate event decoding. Competition is off by default; the standard search and decoding functions never invoke it.

Matrices must be nonnegative and symmetric, and coordinates finite. Both search and event conversion support partial mappings and unequal endpoint sizes. Counts include differences touching mapped atoms and their boundaries, while wholly unmatched bonds are ignored. A completed conversion does not imply all endpoint atoms were mapped.

A candidate is a **signed bond-event class**, with one explicit witness, its event edges, and supporting family IDs. Use `endpoint_events` for bond lists in the original reactant/product indices. A deduplicated final fragment branch (`aam.final_catalogue()`) retains compressed mapping choices. Neither API enumerates every atom bijection.

`decoded.complete` covers conversion of the retained fragment contexts, including preservation of mapped atoms and absence of internal counted events. A capped rematch remains incomplete. `decoded.minimum_candidates` selects the lowest direct branch counts. This is not an exhaustive old-family or global-search claim. The decoder uses the existing fragment matcher and no solver enumeration. `decoded.conversion` retains the rematched compressed graphs and recombination provenance. Hard wall/memory limits belong to the caller's process or shared executor.

## Ask about symmetry after choosing a witness

```python
candidate = decoded.candidates[0]
evidence = decoded.symmetry(candidate)
from mappa.mappa_workflow.mapping_analytics.queries import query_assignments
answer = query_assignments(decoded, candidate, {source_atom: desired_target}, same_events=True)
if answer.status == "allowed":
    alternative_witness = answer.mapping
    certificate = answer.actions
```

This is a conditional existence query over the converted family union. Unspecified atoms may move together. Supply the entire modified mapping to test exactly one swap while fixing everything else. `same_events=True` preserves the selected signed-event class; `False` asks only about the converted matching relation. Answers are `allowed`, `forbidden` within these saved families, or `unknown`. A timeout is never a negative proof.

The evidence retains ordered target generators/pools, required fragment edges, matching policy, and archive/path provenance. Individual orbit memberships are not independent shuffle permissions. Intrinsic fragment automorphisms and event-preserving equivalences are different objects. Final event symmetry uses exact response labels at the selected thresholds, rather than a second approximate tolerance. Changing event thresholds never changes the saved search. Advanced callers can use `convert_fragments` to retain converted branches before event grouping, and `query_path` for explicitly requested mapping queries.

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
    for candidate in directed_decoding.candidates:
        input_mapping = run.to_input_mapping(candidate.mapping)
```

`direction` accepts `forward`, `reverse`, `both`, `smaller_first`, and `larger_first`. Size counts explicit atoms, including H. On ties smaller-first uses input R→P; larger-first uses the reverse. The example converts each direction separately; partial mappings use the same bond-difference rule. Both directions run sequentially under the same worker ceiling. Each direction has its own checkpoint directory when `intermediate_dir` is supplied. Anchors are always supplied in input R→P indexing and reversed automatically. Anchors constrain search, not just the displayed witness. The `reference` and `reused_native` backends accept anchors; the optional `shared_policies` backend currently rejects anchored searches.

Keep compressed results in their search direction. Invert only concrete witnesses. Event IDs are direction-local: merging already reduced forward/reverse event classes is **not** an exact common-frame bidirectional decoder. That joint event quotient is not implemented; use the two directed results for audit. `plan_aam_search` remains the lightweight smaller-first planning API.

## Conditional fragment and fixed-query isomorphism

`match_fragment(source_graph, target_graph, seed=..., context=FragmentMatchContext(locked_mapping=...), config=FragmentMatchConfig(...))` grows a seeded saturated fragment subject to existing assignments. Its result includes compressed placements and an explicit cap flag. It does not promise a globally largest common subgraph.

`match_weighted_subgraph(query, target, anchor_map=..., iso_tol=..., symmetry_wbo_tol=..., graph_floor=..., max_branches=..., node_policy=...)` matches a fixed query, with complete validated query placements and unsuccessful/capped search evidence. These two operations serve different purposes. Both accept `WeightedGraph` or prepared graphs from `build_graph`.

## Adjustable parameters

| API/configuration | Controls |
|---|---|
| `AAMSearchConfig` | `iso_tolerance=1.0`; `graph_floor=0.2`; `cut_floor=0.2`; `sweep_cuts=True`; `seed_count=1`; `random_seed=42`; `seed_selection="random"` or `"distance"`; `branch_limit=100`; `anchors=()`; `task_chunksize=1` |
| `search_aam` / directional helper | `workers`; `execution="reference"`, `"reused_native"`, or `"shared_policies"`; `intermediate_dir`; `resume`; `archive_format="json"` or `"checkpoint"` |
| `EventDecodeConfig` | `threshold=0.5`; `metal_threshold=0.3` (`None` disables the special metal rule); `max_events=None`; `event_slack=None`; `branch_limit=100`; `random_seed=42` |
| `query_assignments` | `assignments`; `same_events=True`; `seconds=10` |
| `FragmentMatchConfig` | `graph_floor`; `iso_tolerance`; `minimum_size`; `branch_limit`; `node_policy`; `allow_mapped_seed`; `orbit_dedup` |
| `FragmentMatchContext` | `locked_mapping`; `islands`; `deferred_edges`; optional reusable `source_orbits`, `target_orbits`, `growth_replay` |
| `CompetitionConfig` | `operation_budget=128`; `seconds=270`; `parent_limit=8`; `depth_limit=2`; `queue_limit=512`; `dependent_component_limit=8` |

The library and this notebook use the branch-limit default of 100. A cap bounds a growth call's live alternatives, not the number of saved histories or all permutations. The default pipeline uses one seed and cut sweep; no-sweep is an explicit ablation (`sweep_cuts=False`). The random seed controls reproducible, independent per-cut streams and is separate from the number of seed orderings.

Raw search does not classify or repair events. `decode_events` owns event conversion and counting through `EventDecodeConfig`. `AAMSearchConfig.event_threshold` and `metal_event_threshold` remain inputs to competition ordering; they do not configure decoding. TS analysis requires a `ReactionContext` supplied by the caller. The old mechanism-grouping API, its compiler and symmetry-repair options are removed. Full-result checkpoints containing those retired options must be read with their producing revision; current readers reject them rather than discard settings. Saved raw/finalized graph-only checkpoints remain readable.

## Experimental fragment competition (optional; default off)

```python
from mappa.mappa_core.core.search.branching.competition import CompetitionConfig, compete_fragments
augmented = compete_fragments(aam, CompetitionConfig(operation_budget=128))
decoded = decode_events(augmented.final_catalogue())
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
`search_aam` and `decode_events` explicitly. Coordinate orientation is not
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
from mappa.mappa_core.post_processing.decoder import EventDecodeConfig, decode_events
from mappa.mappa_core.utils.io.records import decoded_events_record

result = decode_events(aam, EventDecodeConfig(event_slack=0))
print(result.observed_minimum, result.complete)
for candidate in result.candidates:
    print(candidate.total, candidate.mapping, candidate.endpoint_events)
record = decoded_events_record(result)
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
