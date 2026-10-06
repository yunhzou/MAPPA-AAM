# Importable AAM API

Start with [AAM_SIMPLE.ipynb](AAM_SIMPLE.ipynb): a self-contained six-atom reaction with embedded coordinates, bond orders, and display helpers. It inspects the raw AAM, unique event candidates, certified shuffles, anchors, both search directions, conditional matching, serialization, and one growth animation per sweep. It uses one worker and can run outside the checkout after installation. Saved tables and static molecule drawings render on GitHub; interactive views require a trusted Jupyter notebook. The longer [TUTORIAL.ipynb](TUTORIAL.ipynb) covers R/P geometry and TS analysis.

Install `python -m pip install -e ".[notebook]"` for the notebook, or `.[postprocessing]` for symbolic event decoding without notebook packages. A C++17 compiler builds the bookkeeping extensions. The optional fast growth engine requires `python native/build_engine.py`; see [native installation](../native/README.md). Ordinary `search_aam` can use the Python reference growth implementation. `compete_fragments` and `execution="reused_native"` / `"shared_policies"` require that optional engine and fail explicitly if unavailable. No xTB is needed when coordinates and bond-order matrices are supplied.

## Default pipeline; separate post-processing

```python
from mappa import AAMProblem, MolecularEndpoint, AAMSearchConfig, search_aam
from mappa.postprocessing import EventDecodeConfig, decode_events

problem = AAMProblem(
    MolecularEndpoint(elements_R, xyz_R, wbo_R),
    MolecularEndpoint(elements_P, xyz_P, wbo_P))
aam = search_aam(problem, AAMSearchConfig(iso_tolerance=1.0), workers=1)
decoded = decode_events(aam,
                        EventDecodeConfig(threshold=0.5, metal_threshold=0.3))
candidates = decoded.candidates
witness = candidates[0].mapping if candidates else None
```

Default MAPPA uses one-seed cut sweep followed by separate event decoding. Competition is off by default; the standard search and decoding functions never invoke it.

Matrices must be symmetric and coordinates finite. The event decoder additionally requires balanced elements, complete bijections, nonnegative and exactly symmetric raw matrices. Search supports unbalanced/partial mappings; inspect `aam.graph.paths()` and stop reasons for those cases. Do not silently treat a partial mapping as a complete event candidate.

A candidate is a **signed bond-event class**, with one explicit witness, its event edges, and supporting family IDs. Event edges use source atom indices. This is different from a unique mapping modulo exact molecular symmetry (`PatternEquivalence` / `extract_path_patterns` in `mappa.pattern_collection`) and from a deduplicated final fragment branch (`aam.final_catalogue()`). None of these APIs enumerates every atom bijection.

`decoded.complete` only certifies all retained complete families in the requested window; it is false when a family times out or a partial path was skipped. It does not certify exhaustive search. `decoded.minimum_candidates` selects the lowest **observed** event count. With no candidates, inspect family reports and search stops. Budgets are soft; production callers needing a hard wall limit should isolate the call in a process. The default decoder has no class-count cap but a ten-second budget per family.

## Ask about symmetry after choosing a witness

```python
candidate = decoded.candidates[0]
evidence = decoded.symmetry(candidate)
answer = decoded.query(candidate, {source_atom: desired_target}, same_events=True)
if answer.status == "allowed":
    alternative_witness = answer.mapping
    certificate = answer.actions
```

This is a conditional existence query over the saved family union. Unspecified atoms may move together. Supply the entire modified mapping to test exactly one swap while fixing everything else. `same_events=True` preserves the selected signed-event class; `False` asks only about the broader saved matching relation. Answers are `allowed`, `forbidden` within these saved families, or `unknown`. A timeout is never a negative proof.

The evidence retains ordered target generators/pools, required fragment edges, matching policy, and archive/path provenance. Individual orbit memberships are not independent shuffle permissions. Intrinsic fragment automorphisms and event-preserving equivalences are different objects. Final event symmetry uses exact response labels at the selected thresholds, rather than a second approximate tolerance. Changing event thresholds never changes the saved search. Advanced callers can use `SignedEventIndex`, `extract_path_events`, and `query_path` directly.

## Inspect growth directly from a result

```python
from html import escape
from IPython.display import HTML, display
from mappa.viewers import aam_growth_html

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

`direction` accepts `forward`, `reverse`, `both`, `smaller_first`, and `larger_first`. Size counts explicit atoms, including H. On ties smaller-first uses input R→P; larger-first uses the reverse. The example decodes each direction separately and requires balanced endpoints for event decoding. Both directions run sequentially under the same worker ceiling. Each direction has its own checkpoint directory when `intermediate_dir` is supplied. Anchors are always supplied in input R→P indexing and reversed automatically. Anchors constrain search, not just the displayed witness. The `reference` and `reused_native` backends accept anchors; the optional `shared_policies` backend currently rejects anchored searches.

Keep compressed results in their search direction. Invert only concrete witnesses. Event IDs are direction-local: merging already reduced forward/reverse event classes is **not** an exact common-frame bidirectional decoder. That joint event quotient is not implemented; use the two directed results for audit. `plan_aam_search` remains the lightweight smaller-first planning API.

## Conditional fragment and fixed-query isomorphism

`match_fragment(source_graph, target_graph, seed=..., context=FragmentMatchContext(locked_mapping=...), config=FragmentMatchConfig(...))` grows a seeded saturated fragment subject to existing assignments. Its result includes compressed placements and an explicit cap flag. It does not promise a globally largest common subgraph.

`match_weighted_subgraph(query, target, anchor_map=..., iso_tol=..., symmetry_wbo_tol=..., graph_floor=..., max_branches=..., node_policy=...)` matches a fixed query, with complete validated query placements and unsuccessful/capped search evidence. These two operations serve different purposes. Both accept `WeightedGraph` or prepared graphs from `build_graph`.

## Adjustable parameters

| API/configuration | Controls |
|---|---|
| `AAMSearchConfig` | `iso_tolerance=1.0`; `graph_floor=0.2`; `cut_floor=0.2`; `sweep_cuts=True`; `seed_count=1`; `random_seed=42`; `seed_selection="random"` or `"distance"`; `branch_limit=100`; `anchors=()`; `task_chunksize=1` |
| `search_aam` / directional helper | `workers`; `execution="reference"`, `"reused_native"`, or `"shared_policies"`; `intermediate_dir`; `resume`; `archive_format="json"` or `"checkpoint"` |
| `EventDecodeConfig` | `threshold=0.5`; `metal_threshold=0.3` (`None` disables the special metal rule); `max_events=None`; `seconds_per_family=10`; `max_patterns_per_family=None` |
| `decoded.query` | `assignments`; `same_events=True`; `seconds=10` |
| `FragmentMatchConfig` | `graph_floor`; `iso_tolerance`; `minimum_size`; `branch_limit`; `node_policy`; `allow_mapped_seed`; `orbit_dedup` |
| `FragmentMatchContext` | `locked_mapping`; `islands`; `deferred_edges`; optional reusable `source_orbits`, `target_orbits`, `growth_replay` |
| `CompetitionConfig` | `operation_budget=128`; `seconds=270`; `parent_limit=8`; `depth_limit=2`; `queue_limit=512`; `dependent_component_limit=8` |

The library and this notebook use the branch-limit default of 100. A cap bounds a growth call's live alternatives, not the number of saved histories or all permutations. The default pipeline uses one seed and cut sweep; no-sweep is an explicit ablation (`sweep_cuts=False`). The random seed controls reproducible, independent per-cut streams and is separate from the number of seed orderings.

For compatibility, `AAMSearchConfig` also retains `event_threshold`, `metal_event_threshold`, `symmetry_repair`, `symmetry_repair_min_changes`, and `symmetry_repair_max_evaluations`. The raw search does not classify or repair events. These fields configure the older mechanism/geometry pipeline; competition also uses its event thresholds for parent ordering. The separate signed-event decoder uses **only its own `EventDecodeConfig`**, not these legacy event fields. `group_mechanisms`, `compile_mapping_families`, `compile_mechanism_families`, `select_rp_mappings`, and TS routines remain available and unchanged.

## Experimental fragment competition (optional; default off)

```python
from mappa.competition import CompetitionConfig, compete_fragments
augmented = compete_fragments(aam, CompetitionConfig(operation_budget=128))
decoded = decode_events(augmented.final_catalogue())
```

This exposes the bounded local B-priority takeover policy used in the coordinate campaign. It returns the untouched baseline plus independently validated repair families, counters, and proposal provenance. It performs no decoding or reference comparison. Growth/completion use the supplied AAM branch limit and tolerances; hard user anchors remain fixed. Parent selection uses the AAM configuration's observed event counts. Original scripts remain as historical campaign artifacts. Competition requires an explicit call and the native growth engine. It is excluded from the published method and its reported results.

**TODO (experimental):** establish a controlled coverage/cost benefit and improve proposal scheduling and deduplication before considering this extension for the default pipeline. The historical experiments are retained under `reports/`.

## Chirality and audit scope

**TODO: a general stereochemistry-preserving AAM search/decoder contract.** Element/WBO inputs do not encode CIP labels, E/Z geometry, formal charge, or isotope labels as first-class constraints. Do not claim raw AAM preserves them. The existing `select_rp_mappings` stage applies coordinate/index-chirality and fixed-mapping RMSD selection to exact analytical families. That downstream check does not supply a general stereochemical guarantee for the event notebook. Fixed-query matching can accept richer node policies, but the full AAM endpoint contract still needs that extension.

The September 2026 Python audit added importable direction orchestration, signed-event decoding with witness-linked symmetry queries, explicit random/cut controls, optional dependency groups, and the competition adapter. It leaves the computational stages separate. Focused tests verify anchors in both orientations, threshold independence, partial/timeout status, same-class versus broader-family queries, checkpoint identity, and baseline retention under competition. This is software validation using small examples, not a new benchmark campaign.
