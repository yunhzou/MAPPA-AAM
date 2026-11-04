# Importable AAM API

Start with [AAM_SIMPLE.ipynb](AAM_SIMPLE.ipynb): a self-contained six-atom reaction with embedded coordinates, bond orders, and display helpers. It inspects the raw AAM, unique event candidates, certified shuffles, anchors, both search directions, conditional matching, serialization, and one growth animation per sweep. It uses one worker and can run outside the checkout after installation. Saved tables and static molecule drawings render on GitHub; interactive views require a trusted Jupyter notebook. The longer [TUTORIAL.ipynb](TUTORIAL.ipynb) covers R/P geometry and TS analysis.

Install `python -m pip install -e ".[notebook]"` for the notebook, or `.[postprocessing]` for local event spectra and symbolic mapping queries without notebook packages. A C++17 compiler builds the bookkeeping extensions. Normal search requires the fast growth engine: run `python src/mappa/mappa_core/native/build_engine.py`; see [native installation](../src/mappa/mappa_core/native/README.md). Missing native execution is an error. Python comparison runs explicitly select `MAPPA_NATIVE=0` and `search_aam(..., execution="reference")`. No xTB is needed when coordinates and bond-order matrices are supplied.

## Default pipeline

```python
from mappa import AAMProblem, MolecularEndpoint, AAMSearchConfig, search_aam, deduplicate, decode_events
from mappa.mappa_core.utils.io.records import decoded_events_record

problem = AAMProblem(MolecularEndpoint(elements_R, xyz_R, wbo_R),
                     MolecularEndpoint(elements_P, xyz_P, wbo_P))
aam = search_aam(problem, AAMSearchConfig(iso_tolerance=1.0))
catalogue = deduplicate(aam)
decoded = decode_events(catalogue, threshold=0.5, metal_threshold=0.3)
family_id = decoded.retained_family_ids[0]
partition = decoded.window.retained_partitions[family_id][0]
bins = decoded.attainable_counts(family_id, partition=partition)
record = decoded_events_record(decoded)
```

Decoding uses the `independent_fragments_with_sweep` approximation for both
selection and reported counts. Ordinary boundaries contribute fixed breaking and
forming events. Local fragment event spectra combine independently; the saved
sweep boundary retains its coupled correction. These are approximate model counts,
not a physical minimum guarantee. The original mapping catalogue and joint
constraints remain available for reference verification and explicit realizations.

The default window is the model minimum plus six. Change the allowance or keep
all model bins with an explicit `None`:

```python
decoded = decode_events(catalogue, event_slack=2)
print(decoded.minimum_events, decoded.event_window)
print(decoded.retained_family_ids)
unfiltered = decode_events(catalogue, event_slack=None)
```

`attainable_counts(family_id, partition=...)` returns the retained bins for one
partition/cut occurrence in the model. A numeric bin need not have a complete
joint mapping witness: independent composition may combine incompatible local
or placement choices. Counts between the model minimum and maximum need not all
occur. `realize`, `contributions`, `evaluate` and `counts` validate an explicitly
requested joint mapping and reject it when its model count exceeds the window.
A retained occurrence's default witness can be outside that window.

Records use `mappa.decoded_events/v10`, preserve the full catalogue, and label
`approximate=True`, `bin_witnesses_guaranteed=False`, `minimum_proven=False` and the separate
`model_minimum_proven` flag. The CLI defaults to `--event-slack 6`; pools accept
`EventDecodeConfig(event_slack=2)`. Supply the policy through either a config or
the keyword, not both. Historical pickled decoder configs without `event_slack`
require their producing revision.

For explicit matching/deduplication/decoding passes over many reactions, see the
[persistent local and Slurm pool](shared-execution.md). Its RAM and bundle storage
share the same scientific operations. Direct calls above need no pool.

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
    directed_decoding = decode_events(run.aam, event_slack=None)
    for fid in range(len(directed_decoding.catalogue.families)):
        input_mapping = run.to_input_mapping(directed_decoding.realize(fid))
```

`direction` accepts `forward`, `reverse`, `both`, `smaller_first`, and `larger_first`. Size counts explicit atoms, including H. On ties smaller-first uses input R→P; larger-first uses the reverse. The example converts each direction separately; partial mappings use the same bond-difference rule. Both directions run sequentially under the same worker ceiling. Each direction has its own checkpoint directory when `intermediate_dir` is supplied. Anchors are always supplied in input R→P indexing and reversed automatically. Anchors constrain search, not just the displayed witness. The `reference` and `reused_native` backends accept anchors; the optional `shared_policies` backend currently rejects anchored searches.

Keep compressed results in their search direction. Invert only concrete witnesses.
Each direction has its own approximate model minimum and bins; matching counts
do not establish that two directions describe the same event pattern.
`plan_aam_search` remains the lightweight smaller-first planning API.

## Conditional fragment and fixed-query relaxed matching

`match_fragment(source_graph, target_graph, seed=..., context=FragmentMatchContext(locked_mapping=...), config=FragmentMatchConfig(...))` grows a seeded saturated fragment subject to existing assignments. Each extension constrains only the selected growth edge. Returned `Fragment` records expose those edges as `accepted_edges`; other bonds between included atoms remain available in the original molecular graphs for event interpretation. Its result includes compressed placements and an explicit cap flag. It does not promise a globally largest common subgraph.

```python
from mappa import FragmentMatchConfig, match_fragment

result = match_fragment(source_graph, target_graph, seed=0,
    config=FragmentMatchConfig())
```

Fragment growth finds one valid placement per extension, allowing assignments
inside the growing fragment to change while earlier fragment locks stay fixed.
It first tries extending the witness, then searches the accepted fragment if
reassignment is needed. At saturation, target atom sets related by a lossless
global placement group share one canonical `Fragment` and continuation. The
group preserves the full target WBO matrix and earlier locks; its generators
are saved in `result.placement_generators`. Other target sets remain separate.
Each fragment retains all discovery assignments as joint constraints, plus its
anchor-aware source and target relaxed self relationships. The placement group
moves the complete later continuation together and does not replace either
fragment-local relationship. The branch cap counts canonical placement
continuations at every hierarchy level, not represented atom sets or internal
permutations. Python and native
execution use this same algorithm. There is no growth mode to select.
Completed-fragment results have bounded reuse; intermediate growth candidates
are not cached.

```python
from mappa import AAMSearchConfig, search_aam

config = AAMSearchConfig()
aam = search_aam(problem, config)  # reused_native execution
```

Pass the same config to `pool.search(cases, config)` in a work pool. Each prepared
session uses a 64 MiB fragment-result retention budget, evicts old entries, and
returns oversized results without caching them. This is a cache budget, not a
bound on total search memory. The `reference` execution backend runs the same
growth algorithm without fragment-result caching. Native fragment requests
reject unsupported trace events and custom node policies. Select `MAPPA_NATIVE=0` explicitly when these Python capabilities are needed.

`match_weighted_subgraph(query, target, anchor_map=..., iso_tol=..., graph_floor=..., max_branches=..., node_policy=...)` matches a fixed query, with complete validated query placements and unsuccessful/capped search evidence. These two operations serve different purposes. Both accept `WeightedGraph` or prepared graphs from `build_graph`.

## Adjustable parameters

| API/configuration | Controls |
|---|---|
| `AAMSearchConfig` | `iso_tolerance=1.0`; `self_tolerance=1.0`; `graph_floor=0.2`; `cut_floor=0.2`; `sweep_cuts=True`; `seed_count=1`; `random_seed=42`; `seed_selection="random"` or `"distance"`; `branch_limit=100`; `anchors=()`; `task_chunksize=1` |
| `search_aam` / directional helper | `workers`; `execution="reused_native"` (default), `"reference"`, or `"shared_policies"`; `intermediate_dir`; `resume`; `archive_format="json"` or `"checkpoint"` |
| `EventDecodeConfig` | `threshold=0.5`; `metal_threshold=0.3`; `event_slack=6` (`None` retains all approximate model bins) |
| `RefinementConfig` | `threshold=0.5`; `metal_threshold=0.3` (`None` disables the special metal rule); `max_events=None`; `event_slack=None`; `branch_limit=100`; `random_seed=42` |
| `query_assignments` | `assignments`; `same_events=True`; `seconds=10` |
| `FragmentMatchConfig` | `graph_floor`; `iso_tolerance`; `self_tolerance=1.0`; `minimum_size`; `branch_limit`; `node_policy`; `allow_mapped_seed`; `random_seed`; `bond_tolerance` |
| `FragmentMatchContext` | `locked_mapping`; `islands`; `deferred_edges`; optional reusable `fragment_cache`; `placement_fixed_targets=()` |
| `CompetitionConfig` | `operation_budget=128`; `seconds=270`; `parent_limit=8`; `depth_limit=2`; `queue_limit=512`; `dependent_component_limit=8` |

Discovery tolerance constrains matching between endpoints. `self_tolerance`
independently constrains rearrangements within each original fragment; it is not
the event threshold.

`placement_fixed_targets` is an expert context constraint for callers whose
later operation distinguishes particular target labels. Those labels are fixed
only by global placement compression: they remain available for matching and
are not occupied or anchored. Fragment competition fixes every target label
because its takeover plan distinguishes the owners displaced by each offer.

The library and this guide use the branch-limit default of 100. The cap bounds canonical placement continuations and combined live branches at each hierarchy level; it does not count represented symmetry alternatives. After each completed seed step, pruning proceeds bottom-up: remove the largest child subtree until each node fits, and use its reduced count when checking its parent. Equal-size ties preserve earlier decisions. The completed-fragment limit still rejects an oversized placement result. The default pipeline uses one seed and cut sweep; no-sweep is an explicit ablation (`sweep_cuts=False`). The random seed controls reproducible, independent per-cut streams and is separate from the number of seed orderings.

The returned search graph retains the complete fragment paths leading to its
results. Discarded paths release their payloads. Their cap evidence remains in
`graph.stops`: a cap summary has `state=None`, an explicit `context`, and an
`occurrences` count. Its `count` and `limit` describe the pruning decision;
for combined hierarchy pruning, each occurrence removes one child subtree.
Use `graph.stop_context(stop)` for either a summary or a stop attached
to a retained state. Historical checkpoints keep their original attached stops.

Raw search does not classify or repair events. `decode_events` interprets original mappings; `refine_mapping` is separate experimental hierarchical matching. `AAMSearchConfig.event_threshold` and `metal_event_threshold` remain inputs to competition ordering; they do not configure decoding. TS analysis requires a `ReactionContext` supplied by the caller. The old mechanism-grouping API, its compiler and symmetry-repair options are removed. Full-result checkpoints containing those retired options must be read with their producing revision; current readers reject them rather than discard settings. Earlier fragment schemas require their producing revision; current readers reject them instead of reinterpreting their allowed mappings.

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

## Experimental event refinement

```python
from mappa.mappa_core.post_processing.hierarchical_mapping import RefinementConfig, refine_mapping
from mappa.mappa_core.utils.io.records import refined_mapping_record

result = refine_mapping(aam, RefinementConfig(event_slack=0))
print(result.observed_minimum, result.complete)
for candidate in result.candidates:
    print(candidate.total, candidate.mapping, candidate.endpoint_events)
record = refined_mapping_record(result)
```

`aam` can be an `AAMResult` or a `FinalBranchCatalogue` with recorded fragment pairs. A flattened family list alone does not contain the required fragment atom sets and is rejected. `convert_fragments(aam)` exposes the same reusable conversion without event grouping. This is the experimental refinement API. The work pool uses the original-space decoder described above.

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

`search_core_assignments` projects the requested atoms through the saved ordered
fragment discovery constraints and both-sided relaxed self relationships. It uses the optional symbolic-query
dependency from `.[postprocessing]`, enumerates only distinct requested tuples,
and enforces `assignment_limit` without expanding complete atom mappings.

Trusted checkpoint and bundle readers explicitly relocate saved search class
names. Rebuild the native engine after updating: the unused experimental
`FragmentChoices` export has been removed. Shared seed search and optional
fragment competition remain supported.

## Interpret events without rematching

```python
from mappa import decode_events

events = decode_events(aam, event_slack=None)  # All approximate model bins.
counts = events.contributions(0)['counts']
actions = events.family_actions[0]
```

The original compressed mapping relation remains attached. Each action has a
`response_invariant` certificate and per-generator changed-pair counts. An
uncertified action remains compressed; it is not enumerated.

To request a realization, provide one choice per ordered saved action to
`events.realize(family_id, choices, partition=0)`. For a fragment relation, `None`
keeps its discovery witness with identity self rearrangements. A dictionary can
supply `discovery` (original source-to-target indices), `source_word` and
`target_word` (signed 1-based generator numbers), or `source_permutation` and
`target_permutation` in the respective fragment-local index frames.

For a global placement action, supply a word of signed 1-based generator
numbers; `None` or `()` selects identity. The family program applies each local
relation before its placement action, with later fragment decisions applied
first. This moves later assignments and their shared anchors together. A saved
path exposes the same choices by transition ID:

```python
realized = path.realize(
    fragment_choices={transition_id: fragment_choice},
    placement_choices={transition_id: (1,)},
)
```

Discovery and both self rearrangements are checked against their original
constraints. The resulting two-sided mapping may exceed discovery tolerance;
event responses still use original WBOs. `events.counts(family_id, choices)`
streams approximate counts for explicit valid choices in bounded batches.
Without `categories=True`, it returns decrease/increase counts; with it, the
columns are breaking, forming, strengthening and weakening.

```python
tables = events.tables(0, partition=0)
tables.internal                 # Fragment ID -> internal source-pair terms.
tables.reactant_boundary_pairs  # Fixed breaking contributions.
tables.product_boundary_pairs   # Fixed forming contributions.
tables.sweep_cut                # One coupled boundary cut, or None.
result = events.contributions(0, choices=None, partition=0, include_events=True)
rows = events.counts(0, requested_choices, partition=0, categories=True)
report = events.evaluate(0, choices=None, partition=0)
```

`evaluate` reports the mapping and approximate contributions with an explicit
model label; it does not produce an exact canonical event ID. Internal terms
retain original action IDs and correlated endpoint choices. Ordinary boundary
pairs use the saved representative partition. Before composing a sweep, remove
its counted source contribution and reserve one counted target-boundary event,
when available. The coupled term restores any unused reservation and adds the
actual original-WBO event: a contribution from zero to two. An uncounted boundary
event is never subtracted, including when matching uses transformed weights.
A cut internal to one fragment is counted there and receives no boundary
correction.

Each `partition=` index identifies one partition/cut occurrence. Different
histories keep their own cuts, even when they share a family or partition. Their
approximate counts can differ; do not combine their totals or union their cuts.
`events.inventory.cuts[family_id][partition]` exposes this context.

Bond presence uses WBO greater than the pair threshold. An event additionally
requires an absolute WBO change greater than that threshold; equality is
unchanged. Categories are exclusive. Partial mappings include boundaries touching
mapped atoms and omit wholly unmatched bonds.

If matching used transformed weights, pass
`decode_events(catalogue, event_problem=original_problem)`. Elements, atom order
and coordinates must agree. Matching constraints use `catalogue.problem`; event
responses use `events.event_problem`. `EventDecodeSession` accepts the same
keyword. The v10 record always identifies both endpoint hashes and adds
`event_problem` when its weights differ from the matching `problem`.

The embedded `mappa.final_fragment_branches/v3` catalogue stores shared fragment
and placement records with ordered action references. Saved v10 tables retain
internal responses, fixed boundary pairs, sweep context and occurrence indices.
Derived per-atom lookup maps remain runtime data. Binary export preserves sharing
without retaining every temporary family table.

## Local event preparation

The decoder first collects every unique saved fragment across the completed
reaction/direction catalogue. `events.inventory.fragments` contains those
scientific records, while `events.inventory.uses[family_id]` contains the partition
piece references for each occurrence. Preparation is reused when the local
mapping, relevant ordered actions and constraints agree.

Local spectra retain correlations within each fragment. Their independent
composition is the declared approximation; only the sweep endpoints are coupled
across pieces. `tables.components` describes the local contribution groups.
Preparation follows relevant source/product bond neighborhoods, avoiding a
repeated all-atom-pair scan per family. The runtime table statistics expose the
work performed.

`convert_fragments` and `refine_mapping` remain separate experimental rematching
operations. They do not replace ordinary decoding.
