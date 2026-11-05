# Importable AAM API

Start with [AAM_SIMPLE.ipynb](AAM_SIMPLE.ipynb): a self-contained six-atom reaction with embedded coordinates, bond orders, and display helpers. It inspects the raw AAM, unique event candidates, certified shuffles, anchors, both search directions, conditional matching, serialization, and one growth animation per sweep. It uses one worker and can run outside the checkout after installation. Saved tables and static molecule drawings render on GitHub; interactive views require a trusted Jupyter notebook. The longer [TUTORIAL.ipynb](TUTORIAL.ipynb) covers R/P geometry and TS analysis.

Install `python -m pip install -e ".[notebook]"` for the notebook, or `.` for the mapping and decoding API without notebook packages. A C++17 compiler builds the bookkeeping extensions. Normal search requires the fast growth engine: run `python src/mappa/mappa_core/native/build_engine.py`; see [native installation](../src/mappa/mappa_core/native/README.md). Missing native execution is an error. Python comparison runs explicitly select `MAPPA_NATIVE=0` and `search_aam(..., execution="reference")`. No xTB is needed when coordinates and bond-order matrices are supplied.

## Default pipeline

```python
from mappa import AAMProblem, MolecularEndpoint, AAMSearchConfig, search_aam, deduplicate, decode_events
from mappa.mappa_core.utils.io.records import decoded_events_record

problem = AAMProblem(MolecularEndpoint(elements_R, xyz_R, wbo_R),
                     MolecularEndpoint(elements_P, xyz_P, wbo_P))
aam = search_aam(problem, AAMSearchConfig(iso_tolerance=1.0))
catalogue = deduplicate(aam)
decoded = decode_events(catalogue, threshold=0.5, metal_threshold=0.3, event_slack=6)
print(decoded.minimum_events, decoded.event_window)
witness = decoded.realize(0)
model_counts = decoded.contributions(0)['counts']
physical_events = decoded.evaluate(0)
record = decoded_events_record(decoded)
```

Reduction shares each final reactant/product atom-set pair once within the fixed
reaction, policy and explicit anchors. Growth routes and temporary locks remain
in the search archive. Branches retain their compositions, ordered global
placements and sweep provenance. Reduction performs no symmetry or event search.

Decoding inventories those unique completed fragments, improves one starting mapping
within `fragment_seconds` (default five seconds), and shares its binary-symmetry
event response across branch uses. The original matching catalogue remains
intact. The decoded response covers the selected start's symmetry orbit; it need
not cover every other discovery on the same atom sets.

`minimum_events` is an attained model count. Check `minimum_proven` and
`minimum_lower_bound` before treating it as exact. The default window ends at
`minimum_events + 6`. With attained compatible local minima, the initial sweep
window is at most the smallest composition lower bound plus two. Otherwise the
attained upper bound is used, so unfinished local work cannot discard a better
possible count.

Ordinary boundaries use the declared independence approximation; the single
sweep correction preserves its endpoint/internal-event dependence. `evaluate`
reports physical original-WBO events for an explicitly requested mapping.
Physical counts may differ from the model. `realize` validates explicit proofs
against the complete catalogue, so inspection can also cover discoveries outside
the selected decoded orbit. Neither operation asserts event-selection membership.

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

## Bounded fragment witness improvement

This core operation accepts an actual saved `Fragment`. It improves one valid
internal bond-event witness while preserving the complete saved relation. It is
used by `decode_events` once per unique saved fragment.

```python
from mappa.mappa_core.core.symmetry.optimization import optimize_fragment

answer = optimize_fragment(fragment, time_limit=5.0)
print(answer.witness, answer.event_count, answer.lower_bound, answer.complete)
print(answer.elapsed_seconds, answer.reason)
```

The deadline covers preparation, starting-map improvement and shared binary
symmetry optimization. Easy fragments stop early. `complete=True` certifies this
fragment's minimum; otherwise `event_count` is the actual count of a best-found
mapping, with `lower_bound` reported separately. The returned discovery witness
and source/target operations retain its membership proof. Original event weights
can be supplied through `event_problem`; the ordinary and metal thresholds
default to 0.5 and 0.3. This does not compute boundary events or all event bins.

Deadline checks are cooperative. An opaque native call can overrun until it
returns; a strict wall-clock cutoff requires the execution layer's watchdog.

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
| `AAMSearchConfig` | `iso_tolerance=1.0`; `graph_floor=0.2`; `cut_floor=0.2`; `sweep_cuts=True`; `seed_count=1`; `random_seed=42`; `seed_selection="random"` or `"distance"`; `branch_limit=100`; `anchors=()`; `task_chunksize=1` |
| `search_aam` / directional helper | `workers`; `execution="reused_native"` (default), `"reference"`, or `"shared_policies"`; `intermediate_dir`; `resume`; `archive_format="json"` or `"checkpoint"` |
| `decode_events` / `EventDecodeConfig` | `threshold=0.5`; `metal_threshold=0.3`; `event_slack=6` above the best attained model count; `fragment_seconds=5.0` for starting-map improvement |
| `RefinementConfig` | `threshold=0.5`; `metal_threshold=0.3` (`None` disables the special metal rule); `max_events=None`; `event_slack=None`; `branch_limit=100`; `random_seed=42` |
| `query_assignments` | `assignments`; `same_events=True`; `seconds=10` |
| `FragmentMatchConfig` | `graph_floor`; `iso_tolerance`; `minimum_size`; `branch_limit`; `node_policy`; `allow_mapped_seed`; `random_seed`; `bond_tolerance` |
| `FragmentMatchContext` | `locked_mapping`; `islands`; `deferred_edges`; optional reusable `fragment_cache`; `placement_fixed_targets=()` |
| `CompetitionConfig` | `operation_budget=128`; `seconds=270`; `parent_limit=8`; `depth_limit=2`; `queue_limit=512`; `dependent_component_limit=8` |

Discovery tolerance constrains matching between endpoints. Fragment symmetry
preserves binary connectivity, elements and anchors. Original WBO differences
are evaluated as events; there is no separate self-tolerance filter.

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

`result.symmetry(candidate)` returns compressed choices from the converted catalogue. `query_assignments(result, ...)` is a separate, explicitly requested constraint query over those choices. Conversion does not call this query. The solver dependency is installed with MAPPA.

Historical benchmark proofs apply only to their recorded source revisions; they do not certify converted-branch outputs.

### Matching implementation ownership

`mappa.match_fragment` and its public configuration/result types remain available.
Their implementation is under `mappa.mappa_core.core.matching`; shared symmetry
operations belong to `mappa.mappa_core.core.symmetry`. The checkpoint reader
explicitly relocates unchanged saved class names; it rejects incompatible
scientific schemas rather than retaining a second implementation.

### Search implementation ownership

Public `mappa.search_aam`, `search_aam_checkpoints`, direction planning, partial
mapping and subgraph APIs retain their signatures. Implementations live under
`mappa.mappa_core.core.search`: `api`, `planning`, `session` and `graph`, with branch construction under
`branching/` and partial/subgraph objectives under `queries/`. The same `queries/`
owner compiles saved paths for `query_path` and assignment projection; workflow
analytics consumes it without reversing the dependency into core. Local workers use
`mappa.mappa_core.core.scheduler.search`; the shared task queue retains its current API.

`search_core_assignments` projects the requested atoms through the saved ordered
fragment discovery constraints and both sides' binary symmetries. It uses the
symbolic-query dependency installed with MAPPA, enumerates only distinct requested tuples,
and enforces `assignment_limit` without expanding complete atom mappings.

Trusted checkpoint and bundle readers explicitly relocate saved search class
names. Shared seed search and optional fragment competition remain supported.

## Query shared fragment and branch responses

```python
from mappa import decode_events

events = decode_events(aam, threshold=0.5, metal_threshold=0.3, event_slack=6)
events.minimum_events                  # Minimum of the approximate event model
events.event_window                    # Inclusive upper count: minimum + 6
counts = events.contributions(0)['counts']
physical = events.evaluate(0)           # Separate physical check of this mapping
```

The decoder keeps the original catalogue and its mapping constraints. It shares
each unique fragment's conditional responses across compositions, adds fixed
ordinary boundary counts, and resolves the residual for the single sweep cut.
That residual contributes zero, one or two events and is counted once. A cut or
fragment separation alone does not imply an event. The model treats ordinary
fragment boundaries independently; it does not solve correlations across the
whole branch to recover a physical minimum.

Use `events.select(family_id, lower, upper)` to retain a count interval without
expanding mappings. Omit `upper` for an exact count, or omit both bounds for the
computed window. A targeted lookup returns `found`, `absent` or `unknown` under
its deadline; an unknown answer is never absence.

```python
selection = events.select(0, 2)                 # Exact model event count
answer = selection.find(seconds=5.0)           # One joint witness, if found
if answer.status == "found":
    mapping = dict(answer.witness)
    assert events.contributions(0, answer.choices)["total"] == 2

fragment = events.fragments[0]                 # Shared actual Fragment response
local = fragment.response.select(0, 2)
answer = local.find(assignments={source_atom: target_atom}, seconds=5.0)
```

The stored selection retains every allowed choice satisfying its predicate;
`find()` returns one witness rather than replacing the relation by that witness.
Minimization and target queries reuse the same core symmetry operations. No
complete attainable-count or bijection table is constructed.

Numeric bounds and local minima are computed lazily and reused. Cheap boundary
bounds exclude compositions before local preparation when possible. The decoder
does not eagerly build every count bin, enumerate all symmetry bijections or
repeat fragment preparation per branch. Alternative partitions are alternative
views, so their contributions must not be added together.

`contributions` and `counts` report approximate model counts. `evaluate` is an
explicit original-WBO check of a valid concrete mapping, outside model-window
construction. Its physical count can differ; it does not re-bin the model result
or alter `minimum_events`. No global physical minimum or exhaustive matching
coverage is claimed.

### Requested mappings and physical checks

To request a realization, provide one choice per ordered saved action to
`events.realize(family_id, choices)`. For a fragment relation, `None` keeps its
discovery witness with identity self rearrangements. A choice dictionary can supply `discovery`
(original source-to-target indices), `source_word` and `target_word` (signed
1-based generator numbers). Alternatively, supply `source_permutation` and
`target_permutation` in the respective fragment-local index frames.

For a global placement action, supply a word of signed 1-based generator
numbers; `None` or `()` selects identity. The completed family applies its
disjoint local fragment choices, then the saved joint placement sequence.
Raw search paths retain chronological growth decisions and expose their original
choices by transition ID:

```python
realized = path.realize(
    fragment_choices={transition_id: fragment_choice},
    placement_choices={transition_id: (1,)},
)
```

The decoder validates element compatibility, explicit anchors and each self rearrangement
against its own original fragment. It applies target rearrangement composed with
discovery composed with inverse source rearrangement, once. The resulting mapping
may exceed discovery tolerance; explicit physical evaluation uses original WBOs.
The event-count approximation does not relax these realization checks.

```python
mapping = events.realize(0, choices)
model = events.contributions(0, choices)
rows = events.counts(0, requested_choices, categories=True)
physical = events.evaluate(0, choices)
```

Bond presence uses the same pair threshold as event counting: WBO > threshold.
An event additionally needs absolute change > threshold. Equality produces no
event. Breaking, forming, strengthening and weakening are exclusive categories.
Without `categories=True`, `counts` retains decrease/increase counts.

### Shared responses and persistence

`events.inventory.fragments` contains the unique saved Fragment objects.
`events.inventory.uses[family_id]` identifies the disjoint pieces owning internal
bonds in each final partition. Inventory construction precedes response preparation.
Equal atom-set pairs share a local response under the same endpoint problem,
policy and explicit anchors. Global placement order and sweep context remain
branch-level information; they do not create more local fragments.

`events.tables(family_id)` exposes the conditional fragment and boundary data.
Responses preserve local correlations and action order. The branch composition
model uses fixed ordinary boundaries and the separate sweep correction;
ordinary boundaries are not reevaluated for every combination of local choices.
Mapped/unmapped boundaries follow this same convention. Physical evaluation
counts each original bond touching mapped atoms once and omits bonds wholly
inside unmatched material.

`decoded_events_record(events)` uses `mappa.decoded_events/v12` and retains the original catalogue, shared conditional
responses, model minimum/window and upper-window predicate. It does not export
an eager list of all count bins or repeat derived per-atom maps for every family.
`decoded_events_from_record(record)` restores the selected starts, responses and
window decisions without rerunning optimization or symmetry preparation.
Normal decoding performs no whole-reaction event-pattern deduplication. For
optional pattern comparison, `evaluate` includes a canonical event ID for the
requested physical realization.

If matching used transformed weights, pass the original endpoint problem as
`decode_events(catalogue, event_problem=original_problem)`. Elements, atom order
and coordinates must agree. Matching constraints still use `catalogue.problem`;
event responses use `events.event_problem`. `EventDecodeSession` accepts the same
keyword. When weights differ, saved records retain `problem` as the matching
endpoints and additionally identify the original `event_problem`,
`matching_problem_sha256` and `event_problem_sha256`.

`convert_fragments` and `refine_mapping` retain their existing hierarchical
rematching behavior. They are separate from this model and are never invoked by
normal decoding.
