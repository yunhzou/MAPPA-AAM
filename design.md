# MAPPA design

This document defines the approved scientific design and module boundaries.
The goal is fast, accurate mapping with one-witness fragment growth, direct
compressed completion, and separate matching, reduction and decoding stages.
Implementation status, measurements and repair history belong in `reports/`;
usage details belong in `docs/`.

## Scientific pipeline

```mermaid
flowchart LR
    Input[Endpoint graphs and policies] --> Matching[Matching across seeds and cuts]
    Matching --> Finalize[Finalize saved symmetry]
    Finalize --> Reduce[Group branches and deduplicate families across seeds and cuts]
    Reduce --> Catalogue[Catalogue for one reaction and direction]
    Catalogue --> Inventory[Collect all unique fragments and branch references]
    Inventory --> Local[Decode local fragment and boundary responses]
    Local --> Compose[Compose requested branch event results]
```

A cut is a search condition and provenance field, not the unit of decoding.
A representative mapping is one witness of a compressed family. Scheduling must
preserve the scientific stages and their meaning.
Reference verification is external evaluation logic that queries the completed
catalogue; it is not a stage of the MAPPA mapping pipeline.

## Inputs and matching rules

`AAMProblem` contains reactant and product `MolecularEndpoint` records: elements,
coordinates and bond-order matrices. Matching uses those inputs and declared
policies. Preparation must state which chemical attributes are represented;
elements and WBOs alone do not establish stereochemistry, charge or isotope
coverage.

Each direction retains its own seed/cut plan and atom-index frame. Concrete
reverse-direction witnesses may be inverted into input indexing; compressed
families and event identities retain their search-frame meaning.

Growth chooses a boundary atom uniformly, then an incident frontier edge
uniformly, without WBO ranking. Only the selected growth edge adds a matching
constraint. Other molecular bonds between included atoms are not implicitly
required. Accepted edges are saved as `preserved_bonds`; complete original
endpoint bonds remain available for event interpretation. Joining a completed
fragment retains its locked assignments and earlier constraints.

`SearchSession` owns prepared matching state. Splitting work must preserve seed
indices, seed order, random settings and the full search scope. An indivisible
cut backend cannot be split into independently reseeded searches.

## Relaxed relationships and shared symmetry core

A **relaxed global automorph** describes permitted rearrangements of one whole
endpoint. An **anchor-aware relaxed fragment self-automorph** describes permitted
rearrangements inside a selected fragment, with its anchors fixed. A **relaxed
isomorph** describes a source-to-target mapping. Discovery tolerance, relaxed-self
tolerance and event thresholds have distinct responsibilities.

Construct a relaxed self-relationship from a hyper-relaxed binary graph: bonds
above the declared presence threshold share one structural label, elements are
preserved, and each anchor is fixed. Nauty supplies compact correlated structural
generators. Retain those broad generators together with joint original-WBO
constraints; every internal pair must remain within the declared self tolerance.
Compare a final rearrangement against the original structure. The relaxed
relationship is symmetric but need not be transitive, so individually permitted
operations cannot be composed without checking the resulting arrangement.

Global construction uses the whole endpoint; fragment construction uses the
selected fragment's internal graph. External attachments do not become local
constraints. The binary graph preserves connectivity at its declared threshold;
weight comparisons also retain internal subthreshold contacts. Accepted growth
edges remain separate discovery constraints. Neither the binary graph nor its
generators alone establish the complete relaxed source-to-target relation.

`core/symmetry/` owns ordered shuffle programs, joint endpoint-pair projection,
original-weight bond responses and relaxed self-relationship construction.
It has no dependency on matching, catalogue records or postprocessing. Global
and fragment-local construction use this same owner. Decoder-specific inventory,
boundary bookkeeping and branch-event assembly consume its shared responses.
Pair projections are local response domains; choices across different pairs
remain correlated through the original joint action. No complete permutation
group is materialized.

```mermaid
flowchart TD
    Structural[Binary structure, elements and anchors] --> Relaxed[Relaxed self-relationship]
    Weights[Original WBOs and declared tolerance] --> Relaxed
    Core[Core correlated shuffles and bond responses] --> Relaxed
    Relaxed --> Global[Global placement acceleration]
    Relaxed --> Local[Fragment-local compressed relationships]
    Core --> Decoder[Local and boundary event decoding]
    Local --> Decoder
```

## Fragment growth and completion

Fragment matching has two operations: find one witness during growth, then
construct the complete compressed placements when growth saturates.

### Fragment records and ownership

`FragmentGrowth` owns the changing search state. `Fragment` is the completed
source-to-target piece consumed directly by branches and decoding. There is no
additional `FragmentMatch` wrapper or separate one-sided fragment record.

| Record | Logical fields |
|---|---|
| `FragmentGrowth` | `source_atoms`, `accepted_edges`, `deferred_edges`, `frontier`, `witness` |
| `Fragment` | `source_atoms`, `target_atoms`, `accepted_edges`, `source_self_automorph`, `target_self_automorph`, `relation` |

The growth subgraph is exactly `source_atoms` with `accepted_edges`; no duplicate
subgraph representation is retained. Deferred edges are rejected extension
proposals, and the frontier contains remaining proposals. The current witness
may change its target atom set during growth. The parent branch supplies fixed
assignments and the anchor/lock context. Original molecular bonds remain in the
reaction data, separate from the accepted growth subgraph.

At saturation, source self-automorph preparation is shared across placements.
Each distinct target atom set produces one completed `Fragment`, with its
anchor-aware relaxed target self-automorph and complete compressed relaxed
source-to-target relation. That relation owns its witness, anchor pairs and joint
mapping constraints. Accepted edges retain the discovery evidence separately
from the constraints governing retained relaxed alternatives.

When a fragment locks, deferred edges needed for further search remain with the
branch continuation. They neither constrain fragment self-automorphs nor imply
bond events. Branch state also owns cut provenance and inter-fragment boundary
context. Sharing a completed fragment does not merge these continuation states.
Cache identity and checkpoints preserve their relevant context without retaining
discarded growth histories.

Original atom lists, anchors and prepared self-relationships are stored once and
referenced where shared; the table describes logical access, not duplicate
payloads. Relaxed global automorph belongs to prepared reaction search state.
Scheduling, cache transport and storage policy remain outside these records.

```mermaid
flowchart LR
    Parent[Parent branch locks and anchors] --> Growth[FragmentGrowth]
    Growth -->|Saturation: distinct target sets| Fragment[Completed Fragment]
    Growth -->|Deferred edges| Continuation[Branch continuation]
    Fragment -->|Shared reference| Continuation
    Fragment -->|Final aggregation and inventory| Decoder[Fragment and boundary decoding]
```

### One-witness existence

Each proposal asks whether the enlarged fragment has an injective,
policy-compatible placement. Extend the current witness when possible; otherwise
search again with earlier assignments within the growing fragment free to change.
Assignments of previously completed fragments remain fixed.

The same exact solver finds a witness or proves absence. Degree, compatibility
and injectivity checks prune impossible choices without excluding valid mappings.
Only proven absence rejects an edge; a timeout or resource limit is not evidence
that no match exists. Accepted proposals replace the witness. Rejected proposals
are deferred while the remaining frontier is explored.

### Direct compressed completion

At saturation, completion uses the full accepted constraints and the current
witness to discover each distinct target atom set once within its fixed source
fragment, parent and anchor context. Relaxed global automorph information can
share discovery work but must not collapse distinct target sets. Subsequent
placements must differ in at least one target atom; overlapping sets remain
eligible. Constraints excluding an already discovered set participate in search,
with pruning conditioned on those constraints.

Each candidate retains its complete compressed internal mapping relationship,
using anchor-aware relaxed fragment self-automorph information. Global symmetry
does not substitute for the fragment's own relationship. Free pools are a special
case derived from that relationship, not an independent compression rule.
One witness and one structural action are not assumed to cover every permitted
mapping on a target set. The representation must establish complete coverage
while preserving correlations and avoiding full-bijection expansion.

The discovery witness satisfies the source-to-target `iso_tol`. Saved relaxed
alternatives retain their own declared constraints and may exceed that discovery
tolerance; original-WBO decoding evaluates their events. Records must distinguish
discovery evidence from constraints on the retained relationship.

A compressed match preserves the allowed injective mapping relation through its
witness, pools, fixed roles and ordered correlated actions. Equal atom sets,
equal orbit labels or equal witnesses alone do not prove that two relations are
identical. Independent atom choices must not introduce forbidden combinations.

The replacement preserves valid mappings, accepted/deferred bonds, locks and
legitimate compressed branch choices. Duplicate records and intermediate merge
history are not scientific requirements. Arbitrarily regrouping distinct choices
is also unacceptable: it can change which searches encounter a branch cap.

For each parent branch, each retained completed-fragment match produces one child
continuation. Symmetric bijections inside that match are not separate branches.
The existing cap applies at every hierarchy level, including the combined live
leaves across parents. Completion, stalls and caps remain explicit; capped search
is not an exhaustive absence proof.

```mermaid
flowchart TD
    Frontier{Remaining frontier?} -->|Yes| Proposal[Propose one growth edge]
    Proposal --> Query[Find one witness with completed locks fixed]
    Query -->|Found| Accept[Accept edge and replace witness]
    Query -->|Absent| Defer[Defer edge]
    Accept --> Frontier
    Defer --> Frontier
    Frontier -->|No| Complete[Construct each distinct compressed match once]
    Complete --> Branch[Existing hierarchy branching and caps]
```

| Core module | Responsibility |
|---|---|
| `growth/engine.py` and native growth loop | Select proposals and record accepted/deferred constraints. |
| `growth/existence.py` and native `fragment_query.h` | Find one witness or prove absence, reusing the current witness where possible. |
| `growth/completion.py` and native `fragment_completion.h` | Construct distinct compressed matches from the completed constraints. |
| Shared symmetry core and matching representation | Supply relaxed global and fragment-local relationships and compact result records. |
| Reaction search | Lock completed placements, form child branches and enforce hierarchy caps. |

There is one growth algorithm, with no legacy mode or error-triggered fallback.
Candidate-growth replay, generate-then-deduplicate completion, candidate-prefix
checkpoints, candidate-deduplication metadata and symmetry-certificate retention
are outside this design. Cross-seed/cross-cut family reduction is a separate,
required stage.

Python and native execution obey the same scientific contract. Native execution
is the normal path; explicit Python execution supports diagnostics and custom
policies. Unsupported capabilities and execution failures are reported, never
silently routed to another algorithm. Historical records keep their recorded
scientific meaning; loading cannot reinterpret decoding as rematching.

## Completed-fragment reuse

A bounded cache surrounds fragment matching. Reuse requires the same endpoint
problem and policies, with compatible topology, occupancy, locks, joined islands
and deferred boundaries. Cache eviction, disabled retention or an oversized
entry may increase computation but cannot change matching or caps.

The cache contains completed compressed results and their dependency evidence.
It does not retain expanded bijections, candidate prefixes or symmetry
certificates. Its budget bounds retained entries, not total process memory.
Switching reactions releases prepared reaction-specific state.

### Shared reaction cache snapshots

Shared snapshots belong to the execution layer; the matching algorithm exports
and imports cache data without scheduling workers or performing shared-file I/O.
Each reaction/configuration has an immutable first-snapshot **baseline** and a
replaceable **latest** snapshot, with separate directional content. Workers check
the latest published version at task boundaries and reuse their resident version
when current. A failed latest load may use the baseline after integrity and
context checks. This approved data-recovery path is not an algorithm fallback.
If both stored versions are unusable, the error remains explicit.

One coordinator owns publication. Only an accepted publisher writes a complete
staged snapshot; the coordinator replaces the latest snapshot atomically. Readers
retain their opened version while replacement occurs. Stale publishers cannot
overwrite newer versions. Only changed retained entries are published, with
updates coalesced to limit storage overhead. Scientific results are published
independently of cache snapshots.

Snapshot transport validates checksums, endpoint identity, policies and code
compatibility, and reapplies the resident memory budget. Cache I/O, hits and
recovery are measured. Persistent snapshot count is bounded per reaction, never
per seed, sweep or cut.

## Search records and family reduction

The retained representation has three distinct levels:

| Level | Identity and responsibility |
|---|---|
| `AAMSearchGraph` | Contexts, states, transitions, roots and stop evidence, including compressed placements and seed/cut histories. |
| `AAMResult.branches` | Paths grouped by the unordered collection of final reactant/product fragment atom-set pairs. |
| `FinalBranchCatalogue` | Saved mapping families, their alternative fragment partitions and provenance. |

Native transport returns the typed search graph with immutable shared state
blocks. Binary checkpoints and pool handoffs preserve that sharing. Explicit
record/JSON export may expand views and need substantially more memory. Storage
sharing does not merge scientific alternatives or enumerate fragment bijections.

Matching retains the full ancestry of live continuations and returned results,
including every incoming path at a join. Once a path has no surviving descendant,
its mapping, fragment and boundary payloads are released during ordinary search,
before transport or publication. Cap evidence records the context, seed, step,
stage, attempted count, limit and occurrence count without retaining rejected
child states. A cap summary is not a returned mapping. Shared-policy searches
also retain decisions still reusable by their active policies, and release that
working state when those policies finish. Historical archives keep their original
recorded contents when loaded.

Graph combination concatenates and reindexes histories; it does not deduplicate
final families. Public branch grouping retains every path. Seed order, growth
history, island labels, anchors and witnesses do not split an atom-set branch;
different labelled atom sets or fragment pairings do. Live continuation contexts
cannot be discarded merely because their final atom-set groups coincide.

One catalogue per endpoint problem and direction reduces all supplied seeds and
cuts. Its family key includes matching policies, representative mapping, required
preserved edges and the canonical ordered pool/group action program. Only proven
commuting actions may be reordered. Distinct constraints or correlations remain
distinct even when witnesses coincide. This is conservative representation-based
reduction, not a general semantic normal form for all equivalent mapping sets.

The catalogue preserves source histories, cut provenance and alternative fragment
partitions. Raw paths, public atom-set branches and catalogue families are separate
counts. Finalized search graphs may be ingested one at a time; the accumulated
catalogue and provenance still consume memory.
Restoring a catalogue shares equal immutable mapping, constraint and action
records while preserving every family's values, action order and provenance.

## Decoding: unique fragments first

For one reaction and direction, first collect **all unique saved fragments** from
the completed catalogue, together with every branch occurrence, partition,
constraint, action and boundary reference. Complete this inventory before
computing event responses. Structural sharing does not discard restrictions or
make correlated choices independent.

```mermaid
flowchart LR
    Catalogue[Completed catalogue] --> Inventory[Unique fragments and branch references]
    Inventory --> Internal[Shared internal fragment responses]
    Inventory --> Boundary[Shared inter-fragment and unmatched-boundary responses]
    Internal --> Compose[Compose requested valid branch choices]
    Boundary --> Compose
    Catalogue --> Constraints[Original joint constraints and ordered actions]
    Constraints --> Compose
```

`decode_events` interprets the retained mapping relation. It does not rematch,
replace choices, enumerate complete bijections or select a minimum-event mapping.
Local responses are shared only when their fragment, witness, relevant ordered
actions and cut/boundary context agree. A shared action remains one joint choice
across all affected fragments. Required edges, target occupancy and action
dependencies remain attached to their original families.

Preparation is shared with the responses: repeated fragment occurrences reuse
their prepared symmetry and reachable atom sets, rather than deriving them again
before each lookup. Branches bind shared local choices to their original action
positions and retain their own boundary and constraint relationships. Shared
preparation must preserve action order and remain bounded; it does not retain
intermediate symmetry histories or complete working tables for every branch.

Local tables cover internal bonds, bonds between fragments and mapped/unmapped
boundaries. Evaluation validates a requested joint choice, then sums each bond's
contribution once. Alternative partitions decompose the same family and must not
have their totals added together. Constant responses may avoid pair expansion;
they do not permit independent atom choices or removal of correlated actions.
Preparation follows reachable endpoint pairs rather than enumerating full group
elements, though pair-table size can still grow with the reachable domains.

Events use **original endpoint WBOs**, not cut search graphs. A cut or fragment
boundary does not itself imply an event. The ordinary threshold defaults to 0.5
and the metal threshold to 0.3; an event requires a strictly larger WBO change.
Breaking, forming, strengthening and weakening are exclusive categories.
Matching tolerance is separate. Count changes touching mapped atoms, including
mapped/unmapped boundaries, and omit bonds wholly inside unmatched material.
No correspondence is invented for unmatched atoms.

When matching uses transformed bond weights, `event_problem` supplies the
original scoring endpoints to decoding. Both problems must retain the same
indexed atoms and coordinates. Required edges remain constraints on the matching
problem; event responses use the scoring problem. Saved results identify and
retain both endpoint problems when their weights differ.

`EventDecodeSession` owns reusable endpoint preparation and bounded caches.
`EventInterpretation` owns the inventory, shared local responses, dependencies
and family/partition views. These reaction-specific results are released after
publication. Binary export preserves shared tables; explicit JSON may expand
references.

Export processes family/partition working tables one at a time and shares equal
saved metadata without merging families or their choices. Arrays used to evaluate
joint choices are prepared only when evaluation is requested. The retained
catalogue and exported records still consume memory; transient-table release is
not a fixed bound on total decoding memory.
Derived per-atom domain and action lookup maps belong to working tables; saved
output retains the original catalogue and local responses without repeating these
maps for every family or partition.

Normal decoding/export retains the catalogue and local composition without a
whole-reaction event-pattern deduplication pass or automatic representative
selection. Explicit whole-reaction queries may compare event patterns. Different
fragment compositions can yield identical bond changes, so decoding does not
claim exhaustive unique event patterns or a proven global minimum. Experimental
rematching and minimum-event selection are separate operations.

## Local and distributed execution

Simple local calls remain `search_aam`, `deduplicate` and `decode_events`.
`mappa.execution.WorkPool` exposes the same separate stages through typed
`MatchingResult`, `CatalogueResult` and `DecodingResult` handoffs. Execution
strategy must not change scientific settings, task identities or results.

```mermaid
flowchart TB
    User[Direct calls or explicit stage submissions] --> Core[Matching, reduction and decoding]
    User --> Pool[Persistent WorkPool]
    Pool --> Local[Local processes]
    Pool --> Slurm[One srun inside an existing allocation]
    Local --> Adapter[Typed stage adapters]
    Slurm --> Adapter
    Adapter --> Core
    Adapter --> Storage[MemoryStore or aggregated BundleStore]
```

The pool starts persistent workers once, validates source/engine identity and
waits for every requested worker to be ready. Startup failure is explicit.
Slurm does not submit a job per scientific task. Dependencies, scheduling,
leases, budgets, transport and storage belong entirely to the execution layer.
Workers run one compute process each without nested pools and reuse prepared
state within a reaction. Affinity can prefer the same reaction while allowing
other ready work to saturate the pool; it cannot alter search order or policies.

`MemoryStore` communicates intermediates in RAM. `BundleStore` persists immutable
ranges in aggregated worker bundles, with explicit reaction indexes. Persistent
file count grows with workers and byte volume, never with seed, sweep or cut
count. Small payloads may remain in bounded memory spools; larger ones use local
scratch. Durable publication retains integrity checks and shared state.

Massive benchmarks finish and audit matching and cross-seed/cross-cut reduction
before starting decoding. This barrier is scheduling policy, not a requirement
to load the entire dataset into RAM. The decoding job is one reaction/direction
catalogue, with original endpoints, policies, integrity metadata and a manifest
of completed, capped, failed or missing search work. It loads that catalogue,
collects its unique fragments, decodes, publishes and releases the result.
A cut or isolated fragment is not an independent decoding job without a design
that preserves shared constraints and actions.

The coordinator maintains a bounded window of active reactions and an aggregated
campaign ledger. Retired task specifications and attempt histories remain
archived. Separate allocations may own disjoint reaction batches; completion of
one allocation does not imply completion of the matching pass. A fresh benchmark
starts with empty results and caches and imports immutable endpoint/reference
inputs. Re-decoding saved catalogues under a new event policy does not rerun
matching.

Recovery preserves completed outputs, original task identities, attempt history
and consumed compute budgets. Case time counts overlapping computation once;
publication and idle time do not replenish it. Extra retry time or a changed
resource allowance must be explicit. Failure, resource limits, lost leases and
incomplete scope remain distinguishable from scientific absence. Coordinator
loss requires explicit recovery.

Cross-node execution requires compatible code, inputs, artifact paths and clocks.
Storage-root relocation preserves offsets and checksums, including nested
references; it cannot change saved scientific content. Checkpoint migration stops
source writers and gives destination workers disjoint ownership.

## Verification and performance

Reference verification queries the retained mapping relation under an explicit
atom subset and endpoint-equivalence policy. One witness, element sets or event
counts cannot establish coverage of every compressed choice. Heavy-atom
references leave hydrogen correspondence unspecified. Unknown, timeout, capped
search and missing work remain distinct from proven absence.

Matcher acceptance requires independent mapping-relation checks, meaningful
branch and cap checks, and Python/native agreement. Critical cases include C6
into C6 producing one compressed match, C6 into two separate C6 rings producing
two, and a proven impossible query returning zero without unnecessary search.
Broader validation covers locks, tolerances, accepted/deferred edges, symmetry
correlations, finalization, reduction, decoding and reference membership.

Performance is an acceptance condition alongside accuracy. Compare identical
scientific scope, measuring matching, completion, publication, scheduling and
memory separately. Resource measurements include the process tree and sampling
coverage; an unsampled short task is not evidence of zero memory use. A cache
budget does not bound total resident memory, and local tests do not establish
cluster scaling. Use the checks in [README.md](README.md); retain independent
scientific oracles and report the actual validation scope in `reports/`.

## Module ownership

Paths below are relative to `src/mappa/`.

| Owner | Location |
|---|---|
| Shared endpoint, configuration and result records | `mappa_core/domain.py` |
| Fragment growth, existence and completion | `mappa_core/core/matching/` |
| Relaxed relationships, correlated shuffles and bond responses | `mappa_core/core/symmetry/` |
| Seed/cut planning, locking and hierarchy branching | `mappa_core/core/search/` |
| Branch grouping and family reduction | `mappa_core/post_processing/final_branches.py` |
| Unique-fragment inventory and event decoding | `mappa_core/post_processing/decoder/` |
| Persistent local/Slurm pool API | `execution/` |
| Scheduling, task handoffs, recovery and storage | `mappa_core/core/scheduler/` |
| Mapping and event-membership queries | `mappa_workflow/mapping_analytics/` |
| Native scientific implementation and transport | `mappa_core/native/` |
| Versioned persistence | `mappa_core/utils/io/` |
| Explicit reaction-to-TS analysis | `mappa_workflow/ts_analytics/` |
| Saved-result and growth views | `presentation/` |

Schedulers, CLIs, viewers and benchmark adapters consume these owners rather than
implementing competing algorithms. Matching does not depend on presentation or
benchmark orchestration. Display alignment does not change scientific atom
correspondence. Superseded implementations belong in Git, not fallback paths.

Usage guides: [Python API](docs/PYTHON_API.md),
[shared execution](docs/shared-execution.md), and [viewers](docs/VIEWERS.md).
