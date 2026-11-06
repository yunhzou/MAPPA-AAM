# MAPPA design

This document defines the approved scientific design and module boundaries.
The goal is fast, accurate mapping with one-witness fragment growth, direct
compressed completion, and separate matching, reduction and decoding stages.
Implementation status, measurements and repair history belong in `reports/`;
usage details belong in `docs/`.

## Design principles

- Keep scientific operations separate: growth, completion, reduction and decoding
  have explicit inputs and outputs. External reference verification evaluates
  saved results; it does not participate in mapping.
- Store each scientific relationship once. Share immutable fragments and
  correlated symmetry actions; do not reconstruct competing representations or
  expand complete bijection tables for transport or presentation.
- Put reusable mathematics in core. Workflows, execution, viewers and benchmark
  adapters consume core operations; core never imports their orchestration.
- Keep one implementation owner per responsibility. Remove superseded code,
  modes and repair scripts after replacement. Git preserves historical source.
- Preserve explicit incompleteness. Resource limits, bounded optimization and
  the ordinary-boundary approximation cannot become exact scientific claims.
- Optimize without changing the declared mapping relation. Validate scientific
  changes against independent small examples and measure comparable work.

## Scientific pipeline

```mermaid
flowchart LR
    Input[Endpoint graphs and policies] --> Matching[Matching across seeds and cuts]
    Matching --> Reduce[Group branches and deduplicate families across seeds and cuts]
    Reduce --> Catalogue[Catalogue for one reaction and direction]
    Catalogue --> Inventory[Collect all unique fragments and branch references]
    Inventory --> Local[Shared conditional fragment responses]
    Local --> Compose[Compose approximate event counts and window]
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
required. Accepted edges are saved in each fragment relation; complete original
endpoint bonds remain available for event interpretation. Joining a completed
fragment retains its locked assignments and earlier constraints.

`SearchSession` owns prepared matching state. Splitting work must preserve seed
indices, seed order, random settings and the full search scope. An indivisible
cut backend cannot be split into independently reseeded searches.

## Binary symmetry and shared event core

A **global symmetry** describes rearrangements of one whole endpoint. A
**fragment symmetry** describes rearrangements of the fragment's own internal
graph, with its anchors fixed. A **relaxed isomorph** describes a source-to-target
mapping. Discovery tolerance controls growth; event thresholds classify original
bond-weight changes. There is no separate self-tolerance filter.

Construct each symmetry from a binary graph: bonds
above the declared presence threshold share one structural label, elements are
preserved, and each anchor is fixed. Nauty supplies compact correlated structural
generators. Every composition of those generators is permitted. Original WBOs
remain attached for discovery constraints and event counting; numerical weight
differences do not restrict a binary self-rearrangement.

Global construction uses the whole endpoint; fragment construction uses the
selected fragment's internal graph. External attachments do not become local
constraints. The binary graph preserves connectivity at its declared threshold;
weight comparisons also retain internal subthreshold contacts. Accepted growth
edges remain separate discovery constraints. Neither the binary graph nor its
generators alone establish the complete relaxed source-to-target relation.

`core/symmetry/` owns binary symmetry construction, ordered shuffle programs,
joint endpoint-pair projection and original-weight event computation.
It has no dependency on matching, catalogue records or postprocessing. Global
and fragment-local construction use this same owner. Decoder-specific inventory,
boundary bookkeeping and branch-event assembly consume its shared responses.
Pair projections are local response domains; choices across different pairs
remain correlated through the original joint action. Minimum-event computation
separates independent symmetry components only when event terms do not join
them, and skips operations only when they preserve those responses. It retains
one actual witness and explicit completion status, with bounded search and
caches. No complete permutation group or Cartesian mapping table is materialized.
Raw search retains discovery alternatives. Completed atom-set fragments retain
their potential correspondences, while decoding selects one starting orbit.
Event optimization may deliberately examine only a bounded subset of starting
correspondences. Its result is best found unless a valid lower bound is attained
or complete coverage is certified. This does not narrow the saved mapping
relation. An orbit minimum is not a minimum over a larger relation unless coverage
is established. Exhausted search limits remain incomplete, never exact minima.

### Bounded starting-correspondence improvement

One candidate identifies one source/target atom-set placement. Different starting
correspondences on that placement need not be connected by binary symmetry.
Search for a useful starting correspondence within a fixed budget; do not turn
growth routes into additional completed fragments or enumerate all starting classes.

Completed fragments are identified by their reactant/product atom-set pair,
within the fixed endpoint problem and policy. Growth witnesses and accepted
edges do not create additional completed fragments. Raw search archives retain
that evidence for replay and search validation. Explicit endpoint anchors remain
constraints; temporary locks inherited during growth do not become permanent
anchors of the completed fragment.

Use a saved element-compatible correspondence as the incumbent. First obtain a lower bound by
matching the original pair-weight multisets within element-pair types, including
zero weights and subthreshold contacts, while relaxing atom consistency. If this
bound does not settle the saved witness, also relax each atom's incident bond
assignments and match the atoms injectively. Every bond event touches two atoms,
so half that minimum incident cost, rounded up, is another valid lower bound.
Use the stronger bound and stop immediately when it meets the witness. Improve
the incumbent through the shared binary-symmetry event core. If needed, search
the completed atom sets for a better element-compatible, injective starting
correspondence, respecting explicit endpoint anchors, and apply that same event
core to it. Original growth edges do not restrict this search. Partial
original-weight event bounds account for injective occupancy of the remaining
atoms through a linear assignment; bonds between two unassigned atoms contribute
zero to that bound. These bounds guide the search; they do not prove absence of
better unexplored symmetry choices.

The default is one five-second wall-clock search deadline per unique fragment,
shared by preparation, discovery improvement and symmetry optimization. It is a
maximum allowance, not a target runtime. Stop immediately on an attained lower
bound. Otherwise retain one best actual witness and explicit incomplete optimality.
Time-limited answers may vary
with hardware; record the allowance, elapsed time, bound and work counts.

Core checks deadlines between operations and retains a valid incumbent. An
opaque native operation can finish after its deadline; hard preemption belongs
to the execution layer's external watchdog. No scheduler or per-fragment worker
creation belongs in the scientific core. Bounded witness optimization does not
establish all event bins, endpoint-conditioned responses or a whole-reaction
minimum. Decoder assembly must preserve that distinction when consuming it.

```mermaid
flowchart TD
    Structural[Binary structure, elements and anchors] --> Symmetry[Compact symmetry generators]
    Symmetry --> Global[Certified global placement acceleration]
    Symmetry --> Local[Fragment-local relationships]
    Weights[Original WBOs and event thresholds] --> Core[Shared symmetry event computation]
    Local --> Core
    Core --> Decoder[Local and boundary event decoding]
    Local --> Decoder
```

## Fragment growth and completion

Fragment matching has two operations: find one witness during growth, then
construct the complete compressed placements when growth saturates.

### Fragment records and ownership

`FragmentGrowth` owns the changing search state. `Fragment` stores one
source-to-target placement. Raw search records also carry its growth proof;
reduction creates the shared final atom-set record consumed by decoding.
There is no additional candidate or
`FragmentMatch` wrapper, or separate one-sided fragment record.

| Record | Logical fields |
|---|---|
| `FragmentGrowth` | `source_atoms`, `accepted_edges`, `deferred_edges`, `frontier`, `witness` |
| `Fragment` | `source_atoms`, `target_atoms`, initial `witness`, explicit anchors, `source_symmetry`, `target_symmetry`, `relation`; raw search fragments additionally retain `accepted_edges` |

The growth subgraph is exactly `source_atoms` with `accepted_edges`; no duplicate
subgraph representation is retained. Deferred edges are rejected extension
proposals, and the frontier contains remaining proposals. The current witness
may change its target atom set during growth. The parent branch supplies fixed
assignments and the anchor/lock context. Original molecular bonds remain in the
reaction data, separate from the accepted growth subgraph.

At saturation, source self-automorph preparation is shared across placements.
Each retained representative target atom set produces one search `Fragment`, with its
anchor-aware binary target symmetry and complete compressed relaxed
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
payloads. Global relationships used to accelerate placement discovery belong to prepared
reaction search state. A completed matching operation also supplies a shared
placement relationship to its branch transitions. This relationship connects
equivalent target atom sets; it is separate from each fragment's local relation.
Scheduling, cache transport and storage policy remain outside these records.

One matching-operation result carries the fragment collection, shared deferred
continuation and completion/cap status. Each child transition references its
chosen fragment directly in the existing hierarchy. Final reduction and decoding
retain those scientific records rather than reconstructing a second fragment
representation from dictionaries.

```mermaid
flowchart LR
    Parent[Parent branch locks and anchors] --> Growth[FragmentGrowth]
    Growth -->|Saturation: representative target sets| Fragment[Completed Fragment]
    Growth -->|Deferred edges| Continuation[Branch continuation]
    Fragment -->|Shared reference| Continuation
    Global[Certified global placement relationship] --> Continuation
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
witness to discover representative target atom sets within its fixed source
fragment, parent and anchor context. Distinct sets may share one continuation
only when a saved global relationship proves their remaining search problems
equivalent. Every represented placement remains recoverable through that
relationship. Overlapping sets remain eligible.

Placement compression uses a certified subgroup of the global binary symmetry.
It preserves the full target WBO matrix, active connectivity and
all source-to-target node compatibility, and fixes every prior locked target
pointwise. Thus every later matching constraint is preserved under composition.
Callers whose continuation depends on additional target labels declare those
labels fixed for placement compression; this does not occupy them for matching.
Binary connectivity alone does not establish this stronger condition.
Uncertified rearrangements never justify dropping a continuation.

Discover one representative per target-set orbit. Compact orbit exclusions and
conditioned symmetry pruning participate in search; do not enumerate all orbit
images or internal bijections and then merge them. The branch cap counts these
representative continuations at every hierarchy level. A representative count
is not the number of labelled target placements it represents.

Search continues with one representative assignment and ordinary occupancy. The
saved placement relationship subsequently acts on the entire continuation,
moving later fragments together. It never independently moves one fragment into
another's occupied target. Each decision applies its local relation, then its
placement action; complete programs apply decisions in reverse chronological
order. Retain that order through caches, checkpoints, reduction, mapping queries
and decoding, including unmatched and cut-boundary events.

Each candidate retains its complete compressed internal mapping relationship,
using anchor-aware fragment binary symmetry. Global symmetry
does not substitute for the fragment's own relationship. Free pools are a special
case derived from that relationship, not an independent compression rule.
One witness and one structural action are not assumed to cover every permitted
mapping on a target set. The representation must establish complete coverage
while preserving correlations and avoiding full-bijection expansion.

The discovery witness satisfies the source-to-target `iso_tol`. Binary symmetry
alternatives may exceed that discovery tolerance; original-WBO decoding evaluates
their events. Records distinguish discovery evidence from the retained symmetry.

For one fixed source/target atom-set pair, let D contain every mapping satisfying
the accepted discovery constraints. Retain both source and target permitted
self-rearrangements: each final mapping is beta composed with d composed with
alpha inverse, where d belongs to D and alpha/beta each preserve their own
fragment's binary graph, elements and anchors. An accepted discovery bond can become
an event after rearrangement. Decoding counts and classifies those changes under
the original endpoint WBOs. Its bounded starting-correspondence selection and
event window restrict the decoded response, not the saved matching relation.

A completed fragment is one atom-set pair. Its potential starting mappings are
element-compatible bijections respecting explicit endpoint anchors. Decoding
retains the bounded search's selected starting correspondence and its correlated
binary symmetry response. Potential catalogue membership and membership in a
selected decoded response are different measures; report them separately.
Independent atom choices must not introduce forbidden combinations.

Matching preserves valid mappings, accepted/deferred bonds, locks and
legitimate compressed branch choices. Duplicate records and intermediate merge
history are not scientific requirements. Regrouping choices requires the
placement-equivalence proof above; similar events, witnesses or atom labels
alone do not suffice.

For each parent branch, each retained representative completed-fragment match
produces one child continuation with its saved placement relationship. Symmetric
bijections inside that match are not separate branches.
After each completed search step, apply the cap bottom-up at every hierarchy
level. Each node counts distinct surviving continuations below it. If over the
cap, remove its largest child subtree, recount and repeat until within the cap;
ties preserve the earlier decision. Ancestors use the reduced descendant counts,
so pruning below can spare their other subtrees. A shared continuation remains
reachable through its retained parents, without retaining a rejected decision
path. Seed policies apply this rule independently even when computation is shared.
The completed-fragment limit remains an early rejection of a placement result
that exceeds the cap; it does not truncate intermediate one-witness growth.
Completion, stalls and caps remain explicit; capped search is not an exhaustive
absence proof.

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
Entries share equal immutable fragment self-preparations within one reaction.
The budget counts each shared preparation once, including its ownership index;
evicting its last entry also releases that index record.
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
working state when those policies finish. Incompatible historical fragment schemas are rejected; use their producing
revision to inspect them without changing their scientific meaning.

Graph combination concatenates and reindexes histories; it does not deduplicate
final families. Public branch grouping retains every path. Seed order, growth
history, island labels, anchors and witnesses do not split an atom-set branch;
different representative atom sets or fragment pairings do. Saved placement
relationships retain the additional labelled sets without expanding branch
records. Live continuation contexts
cannot be discarded merely because their final atom-set groups coincide.

One catalogue per endpoint problem and direction reduces all supplied seeds and
cuts. Intern each final reactant/product atom-set pair once under the fixed
policy and explicit anchors. Deduplicate branch compositions by their fragment
pairs. Intermediate overlapping growth fragments and accepted-edge histories
do not become additional completed fragments or family identities.

Retain ordered global placement alternatives as branch context: they move whole
continuations and must not permit independent overlapping occupancy. Such
alternatives can share one atom-set branch and its local fragment responses.
An initial placement that preserves every final target atom set and explicit
anchor is already included in the local binary symmetries and can be omitted.
After any placement moves an atom set, preserve the remaining action order.
Reduction performs structural key lookups and provenance aggregation. It does
not solve symmetry, search for starting maps or count events. Binary symmetry
preparation is deferred until a consumer needs it.

Project completed states to shared fragment IDs in one native operation, then
deduplicate their compositions using those IDs. Reuse source partitions and
prepared fragment actions. Consume shared immutable state blocks directly;
temporary conversion reuse is scoped to one input graph. Ordered placement
programs are collected once per shared search state. Keep growth provenance as
parent links and terminal references, with path counts computed on the shared
graph; do not unfold every
growth path during reduction, decoding or publication. Explicit history replay
may expand requested paths. Sweep setup reads the saved cut/partition metadata
directly from those references.

The catalogue preserves source histories, cut provenance and alternative fragment
partitions. Raw paths, public atom-set branches and catalogue families are separate
counts. Completed search graphs may be ingested one at a time; the accumulated
catalogue and provenance still consume memory.
Restoring a catalogue shares equal immutable mapping, constraint and action
records while preserving every family's values, action order and provenance.
Trusted binary transport restores already validated scientific fields without
repeating construction checks. Reduction reuses immutable fragment hashes and
prepared transition actions; each saved path still contributes its provenance.

## Decoding: shared fragment responses and branch queries

First inventory every unique completed Fragment and every branch, partition and
sweep reference. The matching catalogue records potential atom-set placements.
For decoding, select one starting correspondence per unique Fragment using the
bounded improvement above, then retain both sides' element- and anchor-preserving
binary symmetries. For an aligned equal-connectivity pair, one side suffices.
Other starting correspondences need not lie in this selected
symmetry orbit: this is an explicit coverage restriction, not exhaustive decoding
of every potential correspondence on the atom-set pair.

The shared core response stores compact correlated generators, a fixed event
contribution and sparse variable original-WBO event terms. Independent symmetry
components are separated only when no variable event or requested endpoint joins
them. Operations with identical event effects can share computation. Neither a
complete permutation list nor a table of every attainable event count is stored.
A query retains the response and a count interval; finding a member returns a
valid witness or explicit absence/unknown within its compute allowance. Endpoint
conditions preserve occupancy and their correlation with internal event counts.

```mermaid
flowchart LR
    Catalogue[Matching catalogue] --> Inventory[Unique fragments and branch references]
    Inventory --> Start[Bounded starting-map improvement]
    Start --> Response[Shared binary-symmetry event responses]
    Response --> Branch[Ordered branch composition and sweep correction]
    Branch --> Window[Retained branches with a queryable event window]
```

Branches remain unique compositions referencing shared final fragment responses.
Final fragments own disjoint internal bond contributions. Compose their local
choices and retain the saved global placement order, binding a shared placement
to one joint choice. Never replace an ordered placement sequence by unrestricted
group closure or allow two fragments to occupy the same target atoms.

Assembly binds the prepared cost projections to their original shared choices;
it does not expand them back into raw bond terms. Weight-preserving placement
actions retain those internal costs while moving the endpoints. The sweep adds
only its endpoint-pair correction, with endpoint dependencies traced through the
saved action order. When internal costs are constant, resolve the sweep directly
over jointly reachable endpoint pairs, retaining one valid ordered choice as its
witness. This must not enumerate complete symmetry assignments. Variable internal
costs retain their dependence on the endpoints. Reuse local event-invariant compression only when it also
preserves every affected endpoint and other contribution. Accept an attained
lower bound immediately rather than searching again for the same minimum.

The event model uses the accepted **independent ordinary-boundary approximation**.
Add local internal events to a fixed ordinary boundary contribution, including
mapped/unmapped boundaries. Do not recompute ordinary boundaries for each choice.
The single sweep correction retains its joint endpoint/internal-event dependence
and contributes a residual from zero to two events. Count it once. Alternative
partitions and sweep histories are alternative compositions, never additive ones.

For window selection, let L be the smallest composition lower bound before the
sweep residual. When attained local minima compose, an attained branch upper
bound U is at most L+2: resolve only compositions that can improve U, then retain
those potentially within the best attained count plus the requested slack n
(default six). If a local search is unfinished or local minima conflict, retain
its honest lower bound and use an attained U instead; a loose lower bound plus
two is not a valid ceiling. Report whether the model minimum is proved. An
unfinished search never becomes absence or an exact minimum. Branch queries
preserve the surviving choices rather than exporting one chosen bijection.

Events use **original endpoint WBOs**, not cut search graphs. The ordinary
threshold defaults to 0.5 and the metal threshold to 0.3; an event requires a
strictly larger weight change. Breaking, forming, strengthening and weakening
are exclusive categories. Growth tolerance is separate. Count changes touching
mapped atoms, including mapped/unmapped boundaries; omit bonds wholly inside
unmatched material. Invent no correspondence for unmatched atoms.

When matching uses transformed weights, `event_problem` supplies original weights
for the same indexed atoms and coordinates. Discovery permissions remain tied to
the matching problem. Output retains both problems when their weights differ.

`EventDecodeSession` shares fragment preparation across uses. `EventInterpretation`
retains branches, responses, minimum/proof status and window. Explicit catalogue
proof evaluation remains available for inspection, including discoveries outside
the selected decoded orbit; it does not assert membership in a stored event
selection. Physical whole-reaction evaluation may differ from the approximate
boundary model and does not silently reassign model bins.

Saved results retain the unchanged matching catalogue, each selected starting
correspondence, shared response and window decisions. Loading reconstructs query
objects without rerunning starting-map optimization or symmetry preparation.
Working branch tables are released after publication. Persist only compact
relations and sparse costs; their size is still not a bound on total reaction
memory. Scheduling and hard resource watchdogs remain outside the scientific core.

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
counts cannot establish coverage of every compressed choice. Completed catalogue
membership means the reference fits a retained atom-set composition and its
joint placements; it does not establish membership in the selected decoded
symmetry orbit or event window. Report those metrics separately when comparing
reduction semantics. Heavy-atom
references leave hydrogen correspondence unspecified. Unknown, timeout, capped
search and missing work remain distinct from proven absence.

For a completed atom-set composition, verification first tries to construct a
full injective correspondence containing the annotated reference, including
unannotated atoms. Acceptance requires validation through the saved relation.
References needing global placements or endpoint equivalence use the joint query.
That query first rejects references outside conservative per-atom bounds of
the complete ordered relation, including both local self-relationships,
placements and endpoint equivalence. These bounds cannot establish acceptance;
survivors retain every joint constraint. Reuse shared immutable action preparation
and stop work once absence is proven. Construct physical output expressions only
for callers that need a mapping or an event query.
Membership queries may absorb a self action only when it is certified to preserve
every discovery constraint. The returned mapping still carries a valid discovery
and symmetry witness. Reduced finite constraints are transferred to the native
solver in one batch; this changes query execution, not the retained mapping
relation or the reference-equivalence policy.

Matcher acceptance requires independent mapping-relation checks, meaningful
branch and cap checks, and Python/native agreement. Critical cases include C6
into C6 producing one compressed match, and C6 into two rings retaining both
placements: two representatives when attachments distinguish them, or one with a
certified relationship when the whole target permits their exchange. A proven
impossible query must return zero without unnecessary search.
Broader validation covers locks, tolerances, accepted/deferred edges, symmetry
correlations, transport, reduction, decoding and reference membership.

Performance is an acceptance condition alongside accuracy. Compare identical
scientific scope, measuring matching, completion, publication, scheduling and
memory separately. Resource measurements include the process tree and sampling
coverage; an unsampled short task is not evidence of zero memory use. A cache
budget does not bound total resident memory, and local tests do not establish
cluster scaling. Use the checks in [README.md](README.md); retain independent
scientific oracles and report the actual validation scope in `reports/`.

## Module ownership

Paths below are relative to `src/mappa/`.

Reduction and decoding use the existing C++ backend for bulk scientific
operations: assignment-domain updates and event bounds, compact permutation
operations, shared response binding and composition bookkeeping. Python owns
the public records and stage APIs; native operations consume shared integer
IDs and prepared arrays rather than calling Python inside atom-level loops.
The shared symmetry core remains the single owner of the mathematics consumed
by both fragment preparation and decoding. This execution boundary changes no
mapping permissions, event model, action order, cache limits or completeness
rules. Work buffers are bounded or released with their stage; scheduling and
storage remain outside the native scientific operations.

```mermaid
flowchart LR
    API[Python stage APIs and records] --> Native[Native reduction and event kernels]
    Native --> Shared[Shared fragment IDs, arrays and compact generators]
    API --> Storage[Aggregated persistence]
```

| Owner | Location |
|---|---|
| Shared endpoint, configuration and result records | `mappa_core/domain.py` |
| Fragment growth, existence and completion | `mappa_core/core/matching/` |
| Relaxed relationships, correlated shuffles and bond responses | `mappa_core/core/symmetry/` |
| Seed/cut planning, locking and hierarchy branching | `mappa_core/core/search/` |
| Saved-path relation compilation and assignment queries | `mappa_core/core/search/queries/` |
| Branch grouping and family reduction | `mappa_core/post_processing/final_branches.py` |
| Unique-fragment inventory and event decoding | `mappa_core/post_processing/decoder/` |
| Persistent local/Slurm pool API | `execution/` |
| Scheduling, task handoffs, recovery and storage | `mappa_core/core/scheduler/` |
| Mapping analytics and external event/reference evaluation | `mappa_workflow/mapping_analytics/` |
| Native scientific implementation and transport | `mappa_core/native/` |
| Versioned persistence | `mappa_core/utils/io/` |
| Explicit reaction-to-TS analysis | `mappa_workflow/ts_analytics/` |
| Saved-result and growth views | `presentation/` |

```mermaid
flowchart BT
    Core[Core matching, symmetry and saved-path queries]
    Post[Reduction and decoding] --> Core
    Workflow[Mapping and TS analytics] --> Core
    Workflow --> Post
    Execution[Persistent execution and storage] --> Core
    Execution --> Post
    Presentation[Viewers and CLI] --> Workflow
    Presentation --> Post
    Bench[Benchmark adapters] --> Execution
    Bench --> Workflow
```

Arrows mean dependency. Schedulers, CLIs, viewers and benchmark adapters consume
these owners rather than implementing competing algorithms. Presentation owns
rendering, offline packaging and navigation; tools only adapt inputs and write
outputs. Display alignment does not change scientific atom correspondence.

This document is the current design contract. `docs/` describes supported usage;
`bench/` owns maintained experiment runners; `reports/` preserves measurements
and provenance, with historical scripts linked to their producing Git revision.
Historical results do not establish current-method accuracy. Explicit persistence
relocations preserve unchanged records; incompatible scientific schemas are
rejected. Superseded implementations belong in Git, not fallback paths.

Usage guides: [Python API](docs/PYTHON_API.md),
[shared execution](docs/shared-execution.md), and [viewers](docs/VIEWERS.md).
