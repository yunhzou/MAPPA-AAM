# MAPPA architecture

This is the authoritative design for the current code. Usage belongs in
[the Python guide](docs/PYTHON_API.md), [execution guide](docs/shared-execution.md),
and [viewer guide](docs/VIEWERS.md). Measurements and validation records belong
in `reports/`; historical implementations and superseded designs remain in Git.

## Workflow and ownership

```mermaid
flowchart LR
    Input[Endpoint graphs and settings] --> Search[AAM search]
    Search --> Result[Search graph and grouped branches]
    Result --> Conversion[Fragment rematching and recombination]
    Conversion --> Events[Bond-event counts and classes]
    Search --> Matcher[Shared fragment matcher]
    Conversion --> Matcher
    Detection[Precursor fragment detection] --> Matcher
    Detection --> Assembly[Retrosynthesis assembly]
    Result --> View[Viewers and saved artifacts]
    Events --> View
    Events --> Reaction[Explicitly selected reaction context]
    Reaction --> TS[Transition-state analysis]
```

| Owner | Responsibility | Current implementation |
|---|---|---|
| Graph/chemistry primitives | Endpoint graphs, bond rules, input conversion | `frag.py`, `chemistry_computations/`, `smiles.py` |
| Fragment matcher | Conditional injective matches with compressed correlated choices | `fragment.py`, `matcher/`, `growth/`, native engine |
| AAM search | Fragment locking, completed-branch frontier, seeds and cuts | `aam.py`, `alignment/branch.py`, `alignment/sweep.py`, `native_search.py` |
| Representation | Shared decisions, paths, atom-set branches, saved mapping families | `domain.py`, `search_graph.py`, `alignment/post_aam.py`, `final_branches.py` |
| Symmetry | Conditioned group reconstruction and action reduction | `search_symmetry.py`, `conditioned_symmetry.py`, `action_reduction.py` |
| Decoding | Anchored rematching, recombination, direct counts and event classes | `fragment_conversion.py`, `postprocessing.py`, `event_patterns.py` |
| Explicit queries | Membership of a supplied mapping in saved correlated choices | `family_query.py`, `reference_constraints.py` |
| Execution | Task dependencies, workers, budgets and verified publication | `execution/` |
| Retrosynthesis | Detection, coverage ownership, assembly and ranking | `fragment_matching/`, `retrosynthesis/` |
| TS analysis | Supplied reaction context, endpoint-to-TS cores and scoring | `ts.py`, `core_aam.py`, `scoring.py`, `bond_diagnostics.py`, `modes.py` |
| Artifacts/presentation | Versioned persistence, molecular views and growth replay | `artifacts.py`, `viewers.py`, `search_trajectory.py`, `growth_display.py`, `static/` |

CLIs, benchmark adapters and viewers consume these owners. They must not contain
another matcher or decoder. Scheduling changes do not change scientific settings.

## Data architecture

```text
AAMProblem
  reactant, product: MolecularEndpoint
    elements, coordinates, bond-order matrix
AAMResult
  problem, config, metrics
  graph: AAMSearchGraph
    contexts: seed/cut/settings/anchor provenance
    roots: context -> initial state
    states: cumulative mapping, island ownership, continuation position
    transitions: source/target state references and committed fragment evidence
      placement: FragmentMatch
        witness, pools, fixed assignments, automorphism domains/generators
    stops: completion, stall or explicit limit
  branches: SearchBranch[]
    unordered R/P fragment atom-set pairs
    paths: all retained alternatives referencing the graph
```

The graph shares prefixes and exact reconvergence. A path references committed
fragment decisions; it is not a list of enumerated atom bijections. A viewer may
unfold the graph into a tree without changing the stored relation.

`FragmentMatch`, `SymmetryDomain`, `PermutationGroup` and `AAMHierarchy` are
existing typed views of fragment choices. Graph parsing caches these views and
shares exact permutations. A representative is one valid witness, not the entire
represented mapping set. Encounter multiplicity is not a mapping cardinality.

`FinalBranchCatalogue` is a downstream view of saved paths. `FinalFamily` retains
a representative, required bonds, ordered pool/group actions, matching policy
and provenance. Its policy-qualified storage buckets are not the public AAM
branch identity. Reading a catalogue must preserve correlations and coverage.

## Matching and growth invariants

For a fixed fragment, locked context and retained source bonds, the compressed
candidate union must represent every permitted injective match under that
context. This is representation fidelity; it does not assert that greedy search
or a capped frontier discovers every possible fragment.

A candidate stores a witness, set-valued pools, fixed roles and exact correlated
symmetry actions. Equal per-atom orbit labels do not establish equality of whole
mapping relations. Compression requires a relation-preserving certificate and
transporter. Moving a fragment transports dependent assignments consistently;
independent Cartesian products of atom destinations are invalid.

Growth selects boundary atoms uniformly without a bond-order preference. It
extends the represented relation and locks completed fragment placements. The
branch cap applies to the completed-fragment frontier/subtrees, not intermediate
one-to-many iso matches. A cap is a recorded loss of search completeness, not
an equivalence rule or an algorithmic absence result.

Live-state sharing requires compatible continuation contexts and equal state.
It retains incoming histories. Do not merge live states solely because final
atom sets or endpoint automorphisms agree. Separate cuts and seed runs retain
provenance; their contexts are not interchangeable continuation states.

`match_fragment` is the reusable boundary. The Python and native implementations
must preserve the same mapping relation for supported inputs. Explicit native
execution requires the built engine. Node policy, pair-dependent bond tolerance
and trace requirements determine applicable implementations. No failure may
silently select another scientific policy or relax a constraint.

## Final branch identity

A fragment pair is identified by its R atom set and P atom set. A public AAM
branch is the unordered collection of these pairs. Atom ordering, island labels,
seed, history, anchors and representative mappings do not split that identity.
Different labelled atom sets or different pairings remain different branches.

`AAMResult` groups before returning and retains every path/context alternative
inside each group. `retained_branch_count` counts these groups; raw terminal/path
counts remain separate. The first representative or hierarchy is only a witness.
Grouping must not discard other allowed mappings or manufacture new ones.

## Event conversion and counting

AAM matching tolerance and bond-event thresholds serve different purposes.
`convert_fragments` rematches inside the saved R/P fragment atom sets at the event
thresholds. Only directly fixed assignments become anchors. A role moved by a
joint symmetry remains free. The matcher rediscovers free choices; it does not
pin the old representative or require the old free-choice program.

Distinct pair/anchor contexts share rematching work. Identical full conversion
inputs share recombination work and retain every source-family reference. This
work cache is not branch identity. Rematching can split a fragment; compatible
results are recombined into converted compressed branches.

The ordinary event threshold is 0.5 and the metal threshold is 0.3 by default.
An unchanged bond differs by strictly less than its threshold; equality is an
event. Count breaking/weakening and forming/strengthening differences touching
mapped atoms, including mapped/unmapped boundaries. Ignore bonds wholly between
unmatched atoms. Do not invent atoms or correspondences.

`decode_events` counts all retained converted branches, then constructs event
classes within the requested window. The minimum is the lowest count among
these branches. `event_slack=0` selects the minimum; minimum+n and absolute bounds
filter retained results. These are not hidden-bijection optimization bounds.
Identical literal mappings share count work; partition/provenance checks remain
attached to each branch.

Completeness requires uncapped completed rematches, retained mapped atoms and no
counted events internal to a converted fragment. Incomplete internal-event cases
must remain flagged. Neither completion nor `minimum_proven` establishes a global
minimum across unexplored AAM branches or exhaustive old-family enumeration.
Event-class identity uses the correlated decorated-graph relation, not a union
of per-atom orbit IDs.

A user-requested shuffle or reference-membership query is separate from decoding.
It may use the optional solver on saved choices. The public conversion/counting
path does not call solver enumeration or use a solver decoder as a fallback.

## Retrosynthesis

Detection calls the shared matcher; assembly consumes saved fragment candidates.
Explicit hydrogens participate in both detection and coverage. Augmented residual
copies are competitors, not hard anchors, and projected fragments must respect
the original cut attachment. The beta big-block/gap-first workflow is an explicit
recommendation policy, not automatic failure recovery for full-bank detection.

`fragment_matching/` owns candidate generation, augmentation, projection,
serialization and optional RDKit adaptation. `retrosynthesis/` owns candidate
indexes, bounded coverage search, assembly validation and ranking. Catalog I/O
and launchers call these modules rather than duplicating their algorithms.

An assembled target atom has at most one owner. Repeated precursor copies retain
separate occurrence identities. Coverage, exact ownership and geometry are
separate checks; geometric assembly is not proof of chemical feasibility.
Chirality-insensitive screening cannot reject solely on stereochemistry;
explicit chirality-sensitive refinement may rerank or reject. Preserve cap and
incomplete diagnostics and write expensive results before discarding them.

## Transition states and presentation

TS analysis requires an explicit `ReactionContext`: selected mappings and bond
events. `core_aam.py` serves endpoint-to-TS core assignment, not a replacement
full-reaction search. Raw AAM and event conversion do not select 3D orientation or
claim a general stereochemistry guarantee from element/WBO inputs.

Viewers consume recorded mappings, fragments and symmetry. Display alignment
must not alter scientific correspondences. Saved traces show actual search
choices; compressed alternatives must not be presented as independently
selectable atom pairs. Shared reaction/catalog renderers own presentation.

## Execution, storage and evidence

One operation runs on a single CPU, local workers or cluster workers through
`execution/`. Decompose independent case/direction/seed/cut work using original
seed identities; do not create a second cluster implementation. Shared queues
own leases, dependencies and original case deadlines. Scheduling cannot reset
scientific budgets or silently retry completed/failed work.

Use node-local scratch and SQLite writes. Publish verified bounded bundles and
consistent queue backups; do not copy a live WAL database's main file alone or
create a shared file for every tiny intermediate. Preserve source/input hashes,
settings, timings, memory, attempts, output provenance and incomplete flags.
Execution success is not scientific completeness. Resource limits, unknowns,
missing outputs and unattempted tasks are never algorithmic negatives.

## Repository boundaries and remaining cleanup

`src/` and `native/` contain reusable algorithms; `tests/` contains their checks.
`docs/` contains current usage guides and executable examples. `bench/` owns
benchmark adapters; `reports/` owns measurements and provenance. `tools/` and
`hpc/` provide user/launch adapters; `manuscript/` owns paper artifacts. Historical
code is recovered from its recorded Git revision, not kept as runnable fallback
copies in the current tree. Task status belongs in the MAPPA_AAM tracker.

This document does not claim all consolidation is finished. Adaptive/shared
schedulers, explicit analytical/mechanism APIs and historical benchmark scorers
still have callers. Their removal requires migrating those consumers together;
changing directory names does not establish equivalent behavior. New code must
use the existing owners, remove superseded paths, and update this one design.
