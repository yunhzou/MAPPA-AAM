# MAPPA source refactoring plan

Status: matching and reaction-search stages functionally validated; scheduler and later stages pending.

Recorded: 2026-10-28. Inspected source: `19c73fc`.

This document describes the proposed refactor. `design.md` describes the current
implementation until each corresponding change is implemented and documented.

## Central boundary

Fragment iso matching, fragment growth and their compressed symmetry
representation belong together. Reaction search uses that engine to construct
branches. Execution scheduling is a third responsibility.

- **Matching:** determine and extend allowed fragment matches while preserving
  correlated compressed choices.
- **Search:** choose seeds and cuts, commit/lock matched fragments, construct
  branches, apply the branch cap and record search history.
- **Scheduler:** execute independent tasks locally or across workers while
  preserving dependencies, budgets and provenance.

Matching must not depend on the overall reaction-search scheduler. Workflows
consume core functionality; the core must not depend on workflows.

## Target structure

```text
src/mappa/
├── mappa_core/
│   ├── core/
│   │   ├── matching/
│   │   ├── search/
│   │   └── scheduler/
│   ├── post_processing/
│   │   ├── final_branches.py
│   │   └── decoder/
│   ├── native/
│   ├── utils/
│   │   └── io/
│   └── domain.py
├── mappa_workflow/
│   ├── mapping_analytics/
│   └── ts_analytics/
├── presentation/
├── __init__.py
└── cli.py
```

Folders define responsibilities, not quotas for new files. Merge small private
helpers into their consumers. Keep substantial, distinct responsibilities in
separate files. Do not create a new file for every function or class.

## 1. Consolidate the matching engine

Bring `fragment.py`, the shared graph preparation from `subgraph.py`,
`matcher/`, `growth/` and matching-related
symmetry code into `mappa_core/core/matching/`.

This owner contains fragment iso matching, growth, compressed candidates,
assignment blocks, automorphisms, extension, compatibility checks, equivalent
candidate compression, node policies and bond-tolerance rules.

Required refactoring:

- Give fragment growth and fixed-fragment matching one shared graph-preparation
  helper. Currently `fragment.py` imports a private helper from `subgraph.py`.
- Consolidate repeated construction/parsing of fragment records.
- Keep conditioned symmetry and its accelerated workspace together, with shared
  automorphism helpers owned here. Graph traversal that attaches symmetry to
  saved search transitions belongs with search and calls this functionality.
- Preserve distinct live and saved representations where their lifetimes differ;
  do not force them into one mutable object merely to reduce class count.

The matching interface returns compressed fragment placements. It does not
manage the overall reaction's branch frontier.

Implemented as one matching package with dependency levels:
`representation/` holds candidate state and fragment contracts; `symmetry/`
holds automorphism and canonical-equivalence operations; `growth/` extends and
compresses candidates using those lower layers. `api.py` is the public entry;
`policy.py` defines shared matching rules. Internal helpers are no longer
re-exported from the package facade. Node and bond policies share `policy.py`;
fragment configurations/results share `records.py`; symmetry workspaces share
`symmetry.py`. Existing matching and growth function bodies are retained.
`subgraph.py` calls full reaction search, so only its graph preparation moved.
Saved class names resolve through an explicit relocation table in the checkpoint
reader; superseded source modules are removed.

## 2. Consolidate reaction search

Bring the search parts of `aam.py`, `alignment/branch.py`, seed/cut/direction
planning, `search_graph.py`, `native_search.py` and `shared_seed_search.py` under
`mappa_core/core/search/`.

This owner chooses seeds and cuts, calls matching, locks fragments, constructs
branches, applies branch limits and records search history. Native branch
scheduling is an implementation of this search logic; it is distinct from the
worker/task scheduler.

Required refactoring:

- Separate scientific search operations from process-pool management in `aam.py`.
- Extract the live condition/result classes from `adaptive_search.py` beside
  their search consumers. Remove the unused experimental search engine and its
  dedicated `fragment_choices.py` wrapper after checking remaining consumers.
- Retain the shared seed-search behavior and its represented alternatives.
- Keep reusable partial/core-atom mapping in search; move TS interpretation out.
- Preserve optional fragment competition/backbiting in `competition.py` pending
  a separate feature decision. Optional functionality is not automatically dead.

Implemented under `mappa_core/core/search/`: planning, branch construction,
saved graph, per-cut session, shared seed policies, reuse/native adapters, partial
mapping, subgraph search and optional competition. The unused adaptive closure
engine, its wrapper/native adapter and unused historical branch/symmetry helpers
are removed. Live condition/result classes remain with shared seed search.
Local pool loops and worker state were extracted to `execution/local_search.py`;
the next stage owns scheduler consolidation. Existing task/budget policy is unchanged.

## 3. Consolidate execution scheduling

Bring `execution/` and process-pool management embedded in `aam.py` under
`mappa_core/core/scheduler/`.

Use the same scientific task operation for single-process, local-parallel and
cluster execution. Consolidate task definitions and execution management rather
than maintaining separate local and cluster pipelines.

Preserve seed identities, task dependencies, shared directional budgets,
checkpoint behavior, memory accounting and publication guarantees. Do not create
nested worker pools or reset budgets when moving work between workers.

Keep work-reuse caches beside the algorithm that knows whether cached work is
still valid. The task scheduler does not own all matching/search caches.
Cluster-specific submission scripts remain under `bench/`.

## 4. Clean postprocessing and decoding

Under `mappa_core/post_processing/`:

- Merge `action_reduction.py` into `final_branches.py` as private helpers.
- Keep final branch grouping before the AAM result is returned. Folder placement
  does not change execution order or branch identity.
- Put anchored rematching, compatible recombination, direct bond-event counting
  and event-class construction under `decoder/`.
- Separate decoder configuration/result definitions from conversion. Currently
  `fragment_conversion.py` imports configuration from `postprocessing.py`, which
  imports conversion in return.
- Move explicit mapping queries out of `DecodedEvents.query()` and move symbolic
  event-query functions out of `event_patterns.py` into mapping analytics.

Normal decoding must not import or invoke the solver-query implementation.
Any public query API change must be documented and its callers updated; do not
leave a forwarding method that makes core depend on workflows.

## 5. Separate workflows, shared data and I/O

### Mapping analytics

`mappa_workflow/mapping_analytics/` owns `family_query.py`,
`reference_constraints.py`, explicit shuffle/mapping queries and symbolic
same-event queries. It consumes saved core results without modifying matching.

### TS analytics

`mappa_workflow/ts_analytics/` owns TS processing, scoring, bond diagnostics,
TS-specific mode analysis and TS-specific configuration/result classes.
Reusable partial mapping remains in core. Generic mode-file parsing belongs
with input adapters so core I/O does not import the TS workflow.

### Shared data and persistence

Keep shared molecular/reaction data in `mappa_core/domain.py`. Put
algorithm-specific records beside their owners. Preserve the existing meanings
of `_SymCand`, `_SymBlock`, `FragmentPlacement`, `FragmentMatch`, `SymmetryDomain`,
`AAMHierarchy`, `_Branch`, `AAMSearchGraph`, `SearchPath`, `SearchBranch` and
`AAMResult`.

Consolidate endpoint serialization, checkpoints, SMILES, XYZ and xTB adapters
under `utils/io/`. Keep generic serialization separate from workflow-specific
result formatting. Core I/O must not import TS analytics or viewers.

### Presentation

Move existing viewers and static assets under `presentation/`. Separate HTML
bundle generation from data persistence. Preserve behavior; animation/viewer
redesign and chirality feature work remain deferred. Keep `skills/` unchanged.

Alignment ownership cleanup: shared hierarchy/permutation types now live in
`mappa_core/domain.py`, and coordinate alignment/interpolation in
`presentation/interpolation.py`. The old `alignment/` package and two obsolete
interpolation aliases are removed. Core data and presentation no longer share
an artificial owner. Remaining root-domain, I/O and viewer consolidation is
still pending; this step preserves data semantics and interpolation behavior.

## Files to dismantle rather than simply relocate

| Current file | Distribute its responsibilities to |
|---|---|
| `alignment/sweep.py` | Search cut planning, matching symmetry helpers and partial-assignment handling |
| `frag.py` | Shared graph primitives/preparation, event counting and input adapters |
| `artifacts.py` | Core persistence, workflow result formatting and presentation |
| `postprocessing.py` | Decoder configuration/results/counting and mapping analytics queries |
| `event_patterns.py` | Direct event counting/classes and symbolic mapping-query constraints |
| `adaptive_search.py` | Live search condition/result types; retire unused experimental engine |

Delete emptied modules and obsolete exports. Older helpers in
`alignment/branch.py` and older symmetry-completion functions in
`alignment/sweep.py` had no production callers in the inspected tree; confirm
remaining consumers before removing them and their dedicated tests. Static
absence of a caller is evidence to inspect, not proof about external users.

## Invariants and scope

- Preserve compressed mapping relations, correlated symmetry, direct fixed
  assignments, growth policy, branch-cap behavior and event rules.
- Preserve branch identity: an unordered collection of R/P fragment atom-set
  pairs, retaining the alternatives within each grouped branch.
- Preserve scientific execution order, resource/completeness flags and
  provenance. Do not turn incomplete outcomes into negatives.
- Do not replace the algorithm to obtain a cleaner directory layout.
- Remove superseded implementations completely. No old-path forwarding files,
  fallback engines or duplicate local/cluster implementations.
- Keep a deliberate public export surface in `mappa.__init__`.
- Check saved checkpoint class paths before moving types. Use an explicit,
  documented data-format relocation/conversion policy, not silent fallback
  loading or a second runtime implementation. Existing saved AAM evidence must
  not be lost or reinterpreted.
- Keep native source and build configuration synchronized with Python ownership.
  Preserve third-party licenses. Remove native exports only after checking
  internal native consumers as well as Python callers.

## Implementation order and completion criteria

1. Consolidate matching and its data boundaries.
2. Consolidate search and separate execution management.
3. Unify the task scheduler around those search operations.
4. Separate postprocessing/decoder and workflow analytics.
5. Finish shared I/O, presentation packaging and build/import cleanup.

Each step should leave a usable codebase. Update `design.md`, current callers,
public API documentation and the existing MAPPA source-cleanup ticket as each
step is completed. Do not replace the current design with an unimplemented one.

After each module is refactored, run its functional tests before declaring that
stage complete. This supersedes the earlier request to limit validation to
syntax/import/build checks. Matching requires the existing exhaustive compression
oracles and Python/native comparison tests. Search requires branch/cap, partial
mapping, seed/cut, parallel execution and checkpoint-resume tests, plus the public
workflow smoke checks. Fix failures and record the exact tested revision, counts,
skips and scope. Do not rerun dataset benchmarks or restart monitoring for this
validation. Preserve the core algorithm while fixing refactoring regressions.

## Completed-stage functional validation

Matching and reaction search have passed 492 selected functional tests across
39 files on NRT, with no skips. Production source is `fa67d7a`; test corrections
are `657f742` and `a8aee6c`. The main run passed 491 tests; the one failing test
passed after correcting its moved patch target. An earlier collection failure
was fixed by retaining a shared molecule fixture in the test-input module.
Production source did not change during these corrections.

Coverage includes exhaustive small-graph compression, Python/native matching,
ordinary/metal tolerances, branch caps, partial mapping, shared seed policies,
cut reuse, parallel execution, checkpoint resume, saved hierarchies and public
workflow integration. This is module-level functional validation, not a full
repository test run or dataset benchmark. NRT jobs: `68515` and `68563`.
Logs and exact test identities are in Desktop
`reports/source_audit_20261028/functional-validation/summary.json` and ticket 25.

## Pickup information

- Canonical publication checkout at planning time:
  `/private/tmp/mappa-exact-fragment-publish`.
- Repository: `git@github.com:yunhzou/MAPPA-AAM.git`; published branch: `main`.
- Desktop checkout `/Users/yunhengz/Desktop/AAM Writing` contains other work;
  do not reset it or overwrite its changes.
- Tracker: `/Users/yunhengz/reactivo/yunheng-task-tracker`, branch `MAPPA_AAM`,
  existing ticket `25_source_ownership_cleanup.md`.
- Re-read `AGENTS.md`, this plan, current `design.md` and Git status before
  starting. Inspect current code rather than assuming the planning revision is
  still current.
- Implementation was authorized after the plan. Complete one major element at
  a time. Preserve core logic; remove only established redundancies.
