# Golden: fresh one-seed and ten-seed measurements

This paired full-dataset study measures revision
[`27b18ab`](https://github.com/yunhzou/MAPPA-AAM/commit/27b18abad691648cf2a87d2e5ae43db117d80b66),
including the native reduction and event kernels and the incident-bond event bound. Each arm contains all
**1,851 reactions / 3,702 directions**. Matching, reduction, decoding, external
reference verification and reaction caches start empty in both arms. Only
immutable prepared endpoint and reference inputs are reused.

**Status: both matching, reduction and reference-verification arms are complete.
Fresh one-seed decoding at `a5f1f65` completed all 3,702 directions. Fresh ten-seed
decoding stopped admission after compute timeouts and has drained: 732 successful
outputs, six timeouts and 2,964 unattempted directions. Ten-seed decoding remains
incomplete.** Earlier memory failures and all consumed budgets remain preserved.
The window remains model minimum through **minimum + 5**.

## Protocol

| Setting | Both arms, except where stated |
|---|---|
| Seeds | 1 versus 10; random seed 42 |
| Discovery tolerance / branch cap | 2 / 100 |
| Search scope | Both directions; uncut plus full single-edge sweep |
| Bond weights | Original prepared matrices, unchanged |
| Graph / cut presence threshold | 0.2 / 0.2 |
| Ordinary / metal event threshold | 0.5 / 0.3 |
| Decoding window | Model minimum through minimum + 5 |
| Starting-map improvement allowance | Five seconds per unique fragment |
| Scheduling | One persistent 32-worker pool per arm; one numerical thread per worker |
| Memory guard / direction compute allowance | 2 GiB / 600 seconds; no automatic retries |
| Completed-fragment cache budget | 64 MiB; shared reaction snapshots enabled |
| Seed batching | All seeds of one cut together: 1 or 10 respectively |
| Admission window / runner allowance | 128 directions / 9,600 seconds |

The same portable launcher and execution policy serve both arms. Each arm must
finish its full matching and reduction scope before decoding. Cache snapshots
and results are isolated between arms. Grouping all seeds of a cut allows the
ten-seed arm to reuse its resident cache; it does not change seed selection or
branch pruning. The final Linux build passed **277 affected tests in 23.96 seconds** before
submission. On macOS, the full suite passed 1,288 tests; its six failures were
stale imports of a test helper removed in the preceding cleanup. After correcting
those tests to use the supported native operation, all 31 symmetry tests passed.
The independent scientific oracles were retained. Full-suite and paired native-kernel validation are documented in the
[preceding implementation study](../golden_seed10_tol2_cap100_20261105/README.md).
All 32 workers registered in each pool, in 3.672 seconds for one seed and 3.489
seconds for ten seeds. Registration does not establish continuous utilization.
Both arms' 5,553 distinct input/reference blobs passed checksum validation.
Exact settings and Linux source/engine identity are in [protocol.json](protocol.json).

## Evaluation and timing definitions

Primary accuracy is reference recovery in **either direction**, over all 1,851
reactions, under Golden's established heavy-atom and endpoint-equivalence policy.
This checks the retained matching catalogue. It does not establish that the
reference occurs in the selected decoded symmetry orbit or reported event
window. Unrecovered, unknown and failed outcomes remain separate.

The decoder uses binary fragment symmetry with original-weight event counting,
bounded starting-map improvement, the ordinary-boundary approximation and the
joint sweep correction. Minimum proof status and incomplete work remain explicit.

Archived task records retain compute CPU and wall time separately from output
publication, per-seed timing and cache counters, cut scope, sampled process-tree
memory and resource failures. Allocation startup, observed task overlap, stage
elapsed time and scheduler accounting are separate measurements. Summed task
wall time is not parallel elapsed time; process-lifetime peak memory is not a
fresh per-task peak.

The final comparison will report full-scope completion and reaction recovery,
paired gains/losses, separate matching and reduction totals, decoding totals,
publication and verification costs, timing distributions including maxima, and
cache/scaling measurements. Failed or missing work must not disappear from the
denominator. Timing comparisons will identify their common successful scope.
The one-versus-ten seed comparison measures the whole configured system; it is
not by itself a cache-on versus cache-off causal experiment.

Operational configuration, scheduler identifiers and raw logs remain in the
Git-ignored private campaign evidence. Scientific input hashes, version identity
and sanitized measurements accompany the report.

## Initial attempt: stopped before accuracy evaluation

The initial attempt used revision `040f817`. Both pools drained without resource
failures or retries. One seed completed
3,702 matching directions, 3,702 reductions and 1,334 decodes; ten seeds completed
2,239 matching directions. Neither arm reached reference verification. These
partial runs are not full-dataset accuracy or timing results.

Case 890 forward took 137.726 compute CPU seconds in this attempt, compared with
3.025 seconds in the historical one-seed study. Those revisions differ in
completed-fragment semantics as well as implementation, so this difference does
not isolate the native port. A separate profile on the same saved matching
confirmed 137.436 CPU seconds in decoding, with only 0.384 seconds in record
construction. Of 180 unique fragments, 27 exhausted their five-second
starting-correspondence allowance; two more finished without proving optimality.
Catalogue and decoded-window reference checks both recovered the reference;
the model minimum was ten events.

A separate bound experiment used each atom's incident weight profiles, followed
by an injective atom assignment. Half the relaxed incident cost, rounded up,
is a valid event lower bound because each actual bond event touches two atoms.
It met the saved witness on all 29 unproved fragments, taking 0.569 CPU seconds
in the standalone Python experiment. This diagnoses avoidable search under the
weaker pair-weight bound; it is not a full-decoder timing or a completed repair
benchmark. The bound was subsequently integrated in `27b18ab` without changing allowed
mappings, event thresholds or the five-second allowance. Independent exhaustive
small examples validate that it never exceeds the actual optimum. On saved
matching for case 890 forward, current decoding plus record construction takes
1.358 CPU seconds, with both catalogue and decoded-window recovery retained.
Other selected cases retain a slow bounded-search tail; the complete rerun will
measure that cost without changing the search policy.

The ten selected slow directions all completed their point validation, retaining
both catalogue and selected-window reference recovery. Combined baseline stage
compute was 1,521.666 CPU seconds; revised decoding plus record construction was
663.782 CPU seconds. No selected case became slower, but four still exceeded
125 seconds. These selected cases do not establish a full-dataset speedup.
Raw stage timings, output counts and memory measurements are in
[bound-validation.json](bound-validation.json).

## Completed one-seed result

All 3,702 directions completed matching, reduction, decoding and external
verification: 14,808 terminal stage records, zero failed tasks and zero retries.
Reference recovery is **1,841/1,851 reactions (99.46%)**. Directional checks comprise
3,650 recovered, 50 not recovered and two unknown; both unknown directions have
a recovered opposite direction. Ten reactions are not recovered.

Pool elapsed time was **1,147.929 seconds (19m08s)** with 32 workers. Stage compute
CPU was 1,427.875 seconds matching, 794.006 seconds reduction, 16,135.380 seconds
decoding and 629.800 seconds external verification. Publication CPU was separate:
485.145, 90.170, 199.986 and 3.689 seconds respectively. The full planned cut/seed
scope, stage barrier, archived output references and nonoverlapping worker
ownership passed independent checks. Raw results are in
[seeds1-summary.json](seeds1-summary.json) and
[seeds1-directions.json.gz](seeds1-directions.json.gz).

## Ten-seed result and decoder storage validation

Matching and reduction completed all 3,702 directions. External reference
verification subsequently completed all 3,702 checks, recovering **1,847/1,851
reactions (99.78%)**, with 3,687 recovered and 15 unrecovered directions and no
unknowns. Existing scientific rows remained hash-identical through the separate
verification pass. The ten-seed decoder stopped with 259 successful outputs,
nine memory-limit failures and 3,434 directions not attempted. These failures
are resource exhaustion, not mapping absence; this is not a completed full run.
The audit is retained in [seeds10-summary.json](seeds10-summary.json) and
[seeds10-directions.json.gz](seeds10-directions.json.gz).

Replaying case 1654 reverse under the same 2 GiB guard reproduced the failure.
Instrumentation found 1,500,577 retained bond-term bindings: branch-specific
action positions caused copies of otherwise shared fragment terms. Saving
results also copied those expanded records. Revision `ac4583c` replaces those
copies with one shared term sequence and one action-position map per branch
use, preserving the event model and mapping permissions. The persisted record
is version 13, with an explicit group binding. All **1,295 local tests passed
in 117.96 seconds**, including independent event-choice checks, serialization
round trips and a test forbidding bond-term copies during export. Case 1654
reverse then completed in 37.815 compute CPU seconds plus 1.571 publication CPU
seconds, with a process-lifetime peak of 1.613 GiB. Across all nine prior failures,
four completed and five still exceeded the guard. This was a partial repair.

Revision `6fdb42e` removes a reference cycle retaining decoder working state and
releases bond-construction indexes before output construction. Its complete
local suite passed **1,297 tests in 111.46 seconds**. Six of the nine original
failures then completed; 862 in both directions and 1264 forward still exceeded
the guard. Revision `d34ffc7` additionally releases derived responses belonging
only to contexts excluded by the completed event window, preserving selected
responses and all saved starting-map solutions. Its **94 affected tests passed
in 1.35 seconds**, including a test that discards an excluded context while
reusing the selected response and preserving later explicit inspection.
Seven of the nine original failures completed after this change; 862 reverse
failed during computation and 1264 forward during publication. Releasing the
completed worker session in `2aa91ee` did not resolve the latter. These remained
partial repairs, not evidence of a completed full decoding run. Original failures, artifacts and consumed budgets are
preserved; these are separate new-revision diagnostic passes at the same guard.


### Shared core response validation

A bounded profile of 862 reverse found that 1,104,872 of 1,108,858 contextual
bond-term uses (99.64%) belonged to contexts already represented by shared core
fragment costs. This is a sampled intermediate state, not a complete decode.
Revision `a5f1f65` removes those duplicate detailed lookups from normal assembly
and persistence. It saves the selected fragment responses and branch references
directly in schema 14. Explicit bond inspection builds details only when
requested; genuinely joint contexts retain the terms required by their ordered
cost. Mapping permissions, event thresholds, sweep coupling, the event model
and the five-second fragment allowance are unchanged.

All **1,300 local tests passed in 118.53 seconds**, and **109 Linux decoder tests
passed in 6.05 seconds**. Independent scalar event checks, ordered placement
and sweep cases, save/load queries and a test forbidding duplicate bond-table
construction during shared-response decoding passed. All nine former memory failures completed at the unchanged **2 GiB** guard,
including output publication, with one attempt each. Their compute CPU ranged
from 10.155 to 272.652 seconds; the largest sampled process-tree peak was
1.904 GiB. The saved matching catalogues and original failed ledgers remained
unchanged. [Resource validation](decoder-storage-validation.json) preserves the
per-direction measurements and preceding unsuccessful repairs.

### Updated decoding measurement

The separate full one-seed decoding pass at `a5f1f65` completed **3,702/3,702
directions**, with zero failures or retries. All saved-catalogue references and
output checksums were audited. It used 13,817.891 compute CPU seconds and 89.859
publication CPU seconds; elapsed pool time was **501.986 seconds (8m22s)**.
Median directional decoding CPU was 0.492 seconds, p99 was 85.741 seconds and
the maximum was **322.575 seconds for 873 forward**. Sampled/lifetime peak memory
was 0.400 GiB. Compared with the previous pass on identical saved catalogues,
compute CPU decreased 14.4%; this comparison isolates the subsequent shared-response
cleanup, not the earlier native port or atom-set reduction change.
[Audited direction and reaction timings](seeds1-decoding-timings.json.gz) retain
full scope and separate publication measurements.

The ten-seed pass stopped admission after its first compute timeout, then let
already admitted work drain. It finished with **732 successes and six 600-second
compute timeouts**: 817 reverse, 869 reverse, 411 forward, 873 forward, 902 reverse
and 545 forward. The other 2,964 directions were not attempted. There were no
retries or restored budgets. These six directions previously completed in
22.877, 29.123, 38.570, 121.074, 27.658 and 24.172 CPU seconds, respectively,
in the `b51661d` ten-seed study. The measured regression crosses both a completed-
fragment representation change and the native port; those causes require
separate controlled comparisons.

Both new decoding passes used 32 workers, the unchanged 2 GiB guard and minimum
+ 5 window. They reuse audited matching catalogues; matching, reduction and
verification retain their original measured stages. These are staged measurements,
not a new uninterrupted end-to-end wall time. No complete ten-seed speedup or
regression-free release is claimed.

## Controlled regression diagnosis

The earlier fast ten-seed results used `b51661d`, before atom-set reduction
(`b26e2ee`) and the C++ port (`040f817`). Comparing those results with the new
benchmark therefore cannot isolate the native port. Controlled tests now locate
one concrete failure mechanism in the completed-fragment representation change.

On **256 small cases**, the native and immediately preceding Python assignment
search produced identical sequences of improving mappings. Both agreed with
independently enumerated minima. The 173 affected optimizer, native group,
fragment-response and decoder-window tests passed in 10.09 seconds. These checks
do not establish equality of every hardware-dependent, time-limited result.

Four representative fragments from 817 reverse, 869 reverse, 411 forward and
873 forward were then tested on the same compute host, with the same source and
target weights, initial mapping and anchors. The historical catalogue was read
using its producing revision's reader; no schema was silently reinterpreted.
The current and old Python optimizers were compared on identical fragment inputs.
No production source, scientific ledger or failed-stage budget was changed.

| Case / direction | Atoms | Native, original constraints (s CPU) | Native, broader domains only (s CPU) | Broader domains and no growth edges (s wall) |
|---|---:|---:|---:|---:|
| 817 reverse | 57 | 0.0119 | 0.0139 | 5.0007 |
| 869 reverse | 42 | 0.0073 | 0.0087 | 5.0005 |
| 411 forward | 54 | 0.0091 | 0.0126 | 5.0007 |
| 873 forward | 70 | 0.0156 | 0.0232 | 5.0011 |

Removing growth edges alone also remained fast in the separate control. Removing
those constraints **and** replacing prior per-atom partner sets with every
same-element target reproduced the expensive search. In the 54-atom example,
prior partner sets contained one to three targets; most new sets contain 26 or
27 targets. Source/target WBOs, saved witnesses and anchors were unchanged.
The historical Python optimizer also uses its five-second allowance on these
new broad inputs, while the native optimizer is faster on the original inputs.
Thus this demonstrated mechanism comes from the earlier semantic broadening,
not slower C++ execution. Its whole-dataset contribution has not been separately
measured. Initial fragment selection was targeted, not a random speed sample.

The broadened search can change the selected response: one 869 reverse fragment
returned six events instead of ten. This is not evidence of an accuracy loss,
but it disproves treating the representation change as output-identical cleanup.
Deduplication identity and starting-map permissions are separate decisions.

The **minimum + 5** setting is unchanged in both benchmark configurations.
Starting-map preparation precedes the final event-window selection, so the
requested window does not by itself prevent this extra search. The uncommitted
experimental residual bound was removed from the active checkout and its native
module rebuilt. No replacement search policy has been implemented or deployed.
The regression is diagnosed; the full ten-seed study remains stopped pending a
repair that respects the agreed scientific contract. Full controls, work counts,
incumbents and limits are in [logic-validation.json](logic-validation.json).

### Starting-map eligibility

The clarified contract permits extra starting-map search only when growth
connectivity differs from the completed fragment on at least one side and the
current witness has not already attained an established minimum. Permitted local
symmetry rearrangements remain valid; growth constraints do not veto their outputs.
The connectivity gate is absent from `optimize_fragment`: after symmetry
optimization it launches starting-map search whenever the incumbent exceeds its
current lower bound. This is not equivalent to the clarified rule.

A direct six-carbon-ring regression example uses identical complete binary ring
connectivity on both endpoints, with edge weights `(2,2,1,1,1,1)` and
`(2,1,2,1,1,1)`. Growth includes all six ring edges. Shared symmetry optimization
attains two events; independent enumeration of all 720 atom assignments confirms
that minimum. The growth record stops without another starting search. Its
completed atom-set representation invokes starting search, visits **77 partial
assignment nodes**, finds no improvement and reports a loose lower bound of one.
This small case proves unnecessary search on an ordinary fragment; its sub-
millisecond native search is not itself a large-runtime benchmark.

At the diagnostic revision, the restriction was recorded in `design.md` and the
implementation remained pending. The min + 5 window was unchanged.

The reducer contract is now explicit: atom-set deduplication groups existing
fragment records without changing their scientific payloads. Different original
information remains available through shared references under the same identity.
The then-current `Fragment.from_atom_sets` reconstruction violated this
preservation rule by replacing discovery domains and omitting accepted growth
edges. The following repair removes that reconstruction.

### Original-fragment preservation repair

This measured approach at `9b6ee47` was superseded by the representative-selection
correction below. Keeping every original record violated the intended one-record
rule even though those records shared an atom-set branch.

Reduction now retains the original immutable fragment records and their ordered
local/placement actions. Final atom-set compositions remain the branch key;
different original records are retained as alternatives within that branch.
The native history collector shares identical decision programs and parent links
without expanding duplicate paths. The replacement-fragment constructor and
the placement-only projection were deleted. Catalogue schema v6 distinguishes
these records from the prior reconstructed v5 catalogues: original matching
graphs must be reduced again to recover their information.

The first saved-case test exposed repeated validation of shared minimum choices.
The decoder now prepares each shared local minimum operation once and reuses it
when composing branches, preserving choice agreement and action order. A separate
indexing correction makes the symmetry-coverage check compare molecular bond
indices in the same frame. A six-ring fixture using atom numbers 6–11 now stops
after the same two-event symmetry solution as the fixture using numbers 0–5.
Neither invokes extra starting-map search. The matcher and event model are
unchanged; all measurements retain **minimum + 5**.

Validation: **1,301 tests passed** in the full suite before the final indexing
correction; **292 affected tests passed** afterward. Sixty-four tiny independent
mapping oracles now require exact equality across reduction, rather than merely
allowing the reduced relation to contain extra mappings. Tests also retain
different growth constraints under one atom-set branch, preserve correlated
placements and round trips, and reduce a million duplicate histories without
unfolding their paths. Assertions that previously expected growth edges to be
erased were replaced with preservation checks.

The [saved-input measurements](reducer-preservation-validation.json) cover two
Golden directions at one and ten seeds, using the same saved matching inputs and
a 2 GiB guard on the Mac. All four retain their branch-composition hashes and
recover their reference both in the catalogue and in the selected decoded window.
They do **not** establish full Golden accuracy for this revision.

| Direction / seeds | Recorded native decode CPU (s) | First preservation repair (s) |
|---|---:|---:|
| 1640 forward / 1 | 3.008 | 2.540 |
| 1640 forward / 10 | 15.673 | 21.508 |
| 1763 reverse / 1 | 3.046 | 3.217 |
| 1763 reverse / 10 | 4.864 | 23.046 |

Further benchmark expansion stopped after the regression was observed. Shared
minimum operations and the indexing correction subsequently reduced **1763
reverse / ten seeds to 16.539 CPU seconds** for decoding, **2.618 seconds** for
reduction and **4.660 seconds** for export; sampled process-tree peak was
**1.067 GiB**. The original 4,532 atom-set branches now retain 33,928 original
action/witness records, versus 4,893 reconstructed records in the earlier native
measurement. Different records are not necessarily different chemical mappings;
this repair preserves their information rather than assuming equivalence.

The preservation defect is repaired, but performance is **not uniformly improved**.
The full benchmark remains stopped. Recorded pre-port times and every measured
stage are retained in the evidence, separately from the earlier native times.
These are single observations, with correctness tests concurrent with some runs;
do not treat timing differences as precise speed ratios. No failed campaign
ledger was reset, no matching was repeated, and no cluster retry was launched.

### One original fragment per atom-set pair

Reduction now retains one original record per source/target atom-set pair. Its
witness, domains, accepted growth edges, anchors and binary symmetry stay intact;
later records on the same sets are duplicates. Branches reference that retained
record and use its starting witness. The existing symmetry optimization and
bounded starting-position search handle alignment. The matcher, event core,
five-second fragment allowance and **minimum + 5** window are unchanged.

Native history collection groups decisions by atom sets and placement context,
while retaining one original decision from the supplied reachable scope. Growth
provenance remains compact. Shared action conversion is reused within one input
graph; duplicate witnesses do not create additional decoding work. Catalogue
schema v7 requires re-reduction of the saved matching graphs.

All **1,307 tests passed in 123.28 seconds**. Critical checks cover differing
witnesses and constraints on the same sets, the retained witness's use across
branches, target-frame relabelling, selected-path scope, persistence and sweep
coupling. The triangle-tail fixture discards a second, better growth witness and
still discovers its zero-event alignment through the existing bounded search.
A million histories carrying different discovery records reduce without path
expansion. The 64 tiny mapping oracles retain valid mappings and unchanged atom-set
branches: 60 retain the same mapping union and four retain a smaller union.
This is the declared one-representative policy, not preservation of every raw
discovery history; it must not be presented as exhaustive coverage.

The [local saved-input measurements](unique-fragment-validation.json) cover two
directions at one and ten seeds. All four retain their branch-composition hashes,
catalogue reference recovery and exact annotated recovery within the decoded
window. Matching was reused; no previous implementation was rerun.

| Direction / seeds | Recorded native reduction CPU (s) | Current reduction CPU (s) | Recorded native decoding CPU (s) | Current decoding CPU (s) |
|---|---:|---:|---:|---:|
| 1640 forward / 1 | 0.049 | 0.015 | 3.008 | 1.392 |
| 1640 forward / 10 | 0.214 | 0.114 | 15.673 | 7.173 |
| 1763 reverse / 1 | 0.060 | 0.127 | 3.046 | 2.461 |
| 1763 reverse / 10 | 0.251 | 0.583 | 4.864 | 3.596 |

For 1763 reverse at ten seeds, the 4,532 atom-set branches now reference **731
unique fragments**, rather than 1,099 records. The 5,032 ordered placement contexts
replace 33,928 witness/history records. Decoding decreases from the preceding
repair's 16.539 to **3.596 CPU seconds**, with sampled process-tree peak **0.447
GiB** across reduction, decoding and publication. Both current stages together
are faster than the recorded native measurement in all four configurations;
reduction alone remains slower for 1763. The measurements are single observations,
using one numerical thread and the current 2 GiB guard, and do not establish
full-dataset speed or accuracy. Full one-seed and ten-seed validation remains
pending; no benchmark was restarted and no original ledger or budget was changed.

### Full corrected-reducer validation

Revision `7b24850` runs fresh reduction, decoding and reference verification
from all saved matching results at one and ten seeds. The scientific repair is
`7648431`; the later commit adds only benchmark orchestration. There are 1,851
reactions / 3,702 directions per arm, tolerance 2, branch cap 100, the full cut
sweep and minimum + 5 decoding. Each pool retains the 2 GiB worker guard,
600-second direction budget and five-second fragment allowance.

The maintained runner now supports `--postprocess-only`: transport conversion,
reduction, decoding and verification run in separate measured passes, with a
barrier between passes. It stops further admission on a failed task. Old ledgers,
failures and consumed budgets remain unchanged; these are separate measurements
of the corrected revision with explicitly recorded fresh stage budgets. No
matching is rerun. Local orchestration checks passed 42 tests; Linux affected
checks passed 151 tests in 94.11 seconds.

Both full runs started with 32 registered workers each. The first running
observation recorded all 3,702 one-seed transport tasks complete and reduction
starting; ten-seed transport had 541 completions and 32 simultaneous running
tasks. No task failures were observed. These are partial results.
Reaction-level recovery will be compared against the recorded 1,841/1,851 and
1,847/1,851 results. Saved matching time and fresh postprocessing time are
separate measurements, not an uninterrupted end-to-end wall time.
