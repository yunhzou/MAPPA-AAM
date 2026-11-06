# Golden: fresh one-seed and ten-seed measurements

This paired full-dataset study measures revision
[`27b18ab`](https://github.com/yunhzou/MAPPA-AAM/commit/27b18abad691648cf2a87d2e5ae43db117d80b66),
including the native reduction and event kernels and the incident-bond event bound. Each arm contains all
**1,851 reactions / 3,702 directions**. Matching, reduction, decoding, external
reference verification and reaction caches start empty in both arms. Only
immutable prepared endpoint and reference inputs are reused.

**Status: the one-seed run is complete and independently audited. The ten-seed
run stopped admission after decoding exceeded its unchanged 2 GiB memory guard.
Its workers have drained, and a separate verification-only pass completed all
reference checks. Ten-seed decoding remains incomplete.** Historical
accuracy and selected-case speedups are not results of this new full run.

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
Cluster validation of this last change and full updated decoding measurements
remain pending. Unit tests alone do not establish the resource repair or its
full-dataset performance. Original failures, artifacts and consumed budgets are
preserved; these are separate new-revision diagnostic passes at the same guard.
