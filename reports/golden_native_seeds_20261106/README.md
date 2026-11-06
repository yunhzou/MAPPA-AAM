# Golden: fresh one-seed and ten-seed measurements

This paired full-dataset study measures revision
[`27b18ab`](https://github.com/yunhzou/MAPPA-AAM/commit/27b18abad691648cf2a87d2e5ae43db117d80b66),
including the native reduction and event kernels and the incident-bond event bound. Each arm contains all
**1,851 reactions / 3,702 directions**. Matching, reduction, decoding, external
reference verification and reaction caches start empty in both arms. Only
immutable prepared endpoint and reference inputs are reused.

**Status: both complete runs have restarted with empty results and caches;
final benchmark results remain pending.** Historical
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
