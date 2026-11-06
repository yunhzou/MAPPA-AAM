# Golden: fresh one-seed and ten-seed measurements

## Current accepted version and verification audit

The accepted reducer at `7648431` retains one original fragment per atom-set
pair, including its discovery constraints. Its one-seed postprocessing pass
completed all 3,702 reductions, decodes and catalogue checks: **1,834 / 1,851
reactions recovered in the catalogue**, with 17 unrecovered and no unknowns.
This is **catalogue recovery, not post-decoding accuracy**. The ten-seed pass has
now completed all 3,702 reductions and decodes. Catalogue recovery is
**1,843 / 1,851 reactions**; one reverse-direction catalogue query timed out,
but that reaction has a positive forward result. The independent decoded audit
is reported separately below.
See [one-seed evidence](unique-seeds1-summary.json) and
[completed ten-seed stage evidence](unique-seeds10-summary.json).

Latest one-seed decoded verification confirms **1,832 / 1,851 reactions (98.974%)**:
16 are absent from the saved exact-index decoded response and three remain
unresolved: **852, 1568 and 1606**. The targeted retest below recovered 607 and
1586 with replayed witnesses. This remains a confirmed lower bound on coverage.

### Completed ten-seed postprocessing and decoded audit

All **1,851 reactions / 3,702 directions** have successful matching, reduction
and decoding outputs. This continuation reused every saved matching result,
all 3,702 transport snapshots and 433 completed reductions. Those completed rows,
the original plan and consumed budgets remained unchanged. No matching was
repeated, and no scientific retry allowance was added.

Settings remain discovery tolerance **2**, branch cap **100**, ten seeds, random
seed 42, both directions and the complete single-edge sweep. Decoding retains
the **minimum + 5** window, original weights, event thresholds 0.5/0.3 and the
five-second fragment improvement allowance. Every worker retained the **2 GiB**
guard. The matching revision is `27b18ab`; reduction and decoding resumed their
frozen `7b24850` revision, containing the accepted reducer at `7648431`. Later
verifier additions do not change these mapping operations.

| Evaluation | Confirmed reactions | Absent | Unresolved |
|---|---:|---:|---:|
| One seed: saved decoded responses | 1,832 / 1,851 (98.974%) | 16 | 3 |
| Ten seeds: saved decoded responses | **1,841 / 1,851 (99.460%)** | **8** | **2** |
| Ten seeds: catalogue membership | 1,843 / 1,851 (99.568%) | 8 | 0 |

Decoded recovery requires exact annotated heavy-atom indices in the retained
response and window; catalogue recovery allows endpoint chemical symmetry.
These are distinct metrics. Ten-seed decoded unknowns are **590 and 607**.
The eight absent reactions are **7, 19, 986, 1285, 1358, 1377, 1475 and 1553**.
Unknowns remain inconclusive, not failures of mapping or proven absence.

Ten reactions have newly confirmed decoded coverage: **20, 44, 852, 867, 964,
1034, 1228, 1314, 1568 and 1606**. There are no proven reaction-level losses;
**607**, previously confirmed with one seed, remains unknown in this ten-seed
audit. Thus the confirmed total increases by nine. A bounded membership query
cannot establish monotonic coverage across different decoded catalogues.

The latest verifier at `1591a8c` checked all 3,702 directions, then gave nine
unknown directions in five unresolved reactions an explicit 50-second follow-up,
matching the one-seed targeted allowance. It recovered 575 and 852. Cumulative
direction outcomes are **3,651 recovered, 30 absent and 21 unknown**. All 3,651
positive proofs passed replay, injectivity, reference and independent event-count
checks. There were no audit worker failures or invalid certificates. In 32
positive directions the accepted approximate boundary model differs from the
physical count; every returned physical witness still lies inside its saved
numerical window. This is not a claim of exact physical event bins or global
physical minima.

One separate **catalogue** check, **1351 reverse**, exhausted its original
600-second compute allowance, with about 0.63 GiB sampled peak memory. Its forward
catalogue result already recovers the reaction. The timeout and its consumed
budget remain recorded: 3,701 successful catalogue checks plus one resource
failure, not 3,702 successful checks. The independent decoded audit now treats a
catalogue certificate only as an optional search hint, so a failed catalogue
check cannot block evaluation of a successful saved decode. This benchmark-runner
correction changes neither mapping/decoding nor the decoded acceptance test.
Validation passed **42 affected tests and 26 smoke tests**.

The stage audit checked all 14,808 terminal records and all 14,807 published result
checksums, with zero retries and an empty scientific work queue. All nine earlier
ten-seed decoding failures completed successfully in this corrected pass. Detailed
direction timings and replayable proofs remain with the private campaign;
public evidence contains aggregate measurements and source hashes.

### One-seed versus ten-seed timing

Both studies used 32 compute workers. Matching is the previously measured,
reused pass; reduction and decoding are the completed corrected passes. CPU time
sums compute across workers; active elapsed time measures the parallel stage,
including publication. The paused interval between the two ten-seed reduction
allocations is excluded. These are stage measurements across allocations, not
a newly measured uninterrupted end-to-end run.

| Stage | One-seed compute CPU (s) | Ten-seed compute CPU (s) | One-seed active elapsed (s) | Ten-seed active elapsed (s) |
|---|---:|---:|---:|---:|
| Matching | 1,427.88 | 13,450.55 | 401.67 | 957.41 |
| Reduction | 1,001.06 | 8,094.12 | 56.58 | 307.41 |
| Decoding | 3,474.18 | 17,180.81 | 124.95 | **610.12** |

Ten seeds require **9.42× matching CPU, 8.09× reduction CPU and 4.95× decoding
CPU**, compared with one seed. This compares different amounts of retained work;
it is not a fixed-workload acceleration measurement. Matching cache hit fractions
were **89.94%** and **92.41%** respectively. These full-run ratios do not isolate
the cache's causal speedup; selected cache-disabled controls appear below.

| Per-direction decoding compute time | One seed (s) | Ten seeds (s) |
|---|---:|---:|
| Median | 0.476 | **1.964** |
| 90th percentile | 2.093 | 10.428 |
| 95th percentile | 3.299 | **17.502** |
| 99th percentile | 7.477 | 44.601 |
| Maximum | 52.693 | **316.459** |

The ten-seed maximum is **607 forward**. Sampled decoding peak memory was
**1.831 GiB**, with no memory-limit failures. Publication is separate from compute:
ten-seed reduction used 1,100.26 publication CPU seconds / 1,299.72 summed worker
wall seconds; decoding used 683.88 publication CPU seconds / 1,716.41 summed
worker wall seconds. These worker sums are not additional parallel stage durations.
Snapshot transport used 190.27 elapsed seconds and is not matching or reduction.

External catalogue verification took **636.64 elapsed seconds**, including the
timeout: successful tasks used 3,525.57 compute CPU seconds, and the timed-out
task contributed another 599.70 sampled CPU seconds. The latest full decoded
audit took **218.09 elapsed / 6,265.76 compute CPU seconds**; the targeted
nine-direction follow-up used **66.44 elapsed / 455.56 compute CPU seconds**.
These external evaluation costs are excluded from mapping and decoding timings.
Query deadlines are soft; task loading and native operations can extend elapsed
time beyond the query allowance, within the separate worker guard.
The resumed allocation used 1,453 elapsed seconds / 27,968 total CPU seconds and
ended with a failure exit status because of the retained catalogue timeout.
The full and targeted decoded-audit allocations ended successfully, using
239 seconds / 6,373 CPU seconds and 87 seconds / 461.65 CPU seconds respectively.
Allocation accounting includes startup, supervision and shutdown; it is separate
from the stage measurements above.

Evidence: [completed stage and timing summary](unique-seeds10-summary.json),
[decoded coverage and proof-audit summary](decoded-seeds10-summary.json).
The earlier [433-reduction snapshot](unique-seeds10-partial-summary.json) remains
historical evidence, not the current completion state.

### Historical whole-growth cache audit (`4b951a4`)

The now-retired completed-fragment cache skipped both growth and completion on a valid hit.
It persists across seeds within a cut and across cuts while the reaction session
remains resident. Shared snapshots exchange knowledge at task boundaries. Keys
include the starting atom, settings and occupied target atoms; dependency checks
also validate relevant locks, source topology, joined fragments and deferred
edges. The same final atom sets reached under another request do not alone
establish a reusable answer. Branch traversal, continuation construction, caps
and export still execute after a hit. Snapshot saving belongs to publication,
not the matching CPU reported above.

At revision `4b951a4`, fresh local controls on cases **895, 1261 and 1841**, both
directions, compared one and ten seeds with 64 MiB retention versus zero
retention. All **24 direction runs** completed and all **12 paired comparisons**
agreed on catalogue contents, cut/seed scope and cap evidence. No scientific code
changed. Each arm used one compute worker on an Apple M4 Pro, with a 2 GiB guard.

| Seeds | Cache-on matching CPU (s) | Cache-off matching CPU (s) | Speedup |
|---|---:|---:|---:|
| 1 | 1.039 | 1.796 | 1.73× |
| 10 | 12.477 | 20.804 | 1.67× |

These are sums across six directions, with one observation per arm. They exclude
preparation, reduction, verification, audits and file publication. The control
disables retention, preserving the native algorithm and dependency accounting.
It is not an alternative cache-free algorithm. At ten seeds, **89.26% of
requests hit**, but those hits reused only **47.84% of logical extension steps**.
Request hits therefore overstate the fraction of growth work avoided; extensions
themselves also vary in cost. Case 1841 produced all 67,583 evictions; the other
two cases had none.

An additional direct test repeated the identical initial fragment request 200
times in each of three alternating-order measurements. Warm hits took median
**8.62–13.79 microseconds**, versus **48.49–88.39 microseconds** with zero
retention: **5.63–8.36× faster**, with equal scientific outputs. These direct
calls include Python result conversion; the native branch search consumes hits
without that conversion. They demonstrate cheap hits, not full-search speedup.

In the full ten-seed campaign, per-seed matching CPU summed over all cuts drops
from **1,425.21 s for the first seed to 1,212.61 s for the tenth**, about 15%.
Total requests still grow from 74.28 million at one seed to 752.05 million at ten;
misses grow from 7.48 million to 57.07 million. Reuse is active, but does not remove
the repeated branch work or new conditional fragment requests.

The earlier [three-repetition cache study](../golden_seed10_tol2_cap100_20261105/README.md#completed-cache-control-results)
measured a 1.76× ten-seed matching speedup. Its native growth, completion, cache
and seed scheduler sources are unchanged in this audit. Its capacity controls
found only about 2% aggregate benefit from 64 to 256 MiB on eight cases, and
4.5–6.5% from 64 MiB to 1 GiB on two heavier cases. Those remain selected historical
measurements, not proof that capacity never matters. Current evidence supports
useful absolute cache savings for that revision, not strongly sublinear full-dataset seed scaling.
[Aggregate audit evidence](cache-review.json) records counters and timing scope.

### Replacement: unique fragments answer extension queries

The replacement removes whole-growth replay and its dependency tracking. Matching
accumulates one fragment per source/target atom-set pair. A stored witness can
prove a smaller extension feasible after checking its required edges, tolerance,
locks and available targets. Otherwise the ordinary existence solver runs.
Growth proposals always execute afresh. Cut/trial random streams are distinct
and reproducible; saturation chooses its representatives independently of the
existence witness. The saved search graph supplies its fragment inventory to
reduction. Decoder logic is unchanged.

Selected local controls used cases **895, 1261 and 1841**, both directions, one
and ten seeds, tolerance 2, cap 100 and the full cut sweep on an Apple M4 Pro.
All **24 direction runs** completed. All **12 cache-on/off pairs** had identical
catalogue fingerprints, cut/seed scope, graph counts and cap outcomes. Reference
membership recovered all three reactions in every setting. These controls do
not perform event decoding or establish full-Golden accuracy.

| Seeds | Pool-on matching CPU (s) | Pool-off matching CPU (s) | Speedup | Pool-on reduction CPU (s) |
|---|---:|---:|---:|---:|
| 1 | 0.706 | 1.040 | 1.47× | 0.317 |
| 10 | 7.505 | 11.152 | 1.49× | 2.903 |

Times are sums across six directions, one measurement per arm, excluding
preparation, audits, verification and publication. At ten seeds the six pools
retained **6,116 unique fragments** and supplied **26,415 positive existence
answers**. The largest pool accounted for **9.52 MiB**, with no evictions. These
queries follow the cheaper direct extension of the current witness, so their hit
fraction is not the fraction of all extensions accelerated.

The historical cache-on totals above were 1.039 s and 12.477 s. Changed random
streams and representative ordering change the sampled search workload; those
before/after numbers do not isolate a cache speedup. The replacement's on/off
comparison does. This sample still does not show strongly sublinear seed scaling.
The measured worktree identity, per-stage CPU/wall times and all comparison
outcomes are retained in [the aggregate evidence](cache-review.json).

Final validation: **1,313 tests passed in 131.96 s**, including matching oracles,
Python/native agreement, reduction/decoding, trajectory replay and work-pool
storage. An additional eviction/snapshot-budget test passed separately. The suite
caught and verified the correction of an explicit-path-subset reduction bug:
excluded paths cannot supply the retained fragment. The final source identity
and the timing run's earlier identity are recorded separately in the evidence.

### Full-pipeline validation of `465ee34`

Fresh local runs of the same three cases used **two workers**, the **2 GiB**
guard, tolerance 2, cap 100, both directions, every cut and **minimum + 5**
decoding. One and ten seeds each completed all six directions in all four
stages. The independent saved-decoder audit recovered **3/3 reactions and 6/6
directions in each arm**, replaying concrete witnesses with zero certificate
failures or unknowns.

| Measurement | One seed | Ten seeds |
|---|---:|---:|
| Completed cut searches | 270 | 270 |
| Completed seed/cut searches | 270 | 2,700 |
| Matching worker CPU (s) | 0.802 | 7.482 |
| Reduction worker CPU (s) | 0.423 | 4.886 |
| Decoding worker CPU (s) | 1.355 | 8.024 |
| Catalogue verification worker CPU (s) | 0.331 | 1.806 |
| Complete pipeline wall time (s) | 3.476 | 15.455 |
| Separate decoded-reference audit wall time (s) | 0.593 | 2.939 |
| Slowest direction's decoding CPU (s) | 0.671 | 2.627 |

CPU values sum work across six directions; they are not parallel elapsed time.
Worker stages include artifact loading, and decoding includes record preparation.
Publication is separate. The aggregate evidence retains publication costs,
per-direction measurements and narrower matching-only CPU (0.682/6.981 s).
These are single observations, not a claim of a full-dataset speedup.

All **12 catalogues**, matching scopes, graph counts and cap outcomes matched
the earlier sequential, zero-retention controls exactly. Both workers registered;
two tasks actually overlapped in every stage. Audited timestamps preserved each
global stage barrier. There were **no failures or retries**.

Shared snapshots were actually loaded: 3,287 bytes in the one-seed arm and
79,621 bytes at ten seeds. The ten-seed workers reported 26,244 fragment hits;
combined snapshot load/save time was 0.043 s. Storage retained exactly **six
cache files for three reactions** in each arm. Maximum accounted resident pool
size was 8.72 MiB at ten seeds; process lifetime peak RSS stayed below 407 MiB.
The complete evidence is under `pipeline_validation` in
[the cache report](cache-review.json). No scientific code changed for this test.

Earlier catalogue counts from the broadened-fragment reducer must not serve as
accuracy evidence for the current mapping contract. Case 19 demonstrates why:
the earlier successful fragment had no accepted-edge constraints. Its accepted
reference mapped three required fragment bonds onto nonbonds; the fragment's
binary symmetry could not repair those correspondences. The seven reaction-level
differences from that historical count are not automatically seven regressions.
The remaining cases require individual evidence.

The current verification audit queries only saved decoded starting orbits and
their retained event window. Every positive answer must include a concrete
mapping, replayable local and placement choices, and an independently counted
model event total. Physical whole-reaction events are reported separately.
Exact annotated-index recovery and endpoint-symmetry-normalized catalogue
recovery are different policies and must remain labelled. Resource limits are
unknown outcomes. Original outputs and consumed budgets remain unchanged.

### Initial audit of saved one-seed decoding

All **3,702 directions / 1,851 reactions** were queried at the unchanged scientific
identity. Recovery here requires the **exact annotated heavy-atom indices**;
unannotated hydrogens must admit a complete injective witness. This is stricter
than the catalogue evaluation's endpoint-equivalence policy.

| Outcome | Reactions, either direction | Directions |
|---|---:|---:|
| Verified in the saved decoded window | **1,830 (98.865%)** | **3,591** |
| Exact-index reference absent from saved decoded response | 16 | 80 |
| Inconclusive within the query budget | 5 | 31 |

At the end of these initial passes, the five unresolved reactions were
**607, 852, 1568, 1586 and 1606**. They were neither recovered nor demonstrated
absent. That verified percentage was a lower bound, not an estimate assigning
outcomes to unknowns. Ten-seed post-decoding accuracy had not yet been measured
at that initial one-seed audit; the completed ten-seed evidence appears above.

Every positive query retained a concrete full mapping and its ordered local and
placement choices in the private aggregate evidence bundle. Validation replayed
those choices through the original fragment constraints, checked the saved
decoder starting correspondence, injectivity, elements, anchors and exact
reference coverage, then independently counted events from the original matrices.
All **3,591 positive certificates passed**. All positive decoded directions also
had positive catalogue checks. Input and result hashes remained unchanged; no
matching, reduction, optimization or decoding was rerun by the audit.

There are **90 positive directions** where the approved ordinary-boundary model
overcounts the concrete witness: 78 by one event, seven by two, one by three and
four by four. These are explicitly approximate bins. All returned witnesses'
physical event counts still lie within their saved numerical window. This does
not establish a physical global minimum or exact physical membership of every
reported bin.

Case **1654 reverse** supplies a validated reference witness at three model and
physical events, in its saved 3–8 window. Case **19** is absent in both directions.
Case **1553** is the only catalogue-recovered reaction with a definite negative
exact-index decoded result. Its saved reverse catalogue certificate uses
endpoint-equivalent reference labels and family 147, which is outside all retained
decoded contexts. That certificate's full witness has 20 physical events, versus
the saved 8–13 window. This establishes that the catalogue certificate cannot
justify decoded recovery; it does not isolate reference equivalence, chosen
starting orbit and windowing as independent causes of the difference.

The first external audit used 32 workers and a ten-second query allowance per
direction: **71.766 seconds elapsed / 1,816.574 CPU seconds**, with 49 inconclusive
directions. An explicitly separate pass queried only those 49, prioritizing the
saved catalogue family as a search hint and still requiring the complete decoded
proof. Its 16 workers used another ten-second query allowance: **29.718 seconds
elapsed / 366.339 CPU seconds**, verifying 18 more directions. Neither pass had
worker errors or memory failures; both retained the 2 GiB guard and 60-second
external task guard. Query deadlines are soft, and preparation/native operations
can exceed them; actual task measurements are preserved. No scientific retry
budget changed. The earlier 24-direction pilot took 8.001 seconds and is separate
from these totals.

Validation passed **155 affected tests and 26 smoke tests**; the final ordering and
continuation changes passed their 12 targeted tests. Independent tiny-graph
oracles check all 24 assignments for six ring/anchor/tail/window combinations,
plus joint placement, unmatched atoms, forward/reverse queries, sweep handling,
and deliberately forged certificates. The saved-file worker test confirms that
verification preserves source artifacts and retains replayable proofs.

Evidence: [audit summary](decoded-verification-summary.json),
[direction outcomes and certificate checksums](decoded-verification-directions.json.gz).
Operational configuration and molecular certificates remain private. The reusable
audit is documented in [bench/README.md](../../bench/README.md).

### Reference-conditioned event filter

The external verifier at `b6d8c08` counts events already forced by the reference
before constructing a symmetry query. It rejects only contexts whose lower
bound exceeds the saved window, preserves the ordinary-boundary model and sweep
correction, and reuses the boundary preparation for surviving queries. It does
not optimize a minimum or change matching, reduction, decoding or acceptance.
The event-filter helper changes the source fingerprint; the producing scientific
identity and every saved input/output hash remain unchanged.

The final full audit completed **3,702 directions in 63.515 seconds**, using
**1,555.881 compute CPU seconds**, with 32 workers and the original ten-second
query / 2 GiB worker limits. The entire allocation, including startup and shutdown,
took 74 seconds and 1,619.479 CPU seconds. No worker or certificate failures
occurred. The filter rejected **2,707 of 21,991 examined contexts** before their
symmetry queries. Filter-stage time totalled 28.924 seconds across directions;
this includes the boundary preparation subsequently reused by surviving queries.

Reaction results remain **1,830 recovered, 16 absent and 5 unknown**. This pass
returned 3,588 positive directions, 80 absent and 34 unknown. Three previously
positive directions—737 reverse, 823 forward and 1777 reverse—reached the query
limit. Work stopped for investigation: their existing certificates were replayed
with the current code and independently recounted. All three remain valid and
pass the filter; none became an absence result. Their forced bounds are 11, 17
and 4 against window maxima 15, 21 and 8. The new pass retains these three unknown
query outcomes rather than silently substituting the earlier certificates.
This establishes preserved witnesses, not identical time-bounded query behavior.

The earlier full pass took 71.766 seconds / 1,816.574 CPU seconds, but that
comparison also includes the previously added catalogue-family ordering hint.
For the same 49 difficult directions already tested with that hint, CPU time
changed only from 366.339 to 360.348 seconds: this filter does not resolve the
hardest symmetry queries. A useful negative example is **1228 reverse**: it
still reports absence, skips 209 of 294 contexts, and changes from 1.331 to
0.464 CPU seconds. Case 741 reverse skips 252 of 602 contexts and changes from
1.939 to 1.411 CPU seconds, also retaining absence. These are single-run timing
observations, not a controlled hardware-normalized speed guarantee.

The preliminary filter pass took 67.221 seconds / 1,638.832 CPU seconds. Inspection
identified duplicate boundary preparation, removed before the final pass.
Both measurements and their consumed resources remain recorded separately.
Validation passed **134 affected tests and 26 smoke tests**, including independent
small-graph enumeration, partial references, strict/metal thresholds, unmatched
boundaries and ordered sweep choices. All final positive certificates passed
the existing replay and scalar-count audit.

Evidence: [filter summary and replay checks](decoded-event-filter-summary.json),
[full direction measurements](decoded-event-filter-directions.json.gz).

### Targeted retest of the five unresolved reactions

The user requested another test with the current verifier. The seven unresolved
directions of those five reactions were queried at `b6d8c08`, using a separately
announced **50-second query allowance**, a 90-second task guard and the unchanged
2 GiB worker guard. The three already-negative directions were retained without
rerunning them. No matching, decoding or scientific settings changed.

| Newly recovered reaction | Direction | Compute wall seconds | Model / physical events | Saved window |
|---|---|---:|---:|---|
| 607 | Reverse | 14.592 | 16 / 16 | 16–21 |
| 1586 | Reverse | 17.688 | 11 / 11 | 7–12 |

Both positive answers include a concrete mapping and ordered choices replayed
against the saved starting correspondences, with independent scalar event
counts. All input, decoded-output and original ledger hashes remained unchanged.
Cases **852, 1568 and 1606** remain inconclusive; additional time did not settle
them. Unknowns remain distinct from absence. The query also retains its existing
10,000-node limit; a longer deadline does not remove that limit.

The pass used **54.207 seconds elapsed / 282.781 compute CPU seconds** in an
eight-worker pool, with seven tasks and no worker, memory or certificate failures.
The full allocation took 61 seconds / 288.267 CPU seconds. It returned two
positive and five unknown directions. Combining the retained valid evidence
from all passes gives **1,832 recovered, 16 absent and 3 unresolved reactions**;
direction outcomes are 3,593 recovered, 80 absent and 29 unknown. Recovery is
defined at reaction level, so an unresolved opposite direction does not undo a
validated witness. Every pass's consumed resources remain separate.

Evidence: [targeted retest outcomes, timing and certificate hashes](decoded-five-unresolved-retest.json).

## Earlier measurements of the scientific pipeline

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
