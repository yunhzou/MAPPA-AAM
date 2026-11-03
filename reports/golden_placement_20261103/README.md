# Full Golden benchmark with placement compression

A fresh full Golden rerun was submitted at **2026-11-03 06:19:04 UTC** for
[`1710f69`](https://github.com/yunhzou/MAPPA-AAM/commit/1710f690f5c49b34d11de69782aa54b1f4d79304).
The first campaign exposed a slow symmetry-membership check and was intentionally
stopped before decoding. The repaired implementation has passed **883/883 Mac
tests**, **79/79 Linux tests** and five targeted matching probes. The matching
pools registered all **141 compute workers** at **06:20:04 UTC**. All **3,702
directions have completed matching and family aggregation/deduplication**, called
reduction below. It combines each direction's
seed/cut results into one catalogue while preserving distinct relationships and
provenance. At **08:03:31 UTC**, the frozen pass has **3,655 successful decodes and 47
remaining memory-limit errors**. Reference verification has **3,696 completed
results and six timeouts** out of 3,702. All required direction-stage records are
terminal, with no active task attempts or original finish allocations remaining.
The pass has ended with unresolved failures; it is not a successfully completed
benchmark. Final audited resource and reference comparisons remain pending.

## Scope and settings

The full benchmark covers **1,851 reactions, 3,702 directions and 177,552 matching
cuts**, including the uncut condition. Each reaction runs in both directions.
The input audit found every expected reaction/direction exactly once and verified
all **5,553 unique input/reference blobs** against their checksums. The replacement
run uses fresh results and caches, reusing only immutable inputs and references.

| Setting | Value |
|---|---|
| Seed count | 1 |
| Branch cap | 2,000 representative continuations at each hierarchy level |
| Discovery / relaxed-self tolerance | 1 / 1 |
| Random seed | 42 |
| Graph / cut presence threshold | 0.2 / 0.2 |
| Ordinary / metal event threshold | 0.5 / 0.3 |
| Search scope | Both directions; uncut plus every single-edge cut |
| Input weights | Original prepared endpoint WBOs |
| Initial worker memory guard | 2 GiB for 129 workers; 4 GiB for 12 workers in the batch containing case 1561 |

The user explicitly selected cap 2,000. Placement compression retains equivalent
labelled target sets through shared joint actions; representative branch counts
are distinct from the number of mappings represented. Fragment-local relaxed
relationships remain part of the saved result. The
[integration report](../relaxed_fragments_20261102/README.md) records the approved
design and independent mapping/event checks.

Ten disjoint persistent pools request **141 compute workers** and **171 CPUs
in total**, including execution overhead. They share reaction caches.
All matching and reduction must
finish across the entire campaign before decoding begins. Reference verification
is a separate evaluation of the reduced catalogues. Requested workers, actual
readiness and measured utilization will be reported separately. Operational
configuration and raw scheduler records remain private.

| Stage | Full denominator | Status at 08:03:31 UTC |
|---|---:|---|
| Matching | 3,702 directions / 177,552 cuts | All directions complete; no errors |
| Reduction | 3,702 directions | All complete; two earlier publication failures recovered |
| Decoding | 3,702 directions | 3,655 complete; 47 memory-limit errors remain |
| External reference verification | 3,702 directions / 1,851 reactions | 3,696 complete; six compute timeouts |

Two reductions in this fresh run exceeded their 2 GiB guard during publication:
**1799 forward** reached **2.298 GiB**, and **852 reverse** reached **2.043 GiB**.
Both completed with isolated **4 GiB** retries using their original remaining
budgets. Their original task specifications, failed attempts and all other
completed results were preserved; no fresh retry time or scientific change was
introduced. Final resource totals will include both attempts for each case.
The original dependent finish submissions were replaced only after matching and
reduction were complete and their owners had released the work. New finish-pool
readiness was not inferred from submission.

The first decoding pass produced **3,639 outputs and 63 failures**. Of those
failures, **61 occurred at 2 GiB**. The other two, **1561 forward
and 830 forward**, exceeded **4 GiB during publication**. Planned recovery uses
4 GiB and 8 GiB respectively, retaining original remaining budgets and attempt
histories. Seventeen failed directions have since been retried at 4 GiB after
their original owners released the work: **16 succeeded**, while **852 forward**
exceeded 4 GiB during publication. Remaining errors and consumed budgets are
preserved. The cluster implementation remains frozen during local repair work.

Decoder-stage computation includes catalogue loading and output preparation.
Static review found extra representation cost: distinct numeric weight matrices
expand into Python lists while the original catalogue remains resident, and
shared working indexes persist across the catalogue. These costs have not been
profiled per failed case. A compute-phase memory failure therefore does not, by
itself, identify local event evaluation as the cause.

At **07:47:20 UTC**, the **3,655 successful decodes** used 91,822.40 compute CPU
seconds. Their per-direction compute CPU distribution was **4.24 s median,
70.60 s p90, 122.74 s p95, 257.29 s p99 and 538.52 s maximum**. This success-only
distribution excludes the 47 directions still carrying decode errors. Successful
worker CPU totalled 93,345.92 seconds; the 64 failed attempts had consumed another
24,806.21 sampled CPU seconds. Successful full-benchmark completion remains outstanding.

## Matching repair and local validation

The original full run exposed a bottleneck in checking whether two selected
target atom sets belong to the same certified placement relationship. A bounded
Mac probe of **case 1801, reverse, uncut** did not finish within **30 seconds**,
using **29.80 sampled CPU seconds** and peaking near **57 MiB**. A separate native
profile concentrated in recursive set-membership checks during saturation.

The repair rejects an impossible comparison earlier: every remaining symmetry
action preserves the number of selected atoms within each of its point orbits.
Unequal counts therefore prove that the two sets cannot be related. Equal counts
still require the existing exact correlated check. This necessary condition
removes redundant search without changing the retained mappings, branch cap,
tolerances or decoder.

All five isolated cold Mac probes now finish:

| Case / direction / cut | Matching CPU seconds | Matching elapsed seconds | Families |
|---|---:|---:|---:|
| 1410 / reverse / 40 | 0.200953 | 0.233380 | 1,997 |
| 1606 / reverse / 1 | 0.056183 | 0.057017 | 1,136 |
| 1799 / reverse / uncut | 0.007319 | 0.008104 | 1 |
| 1799 / reverse / 1 | 0.007511 | 0.008424 | 2 |
| 1801 / reverse / uncut | 0.006416 | 0.006879 | 1 |

For a real output-preservation check, the completed **1410 reverse cut 40** from
the original campaign was reduced and decoded using the same consumers as the
repaired probe. Both produce **1,997 families and an identical complete decoded
record**, verified by SHA-256. The saved old search took 362.10 seconds on Linux;
the repaired timings above are from the Mac, so these are not a controlled
cross-version speedup measurement. This cut records branch-cap stops; successful
execution does not establish exhaustive coverage.

The repaired implementation passed **49 focused tests in 7.06 seconds** and all
**883 Mac tests**, comprising 880 existing tests and three new regressions.
Operating-system access was required to rerun 56 process/socket tests after
sandbox permission failures; the other 827 tests had already passed. The measured
Mac source/engine identity is
`c557077902062d7a0e08940544f1f4e0d57b92ce56d9e31b1413ff40acba6f4a`.
The Linux build passed **79/79 tests in 16.80 seconds** before submission; its
pinned source/engine identity is
`62a8c461df57aaab40c69057008443210ef205cfd8d2c0665f090c9c3eed50ab`.
[Point evidence](point-evidence.json) records inputs, stage timings, memory,
cache measurements, cap stops, output equality and evidence hashes. These probes
do not measure full Golden reference recovery.

## Local reduction analysis

Two Mac diagnostics separate family aggregation from loading and publication.
They use seed 1, cap 2,000 and tolerance 1; these are local measurements, separate
from the campaign CPU totals.

| Scope | Matching CPU (s) | Reduction CPU (s) | Publication CPU (s) | Catalogue reload CPU (s) |
|---|---:|---:|---:|---:|
| 1518 forward, all 104 cuts | 5.092 | 2.623 | 0.440 | 1.170 |
| Saved 1410 reverse cut 40, median of three runs | Not rerun | 0.257 | 0.009 | 0.168 |

Case 1518 retains **31,687 paths and families across 31,676 distinct fragment
compositions**. Its 1,814 full fragment keys correspond to 1,701 atom-set pairs;
removing only the witness field leaves 1,749 keys. The outer composition grouping
already places equal final atom-set compositions in one branch. That branch can
retain several compact alternative relationships; grouping the composition does
not require proving those alternatives equivalent. Equal atom sets do not
authorize dropping their different constraints or choices.

A second 1518 probe changes only the cap under the current implementation:

| Cap | Matching CPU (s) | Reduction CPU (s) | Families | Distinct compositions |
|---|---:|---:|---:|---:|
| 100 | 0.869795 | 0.274304 | 3,553 | 3,542 |
| 2,000 | 5.092056 | 2.622537 | 31,687 | 31,676 |

The cap-100 run has 1,407 fragment records and reports 3,553 stalled and 12,150
capped stops. The higher cap expands retained families almost ninefold on this
case, while only 11 family records exceed the number of compositions. This
isolates the cap setting on one case: more retained families increase aggregation
work. It is neither a historical-version comparison nor a reference-accuracy
measurement.

The saved 1410 cut has **1,997 families with 1,997 distinct compositions**.
Its 15,976 fragment occurrences share 58 unique fragment keys. Profiling found
**29,955 prepared-action constructions for 204 retained actions**, identifying
repeated preparation despite existing storage sharing. Loading that cut costs
0.039 CPU seconds. Historical revision `06f7941` already grouped compositions
and families and reused prepared actions per transition; the current path lost
that preparation reuse, a concrete preparation regression. This supports
investigating shared preparation while preserving each family's relationship,
order and provenance.

In three paired 1518 repetitions, a temporary in-process action-reuse diagnostic
reduced median reduction CPU from **2.637507 to 2.254886 seconds (14.51%)**.
Every pair preserved the complete family mappings, edges, actions, policy,
provenance, branch membership and counts exactly. Production source remained
unchanged. This identifies a useful preparation cost, but does not account for
all of the campaign slowdown or establish a dataset-wide speedup.

## Decoder preparation repair

For the frozen **1518 forward catalogue** (104 cuts, seed 1, cap 2,000,
tolerance 1), a Mac diagnostic measured **1.646 CPU seconds for decoder
preparation and 31.089 seconds for lazy assembly/export**, totalling **32.736
seconds**. Assembly/export performs the local and boundary event work as well as
record conversion; it is not just file I/O. These timings exclude catalogue
loading, construction and file publication.

A separate profiled pass found **1,457,602 per-atom preparation calls**,
**39,529,272 set intersections**, 45,169 neighboring-pair calls and 31,687
dependency-component constructions for 31,687 exported family tables. The export
retains **22 shared numeric responses**. This identifies substantial repeated
ordered-program preparation across compositions, although their complete
contexts are not necessarily identical. Final per-program atom rows had no
reuse during that pass; caching those finished rows alone would not eliminate
the first calculation.

The local repair removes that final-row cache and the repeated whole-program
walk for every atom. Each distinct local relationship retains its prepared
reachable sets. Branch assembly binds those shared sets only to affected atoms,
preserving the original ordered joint choices. Overlapping relationships reuse
bounded intermediate sets; each relationship is applied once, without assuming
transitivity. The obsolete checkpoint relocation for the removed cache was also
deleted. Matching, family content and event rules are unchanged.

Both saved catalogues were compared in separate fresh processes, once per
implementation. Historical versions of the two changed modules were loaded in
memory from the recorded commit; the inputs and all other code were identical.
The order was old/fixed for 1518 and fixed/old for 1410.

| Saved catalogue | Families | Old decode CPU (s) | Repaired decode CPU (s) | CPU reduction | Old / repaired peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|
| 1518 forward, all 104 cuts | 31,687 | 32.433454 | 24.137587 | 25.6% | 1,017.59 / 1,012.91 |
| 1410 reverse, cut 40 | 1,997 | 3.441205 | 2.980039 | 13.4% | 120.02 / 116.13 |

Decode CPU includes preparation and complete lazy assembly/export; catalogue
loading and binary publication are measured separately in the
[decoder evidence](decoder-evidence.json). RSS includes import, load, decoding
and publication in each fresh process, measured before validation loads a second
full record. The earlier profiler's process peak included matching and is not
used for this memory comparison. Profiler overhead is excluded from the timings.

All four executions preserved **every field of the complete decoded record**,
verified by both canonical SHA-256 and deep value equality. The catalogues retain
31,687/1,997 family tables and 22/5 shared numeric responses respectively.
The measured repaired source/engine identity is
`e065f2765298deff05b5cb72a5b6521b4d656fef650dad479d1044bdb8fdb6af`.
The evidence records the baseline identity, exact module hashes and output
digests. Validation passed 15 inventory tests, including preparation sharing,
ordered placement/cut-boundary preservation and non-transitive relationships;
64 event checks, 68 related checks, 26 smoke checks and five post-cleanup
checkpoint checks also passed. These groups overlap and are not a unique total.

These two local comparisons establish a CPU improvement with exact saved-output
preservation on their measured scope. The large case's memory is essentially
unchanged; this repair does not demonstrate recovery of the campaign's memory
failures. The cluster implementation and results remain frozen, and these
measurements are not a full Golden rerun or a universal no-regression claim.

## Tolerance-2 recovery of the remaining reactions

The frozen full run recovered **1,831/1,851 reactions (98.92%)**, counting either
direction. Nineteen reactions had two completed negative reference checks;
case 1799 remained unresolved because both checks timed out. These 20 reactions
were tested locally with **one seed, cap 2,000 and matching tolerance 2**.
Only matching tolerance changed: self tolerance remained 1, with the original
inputs, full bidirectional cut sweep, random seed 42 and reference policy.

All 20 reaction outcomes are now determined: **four recovered and 16 did not**.
Cases **1553, 1739 and 1740** were completed misses at tolerance 1; **1799** resolves
a previous timeout. The combined recovery is **1,835/1,851 (99.14%)**. This is the
union of the baseline and a selected follow-up, not the accuracy of a full run
at tolerance 2. A capped miss does not establish mapping absence.

The [saved evidence](recovery-evidence.json) currently contains 39 completed
reference checks, without errors or unknown verdicts. The remaining check is
1799 forward; its completed reverse check already establishes reaction recovery.
Whole-pass verification timing remains pending. No event decoding was run.

The Mac study used six workers with a 4 GiB guard and fresh results and caches.
All **2,033 matching cuts and 40 reductions** completed without failures or retries.
Their task spans were **23.70 seconds** and **65.72 seconds**, respectively;
compute CPU totals were 112.61 and 137.20 seconds, excluding publication.
Peak sampled reduction memory was **3.03 GiB**, during publication.
Native task, cache and saved-fragment checks confirmed discovery tolerance 2
and self tolerance 1; all 685 inspected case-7 fragment witnesses remained valid.

This is separate from the historical
[13/21 recovery study](../golden_seed_scaling_20261102/README.md), which used
**ten seeds**, cap 2,000 and adaptive tolerances under the earlier implementation.
That study recovered eight cases at tolerance 1, one additional case at 1.5 and
four at 2. Its seed count and cohort differ from the current one-seed follow-up.
The proposed cap-100 speed comparison has not been started.

## Stopped initial campaign

The initial campaign used revision
`2658966b5cc1c73e5510ee2d8001b29e7db989ff` and Linux identity
`ab93d412975f9194d794403d7eb2e20716d7a27a464a29ae76a69cef7e10e1bf`.
The user approved repair followed by a fresh run. Cancellation was verified at
**2026-11-03 06:10:52 UTC**, with no campaign jobs or steps remaining.

| Stage | Full denominator | Preserved terminal results |
|---|---:|---:|
| Matching | 3,702 directions / 177,552 cuts | 3,699 directions complete |
| Reduction | 3,702 directions | 3,202 complete |
| Decoding | 3,702 directions | Not started |
| External reference verification | 3,702 directions / 1,851 reactions | Not started |

The saved queue contained **31 running attempts**, interrupted deliberately, and
264 pending tasks. These interruptions are separate from scientific failures;
their final CPU and memory measurements are unavailable. No accuracy or full-run
timing result is derived from this incomplete campaign.

One reduction attempt, **1561 forward**, exceeded its 2 GiB guard while publishing
output: sampled peaks were **1.364 GiB during computation** and **2.141 GiB during
publication**. An isolated 4 GiB retry completed, peaking at **2.040 GiB** and using
**63.19 compute CPU seconds plus 4.819 publication CPU seconds**. Its original
task and failed attempt are preserved. Consumed computation time increased from
99.19 to 162.98 seconds within the original 3,600-second budget; no fresh retry
time was granted. The replacement resource plan accounts for this known need.
[Stopped-run evidence](stopped-baseline-summary.json) retains the stage counts,
ledger hashes, interruption status and both reduction attempts. Its outputs and
caches are not reused by the fresh benchmark.

## Historical comparison and interpretation

The [previous completed Golden run](../golden_retest_20261101/README.md), revision
`096a010961f53944ddd4ddc433bd8d56f282a2f2`, recovered **1,830/1,851 reactions
(98.87%)** and **3,618/3,702 directions**, with no unknown or failed verification.
It used one seed and **cap 100**. The new full run changes both the scientific
implementation and the branch cap, so differences cannot be attributed solely
to placement compression.

The completed matching/reduction audit verified all **7,404 stage rows**, all
**181,254 logical tasks**, and **181,256 distinct attempts**, including the two
failed reductions. All **188,658 published byte ranges (5.22 GB)** passed their
checksums, with no overlapping ownership of a task or worker.

| Stage | Historical cap 100 compute CPU (s) | Current reported compute CPU (s), 07:34 UTC | Ratio |
|---|---:|---:|---:|
| Matching | 5,136.07 | 10,282.05 | 2.00× |
| Reduction | 2,913.60 | 14,405.38 | 4.94× |
| Decoding, recovery unfinished | 8,044.71 | 99,290.16 | 12.34× |
| Verification, 3,693 checks complete | 1,334.67 | 107,790.28 | 80.76× |

Current totals include failed attempts with exact computation measurements.
Compute-phase memory failures may retain only sampled CPU, so the decoding
compute column does not include all failed work. Across 3,716 completed decoding
attempts, **3,652 successes used 92,412.41 exact worker CPU seconds**, and
**64 failed attempts used 24,806.21 sampled CPU seconds**: 117,218.62 observed
worker CPU seconds in total, without double-counting. The nine verification tasks that were live at that measurement
are excluded from its CPU total.

Compute columns exclude publication and idle time. The repaired
single-cut timings do **not** establish a full-dataset speedup: the current run
uses more matching and reduction CPU than the historical run. Final comparison
will include decoding, verification and publication, alongside reference
coverage and the different scientific settings.

The matching audit measured a **97.63% cache hit fraction** across 1.095 billion
fragment calls, with caching enabled for every cut and no recorded cache errors.
The largest retained cache snapshot was **64 MiB**. Loading and saving consumed
50.05 and 1,502.98 summed elapsed seconds across tasks respectively; 12.65 GB of
cumulative cache writes is not the retained disk footprint. The sampled matching
process-tree peak was **341.47 MiB**. These measurements do not imply that the
cache bounds total worker memory.

Verification uses the complete heavy-atom reference relation, including
unmatched atoms, modulo chemically colored endpoint symmetries. Hydrogen
correspondence is unconstrained. Recovery means the annotation is present among
retained alternatives; it does not measure first-output selection accuracy. A
capped miss is not proof of mapping absence. Final reporting will distinguish
completed outputs, failures, unknown verification and missing work, with stage
CPU, publication cost, elapsed time and sampled memory reported separately.

The verifier tries retained families sequentially. Its **nominal 10-second
per-family query budget does not interrupt preparation**; it is not a hard
per-family wall-time bound. The worker retains the original **3,600-second
direction task guard**. The final six timeouts were **1799 in both directions,
1798 forward, 575 reverse, 1606 forward and 1351 reverse**, each retaining its
original budget and with no retries. All six stayed below their 2 GiB memory
guard. These are execution timeouts, separate from a completed query returning
unknown or not recovered. Completed direction results currently comprise
3,615 recovered, 80 not recovered and one unknown. The run does not establish a
universal absence of speed or memory regressions.
