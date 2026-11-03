# Full Golden benchmark with placement compression

A fresh full Golden rerun was submitted at **2026-11-03 06:19:04 UTC** for
[`1710f69`](https://github.com/yunhzou/MAPPA-AAM/commit/1710f690f5c49b34d11de69782aa54b1f4d79304).
The first campaign exposed a slow symmetry-membership check and was intentionally
stopped before decoding. The repaired implementation has passed **883/883 Mac
tests**, **79/79 Linux tests** and five targeted matching probes. The new run
requests **141 compute workers**; submission does not establish worker readiness
or utilization. Full-benchmark accuracy and timing results remain pending.

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
| Worker memory guard | 2 GiB for 129 workers; 4 GiB for 12 workers in the batch containing case 1561 |

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

Verification uses the complete heavy-atom reference relation, including
unmatched atoms, modulo chemically colored endpoint symmetries. Hydrogen
correspondence is unconstrained. Recovery means the annotation is present among
retained alternatives; it does not measure first-output selection accuracy. A
capped miss is not proof of mapping absence. Final reporting will distinguish
completed outputs, failures, unknown verification and missing work, with stage
CPU, publication cost, elapsed time and sampled memory reported separately.
