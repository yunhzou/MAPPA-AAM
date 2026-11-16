# Fresh Golden benchmark: tolerance 2, branch cap 100

A fresh full Golden benchmark was launched on 5 November 2026 using
[`e4e1560`](https://github.com/yunhzou/MAPPA-AAM/commit/e4e15608c9a8bfedee518fa3abb24fd2c5c6f118).
It covers **1,851 reactions / 3,702 directions**. Only immutable endpoint and
reference inputs are reused: matching results, reduced catalogues, decoded
outputs, verification outcomes and reaction caches all start empty.
All **5,553 distinct input/reference blobs** passed checksum verification.

| Setting | Value |
|---|---|
| Discovery tolerance | 2 |
| Branch cap | 100, with recursive subtree pruning |
| Seeds / random seed | 1 / 42 |
| Search scope | Both directions, uncut plus the full single-edge sweep |
| Input weights | Original prepared bond-order matrices |
| Fragment symmetry | Element- and anchor-preserving binary symmetry; no self-tolerance filter |
| Graph / cut presence threshold | 0.2 / 0.2 |
| Ordinary / metal event threshold | 0.5 / 0.3 |
| Decoding window | Model minimum through minimum + 5 |
| Fragment optimization allowance | Default five seconds per unique fragment |
| Compute workers | 128 in one persistent pool |
| Allocated CPUs, including execution overhead | 160 |
| Per-worker memory guard | 2 GiB |
| Direction compute allowance | 600 seconds; no automatic retries |

The pool runs matching, reduction, decoding and external reference verification
as separate stages. All matching and reduction must finish before decoding
starts. Shared reaction caches are enabled within this fresh campaign. Operational
configuration and raw scheduler records remain in the ignored private directory.

The decoder retains shared fragment responses, the independent ordinary-boundary
approximation and the joint sweep correction. Model minima and selected event
windows are distinct from whole-reaction physical minima. External reference
verification queries the complete retained matching catalogue; it does not
measure reference membership inside the reported event window.

Preflight: the cleanup passed 1,133 local tests; the synchronized Linux checkout
passed **109 tests in 28.19 seconds**, covering binary symmetry, saved-path queries,
windowing, persistence and the Golden runner. The existing Linux native engines
were retained after confirming that native sources were unchanged.

## Reference recovery

All **3,702 directions completed matching, reduction and decoding successfully**.
The full run is terminal, including external verification: **1,834/1,851 reactions
recovered (99.08%)**, counting either direction, and **17 not recovered**.

The original external verification produced **3,629 recovered, 69 not recovered, one unknown
and three timed-out directions**. The unknown is 561 reverse. The 600-second
compute timeouts are **1351 reverse, 1799 reverse and 1801 reverse**; each used
approximately 600 sampled CPU seconds and 0.26–0.32 GiB sampled memory. All four
reactions recovered in the forward direction, leaving no unresolved reaction-level
outcomes. These reference-query limits are not mapping or decoding failures.
There were **no memory-limit failures, retries or increased allowances**.

Compared with the older tolerance-1/cap-100 full run, the new run gains cases
**1227, 1553, 1739 and 1740**, with no previously recovered reaction lost.
Compared with the tolerance-1/cap-2,000 full run, it gains confirmed recovery of
**1553, 1739, 1740 and 1799**, but does not recover **1606**. Case 1799 was previously
unresolved rather than a confirmed miss. Tolerance, pruning, symmetry and decoder
design differ between these historical measurements; this is not an isolated
test of a code or parameter change.

The 17 unrecovered reactions are **7, 19, 20, 44, 590, 867, 964, 986, 1034, 1228,
1285, 1314, 1358, 1377, 1475, 1568 and 1606**. These are reference-recovery outcomes
for the retained capped catalogue, not proofs that valid mappings do not exist.

## Time and resource use

All **128 workers registered in 16.953 seconds**. Saved attempt intervals confirm
128 simultaneous matching tasks and 128 simultaneous decoding tasks. Registration
and peak overlap do not establish continuous full utilization.

| Stage | Parallel elapsed seconds | Compute CPU seconds | Publication CPU seconds |
|---|---:|---:|---:|
| Matching | 373.375 | 2,210.395 | 358.795 |
| Reduction | 82.480 | 1,124.317 | 144.562 |
| Decoding | 98.126 | 7,048.248 | 162.461 |
| External reference verification | 603.718 | 4,369.617 successful + 1,799.830 sampled on timeouts | 4.557 |

Matching plus reduction took **455.855 seconds (7m36s)** in the pool. Matching
covered **177,552 cuts**, including uncut searches. Stage elapsed time includes
scheduling, publication and coordinator overhead; worker CPU totals measure
computation and publication separately. These are different measures and should
not be added together.

Compute CPU divided by stage elapsed time and 128 worker slots is **4.6% for
matching, 10.6% for reduction and 56.1% for decoding**. Matching tasks averaged
only 0.0124 compute CPU seconds per cut. This allocation therefore establishes
parallel execution, not efficient saturation during matching. The recorded
worker timings do not separately identify all coordinator and scheduling costs.

The allocation lasted **1,183 seconds (19m43s)**; the pool measured 1,174.788
seconds including startup. The external verifier's three timeouts account for
the long tail, although their reactions were already recovered in the opposite
direction. The allocation ended nonzero because those checks timed out; all
required stage records are terminal. No retry was launched. Allocated capacity
was 189,280 CPU seconds; scheduler-accounted use was 18,184 CPU seconds, including
execution overhead. Worker scientific and publication timings remain separate
from this allocation-wide accounting.

| Per-direction measure | Median | 90th percentile | 95th percentile | 99th percentile | Maximum |
|---|---:|---:|---:|---:|---:|
| Matching + reduction, compute CPU seconds | 0.463 | 1.825 | 2.927 | 8.763 | 33.294 |
| Decoding, compute CPU seconds | 0.841 | 4.345 | 6.820 | 13.680 | 56.036 |
| Decoding, wall seconds | 0.846 | 4.347 | 6.828 | 13.687 | **56.061** |

The slowest decode was **1654 reverse**. Its publication took another 1.362
seconds. Decoding timings include catalogue loading, decoding and record
construction; publication is separate. Matching CPU sums the cut tasks belonging
to each direction. Summed cut wall time is work duration, not direction latency.

Peak sampled decoding memory was **1.620 GiB**, below the unchanged 2 GiB guard.
The completed-fragment cache recorded **65,898,673 hits / 74,281,547 calls
(88.71%)**, with no cache errors or baseline-recovery loads. Aggregate cache
loading and saving took 447.725 and 205.778 worker wall seconds respectively;
these costs are already included in the stage measurements.

## Verification-only acceleration

Revision [`b93591e`](https://github.com/yunhzou/MAPPA-AAM/commit/b93591ee493741466efc71aa81ddc70fb8469586)
accelerates the shared saved-family membership query. It batches finite
constraints into the native solver, uses its finite-domain engine and removes
only self actions certified to be already covered by discovery. Reference
identity, unmatched atoms, joint constraints and returned mapping witnesses retain
their meaning. Matching and decoding were neither changed nor rerun.

All **3,702 saved catalogues** were verified again in one persistent 64-worker
pool, with the same **2 GiB guard and 600-second direction allowance**. This was
a separately recorded measurement, not a retry that overwrote the original
campaign. All original 14,808 stage rows remain hash-identical. The audit checked
every task's catalogue/reference identity, output checksum and single attempt;
the final queue is empty.

All **3,698 previously conclusive answers agree**. The new pass returns
**3,632 recovered and 70 not recovered directions**, with **zero unknowns,
timeouts, memory failures or retries**. Reaction-level recovery remains
**1,834/1,851 (99.08%)**.

| Previously unresolved direction | New result | New compute CPU seconds |
|---|---|---:|
| 561 reverse | Recovered | 3.208 |
| 1351 reverse | Not recovered | 182.224 |
| 1799 reverse | Recovered | 6.134 |
| 1801 reverse | Recovered | 3.315 |

The three former failures each exhausted 600 compute seconds; 561 reverse was
previously an inconclusive family-query result. Non-recovery here means absence
from the saved capped matching catalogue, not absence of a chemically valid map.

| Verification measure | Original | Optimized |
|---|---:|---:|
| CPU seconds, same 3,698 conclusive directions | 4,358.311 | 1,497.043 |
| CPU seconds, all attempts | 6,169.447 | 1,691.924 |
| Compute workers | 128 | 64 |
| Recorded elapsed seconds | 603.718 verification stage | 194.173 whole verification pool |

The common conclusive scope uses **2.91× less CPU**. Total recorded CPU is
**3.65× lower**, including sampled CPU from the original timeouts; the new pass
also establishes four additional conclusive outcomes. The new pool elapsed time
includes 5.629 seconds of startup. Its actual task interval is 186.187 seconds,
with 64 simultaneously running tasks. The allocation completed in 201 seconds.
CPU totals exclude publication, which added 4.631 CPU seconds in the new pass.

The maximum direction verification took **182.224 seconds**, for 1351 reverse; median and 99th
percentile CPU were **0.223 and 2.718 seconds**. Peak sampled worker memory was
**0.279 GiB**. Individual directions are not uniformly faster: three took over
one extra CPU second, with a largest increase of 2.638 seconds (575 reverse).
The cluster comparison uses separate allocations. An identical ten-family
paired test on one Mac took **6.491 → 0.683 CPU seconds (9.51×)** with identical
answers; this small hard-query sample is not the full-dataset speedup.

Validation passed **205 affected tests**, **26 smoke tests**, and **42 Linux
query tests**, including independent finite mapping oracles and reconstruction
of valid witnesses after symmetry simplification. Both Golden and EnzymeMap use
the accelerated query engine; these performance numbers measure Golden only.
The existing JSON and compressed direction records retain the original evidence
and add the follow-up under `verification_acceleration`.

## Completion audit and evidence

The independent read-only audit checked the exact **1,851-reaction / 3,702-direction
scope**, input-index hashes, frozen scientific settings and Linux source/engine
identity. All **14,808 required stage rows** are terminal: **14,805 successful and
three external verification timeouts**. All **188,658 tasks have exactly one
archived attempt**, with no retries or duplicate task ownership. The saved queue
is empty, and attempt timestamps verify every global stage barrier.

Every successful ledger output reference was checked against the archived task
output, or against the completed cut references in the coordinator's matching
handoff. All task-audit and verification blobs passed checksum validation. This
audit reads descriptors and execution evidence; it does not rerun scientific
operations or deserialize every complete decoded payload. The immutable input
blobs were independently checksummed before launch.

[Summary and reaction-level outcomes](summary.json) retain settings, counts,
distributions, baseline gains/losses and failure evidence.
[Per-direction measurements](directions.json.gz) retain separate stage CPU,
publication, memory, attempt counts and reference outcomes. Missing complete-task
timings for timeouts are null, with their sampled CPU and consumed budget recorded
separately. Operational configuration and raw accounting remain private.

Historical reference points are the [tolerance-1/cap-100 full run](../golden_retest_20261101/README.md)
(1,830/1,851 recovered) and the [tolerance-1/cap-2,000 full run](../golden_placement_20261103/README.md)
(1,831/1,851). The [binary-symmetry decoder measurement](../golden_binary_decode_20261105/README.md)
used saved tolerance-2/cap-200 catalogues and did not repeat matching or reference
verification. These differences prevent attributing outcome changes to code
alone.
