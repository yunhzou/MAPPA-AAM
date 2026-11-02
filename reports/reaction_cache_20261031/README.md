# Reaction cache snapshots

Implemented the user-approved shared cache with a permanent first snapshot.
The final distributed pilot preserved all 1,842 matching graphs and recovered all 20
directional references. Shared exchange cost 4.8% more wall time and 7.0% more
worker CPU than RAM reuse alone on this small workload; it does not demonstrate
a throughput improvement. Implementation and validation are complete.

## Loading, using and saving

Each reaction/configuration has two files: `baseline.cache` and `latest.cache`.
Forward and reverse occupy separate sections. The first saved section for each
direction remains unchanged in the baseline. Latest versions atomically replace
their predecessor, so open readers can finish using the previous inode. A failed
latest load tries the validated baseline and records the failure. Both unusable
copies cause an explicit error. Versions do not accumulate into seed/cut files.

Workers receive the latest version with their existing task claim. They keep
prepared sessions and useful cache entries in RAM, load only changed directional
payloads, and do not reread their own publication. Import retains useful local
entries within the same 64 MiB accounted retention budget and skips identical
requests. Every hit still checks occupancy, locks, topology and boundaries.
Most recently successful entries are checked first within a conditional key.
Snapshots contain completed compressed results and their dependencies, without
expanded bijections or the separate symmetry-certificate memo.

One coordinator owns the cache directory and admits one publisher per reaction
version. Only that worker serializes and stages an update. The coordinator
accepts the task before replacing the published file. First snapshots are
immediately eligible; later writes are coalesced to at most one per second per
reaction. Workers continue checking the latest published version at task start.
Complete scientific outputs are saved independently of cache publication.

Workers prefer their current reaction/direction, then an unassigned case, then
any runnable task. This preserves reaction reuse while letting surplus workers
process cuts. Scheduling and file I/O remain outside the matching algorithm.
See [design.md](../../design.md#shared-reaction-cache-snapshots) and the
[execution API](../../docs/shared-execution.md).

## Final validation

The final source revision is `8ba48aa0c4bae1c6d69f766311eb99018a6625c2`.
The Mac full suite passed **733 tests in 56.68 seconds**. Linux cache contract
tests passed **40 tests in 5.63 seconds** before the final pilot. Coverage includes
exact raw/finalized matching, lock/cut invalidation after transport, branch caps,
budget-limited imports, cross-process element IDs, immutable baseline sections,
failed/missing latest loads, simultaneous readers and publication, stale writers,
pool restart, RAM/bundle transport, global deduplication and decoded records.

The distributed pilot completed with exit `0:0` in **1 minute 37 seconds**,
including the native build and Linux tests. One persistent pool used one `srun`
and 16 compute workers with a 4 GiB worker limit. All 16 workers processed tasks; peak task overlap was 16. This measures
task concurrency, not continuous saturation of all CPUs.

The pilot used EnzymeMap record IDs 0–9 in both directions. These **ten records
represent five unique endpoint-graph pairs**, not ten independent chemistries.
Both passes used the full cut sweep, one seed, branch cap 100, tolerance 1,
random seed 42 and existence growth. RAM-only ran first, shared exchange second
in the same pool. All finalized cut graphs were compared exactly; shared results
then underwent global deduplication, decoding and reference verification.

| Final measurement | Worker RAM cache | Shared snapshots |
|---|---:|---:|
| Matching wall time | 5.667 s | 5.942 s |
| Total worker CPU | 44.419 s | 47.542 s |
| Fragment cache hits / 92,132 requests | 55,024 | 56,670 |
| Dependency checks | 77,322 | 80,523 |
| Cache load + save, summed worker time | 0 | 1.564 s |
| Snapshot publications | 0 | 22 |
| Cache errors | 0 | 0 |

Deduplication took 1.483 s, decoding 6.504 s and verification 1.722 s. All
**20/20** directional heavy-atom reference queries recovered the annotation.
Peak observed worker memory across the stages was **221.8 MiB**. The cache
directory ended with ten snapshot files plus one ownership file, **2,598,652
bytes** total, and no pending files. No fallback was needed in this pilot;
fallback recovery was exercised by the contract tests.

Full inputs, configuration, code/engine fingerprint, task metrics, file sizes
and reference witnesses are retained in the [public pilot projection](pilot.json).
The maintained runner is [cache_pilot.py](../../bench/cache_pilot.py) with
[portable allocation guidance](../../bench/README.md). One ordered pair of
short runs cannot establish a general speedup or slowdown on the full benchmark.
The measured cost remains visible rather than claiming zero overhead.

## Large-cache local serialization measurement

The archived 374-atom EnzymeMap 10383 forward / cut 334 cache contained 1,457
entries. Measurements on the M4 Pro Mac are in
[snapshot-profile.json](snapshot-profile.json):

| Operation | Measurement |
|---|---:|
| Serialized cache | 20,827,869 bytes |
| Level-1 compressed cache | 2,071,374 bytes |
| Median export | 7.73 ms |
| Median compression | 37.78 ms |
| Median decompression | 6.96 ms |
| Native import | 6.20 ms |

Five export/compression samples and one fresh import were measured before the
later import-retention and lookup-order improvements. The receiving session had
an empty symmetry-certificate memo. Its complete raw search graph equalled the
original, and all 2,476 fragment requests hit imported completed results. These
are local serialization measurements, not shared-filesystem timings.

## Earlier pilot evidence

All three earlier pilots also preserved all 1,842 graphs and recovered all 20
references. Their results document why local entries are retained on import and
why reaction placement and publication coalescing are part of the final design.

| Public evidence | RAM / shared wall time | RAM / shared worker CPU | Publications |
|---|---:|---:|---:|
| [Initial pilot](pilot-initial.json) | 5.342 / 5.810 s | 48.989 / 51.607 s | 151 |
| [Retain local entries](pilot-retained.json) | 5.299 / 5.841 s | 48.480 / 52.267 s | 150 |
| [Placement and coalescing](pilot-locality.json) | 8.443 / 5.318 s | 41.331 / 46.581 s | 20 |

The placement/coalescing pilot had substantially greater output-publication delay in its RAM-only
pass (18.268 versus 9.478 summed worker seconds); its wall-time difference does
not establish a cache speedup. The final pilot also records dependency checks
after changing lookup order while retaining every validity check.

Raw operational records and full artifacts are retained privately.
