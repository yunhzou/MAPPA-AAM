# Reaction cache snapshots

The user approved sharing completed-fragment caches across workers, with the
first snapshot retained as a permanent fallback. The implementation stores
baseline/latest files per reaction/configuration, with separate directional
sections. Native snapshot import/export preserves compressed results and their
lock/topology/boundary dependencies; scheduling stays in the execution layer.

Workers receive the latest version with the existing task claim. They load only
changed versions, retain RAM caches across seeds/cuts, and do not reread their
own successful publication. Only the coordinator-selected publisher serializes
and writes. Latest updates atomically replace old versions; first directional
entries remain unchanged in the baseline. Failed latest loads use the checked
baseline and record the failure. A cache directory has one active coordinator.

## Local validation

**724 tests passed in 55.08 seconds.** Tests include exact raw/finalized matching,
locks and cut invalidation after transport, branch caps, budget-limited imports,
cross-process element-table differences, immutable first snapshots, failed/missing
latest loads, explicit failure when both copies are unusable, concurrent publisher
admission, readers holding the replaced file open, failed publication cleanup,
pool restart, RAM/bundle transport and complete deduplication/decoded records.

The archived 374-atom EnzymeMap 10383 forward / cut 334 cache contains 1,457
entries. Export/import measurements on the M4 Pro Mac are saved in
[snapshot-profile.json](snapshot-profile.json):

| Operation | Measurement |
|---|---:|
| Serialized cache | 20,827,869 bytes |
| Level-1 compressed cache | 2,071,374 bytes |
| Median export | 7.73 ms |
| Median compression | 37.78 ms |
| Median decompression | 6.96 ms |
| Native import | 6.20 ms |

Five export/compression samples and one import were measured. The receiving
prepared session had an independent empty symmetry-certificate memo. Its complete
raw search graph equalled the original, and all 2,476 fragment requests were
served from imported completed results. No mapping bijections were expanded.
These are local serialization measurements, not shared-filesystem timings.

## DFW pilot

`bench/cache-pilot-dfw.sbatch` allocates two cpu_short nodes, ten CPUs and 48 GiB
per node, and starts one persistent pool of eight compute workers per node
(4 GiB per compute worker). It uses the first ten distinct prepared EnzymeMap
reactions in both directions, the full cut sweep, and the original one-seed,
cap-100, tolerance-1 settings with `growth_policy='existence'`.

`bench/cache_pilot.py` compares native RAM reuse alone against shared snapshot
exchange in the same pool, checks every complete finalized cut graph for exact
equality, then deduplicates across all cuts/seeds before decoding and verification.
The full code/engine identity, settings, source bundle checksums, per-stage times,
cache traffic/errors, memory, workers and task overlap are saved in results.json.
First the batch builds the native Linux engine and runs the cache contract tests.

Pilot measurements will be added after the actual Slurm allocation completes.
