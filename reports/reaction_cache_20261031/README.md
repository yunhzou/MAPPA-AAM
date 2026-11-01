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

The first pilot, job **21132532**, completed successfully with 16 workers on
`cpu1-00102` and `cpu1-00105`. All 1,842 complete finalized cut graphs were exactly
equal. Deduplication and decoding completed for all 20 directions; all 20 exact
heavy-atom reference queries recovered the annotation. There were no cache
errors. The ten input records contain five unique endpoint-graph pairs, so the
cache directory held ten snapshots plus one ownership file (332,113 bytes total).

However, replacing each worker's RAM cache on import cost more than the additional
reuse saved in this small pilot: shared exchange took 5.810 s / 51.607 CPU s,
versus 5.342 s / 48.989 CPU s for RAM reuse alone. Cache load/save summed to
3.272 worker-seconds; hits increased from 53,036 to 54,040 out of 92,132 requests.
The initial result is retained in [dfw-pilot-initial.json](dfw-pilot-initial.json).
This motivated retaining useful local cache entries during import, with exact
request deduplication and the same memory budget. Updates to the other direction
also no longer force loading an unchanged directional payload. Coordinator
publication still selects one snapshot; it does not merge caches or branches.

A follow-up pilot will measure that correction on the same ten records and scope.
The local full suite after retaining imported/local entries passed **726 tests
in 54.19 seconds**. A further targeted check covers acknowledgement when another
direction has published before the original worker's next task.

The second pilot, **21132861**, again completed all 1,842 exact graph comparisons
and all downstream stages. Retaining local entries raised cache hits to 55,273,
but shared exchange still took 5.841 s / 52.267 CPU s versus 5.299 s / 48.480 CPU s
for worker-only reuse. Loading/saving cost 4.061 summed worker-seconds. The
measurements are retained in [dfw-pilot-retained.json](dfw-pilot-retained.json).

The remaining correction is outside the matcher: prefer a worker's current
reaction/direction, then an unassigned case, then any runnable task to keep the
pool saturated. Coalesce rapid latest-snapshot writes to at most once per second
while allowing each first directional snapshot immediately. This preserves
latest-version checks and baseline recovery without writing after every tiny cut.
The final pilot will compare snapshot exchange on/off under this same placement.
The final local implementation passed **729 tests in 52.46 seconds**. Its
execution tests explicitly check that distinct reactions are assigned first,
workers retain their current reaction, surplus workers still receive cuts,
and first snapshots are exempt from the publication interval.
