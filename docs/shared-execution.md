# Separate stages, one persistent pool

Ordinary local use remains three calls, without a scheduler:

```python
from mappa import search_aam, deduplicate, decode_events

matched = search_aam(problem, config, workers=4)
catalogue = deduplicate(matched)
decoded = decode_events(catalogue)
```

`deduplicate` shares completed fragments by reactant/product atom-set pairs and
groups their branch compositions across seeds and cuts. It preserves explicit
anchors, ordered global placements and sweep provenance; intermediate growth
proofs stay in the raw matching archive. Decoding prepares local fragment and boundary responses from
that catalogue. A cut is a matching work unit, not an independent decoding input.

## Many reactions on a computer

Install `.[execution]` for worker memory supervision. Put process-pool calls in
`main()` under `if __name__ == '__main__':` so spawned workers can import the script.
Each case names one reaction and direction in its own atom-index frame.

```python
from mappa.execution import WorkPool

with WorkPool.local(workers=4) as pool:
    matching = pool.search(problems, config).wait_and_validate()
    catalogues = pool.deduplicate(matching).wait_and_validate()
    decoding = pool.decode(catalogues).wait_and_validate()
    result = decoding.load('reaction-1-forward')
```

`problems` is a dictionary of names to `AAMProblem` objects. Pool entry waits until all requested compute workers are ready.
`startup_elapsed_seconds` reports that cost separately; `startup_seconds` sets
the startup timeout (120 seconds by default). Workers stay alive through all
three passes. Every requested matching task must finish before the
reduction pass can start; every catalogue must finish before decoding can start.
Loading a reaction is explicit. The coordinator handles references, not a merged
dataset-wide graph. Matching uses `reused_native` by default; `execution='reference'`
selects the existing reference implementation. Seed scheduling does not support
the indivisible `shared_policies` backend.

`seed_batch` changes task granularity without changing seed count, order or random
policy. A cut with several batches has one dependent task that combines them
in their original seed order. If all its seeds fit in one batch, the worker
publishes the completed cut directly. There is no intermediate scheduling or storage round trip for that cut.
Matching counts retain each seed's index, matcher CPU/wall time and cache
call/hit/dependency-check deltas. These small measurements also accompany the
task audit, so timing analysis need not deserialize its search graph. Task
totals separately include preparation and graph combination; cache transport
and publication remain separately measured by the worker.
Workers prefer their current reaction/direction, then an unassigned one. If only
already assigned reactions have ready tasks, workers take those tasks too, so a
single reaction or a long tail can still use the full pool. This is a scheduling
preference; scientific seed/cut plans and cache-validity rules do not change.

The typed handoffs are exported from `mappa.execution`:

| Result | Contents |
|---|---|
| `MatchingResult` | Input/configuration reference, planned cut/seed scope, completed cut references. |
| `CatalogueResult` | One reaction/direction catalogue deduplicated across the declared cuts/seeds, search configuration, scope and stop counts. |
| `DecodingResult` | Original-space decoded record, event/search configurations, scope and stop counts. |

`wait_and_validate()` checks execution completion. It does not certify exhaustive
matching: branch caps and scientific stop reasons remain in the saved results.
Reference verification is a separate query on the retained mapping relation.

## RAM or aggregated storage

`MemoryStore` is the default. Intermediate objects remain in coordinator RAM;
workers transfer objects through the pool connection. Loading locally returns
the stored object, which callers must treat as immutable. RAM storage retains
published artifacts for the lifetime of the store. It is not bounded-memory
storage for an arbitrarily large benchmark.

For durable handoff or a large dataset:

```python
from mappa.execution import WorkPool, BundleStore, StageResults

store = BundleStore('/shared/run/results')
with WorkPool.local(4, store=store) as pool:
    matching = pool.search(problems, config).wait_and_validate()
    catalogues = pool.deduplicate(matching).wait_and_validate()
    catalogues.save_manifest('/shared/run/catalogues.json')

# A later process can decode without repeating matching or deduplication.
catalogues = StageResults.open('/shared/run/catalogues.json', store=store)
with WorkPool.local(4, store=store) as pool:
    decoding = pool.decode(catalogues).wait_and_validate()
    decoding.save_manifest('/shared/run/decoded.json')
```

Workers append checksum-verified records to their own bundles, rolling over at a
256 MiB target. One larger record can exceed that target. **There is no file per
seed or cut.** Small compressed records use an 8 MiB memory spool; larger ones
spill to node-local scratch. Manifests contain references to immutable byte
ranges and the publishing source fingerprint. Keep the referenced bundle paths
available. These are trusted Python checkpoints, not an interchange format for
untrusted inputs. Stage transitions use the same storage mode; RAM transitions
also require the same `MemoryStore` instance.

BundleStore pools also share completed-fragment caches between workers by
default. At task start a worker loads a newer reaction snapshot if needed;
otherwise it keeps using its prepared RAM cache.
Snapshot import retains useful local entries within the existing memory budget;
identical requests are not inserted twice. An update to the other direction
does not force a reload of an unchanged directional snapshot.
Existing lock, occupancy, topology and boundary checks still decide every cache hit. Only a worker granted
publication by the coordinator writes an update. Its own next task acknowledges
that update without reading the same snapshot back from disk.

The `cache/` directory retains two files per reaction/configuration: the first
cache for each direction in `baseline.cache`, and an atomically replaced
`latest.cache`. Forward and reverse use separate sections within these files.
Failed latest loads try the checked baseline and record the error. Versions do
not accumulate; there are no seed/cut cache files. One ownership file prevents
two independent coordinators from writing the same cache directory. Result
bundles remain independent of cache snapshots.

The native cache keeps its 64 MiB accounted retention budget per prepared
direction. Snapshots retain completed compressed results and dependency checks.
Task metrics report cache hits, bytes, load
and save time, successful publication and baseline recovery. For comparisons,
`pool.search(problems, config, shared_cache=False)` disables disk exchange while
keeping native RAM reuse. `WorkPool.local(..., shared_cache=False)` disables it
for the whole pool. RAM-only pools create no cache files. Snapshot compatibility
includes the actual code/native-engine fingerprint as well as inputs/settings.
First snapshots are published immediately. Later writes are coalesced to at most
one per second per reaction by default; workers continue accumulating entries in
RAM and always check the latest published version at task start. Set
`cache_publish_seconds=0` on the pool to allow every eligible update. This changes
cache I/O frequency only; complete scientific outputs are always saved normally.

The ordinary `search_aam(..., intermediate_dir=...)` checkpoint API also uses
worker bundles and one reaction index. Resume requires that index and the same
input/configuration manifest. Earlier directories containing standalone cut
files require their producing revision. `AAMCheckpointResult.cut_refs` exposes
bundle references rather than one path per cut.

## One Slurm allocation

Prepared EnzymeMap and Golden campaigns use the
[portable benchmark launcher](../bench/README.md), including the global
matching/reduction barrier before decoding. Keep site accounts, partitions,
hostnames, paths and raw scheduler logs in private configuration or ignored
storage outside the published evidence.

For a custom Python campaign, submit one batch job, for example:

```bash
#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=40G
#SBATCH --time=02:00:00
python campaign.py
```

Inside `campaign.py`, use the same three stages:

```python
from mappa.execution import WorkPool, BundleStore

with WorkPool.slurm(workers=16, memory_bytes=2 * 2**30,
                    store=BundleStore('/shared/run/results')) as pool:
    matching = pool.search(problems, config).wait_and_validate()
    catalogues = pool.deduplicate(matching).wait_and_validate()
    decoding = pool.decode(catalogues).wait_and_validate()
    decoding.save_manifest('/shared/run/decoded.json')
```

`WorkPool.slurm` requires an existing allocation. It makes **one startup `srun`**
with one worker service per node. Each service starts the requested number of
persistent compute processes; tasks never invoke Slurm. The coordinator's queue
lives on local scratch. Worker services must reach its hostname over the private
cluster network, share the Python environment and bundle paths, and use matching
Python source and compiled engines. Source/engine fingerprints are checked at
worker startup. RAM transport is also available, with coordinator memory and
network traffic proportional to the retained artifacts.

Memory limits are per compute worker. Leave additional allocation memory and CPU
capacity for the coordinator, supervisors and operating system. Each worker runs
one scientific operation at a time with native thread counts set to one. There
are no nested compute pools. Worker memory is sampled on a clock, normally every
0.2 seconds; zero samples on a short task do not mean zero memory use.

## Completion and recovery

Dependencies become ready after successful publication. Exceptions, timeouts,
memory limits and lost leases remain explicit failures. A failed compute process
is replaced; its task is not silently retried. A case budget counts overlapping
active computation once and pauses when no computation is active. Publication
has a separate 120-second guard.

Bundle-backed pools checkpoint the queue at stage completion and context exit.
Completed stage manifests can be reopened in a new pool. The lower-level
`TaskQueue` provides deliberate retry and recovery from queue backups, retaining
attempt history. The public pool does not automatically resume an interrupted
pass or fail over its coordinator. Work after the latest checkpoint can require
repeating after node loss. Plan allocation walltime accordingly.

After diagnosing a decoding memory failure, `pool.retry_decoding(case, audit,
catalogue=original_reference)` restores that archived task into a new pool. It
preserves the task ID, decoder configuration, previous attempts, memory guard
and consumed case time. Only the remaining original budget is available; an
exhausted budget or changed guard is rejected. The original catalogue reference
must match the failed task. Use a separate `BundleStore` for the recovered output
when retaining a completed campaign as immutable evidence. The
[EnzymeMap recovery command](../bench/README.md#recover-decoding-after-a-repair)
adds original/current code provenance and explicit resume checks around this API.

Local tests exercise RAM/file transport, full scientific-record equivalence,
deadlines, worker replacement, dependency barriers, checkpoint resume and the
actual network worker entry point with an intercepted single `srun` launch.
A saved [two-node pool test](../reports/slurm_pool_20261031/README.md) also verifies
BundleStore handoff, persistent workers and complete scientific equivalence. Its
small-case timing does not establish a throughput improvement or network recovery.
Scheduling stays in `mappa.execution` and `mappa_core/core/scheduler`; matching,
family reduction and event algorithms remain their existing independent owners.
