# Separate stages, one persistent pool

Ordinary local use remains three calls, without a scheduler:

```python
from mappa import search_aam, deduplicate, decode_events

matched = search_aam(problem, config, workers=4)
catalogue = deduplicate(matched)
decoded = decode_events(catalogue)
```

`deduplicate` uses MAPPA's existing family catalogue. It groups results across
seeds and cuts while preserving fragment partitions, ordered actions, constraints
and provenance. Decoding prepares local fragment and boundary responses from
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

`problems` is a dictionary of names to `AAMProblem` objects. Workers stay alive
through all three passes. Every requested matching task must finish before the
reduction pass can start; every catalogue must finish before decoding can start.
Loading a reaction is explicit. The coordinator handles references, not a merged
dataset-wide graph. Matching uses `reused_native` by default; `execution='reference'`
selects the existing reference implementation. Seed scheduling does not support
the indivisible `shared_policies` backend.

`seed_batch` changes task granularity without changing seed count, order or random
policy. A cut with several batches has an explicit finalization dependency. If
all its seeds fit in one batch, the worker finalizes it directly and publishes
once. There is no intermediate scheduling or storage round trip for that cut.

The typed handoffs are exported from `mappa.execution`:

| Result | Contents |
|---|---|
| `MatchingResult` | Input/configuration reference, planned cut/seed scope, finalized cut references. |
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

The ordinary `search_aam(..., intermediate_dir=...)` checkpoint API also uses
worker bundles and one reaction index. Resume requires that index and the same
input/configuration manifest. Earlier directories containing standalone cut
files require their producing revision. `AAMCheckpointResult.cut_refs` exposes
bundle references rather than one path per cut.

## One Slurm allocation

Submit one batch job, for example:

```bash
#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
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

Local tests exercise RAM/file transport, full scientific-record equivalence,
deadlines, worker replacement, dependency barriers, checkpoint resume and the
actual network worker entry point with an intercepted single `srun` launch.
A [real two-node NRT test](../reports/slurm_pool_20261031/README.md) also verifies
BundleStore handoff, persistent workers and complete scientific equivalence. Its
small-case timing does not establish a throughput improvement or network recovery.
Scheduling stays in `mappa.execution` and `mappa_core/core/scheduler`; matching,
family reduction and event algorithms remain their existing independent owners.
