# Shared execution

`mappa.execution` distributes independent work across cases. One scientific
operation runs unchanged on a single CPU, local processes, or remote workers.
Each worker owns one compute process and one prepared case; it reuses that
case's existing bounded native caches. It never starts a nested compute pool.

`SearchSession` is the common search implementation behind the existing public
search APIs and distributed search tasks. Matching, locking, symmetry, branch
admission and random growth are unchanged. `plan_search` splits the original
seed list into batches, verifies the executed seed orders, then finalizes each
cut after all its batches complete. `rebenchmark.shared.plan_direction` adds
independent reference-evaluation tasks using the existing Golden evaluator.
The benchmark adapter preserves direction and partial-mapping policy.

## Local use

Install the `execution` extra for process-tree memory supervision. Place this
example inside `main()` and call it under `if __name__ == '__main__':`, as required
by Python's portable process startup. `problem` is an `AAMProblem` and `config`
is the ordinary `AAMSearchConfig`.

```python
from pathlib import Path
import tempfile
from mappa.execution import TaskQueue, BundleWriter, run
from mappa.execution.search import plan_search, iter_cuts

with tempfile.TemporaryDirectory() as scratch:
    queue = TaskQueue(Path(scratch) / 'queue.sqlite', identity='frozen-manifest-id')
    writer = BundleWriter('/absolute/output/inputs', scratch)
    tasks, cuts = plan_search('case-1-forward', problem, config, writer,
                             seed_batch=10, memory_bytes=4 * 2**30)
    queue.add(tasks)  # Add other cases/directions to this SAME queue.
    state = run(queue, workers=4, identity='frozen-manifest-id', scratch=scratch,
                bundles='/absolute/output/results', memory_bytes=4 * 2**30)
    queue.backup('/absolute/output/queue-backup.sqlite')
    for cut in iter_cuts(queue, cuts, scratch):
        consume(cut['graph'])
    queue.close()
```

Use `workers=1` for one CPU. Memory is **per worker**: four workers at 4 GiB
require at least 16 GiB plus coordinator/supervisor overhead. Task memory is an
admission estimate; workers enforce the specified aggregate child-process-tree
guard. An insufficient allocation leaves tasks pending, not successful.

## Cluster use

Prepare the same queue and inputs once. Put the live queue database on the
coordinator's node-local disk. Give every worker the same frozen code, manifest
identity, shared bundle directory and private authentication key. Slurm only
allocates CPUs and starts the common worker; no Slurm logic enters the matcher.
Use a private cluster network or SSH tunnels for this authenticated Python
transport. It is for trusted internal artifacts/workers, not public clients.

The coordinator and worker commands are:

```text
python -m mappa.execution serve --database /tmp/queue.sqlite \
  --backup /shared/run/queue-backup.sqlite --identity MANIFEST \
  --host COORDINATOR --port PORT --auth-file /shared/private-worker.key

python -m mappa.execution worker --host COORDINATOR --port PORT \
  --auth-file /shared/private-worker.key --identity MANIFEST \
  --workers 4 --memory-mib 4096 --scratch /tmp --bundles /shared/run/results
```

The key must contain at least 32 random bytes and have private permissions.
Workers exit when no ready/running work remains. The coordinator stays available
until its launcher stops it after workers exit, avoiding a shutdown race with
final acknowledgments. The embedded `coordinator()` context supports that same
lifecycle in a Python launcher. The launcher must verify the frozen code/input
manifest on every machine; the queue checks its supplied identity at connection.

## Completion and recovery

A task is accepted once, with a lease token. Dependencies become ready only
after successful publication. Timeouts, memory limits, exceptions, worker exits
and lost leases are separate states. A failed worker is replaced to serve other
work; its failed task is not silently retried. `retry(key)` is an explicit action
after diagnosis, retains prior attempts, and does not reset the case deadline.
Completed tasks cannot be retried or overwritten. A case budget starts on its
first claim and is shared by all tasks bearing that case identity, including
subsequent queue wait. It is a wall deadline, not summed CPU time.

Publication has a separate 120-second guard after computation finishes. Results
are compressed on node-local scratch, then appended as verified immutable ranges
to per-worker bundles (256 MiB target size). One large result may exceed that
size. Shared storage has no file per task, seed or solver query. Queue backups
contain the result indexes and checksums. A crash may leave an unreferenced
partial range; that is never accepted as evidence. Restore the last queue backup
locally, wait for stale leases to expire, diagnose them and explicitly retry.
Work accepted after the last backup can be repeated after coordinator loss;
this implementation does not claim transactional durability across node loss.

Execution completion only means the callable returned and its result was
published. A scientific result can still contain solver unknowns or incomplete
decoding. Reporting must audit those results separately. In particular, the
Golden adapter runs search and reference analysis; it does not claim that
full-event decoding or all benchmark-report metrics are complete.

`execution.decoding.plan_event_decoding` schedules the existing balanced-family
decoder with an explicit absolute bound (or unbounded enumeration). It preserves
incomplete family reports. It rejects relative minimum+n configurations until
the caller supplies the proper bound for the requested prefix. It does not
replace Golden's partial-mapping decoding policy.

Additional stages use `Task(operation='module:function', dependencies=(...))`.
The function receives its payload and a `Workspace`; `workspace.dependency(key)`
loads a verified predecessor result. Other family solvers can use this interface after catalogue construction
and valid bounds for the requested seed prefixes. No decoder semantics or minimum/window policies are changed here.

## Validation and performance

`tests/test_shared_execution.py` checks original graph and seed-prefix equality,
finite balanced decoded-set equality, partial mappings, reference evaluation,
deadlines, actual memory interruption, worker exit, dependency ordering, identity,
leases, explicit recovery and bundle checksums. Computational tests run on NRT.

`python -m rebenchmark.execution_pilot` is a two-node Slurm validation on two
small reactions in both directions. It compares seed batch sizes 1, 2, 5 and 10,
checks every cut graph and seed verdict against the original operations, and
records actual task distribution, CPU, wall, publication and memory. It is an
execution control, not a Golden or EnzymeMap coverage benchmark. Production-scale
throughput and the final EnzymeMap launch remain separate acceptance gates.

Validated on NRT on 2026-10-27/28: commit `146327c` passed 58 checks
(job 7700833) and the two-node pilot (7700834). Across seed batches 1/2/5/10,
120 cut graphs and 1,200 seed evaluations matched the original operations.
Each variant used 11 shared files for 360–630 tasks. These small controls
establish execution equivalence, not production throughput.

The final process-tree monitoring change, `ba12951`, passed the same 58 checks
(job 7703082, 81.95 s). On that node, 30 samples averaged 20.94 ms with the
previous whole-node process scan and 0.070 ms with the worker-tree scan.
This measures monitoring alone. The two-node end-to-end pilot predates this
last monitoring change; no production speedup is inferred from these numbers.
Its repeat and campaign rollout require restored shared-project file capacity.
