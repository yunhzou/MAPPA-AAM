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

`execution.decoding.plan_event_decoding` schedules one saved case/cut catalogue
through the public anchored-fragment converter. It accepts partial mappings and
the ordinary `EventDecodeConfig`, including relative minimum+n filtering. It
retains converted branch records, direct counts and explicit incomplete flags.
There is no separate solver/minimum stage or cluster-specific decoder.

Additional stages use `Task(operation='module:function', dependencies=(...))`.
The function receives its payload and a `Workspace`; `workspace.dependency(key)`
loads a verified predecessor result. The operation owns its scientific policy;
the executor does not change thresholds, mapping coverage or result scope.

## Verification

`tests/test_shared_execution.py` covers graph and seed-prefix equivalence,
partial mappings, deadlines, process memory interruption, dependencies, identity,
leases, explicit recovery and bundles. `tests/test_saved_event_decoding.py` covers
the public conversion adapter. Historical performance records belong in `reports/`.
