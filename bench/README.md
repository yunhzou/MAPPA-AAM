# Benchmarks

Golden and EnzymeMap use one pipeline and the public MAPPA algorithms. Input
preparation is independent of execution: prepare an immutable artifact once,
then name its index and hash in each compatible run configuration.

| Owner | Responsibility |
|---|---|
| `rebenchmark/inputs.py` | Parse datasets, freeze task/input bundles, validate reuse |
| `rebenchmark/operations.py` | Search-dependent reference and event task adapters |
| `rebenchmark/execution.py` | Queue construction, pilot reuse, workers, checkpoints and resource stages |
| `rebenchmark/launch.py` | Slurm submission, allocation handoff, completion routing and status |
| `rebenchmark/report.py`, `timings.py` | Saved-output checks, coverage and bounded-memory timing aggregation |
| `rebenchmark/evidence.py` | Atomic receipts, hashes and verified archives |
| `rebenchmark/run_conversion.sbatch` | The single benchmark allocation wrapper |

The core scheduler owns claims, dependency ordering, directional budgets and
verified result publication. Benchmark modules do not implement another scheduler
or another matching/decoding algorithm. `decoding.py` calls the public decoder
and independently recounts representative events. `golden_evaluation.py` checks
original AAM reference coverage. Coverage is separate from representative event
counts and non-exhaustive correlated event tables; no global minimum is claimed.

## Prepare and run

Use `PYTHONPATH=src:bench` with the production Python environment. Commit source
before submission. All paths in plans must be absolute and accessible on compute
nodes. Cluster account, partition, QoS and stage resources live in the saved plan.

```bash
python -m rebenchmark.launch prepare-inputs /absolute/input-artifact --plan /absolute/settings.json
```

After preparation completes, add `index` and `index_sha256` to the run plan,
pointing to the artifact's `index.json`. Dataset hash, search configuration,
directions, cut/seed grid, reference policy and task resource specifications must
match. Incompatible inputs require explicit new preparation; submission never
silently reparses the dataset. Resource placement (worker count, Slurm QoS and
allocation size) can change without rebuilding compatible inputs.

```bash
python -m rebenchmark.launch submit /absolute/new-run --plan /absolute/run-plan.json
python -m rebenchmark.launch status /absolute/new-run
```

Submission creates a real-input pilot, separate queue preparation and a production
controller. Only that run's verified pilot computations are reused. Old scientific
results are not imported into a fresh run. The controller requests the configured
compute target, including capacity currently busy; status distinguishes requested
workers, running workers and observed CPU use.

The controller stops claims early enough for active work and publication to
finish, waits for worker retirement, hashes a consistent queue and submits its
successor. Continuations retain task identities, attempts and remaining budgets.
Authorized memory escalation selects only outcomes from the specified smaller
memory guard. Both normal execution and the memory stage publish checkpoints and
route terminal results to reporting. Errors, timeouts and missing outputs remain
explicit in the denominator. They are never converted into negative coverage.

All live queues and scratch use verified node-local disk (`/var/tmp` on NRT).
Shared storage receives bounded bundles and consistent backups. Failed queue
owners retain their database/WAL and record the location. Scheduler transport
retries apply only to idempotent queries/cancellation, never blind resubmission.

## Reporting and retained evidence

Reporting reads exact accepted output ranges and checks every planned task key,
specification and parent identity. It processes bounded batches of cases and stores
exact timing samples on local disk. Seed prefixes, uncut/sweep and direction unions
retain unknown outcomes; execution completion is not scientific completeness.
CPU time, compute/publication wall time, peak process-tree RSS and claim-to-accept
concurrency are reported separately. Claim concurrency is not CPU utilization.

Saved archive readers, hierarchical-mapping checks and scientific metric adapters
remain for inspecting existing evidence: `analysis.py`, `saved_analysis.py`,
`check_conversion_outputs.py`, `decoder_metrics.py`, `metric_tasks.py`, and
`decoder_inventory.py`. They share the common evidence/execution APIs; obsolete
campaign, drain and duplicate queue-serving entry points are removed. Historical
campaigns are reproduced from their recorded Git revisions, not fallback runners.
The independent TS allocation wrapper remains `run_typed_ts_array.sbatch`.

For structural changes run the README smoke checks and affected module tests once.
Do not rerun a dataset benchmark to validate code movement.
