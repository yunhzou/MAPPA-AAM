# Real Slurm pool test: EnzymeMap case 446

Job **187503** completed with exit code **0:0** on NRT, partition `cpu`, QoS
`cpu-short`, nodes `cpu-0012` and `cpu-0006`. Tested source: `7548fd2c413756aa7edc0cd8561c309e749d6f67`.
The allocation requested four CPUs per node; two compute workers per node left
capacity for supervision and the coordinator. Slurm accounting shows one batch
submission and one application step (`187503.0`), which stayed alive for all passes.
Total allocation elapsed time was 2:01 including native builds.

## Scientific and execution checks

- All **821 tasks** completed with no failed, lost or retried task.
- The **four compute PIDs** stayed stable throughout the run. The first matching
  pass used only two workers on one node; later matching passes used all four
  across both nodes.
- Complete graphs and decoded records matched direct execution; the one-seed
  decoded record also matched the saved pre-refactor result. Only archive-local
  provenance IDs were normalized for that comparison.
- A three-seed/full-sweep pass with one seed per task exercised explicit seed-batch
  dependencies and matched the complete direct graph and decoded record.
- Reopening the saved catalogue manifest and decoding through the same pool
  reproduced the saved output.
- Exact reference membership recovered **64/64 annotated heavy-atom pairs**.
  Independent dense event counts agreed for **778 fragment partitions**.
- All task artifacts shared **five bundle files**; none was a per-seed/per-cut file.

This run tested BundleStore over the real cross-node connection. RAM transport
has local/network integration coverage but was not exercised in this Slurm job.

## Timing and acceptance

Same reaction, forward direction, full 116-cut sweep, one seed, branch cap 100,
random seed 42, tolerance 1.0. Times include matching, global deduplication and
complete decoded-record preparation/publication, measured sequentially.

| Run | Total seconds |
|---|---:|
| Direct execution, median of three | 4.753 |
| Distributed first pass, including worker startup | 7.729 |
| Distributed second pass | 5.360 |
| Distributed third pass | 4.187 |
| Distributed median of all three | 5.360 |

**Correctness and lifecycle passed. The no-slowdown performance requirement is
not established by this run.** The all-trial median is 12.8% slower than direct
execution for this small case. The two post-startup trials straddle the direct
time; their mean is 4.770 s versus 4.753 s direct. Those two samples are not a
steady-state throughput benchmark and do not justify claiming a speedup.

Task metrics show approximately 0.82–0.88 summed worker-seconds of matching
publication per one-seed pass, plus catalogue loading/publication and coordinator
traffic. These summed worker timings are not additive to wall time because tasks
overlap. No scientific algorithm change was made to hide transport costs.

`summary.json` contains configurations, stage timings, worker identities and
verification results. The original output directory retains the Slurm log,
compressed task metrics, manifests and result bundles. The test driver and batch
script are versioned here; this is one case-study validation, not a benchmark
campaign framework.

## Measured parallelism

This post-run calculation joins coordinator-recorded claim/completion timestamps
with each compute process's CPU timings. Peak overlap includes computation and
publication; CPU seconds divided by the first-claim-to-last-completion interval
independently measures simultaneous CPU consumption. It excludes pool startup
and submission time and is not an end-to-end speedup ratio.

| Matching pass | Workers / nodes used | Peak overlapping tasks | Computation CPU s | Active wall s | Average compute cores |
|---|---:|---:|---:|---:|---:|
| One seed, first | 2 / 1 | 2 | 2.751 | 2.040 | 1.35 |
| One seed, second | 4 / 2 | 4 | 2.618 | 1.530 | 1.71 |
| One seed, third | 4 / 2 | 4 | 2.636 | 1.172 | 2.25 |
| Three seeds, split batches | 4 / 2 | 4 | 9.571 | 3.832 | 2.50 |

The split-seed pass had four tasks in flight for 3.404 seconds of its 3.832-second
active window. This confirms real computation in parallel, with average compute
utilization of about 62% of the four compute workers for that pass. Single-reaction
deduplication and decoding each used one task and one worker. The four-worker
configuration did not make those single-reaction stages parallel.

`parallelism.json` retains the calculations for every pass. The first pass's
limited worker participation was not identified in the original completion/PID
audit; the timestamp analysis above corrects that incomplete characterization.

## Readiness diagnosis: job 189373

A readiness-instrumented rerun on NRT `cpu-short` reproduced the first-pass symptom
without changing scheduling or adding reactions. All timestamps below are from
the coordinator clock, relative to submitting the first 116-cut sweep.

| Event | Seconds |
|---|---:|
| First worker on cpu-0004 ready; 116 runnable cuts | 6.3574 |
| Second worker on cpu-0004 ready; 111 runnable cuts | 6.5015 |
| Last cut in the first sweep completed | 8.3484 |
| Both workers on cpu-0008 ready; no runnable cuts | 9.1156–9.1158 |
| Same sweep submitted after all four were ready | 9.1279 |
| Second sweep's first task claimed | 9.1952 |
| Second sweep's last task completed | 10.3008 |

Each first-node worker claimed its first task about 0.2 ms after reporting ready.
The second-node workers became ready after the entire first sweep had finished.
After all four were ready, the same 116 cuts used all four workers, with task counts
30, 33, 35 and 18. The active matching window was 1.106 s, peak overlap was four,
and summed computation CPU time was 2.464 s. Both retained graphs exactly matched
direct execution; all 232 tasks completed and the job exited 0:0.

The reproduced cause is asynchronous pool startup without a worker-readiness
barrier. `WorkPool.__enter__` returns after launching srun, while services and
compute processes are still starting. The queue dispatches runnable work to the
workers already ready. The first sweep finishes before the late node joins.
There is enough cut parallelism; this observation does not show ready workers
being starved by task dispatch. The original test mixed startup with execution.
A ready-pool API should await all requested workers and report startup separately.
This diagnostic adds telemetry only; it has not yet changed that startup behavior.

`readiness.json` retains worker-ready times, runnable backlog at readiness and all
claim/completion timestamps. `readiness.py` contains the controlled comparison.
