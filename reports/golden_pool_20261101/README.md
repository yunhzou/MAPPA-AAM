# Golden persistent-pool run, 2026-11-01

**Complete:** all 3,702 directions now also have the current **v4 decoded output**.
The decoding-only pass finished in **3 minutes 45 seconds**, including staggered
allocation starts, with no failures and **1.03 GiB peak worker RSS** under the
unchanged 2 GiB limit. All saved output checksums were verified. Matching,
deduplication and reference-verification records were reused unchanged.
[v4 run report](v4-decoding.md) · [combined decoded manifest](v4-decoded.json)

The original run completed all four stages. Its unique-fragment refactor recovered
the last failed decoded output at the unchanged 2 GiB worker guard, with
**1.14 GiB peak worker RSS**. Completed scientific results were preserved.
[Refactor and equivalence evidence](decoder-refactor.md) ·
[Final ledger counts and recovery audit](decoder-refactor-completion.json)

All **3,702 directions completed matching** at 2 GiB per worker, with no matching
failures. The five allocations ran from 05:00:27 to 05:06:40 UTC: **6 minutes
13 seconds**, including their staggered starts. They processed **177,552 cut
tasks**. Peak observed matching-worker RSS was **414.2 MiB** and no cache errors
were recorded. [matching-results.json](matching-results.json) contains the audit
totals and Slurm accounting.

Every 100-worker pool used all its workers and reached 100 overlapping tasks.
Two pools ran concurrently. Average worker CPU consumption per pool was only
**16.7–22.0 cores**, including publication, measured over each pool's task span.
Worker concurrency therefore does not establish full CPU saturation. These
measurements alone do not separate dispatch, shared-storage and task-size costs.

The initial pass saved all **3,702 deduplicated catalogues** and **3,701 decoded outputs**.
Golden **1561 forward** exceeded its guard while publishing its decoded
record. In job 21135190, its compute peak was 1,752,666,112 bytes (1.63 GiB),
and publication reached 2,186,715,136 bytes (2.04 GiB). A fresh-process retry,
job 21135465, reproduced the failure: 1,752,473,600 bytes during computation and
2,169,737,216 bytes during publication. Both guards remained 2,147,483,648 bytes.
This isolates the failure to publication even without earlier tasks in the worker;
it does not identify which serialization allocation caused the increase.
The subsequent [Mac memory diagnostic](decoder-memory.md) traces the retained
family tables, million bond-term records and serializer bookkeeping responsible
for the large footprint. Job **21136551** subsequently completed this direction
with the refactored decoder and saved its output. The original failure and fresh
retry audits remain available; no worker memory limit was increased.

The normal runner stopped batch 1 before verification because of that failure.
Job 21135583 completed its reference queries independently from the saved complete
catalogues, as allowed by the scientific stage design. It does not rerun matching,
alter decoding or mark the failed decoded output complete. Its recorded command
and the fresh-process retry are retained with the run evidence. No worker memory
limit was increased.

Reference verification completed for all 3,702 directions: **3,617 recovered**,
85 not recovered, no unknowns or operational verification failures. Counting a
reaction recovered when either direction succeeds gives **1,830 / 1,851 (98.87%)**.
The remaining 21 annotations were not found in the retained capped search results;
this is not a claim that no valid matching exists. See
[results-summary.json](results-summary.json) and the complete
[per-reaction outcomes](reference-results.json).

The last ordinary downstream allocation ended at 05:16:31 UTC, **16 minutes
4 seconds** after matching started. Independent verification of the blocked
batch ended at 05:21:38 UTC, **21 minutes 11 seconds** after that start, including
investigation and the fresh-process retry. At that point one decoded output was
still missing. The final refactor retry completed it in a 1 minute 45 second
allocation; the intervening investigation and development time must not be
reported as fresh benchmark execution time. All allocations have exited.
The stage audit records sampled RSS separately from process-lifetime peak RSS;
the latter can include earlier stages in a persistent worker.

The user authorized a fresh high-parallelism Golden run with **2 GiB per worker**,
replacing the earlier 4 GiB resource choice for this run. Scientific limits remain
one seed, full bidirectional cut sweeps, branch cap 100, matching tolerance 1,
graph/cut floor 0.2 and random seed 42. Existence growth is explicitly selected.

All **1,851 original Golden records / 3,702 directions** are included. Each
prepared input was checked against the original per-record SHA256 manifest
(`652285008cea1d10a21ae80e476e0f6c0ff17df831f546dee426d1fd0306442b`).
Endpoints are copied unchanged from the trusted archived JSON inputs. No old
matching implementation or scientific result is executed or reused. Input and
reference records occupy one aggregate bundle, with five small batch indexes.

Golden verification retains the original complete heavy-atom relation, including
unmatched atoms, modulo chemically colored endpoint automorphisms. Generators
are prepared from the archived element/charge/isotope/H/stereo colors and bond
colors, then lifted to original atom indices. Reference data enters verification
only. This differs from EnzymeMap's exact-index verification policy; their
coverage percentages must not be conflated.

The existing campaign runner owns all stage orchestration. `bench/golden.py`
only prepares verified inputs and provides the Golden reference query using
MAPPA's existing correlated-family solver. The runner now accepts direct prepared
input references and a manifest-declared verification operation/policy. Core
matching, hierarchy caps, cache logic, deduplication and decoding are unchanged.
Eight affected tests passed in 10.41 seconds, including independent exhaustive
small colored-graph permutations and complete pipeline checks of unmatched-atom
semantics. The Slurm script passed shell syntax checking.

DFW `cpu_short` permits ten nodes per user, with 96 CPUs per node. The initial
two-node / 100-worker batch requests could not start: most spare CPUs belonged to
nodes whose memory was already reserved. A revised batch spreads 100 workers over
four nodes, requesting 28 CPUs and 54 GiB per node, with 25 compute workers per
node. Each worker remains guarded at 2 GiB; each node has 4 GiB above the sum of
worker limits. Five disjoint batches request 500 workers in total, but the
ten-node quota permits at most two revised batches (200 workers) simultaneously,
subject to other jobs. The allocation limit is one hour; admission stops after
50 minutes and drains active work. Requests are not proof of running parallelism.

Case IDs modulo five establish disjoint ownership. Each allocation starts one
persistent pool with one `srun`; tasks never invoke Slurm. Shared reaction caches
retain the permanent baseline and latest versions. Matching is submitted with
`--stop-after matching`; a dependent array resumes the same ledgers only after
the entire matching array succeeds, then performs global per-direction
deduplication, decoding and verification. Resource failures remain explicit and
block the next pass; no automatic memory increase or silent retry is enabled.

Launch and live-progress evidence are recorded alongside this report. Aggregated
inputs/results live under the DFW MAPPA root in `golden_pool_20261101/`.

The original Golden case 0 completed all four stages locally in both directions
in 15.05 seconds with two workers at 2 GiB each. Both reference queries recovered
the annotation; [local-preflight.json](local-preflight.json) records the checks.
The [input manifest](inputs.json) retains every input/reference hash.
