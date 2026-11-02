# EnzymeMap rerun with discarded-history release

This fresh benchmark uses revision `69ac5a91d589cd6b253482488c51e6856f163af8`,
including the validated [history-retention repair](../enzymemap_full_20261102/history-retention-repair.md).
DFW was synchronized through Git and its native engine rebuilt on compute job
`21164734`. All 27 smoke tests passed on that Linux build.

The scope remains **47,974 reactions / 95,948 directions**, in 40 disjoint
batches. Only immutable endpoint/reference inputs were imported. Previous
matching results, ledgers and cache snapshots were not reused. The superseded
search array `21162459` was stopped; its saved outputs remain intact.

Scientific settings are unchanged: one seed, random seed 42, branch cap 100,
tolerance 1 and the full cut sweep. Verification queries exact annotated
heavy-atom indices, without endpoint-symmetry relaxation, using the existing
one-second per-family limit. Shared reaction cache snapshots remain enabled.
Every compute worker retains the **2 GiB** memory guard.

## Submitted work

| Batches | Matching and reduction | Decoding and verification | Workers per pool | Allocation per node |
|---|---|---|---:|---|
| 0–23 | `21164779`, one batch at a time | `21164781` | 155 across five nodes | 33 CPUs, 66 GiB, 31 workers |
| 24–39 | `21164780`, one batch at a time | `21164782` | 100 across four nodes | 27 CPUs, 54 GiB, 25 workers |

The two pools request 255 compute workers and 273 Slurm CPUs in total.
Each node reserves 4 GiB beyond its workers' combined memory limits. The
admission windows are 256 and 128 directions respectively. These are allocation
requests; actual startup and progress are recorded separately below.

The pool shapes fit the observed unallocated CPU and memory capacity while
leaving the unrelated service allocation intact. DFW `cpu_short` limits this
user to ten nodes; that service occupies one. Node selection remains Slurm's
responsibility. Availability and priority can change after inspection.

Both finish arrays require `afterok:21164779:21164780`. Matching and reduction
must succeed across the entire campaign before decoding begins. A failure or
graceful retirement blocks this dependency and requires an explicit continuation;
there are no automatic retries or hidden changes to scientific settings.

DFW source/engine identity:
`3b4b2c193c3e122b846bbf7966397ac8155aab9b61627c9e8acf3690ce2e5bc0`.
Every allocation checks this identity before starting its pool.

## Observed startup

Both matching pools started at **02:04:12 UTC**, 155 seconds after submission.
Slurm rounded the CPU requests to 170 and 112: **282 allocated CPUs across nine
nodes**. All **255 compute workers** registered with the expected code identity;
pool startup took 3.84 and 3.55 seconds.

At **02:08:06 UTC**, the published checkpoints contained **428 completed matching
directions and no failed directions**. The pools reported 129 and 82 running
tasks at their respective checkpoint times. Recorded attempts confirm that all
155 and 100 workers had executed work; retained snapshot histories showed peak
overlap of 137 and 95 tasks respectively. These measurements establish parallel
execution, not continuous full CPU utilization. Archived tasks are absent from
the live queue counts, so those counts are not cumulative campaign totals.

The former memory failure, **reaction 37270**, completed all **665 cuts in each
direction**, with 1,330 successful attempts and no retries. Its recorded peak
process-tree RSS was 0.283 GiB forward and 0.267 GiB reverse; the largest process
lifetime RSS was 0.293 GiB. The audits contain 12,969 memory samples and the
unchanged 2 GiB guard. These measurements cover this reaction, not the full
dataset or future decoding. The two directions recorded 109,875 fragment-cache
hits and 931 task-boundary snapshot loads.

The full benchmark is still running. Decoding and verification have not started;
their global dependency remains in place. This report does not claim a final
coverage score or completion time.

Evidence: [input manifest](manifest.json), [launch commands](launch.json),
[startup and progress](startup-observation.json),
[task overlap and cache samples](parallel-observation.json), and
[all-cut audit for reaction 37270](previous-failure-observation.json).
Remote campaign:
`/lustre/fsw/portfolios/nemotron/projects/nemotron_omni_vision/users/yunhengz/mappa/enzymemap_history_pruned_20261102`.
