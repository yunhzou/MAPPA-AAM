# Fresh EnzymeMap launch

Launched on DFW on 2026-11-01 using scientific revision
`edf2e8ff94e98e40b1db48e62911050c13c26de6`. The obsolete growth selector
and implementation have been removed. The native source/build identity is
recorded in [launch.json](launch.json).

The fresh manifests cover **47,974 reactions / 95,948 directions**. Each
reaction belongs to one batch, with both directions together. Only immutable
endpoint/reference inputs are reused; no previous results or caches were imported.
Settings remain one seed, branch cap 100, tolerance 1, random seed 42 and the
full cut sweep. Shared reaction cache snapshots are enabled, with 2 GiB per worker.

| Batch | Directions | Matching job | Compute workers | Downstream job |
| --- | ---: | --- | ---: | --- |
| 0 | 11,492 | 21140383 | 98 across 7 nodes | 21140390 |
| 1 | 21,114 | 21140386 | 180 across 2 nodes | 21140391 |
| 2 | 21,114 | 21140387 | 180 across 2 nodes | 21140392 |
| 3 | 21,114 | 21140388 | 180 across 2 nodes | 21140393 |
| 4 | 21,114 | 21140389 | 180 across 2 nodes | 21140394 |

All downstream jobs require successful completion of **all five** matching jobs.
They resume saved results and run deduplication, decoding and verification.
Failed matching or a drained allocation does not release this dependency.
Continuation/retry therefore remains necessary if either occurs.

At **08:38:15 UTC**, batch 0 was running with 112 allocated CPUs, including
coordination overhead. The other four matching jobs were pending. These are
requested capacities, not a claim of 818 concurrent compute workers. DFW's
partition limit is ten nodes per user, and an existing service allocation
already occupied one of those nodes. Smaller pending resource requests did not
produce better scheduling estimates; the submitted requests were retained.

The [startup audit](startup-observation.json) observed **98 distinct cuts running
on 98 workers**, no duplicate active task attempts, and 10,025 completed cuts.
Completed tasks recorded 18,608 CPU seconds of computation over a 331-second
observation span; unfinished work is excluded from this CPU total. Fragment
cache reuse accounted for 961,035 hits in 1,712,419 calls, with no cache I/O errors.
Cache saving consumed 855 worker-seconds, compared with 18,624 worker-seconds
of computation. These initial large inputs do not establish a dataset-wide ETA.

**24 cuts exceeded the 2 GiB memory guard** in this observation. Their exact
task identities, measured peaks and compute/publication phases are preserved
in the audit. Other matching work continues. These are failures requiring
resolution, not completed directions; the downstream dependency must remain
blocked until matching succeeds. No memory allowance or algorithm was changed.

Remote campaign directory:
`/lustre/fsw/portfolios/nemotron/projects/nemotron_omni_vision/users/yunhengz/mappa/enzymemap_current_20261101`

Each batch publishes `result/progress.json`, `result/ledger.sqlite` and aggregated
bundles. During the initial admission window, the live queue can be ahead of the
published progress checkpoint. This audit read the running coordinator's queue.
Worker registration totals include replacement processes after memory-limit exits;
they must not be interpreted as current concurrency.
