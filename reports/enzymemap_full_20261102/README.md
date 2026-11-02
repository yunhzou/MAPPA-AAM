# Full EnzymeMap launch after Golden validation

The fresh campaign covers **47,974 reactions / 95,948 directions** in 40 disjoint
batches. It uses revision `0f86f117d5705da31c8bc022419c6b80b3add761`, with the same
scientific code and DFW native engine as the successful
[Golden rerun](../golden_retest_20261101/README.md). No algorithm changes were made.

Only immutable prepared inputs are reused. Previous results and cache snapshots
were not imported. Both directions of each reaction stay together; matching uses
one seed, random seed 42, cap 100, tolerance 1 and the full cut sweep. Shared
reaction caches remain enabled. Verification uses exact annotated heavy-atom
indices with a 1-second per-family query limit and no endpoint-symmetry relaxation.

## Resources actually started

The first job was submitted at **00:47:14 UTC on November 2** and started at
**00:47:45 UTC**, a **31-second wait**. It allocated **192 Slurm CPUs across eight
nodes**, supporting **168 compute workers**: 21 workers, 24 allocated CPUs and
46 GiB per node. Each worker retains its **2 GiB** memory guard. Pool startup took
3.27 seconds, and all 168 workers registered successfully.

The live resource snapshot had spare CPUs but much less unallocated memory.
DFW's `cpu_short` QoS limits this user to ten nodes across jobs, including an
unrelated service allocation. This launch demonstrates prompt admission for the
168-worker shape; it does not establish that 800 or 2,000 workers can start
without delay. Slurm's preliminary test-only estimates were substantially later
than the actual start and are not used as a reliable wait-time prediction.

Initially submitted 180-worker bulk pools required 184 GiB free per node. Those
jobs were verified pending, held and cancelled before doing any work. Their
replacement uses the same eight-node/168-worker shape that actually started.
This avoids waiting for much larger per-node memory reservations. The running
first batch was unchanged. Scheduling can be resized again if capacity changes;
the present bulk array admits one pool at a time.

| Scope | Job | Configuration |
|---|---|---|
| Batch 0 matching and reduction | `21162307_0` | 168 workers; running |
| Batches 1–39 matching and reduction | `21162459` | 168 workers per pool; one pool at a time |
| All decoding and verification | `21162460` | Requires successful completion of both matching/reduction jobs |

## Startup observation and explicit failures

At **00:52:06 UTC**, the published checkpoint reported **323 completed matching
directions**, two failed directions and 160 running cut tasks. The queue's
completed-task count excludes histories already archived, so it is not the
campaign's cumulative cut count. Worker registration is not evidence that every
allocated CPU stays continuously busy.

Both failed directions belong to **reaction 37270**, with 661 atoms per endpoint:

| Direction | Completed cuts | Cuts exceeding 2 GiB | Compute / publication failures |
|---|---:|---:|---:|
| Forward | 644 / 665 | 21 | 18 / 3 |
| Reverse | 658 / 665 | 7 | 6 / 1 |

The failing operation combines matching and symmetry finalization. These audits
do not distinguish the internal compute substep or exclude retained worker state
as a contributor. Four tasks completed computation and crossed the guard while
publishing. No decoder had run. The failures and completed cuts remain recorded;
no memory increase, retry or algorithm change was made. Other matching continues.
Decoding stays blocked until matching and reduction are complete without failed
directions. Drained allocations also require explicit continuation; submission
of the full campaign does not mean it has completed.

Evidence: [manifest](manifest.json), [input checks](preflight.json),
[launch commands and replacement history](launch.json),
[Slurm and worker observation](startup-observation.json), and
[failure counts and audit references](early-failures.json).

Remote campaign:
`/lustre/fsw/portfolios/nemotron/projects/nemotron_omni_vision/users/yunhengz/mappa/enzymemap_full_20261102`.
