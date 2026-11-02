# Full EnzymeMap launch after Golden validation

The fresh campaign covers **47,974 reactions / 95,948 directions** in 40 disjoint
batches. It uses revision `0f86f117d5705da31c8bc022419c6b80b3add761`, with the same
scientific code and native engine as the successful
[Golden rerun](../golden_retest_20261101/README.md). No algorithm changes were made.

Only immutable prepared inputs are reused. Previous results and cache snapshots
were not imported. Both directions of each reaction stay together; matching uses
one seed, random seed 42, cap 100, tolerance 1 and the full cut sweep. Shared
reaction caches remain enabled. Verification uses exact annotated heavy-atom
indices with a 1-second per-family query limit and no endpoint-symmetry relaxation.

## Observed startup

The first pool began after **31 seconds** of waiting. Startup took **3.27 seconds**,
and all **168 compute workers** registered successfully with a **2 GiB** guard.
This is an observed startup, not a general scheduling-delay or capacity claim.

Matching and reduction cover all 40 disjoint batches. Decoding and verification
require successful completion of the entire matching/reduction pass.

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
[Public launch projection](launch.json),
[public worker observation](startup-observation.json), and
[failure counts and audit references](early-failures.json).

Subsequent [Mac point diagnosis](mac-diagnosis.md) reproduced one compute failure
and one publication failure in fresh processes. Repeated search-history storage
and transport explain substantial memory use; the cluster settings remain
unchanged. This local diagnosis does not update the startup progress snapshot.

The subsequent [history-retention repair](history-retention-repair.md) preserves
both sampled cuts' returned mappings and cap evidence while reducing peak memory
to 0.43/0.34 GiB. The report records controlled timings and validation; it does
not claim deployment into the then-running campaign.

Raw operational records and storage locations are retained privately.
