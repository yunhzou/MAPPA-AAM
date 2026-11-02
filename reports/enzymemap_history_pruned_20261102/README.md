# EnzymeMap rerun with discarded-history release

This fresh benchmark uses revision `69ac5a91d589cd6b253482488c51e6856f163af8`,
including the validated [history-retention repair](../enzymemap_full_20261102/history-retention-repair.md).
The Linux native build passed all 27 smoke tests.

The scope is **47,974 reactions / 95,948 directions**, in 40 disjoint batches.
Only immutable endpoint/reference inputs were imported; previous results and
cache snapshots were not reused.

Scientific settings: one seed, random seed 42, branch cap 100, matching tolerance
1 and the full cut sweep. Verification queries exact annotated heavy-atom
indices, without endpoint-symmetry relaxation, with a one-second per-family
limit. Shared reaction caches are enabled. Each compute worker has a **2 GiB**
memory guard. Matching and reduction must complete across the full campaign
before decoding begins. Continuations preserve completed work, task ownership
and consumed budgets.

The checked source/engine identity is
`3b4b2c193c3e122b846bbf7966397ac8155aab9b61627c9e8acf3690ce2e5bc0`.

## Recorded early validation

At the recorded startup check, 428 directions had completed matching with no
failed directions. Execution records confirmed overlapping tasks; concurrency
alone does not establish continuous CPU utilization.

Former memory failure **reaction 37270** completed all **665 cuts in each
direction**, with 1,330 successful attempts and no retries. Recorded peak
process-tree RSS was 0.283 GiB forward and 0.267 GiB reverse; the largest
process-lifetime RSS was 0.293 GiB. The audit contains 12,969 memory samples
under the unchanged 2 GiB guard. These measurements cover this reaction, not
future decoding or the whole dataset. The directions recorded 109,875 fragment
cache hits and 931 task-boundary snapshot loads.

The benchmark remains in progress. This report does not claim final reference
recovery, decoding completion or completion time.

Public scientific projections: [inputs](manifest.json), [run provenance](launch.json),
[progress](startup-observation.json), [parallel execution and cache measurements](parallel-observation.json),
and [all-cut audit for reaction 37270](previous-failure-observation.json).
Private site configuration, original operational records and continuation
instructions are retained outside the published repository. See the
[portable runner](../../bench/README.md) for execution interfaces.

## Latest monitored progress

At **2026-11-01 23:35 Toronto time (2026-11-02 04:35 UTC)**, published ledgers
across all 40 batches recorded:

| Stage | Complete directions / 95,948 | Percent |
|---|---:|---:|
| Matching | 60,060 | 62.60% |
| Family reduction | 57,568 | 60.00% |
| Decoding | 0 | 0% |
| Reference verification | 0 | 0% |

No direction failures or retired retries were recorded. All started batches
reported the required source/engine identity, the same scientific configuration
and the unchanged 2 GiB guard. Matching advanced by 3,441 directions since
04:25:19 UTC, approximately 20,375 directions per hour over that 10-minute,
8-second interval.
This is observed matching throughput, not an end-to-end completion estimate.

The scheduler/log check found 24 successfully completed matching
allocations, two running allocations, no failed or drained allocations and no
errors in the matching logs. Both finish arrays still depend on both matching
arrays; the global barrier remains intact. No intervention was needed.

The two running pools registered 255 workers, all with retained execution
records. Their latest queue snapshots contained 2,949 and 840 completed attempts,
with peak overlaps of 148 and 83 retained attempts during each snapshot's
preceding minute. Both listed no running tasks at their captured instants;
neither contained duplicate running task ownership. Published snapshots can lag,
and retired attempts leave the active queue. Worker registration and task overlap
do not establish continuous CPU utilization. Raw operational snapshots remain
private.
