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
