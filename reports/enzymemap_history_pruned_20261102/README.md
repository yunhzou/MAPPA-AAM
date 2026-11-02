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

At **2026-11-02 00:52 Toronto time (2026-11-02 05:52 UTC)**, published ledgers
across all 40 batches recorded:

| Stage | Complete directions / 95,948 | Percent |
|---|---:|---:|
| Matching | 83,956 | 87.50% |
| Family reduction | 83,956 | 87.50% |
| Decoding | 0 | 0% |
| Reference verification | 0 | 0% |

Thirty-five batches have complete matching and reduction records. The audit
checked all 95,948 direction identities across the 40 disjoint inputs, each
completed stage's exact membership, input hashes, scientific settings and the
required source/engine identity. No scientific direction failures or retired
retries were recorded; the 2 GiB guard remains unchanged. Five batches remain.
This is matching/reduction progress, not final benchmark completion.

One allocation failed during worker startup: no workers became ready and no
scientific tasks started. Its only saved artifacts were the input plan and an
empty cache ownership file. All recorded steps had terminated. Direct inspection
found no remaining job processes on the accessible former nodes; access to one
former node was unavailable after allocation release. Subsequent allocations
started successfully, so the evidence identifies an isolated startup failure,
not a demonstrated algorithm or memory failure. Its underlying scheduler cause
has not been established.

The startup evidence was preserved privately. Only that untouched batch was
resubmitted with the same validated code, inputs, resource guard and budgets;
completed work was preserved. Both decoding prerequisites were repaired to
require the replacement batch and every still-unfinished original matching
batch. All other prerequisites were independently verified complete. At this
observation the replacement was queued; successful worker startup has not yet
been confirmed. The global matching/reduction barrier remains intact.

Raw operational snapshots and exact continuation commands remain private.
