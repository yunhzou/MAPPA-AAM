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

At **2026-11-02 02:18 Toronto time (2026-11-02 07:18 UTC)**, published ledgers
across all 40 batches recorded:

| Stage | Complete directions / 95,948 | Percent |
|---|---:|---:|
| Matching | 95,948 | 100.00% |
| Family reduction | 95,948 | 100.00% |
| Decoding | 33,553 | 34.97% |
| Reference verification | 9,588 | 9.99% |

Matching and reduction are complete for the full declared scope. The audit
checked all direction identities, disjoint batch membership, input hashes and
exact completed-stage membership. No matching or reduction failures occurred.
All started batches retain the required source/engine identity, scientific
settings and 2 GiB guard. The final matching allocation ended at **06:38:01 UTC**;
the first decoding allocation started at **06:38:37 UTC**, confirming the global
barrier was respected.

The 9,588 completed reference checks contain **9,329 recovered**, **259 not
recovered** and **zero unknown** results. These are partial direction-level
results, not final reaction-level benchmark accuracy. A completed verification
task is not necessarily a recovered reference.

Decoding recorded **23 memory-limit failures**. Among the first five failures
inspected in detail, four completed computation and crossed the guard while
publishing their outputs: 13967 reverse, 03861 reverse, 03856 forward and
09159 reverse. Their sampled publication peaks were 2.0115, 2.0905, 2.0906 and
2.1046 GiB respectively. Reaction 37270 forward crossed the guard during
computation after approximately 407 seconds. Its saved metrics do not distinguish
catalogue loading from later decoder operations. The other failures have not yet
received this detailed phase audit. These failures are not evidence of an absent
mapping. No decoding retries, budget resets or memory-guard changes have been made.

The explicit verification-only continuation uses runner revision `a6224c5`;
the scientific source/engine identity is unchanged. Validation comprised
13 affected tests, followed by five focused checks after adding the guard against
unarchived queued work. Before submission, six inactive batches were checked for
complete catalogues, terminal decoding attempts, released ownership and empty
saved queues. The continuation reads their original catalogues without retrying
decoding and preserves their original windows and time budgets. It registered
all **31 workers in 2.115 seconds**, and its first batch completed verification.
All pre-existing matching, reduction and decoding ledger rows—including failure
records and audit references—remain exactly identical across those six batches.
No retired retries were recorded. Verification completion does not mark failed
decoding or the full benchmark complete.

With the user's explicit approval, the first five failed cases' saved catalogues
and necessary endpoint inputs were copied into one ignored local diagnostic
bundle. All ten saved byte-range checksums were verified. These private artifacts
are not committed; local profiling remains pending.

The earlier matching-startup failure was recovered by resubmitting only its
untouched batch. The replacement registered 155 workers in 3.85 seconds and
completed all 2,398 directions with an identical input plan; no completed
scientific work was rerun. Its original startup evidence remains private.

The latest full worker snapshot found no duplicate running task ownership.
Saved snapshots describe task activity, not continuous CPU utilization. Raw
operational records, saved failure metrics and exact continuation commands
remain private.
