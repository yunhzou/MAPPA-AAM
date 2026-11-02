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

At **2026-11-02 03:06 Toronto time (2026-11-02 08:06 UTC)**, published ledgers
across all 40 batches recorded:

| Stage | Complete directions / 95,948 | Percent |
|---|---:|---:|
| Matching | 95,948 | 100.00% |
| Family reduction | 95,948 | 100.00% |
| Decoding | 85,365 | 88.97% |
| Reference verification | 81,927 | 85.39% |

Matching and reduction are complete for the full declared scope. The audit
checked all direction identities, disjoint batch membership, input hashes and
exact completed-stage membership. No matching or reduction failures occurred.
All started batches retain the required source/engine identity, scientific
settings and 2 GiB guard. The final matching allocation ended at **06:38:01 UTC**;
the first decoding allocation started at **06:38:37 UTC**, confirming the global
barrier was respected.

The 81,927 completed reference checks contain **79,688 recovered**, **2,236 not
recovered** and **three unknown** results. These are partial direction-level
results, not final reaction-level benchmark accuracy. A completed verification
task is not necessarily a recovered reference.

Decoding recorded **70 memory-limit failures**. Among the first five failures
inspected in detail, four completed computation and crossed the guard while
publishing their outputs: 13967 reverse, 03861 reverse, 03856 forward and
09159 reverse. Their sampled publication peaks were 2.0115, 2.0905, 2.0906 and
2.1046 GiB respectively. Reaction 37270 forward crossed the guard during
computation after approximately 407 seconds. Its saved metrics do not distinguish
catalogue loading from later decoder operations. The other failures have not yet
received this detailed phase audit. These failures are not evidence of an absent
mapping. No decoding retries, budget resets or memory-guard changes have been made.

Runner revision `9f7d850` keeps external reference verification in the same pool
after terminal decoding failures. Verification reads completed catalogues;
decoding errors remain recorded and still make the overall run unsuccessful.
The scientific source/engine identity is unchanged. This orchestration change
passed all 19 affected campaign tests. The earlier explicit verification-only
continuation (`a6224c5`, 13 affected tests followed by five queue-guard checks)
completed its six inactive batches. A second continuation completed verification
for 11 more inactive batches, covering 26,390 directions without retrying decoding.
Both continuations preserve the original catalogues, decoding failures and
consumed budgets.

Task timing showed prolonged low occupancy near the end of the original large
pools. The 20 remaining unstarted batches, covering 47,974 directions, were
therefore moved to smaller persistent pools. Only their pending allocations were
held and replaced; active scientific work continued. Before replacement, every
batch had complete matching and reduction, no decoding or verification records,
an empty saved queue and released ownership. The replacement retains the exact
input plans, original windows and budgets, and 2 GiB guard.

All 40 batches have now started decoding. The latest snapshots show
**155 running decoding tasks and 34 running verification tasks across eight
pools**, with no duplicate task ownership. The first eight replacement pools
registered their workers in 2.27–3.27 seconds.
All pre-existing scientific ledger rows remain exactly identical across the
20 migrated batches and 17 verification-only batches; all migrated input and plan
hashes are unchanged. No retired retries were recorded. Verification completion
does not mark failed decoding or the full benchmark complete.

With the user's explicit approval, the first five failed cases' saved catalogues
and necessary endpoint inputs were copied into one ignored local diagnostic
bundle. All ten saved byte-range checksums were verified. These private artifacts
are not committed and support the separate local diagnosis.

The earlier matching-startup failure was recovered by resubmitting only its
untouched batch. The replacement registered 155 workers in 3.85 seconds and
completed all 2,398 directions with an identical input plan; no completed
scientific work was rerun. Its original startup evidence remains private.

The latest full worker snapshot found no duplicate running task ownership.
Saved snapshots describe task activity, not continuous CPU utilization. Raw
operational records, saved failure metrics and exact continuation commands
remain private.

## Local decoder memory diagnosis

The five inspected failures were replayed sequentially on the Mac from their
saved Linux catalogues, using the unchanged production decoding operation and
bundle writer. All five completed under a local 16 GiB diagnostic guard. No
matching was rerun, no core code changed, and the campaign guard was unchanged.

RSS measurements are in GiB. The first four columns are samples at operation
boundaries; publication is the maximum of periodic and boundary samples.
“Before record” follows inventory construction and initial symmetry-response
checks. Building the record then evaluates the lazy local tables.

| Reaction/direction | Inventory built | Before record | Record built | Adapter returned | Publication peak |
|---|---:|---:|---:|---:|---:|
| 13967 reverse | 0.685 | 0.696 | 2.140 | 2.128 | 2.617 |
| 37270 forward | 1.074 | 1.105 | 2.567 | 2.534 | 3.023 |
| 03861 reverse | 0.709 | 0.715 | 2.002 | 1.979 | 2.468 |
| 03856 forward | 0.708 | 0.715 | 2.044 | 2.022 | 2.518 |
| 09159 reverse | 0.805 | 0.812 | 1.930 | 1.909 | 2.407 |

Each C Pickler reported **256 MiB** of internal storage after serialization,
while the compressed outputs were only 12.2–26.5 MiB. Publication added roughly
0.49–0.50 GiB above adapter-return RSS. The writer streams protocol-5 pickle
through gzip into an 8 MiB spool that spills to disk; it does not build a complete
serialized byte buffer in RAM. Pickler storage includes its reference memo and
buffers. This measurement does not attribute the entire RSS increase to that
storage or capture every transient allocation.

The five outputs contain 8,309–25,448 family tables but only 16–17 shared response
records. In the four baselines with field counts, their per-table dictionaries
contain **5.49–7.13 million atom-domain entries**, and equally many atom-action
entries. These are entry references, not unique-object counts; structural sharing
has not been quantified. The late table/record growth is therefore a separate
problem from the initial unique-fragment inventory. In particular, 37270 forward
crossed 2 GiB during this stage in a fresh worker, before publication.

The endpoint-pair cache was small in the three baselines where it was counted:
3,228 pairs for each of 03861 reverse and 03856 forward, and 4,252 for 09159
reverse; the largest individual projections contained 8, 8 and 12 pairs.
Those measurements do not support large pair-orbit expansion as the explanation
for these three cases. Counts were not collected for the other two baselines.

A further 13967 reverse replay released the prepared decoder session immediately
after the operation returned, before publishing its unchanged output. The
published bytes were identical: 18,539,804 bytes, SHA-256
`9f6f7ec0cdf236b0c4554205c637ac8d5f6ec10e83deeb4082289b91abdbb14d`.
Publication still peaked at **2.553 GiB**, and session release itself did not
reduce RSS at that boundary (1.792 to 1.817 GiB). This ownership change alone
did not resolve the memory failure.

All timings are diagnostic, not performance comparisons. Separate local recovery
workers overlapped the fifth baseline and the lifetime experiment. Their memory
residency differed, so the small cross-run peak difference is not a measured
causal saving. RSS drops after releasing the output or session also reflect
allocator and operating-system behavior; they are not exact cache-owned byte
counts. All six probes used fresh workers. Persistent-worker overlap between an
old session and the next catalogue remains unmeasured, but was not required to
reproduce these intrinsic peaks. The measurements identify compute and
publication costs; they do not establish a repaired decoder or final accuracy.
