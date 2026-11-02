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

The original pass used source/engine identity
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

The original declared attempt scope is terminal. Its outcomes and resource
measurements are recorded below; the separate decoding repair and recovery
preserve that history.

Public scientific projections: [inputs](manifest.json), [run provenance](launch.json),
[progress](startup-observation.json), [parallel execution and cache measurements](parallel-observation.json),
and [all-cut audit for reaction 37270](previous-failure-observation.json).
Private site configuration, original operational records and continuation
instructions are retained outside the published repository. See the
[portable runner](../../bench/README.md) for execution interfaces.

## Original pass results

The last original-pass allocation ended at **2026-11-02 03:18:25 Toronto time (08:18:25 UTC)**.
The independent final audit finished at 08:25:21 UTC. All 47,974 reactions and
95,948 directions have terminal records for every stage. Matching, reduction and
reference verification completed throughout the declared scope. **84 directions
lacked successful decoding outputs at that audit**. Their later recovery is
recorded separately below.

| Stage | Complete directions / 95,948 | Errors | Successful completion |
|---|---:|---:|---:|
| Matching | 95,948 | 0 | 100.00% |
| Family reduction | 95,948 | 0 | 100.00% |
| Decoding | 95,864 | 84 | 99.91% |
| Reference verification | 95,948 | 0 | 100.00% |

A reaction is recovered when its exact annotated heavy-atom reference is
recovered in **at least one search direction**. By that rule, **47,500 / 47,974
reactions were recovered (99.01%)**. Both directions recovered the reference for
45,792 reactions (95.45%). Neither direction recovered it for 472 reactions;
two additional reactions remain unresolved because a direction returned unknown.
These outcomes describe the retained catalogues under the stated seed, cap and
query settings.

Across directions, verification returned **93,292 recovered (97.23%)**, **2,652
not recovered**, and **four unknown**. Unknown is unresolved, never counted as
recovered or as proven absence. The unknown directions are 23473 forward,
23478 reverse, 23479 forward and 27963 forward. A completed verification task
therefore does not necessarily mean its reference was recovered. Verification
queries the completed catalogues independently of decoding.

All 84 decoding errors were **memory-limit failures**, affecting 67 reactions:
**42 during computation and 42 during publication**. The computation-phase
metrics do not distinguish catalogue loading from later decoding operations.
Publication-phase failures occurred after the decoding operation returned.
The full case list, phase and sampled memory measurements are in
[completion.json](completion.json). The original pass made no decoding retries, budget resets or
memory-guard changes. The separate five-case local diagnosis
below explains the measured memory growth without changing this campaign.

## Time and resource use

The campaign took **6 h 14 min 13 s**, measured from the first allocation start
at 02:04:12 UTC through the final allocation end at 08:18:25 UTC. This includes
startup, idle periods and orchestration changes, but excludes input preparation
and the initial queue wait. Matching and reduction finished at 06:38:01 UTC;
the first decoding allocation started at 06:38:37 UTC, confirming the global
barrier.

| Allocation purpose | Reserved CPU hours | Accounted CPU hours used |
|---|---:|---:|
| Matching and reduction | 1,141.64 | 221.67 |
| Decoding and verification | 400.56 | 214.05 |
| Verification-only continuations | 15.64 | 8.10 |
| **Total** | **1,557.85** | **443.81** |

Reserved CPU time is allocated CPU count multiplied by elapsed time.
Accounted CPU usage is the scheduler's cumulative CPU total, counted once per
allocation. It includes workers, coordinators and supervisors; it is not an
algorithm-only measurement or a separate computation/publication breakdown.
No complete campaign-wide task-phase CPU aggregation was performed. Accounted
CPU usage was **28.49% of reserved CPU time**, so these results do not demonstrate
high overall utilization. Parent totals were cross-checked against their steps
within reporting precision; step totals were not added again.

The resource audit includes 98 used allocations, including the startup recovery
and verification continuations. Twenty cancelled pending allocations had zero
elapsed time and contributed no reserved CPU time. Unrelated jobs are excluded.
The observed execution snapshots demonstrated overlapping work and no duplicate
running ownership; snapshots are not continuous CPU-utilization measurements.

## Completion audit and continuations

Every stage's direction identities exactly match its immutable input index, and
all 40 batches are disjoint. Every verification artifact passed its saved
checksum. The scientific source/engine identity, configuration and **2 GiB guard**
remained unchanged. At the original-pass completion audit, all saved queues were drained, all
coordinator ownership locks were released, and no original-pass allocations
remained active or pending. The separate decoding recovery is recorded below.

Runner revision `9f7d850`, validated by 19 affected campaign tests, allowed
external verification to continue in the same pool after terminal decoding
failures. The failures remain recorded and still cause an unsuccessful overall
run status. Earlier explicit verification-only continuations covered 17 inactive
batches without retrying decoding. Twenty unstarted batches were moved to smaller
persistent pools after confirming empty queues and released ownership. Their
input plans, windows and consumed budgets were preserved.

All pre-existing scientific ledger rows remain identical across the 20 migrated
batches and 17 verification-only batches. There are **zero recorded retry
entries**. The earlier matching-startup failure was recovered by resubmitting
only its untouched batch; no completed scientific work was rerun. Its original
evidence remains private.

[The completion record](completion.json) contains the final counts, input hashes,
reference outcomes, unknown directions, failure list, resource boundaries and
audit checks. Private site configuration and raw operational records are not
included. The first five failed cases' catalogues and necessary endpoint inputs
were copied to one ignored local diagnostic bundle with the user's explicit
approval; all ten byte-range checksums were verified.

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
counts. All six probes used fresh workers. These baseline probes did not measure persistent-worker overlap between an
old session and the next catalogue. That overlap was not required to reproduce
the intrinsic peaks. The baseline measurements identify compute and publication
costs; follow-up repair and validation are recorded below.


## Decoder storage repair and recovery

The first repair (`b8c923e`) keeps the saved mapping relation and the decoded-event
schema. It changes object ownership and storage sharing:

- Export releases each temporary family/partition table after recording it.
  Explicit table queries still retain the tables they request.
- Arrays used to evaluate a requested joint choice are built only for evaluation;
  exporting local responses does not build those unused arrays.
- Equal saved atom-domain and action metadata share storage, including immutable
  tuple values. Families, correlated actions, constraints, cuts and alternative
  partitions remain distinct and in their original order.
- Catalogue loading restores shared immutable mappings, constraints and ordered
  actions instead of allocating equal copies for every family.

A fresh local replay of **37270 forward** completed and published under the
unchanged **2 GiB** guard. Peak sampled RSS fell from **3.023 GiB** in the baseline
to **1.878 GiB**. Its output retained 8,309 family tables, 16 response records and
604,739 term references. The compressed output is 11,014,645 bytes. This one-case
measurement is not a general bound on decoding memory. Persistent-worker
validation and recovery of the full failed scope are separate checks.

For **13967 reverse**, the complete repaired typed result and saved record
matched the retained baseline exactly, including every nested container's type,
values and order: 25,448 families/tables and 16 responses. Structural sharing
changes serialized bytes, so equality was checked on restored values rather
than requiring identical byte streams. Independent event-oracle tests cover
correlated actions, preserved constraints, cuts, partitions and unmatched bonds.

The separate **37270 forward** comparison also passed for the complete typed
result and saved record, including nested types, values and ordering across
8,309 families/tables and 16 responses. It reconstructed the pre-repair decoder
methods from Git and used the same saved catalogue and event thresholds. The
comparison held both results under a 6 GiB diagnostic guard; this was an accuracy
check, not the 2 GiB recovery execution or a speed measurement.

Validation passed 27 smoke checks, 51 decoder checks and 26 catalogue/transport
checks. The execution/campaign suite passed 44 checks before the final ownership
registration change; all 24 campaign checks passed again after that change.
These are the affected suites, not a full-suite or matching benchmark rerun.

The Mac persistent-worker check completed the first four former failures under
2 GiB, then reached the guard on 37270 forward during table construction. Weak
references confirmed that prior catalogues, interpretations, outputs and sessions
were released. An isolated collection check found zero unreachable cycles and
reduced RSS by only 19 MiB. Allocator diagnostics nevertheless showed hundreds
of MiB in unused small-object blocks and pools. Thus cold-worker success does
not establish persistent-worker safety. An earlier session-release experiment also exceeded 2 GiB on the second
consecutive case, so it was not applied to production. Linux validation and
explicit recovery results are required before claiming all failures recovered.

The Linux validation then completed **all five cases in one persistent worker**
under the unchanged 2 GiB guard. All five published outputs matched the repaired
Mac outputs byte for byte. The external sampler observed a maximum of
**1.707 GiB**; the operating system's whole-process high-water measurement was
**1.738 GiB**. The five tasks took 24 minutes, with 1,395.51 seconds of compute
CPU and 8.94 seconds of publication CPU. These are instrumented validation
measurements, not a baseline speed comparison. The same 27 scientific smoke
checks passed on Linux. The observed Mac persistent-process memory limit does not imply a Linux
failure, and these five successful tasks do not establish the outcome of the
remaining recovery scope.

Recovery uses the original 84 failed decoding tasks and saved catalogues. The
maintained runner preserves task identity, attempt history, the 2 GiB guard and
unused portions of the original budgets. Results go into a separate recovery
ledger; original scientific records remain unchanged. One operational
registration and exclusive batch ownership prevent simultaneous or sequential
replays from reusing the same pre-retry budget. No matching, reduction or
reference verification is repeated.

The first recovery round completed **75 of 84 directions** under the unchanged
2 GiB guard. Nine still exceeded the guard: seven during computation and two
during publication. After this round, the combined successful decoding count was
**95,939 / 95,948**; matching, reduction and external verification outcomes were
unchanged.

The round took **16 min 46 s**, with 84 simultaneously claimed tasks on 84
workers. It used 9.62 accounted CPU hours from 25.15 reserved CPU hours (38.25%).
Task overlap is not a claim of sustained CPU saturation. The maximum sampled
RSS among successful tasks was 1.942 GiB. Failed tasks and their memory
measurements remain separately recorded.

The terminal audit verified all 84 input catalogue ranges and all 75 published
outputs. All 121 hashed original campaign files remained unchanged. Every task
retained its original identity, history, settings, guard and consumed budget,
with exactly one new attempt and no additional time grant. All dispatch queues
were drained at the end of this round. At that point, the nine remaining
directions had 2,221–3,008 seconds left from their original budgets.

A further repair (`6974f54`) shares repeated exported fragment-group records. In
the saved 37270 forward result, exact sharing of group dictionaries and component
tuples removed **118.69 MiB** of duplicate containers. The follow-up production
replay completed under 2 GiB with a **1.674 GiB** process peak, compared with
1.878 GiB after the first repair. Pickler storage fell from 64 to 32 MiB. The
complete typed output retained identical values, nested types and ordering;
2,702,578 container pairs were checked. Compute CPU was 367.54 s versus 364.70 s
in the earlier instrumented run, a 0.78% difference, not a controlled speed
comparison.

The additional export change passed 52 decoder tests and 27 smoke checks. The
chained-recovery runner passed 32 EnzymeMap/Golden campaign tests. It accepts the
latest terminal recovery as input, preserves all prior attempts and latest spent
time, and selects only remaining memory failures. Successful outputs stay in
their existing ledgers.

The second round finished at **11:51:45 UTC on 2 November 2026**. It recovered
both remaining publication failures, **06726 forward and 37269 forward**, under
the same 2 GiB guard. Seven directions still exceeded the guard during
computation: 10376, 10379, 10380, 10381, 10382, 10383 and 27963, all forward.
Across both rounds, **77 of 84 original failures are recovered**, bringing
successful decoding to **95,941 / 95,948 directions**. Matching, reduction and
external reference verification outcomes are unchanged.

| Recovered direction | Peak sampled RSS (GiB) | Compute CPU (s) | Publication CPU (s) |
|---|---:|---:|---:|
| 06726 forward | 1.797 | 699.563 | 2.646 |
| 37269 forward | 1.723 | 1,291.782 | 1.910 |

This round took **21 min 53 s**, with all nine workers claiming overlapping
tasks. It used **2.3025 accounted CPU hours from 4.3767 reserved CPU hours
(52.61%)**, including failed work, smoke tests and orchestration. Parent
allocation CPU is counted once; task ownership overlap does not establish
sustained CPU saturation. These runs used different allocations and are not a
controlled speed comparison.

The terminal audit verified the nine input catalogue ranges, both new outputs
and all 75 successful outputs from the first round. The predecessor ledger and
all 121 hashed original campaign files were unchanged. Each direction retained
its two earlier attempts and added exactly one attempt, with no extra time,
memory increase or automatic retry. All recovery queues were drained. The
seven unresolved directions retain **1,665–2,580 seconds** of their original
budgets. Output checksums establish saved-byte integrity; they are separate
from the scientific value comparisons described above.

Detailed repair measurements, comparison scope and recovery status are recorded
in [decoding-recovery.json](decoding-recovery.json).

Further bounded diagnosis found a larger export cost in the remaining cases.
At 12,672 of 32,624 tables for 10383 forward, the saved per-atom domain/action
dictionaries alone occupied **399.73 MiB**. Runtime fragment maps occupied
4.30 MiB. These stopped diagnostic probes produced no scientific outputs and
did not change recovery budgets. Memory after catalogue loading also includes
allocator-retained space; it is not a measurement of live scientific data alone.

The next storage repair removes these derived maps from the saved event record,
using schema **`mappa.decoded_events/v5`**. Runtime table queries retain them.
The original catalogue, ordered actions, constraints, local responses, bindings,
partitions and occupancy boundaries remain saved. No production reader used the
removed export fields; external consumers must recognize the version change.
This deletes redundant export work without changing matching or event rules.

Catalogue restoration also shares equal immutable fragment atom tuples and
fragment-pair tuples across branches, preserving order, provenance and separate
mutable owners. The approved 37270 input contained **112.60 MiB** of duplicate
branch tuple containers; this is a container measurement, not an RSS reduction.
Bounded inspection of 27963 confirmed the same storage issue. Provenance cut
tuples were already shared and were left unchanged.

The changes passed 64 decoder/CLI checks, 14 catalogue checks and 27 smoke checks
(these scopes overlap). The saved local responses were compared with the
independent scalar event oracle across correlated choices, cut boundaries,
partial occupancy and metal thresholds.

The cold production replay of 37270 forward completed under 2 GiB at a
**1.203 GiB** lifetime peak, compared with 1.674 GiB for `6974f54` (**482.56 MiB
lower**). It retained 8,309 tables and 16 responses. Compute CPU was 363.62 s
versus 367.54 s, and publication CPU was 0.24 s versus 0.64 s; these are single
instrumented runs, not a controlled throughput comparison. The saved output
fell from 9,511,318 to 4,093,248 bytes.

The full typed result matched the preceding saved result after changing only
the schema identifier and deleting the two declared fields: **1,174,466 nested
container pairs** had identical values, types and order. The comparison held
both saved outputs under a separate 6 GiB diagnostic guard; production decoding
and recovery retain 2 GiB. This comparison establishes that result's preservation,
not universal equivalence. The seven remaining directions are prepared for a
separate recovery using their latest remaining original budgets.
