# Release discarded matching histories

The two reproduced failures now finish under the original **2 GiB** guard on the
Mac. Three fresh-process runs per cut and version used the exact saved inputs,
seed orders and scientific settings. Versions ran sequentially, alternating after
the initial repaired run, with no concurrent tests. Time includes matching,
symmetry finalization and saving; memory is the largest sampled RSS across the
three runs.

| Forward cut, reaction 37270 | Previous time, median | Repaired time, median | Previous peak RSS | Repaired peak RSS |
|---|---:|---:|---:|---:|
| 432 | 20.32 s | **12.38 s** | 4.24 GiB | **0.43 GiB** |
| 410 | 13.38 s | **9.52 s** | 2.30 GiB | **0.34 GiB** |

Both repaired cuts succeeded in **3/3 runs**. These points use about **85–90% less
peak memory** and **29–39% less total time**. Median saving time fell from 2.33 s
to 0.009 s for cut 432 and from 1.09 s to 0.010 s for cut 410. The old native
binary used an 8 GiB diagnostic ceiling so the comparison could finish.

## What changed

Ordinary matching releases histories that no longer lead to a live branch after
each seed admission that discards branches. It retains every incoming path at a
join. Compaction moves retained records and reindexes their references; it does
not change branch keys, admission order, growth, completion, caps or random draws.
Repeated deferred-edge lists also share storage during native export.

Cap evidence is stored without rejected child payloads. Equal cap descriptors
share one summary with an occurrence count; attempted branch count and limit
remain separate fields. Search metrics, reduction, serialization and the viewer
understand these summaries. Historical archives retain their original contents.
The explicit shared-policy backend keeps reusable decisions while policies are
active, then releases its working caches and discarded history at completion.

| Cut | Stored states, before → after | Stored transitions, before → after | Cap occurrences preserved | Cap summaries |
|---|---:|---:|---:|---:|
| 432 | 120,139 → **501** | 120,138 → **500** | 111,373 | 176 |
| 410 | 67,124 → **345** | 67,123 → **344** | 60,471 | 119 |

Every retained transition in both repaired results reaches a returned terminal.
The fragment cache remains enabled with its original 64 MiB budget. All cache and
symmetry-work counters match the baseline in every repetition.

## Accuracy and validation

For each cut, a complete saved result from each version was compared after
loading. **All 100 ordered surviving paths and all 100 catalogue families per
cut are identical**, including accumulated mapping/fragment state, accepted and
deferred edges, compressed fragment records, conditioned symmetry generators,
ordered family actions and matching policies. Catalogue branch groups and
weighted stop descriptors also agree. Graph IDs and discarded payloads differ
intentionally; no returned scientific choice was lost.

Validation passed:

- **112 focused tests:** exact Python/native agreement, branch caps, join
  ancestry, payload release, catalogue relations, decoded event choices and
  checkpoint round trips.
- **27 smoke checks.**
- **61 integration tests:** shared execution, decoder behavior, partial mappings
  and hierarchy transport. Fourteen initially hit Mac sandbox restrictions on
  process inspection or local sockets; all fourteen passed when rerun with those
  permissions. The other 47 had already passed.

New tests include a rejected subtree that first contributes a valid joined path,
an independent two-mapping/two-event oracle, repeated discarded-payload release,
all-capped results, and weighted cap evidence across JSON, binary and combined
historical contexts. Tests reuse existing modules. The full suite and the full
EnzymeMap accuracy benchmark were not rerun.

## Reproduction and scope

[Raw phase traces, comparisons and source identity](history-retention-repair.json)
contain all twelve measurements. Input provenance is unchanged from the
[initial diagnosis](mac-diagnosis.md). Measurements use the existing
[phase probe](../native_transport_20261031/memory_probe.py) and installed repaired
engine. The previous binary is an ignored benchmark artifact, selected only by
the probe's existing engine override; it is not a production fallback.
Both binaries run through the same current Python probe and record classes;
the previous binary retains its original native history behavior.

The worker memory guard is simulated by the probe's 20 ms RSS sampler; separate
execution tests cover the production guard. Shared cache snapshot I/O was
omitted in both versions, as in the initial diagnosis. No cluster throughput
claim follows from these local points. DFW still had an active batch at validation
time, so its running checkout and jobs were left unchanged.
