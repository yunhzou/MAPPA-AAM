# Relaxed fragment integration

The completed fragment is the candidate: one source/target atom-set pair, its
complete discovery constraints, and both anchor-aware relaxed self relationships.
The hierarchy, catalogue, decoder and reference query share this record. Retired
candidate pools, exact-global finalization and duplicate fragment wrappers are
removed. Search, reduction and decoding remain separate public operations.

The retained local relation is `beta composed with d composed with alpha inverse`.
Discovery `d` satisfies accepted growth constraints; each self rearrangement is
checked against its own original fragment. This is one expansion, not repeated
closure. Decoding uses unique fragments, correlated choices and original WBOs,
including cut, inter-fragment and mapped/unmapped boundaries. Every bond is
counted once. No discovery-tolerance filter removes the expanded event choices.

Global acceleration uses a certified subgroup that preserves current discovery
constraints and self tolerance under composition. One witness supplies an orbit
of distinct target atom sets. The search excludes that complete orbit before
pruning interchangeable assignments. It neither merges different atom sets nor
enumerates their internal bijections.

## Validation

The frozen implementation passed **824/824 repository tests** in 109.20 seconds
on an Apple Silicon Mac, using `python -m pytest -q --tb=short`. This includes
native/Python comparisons, independent mapping and event oracles, hierarchy caps,
cache persistence/exchange, RAM and bundle work pools, checkpoint recovery and
reference queries. Two subsequently added Python/native tests also passed
(0.27 seconds) without production-source changes. They independently verify all
18 target sets for a six-atom path on two connected C6 rings: two whole-ring sets
and 16 crossing the bridge, including every discovery mapping within each set.
Local worker tests ran with localhost sockets and process inspection available.
No cluster benchmark was launched.

Critical checks cover one C6 placement, two externally distinguished C6
placements, overlapping target sets, independently proven absence, continuous
weights, non-transitive rearrangements, anchors and complete same-set discovery
relations. The four-atom relation oracle retains all four discovery mappings and
all six mappings after the approved two-sided expansion.

For EnzymeMap 10383, every returned family witness and partition passed an
independent scalar original-WBO event check: **100/100 families, 100/100
partitions**. All observed witnesses have 59 events; this is not a minimum-event
proof or a reference-recovery measurement.

## Local performance

These are single-cut diagnostics with one seed, cap 100, discovery/self tolerance
1 and random seed 42. The comparison isolates the missing global pruning within
the new relation implementation. It does not assert output equivalence with the
retired fragment algorithm.

| Case | Without global pruning | With global pruning | Result |
|---|---:|---:|---|
| 10383, cut (131,132) | 120 s guard timeout | 5.37 CPU s matching; 6.79 CPU s through export/roundtrip | 100 families, 47 unique fragments; cap evidence retained |
| 37270, cut (65,66) | 2.41 CPU s matching | 0.193 CPU s matching | Same explicit 101 > 100 cap; zero returned mappings |

The measured complete processes peaked at approximately **231 MiB** and
**198 MiB**, respectively. The 10383 result cache remained within its 64 MiB
retention budget. That cache budget does not bound total process memory.

Global preparation has a cost: the 32-arm cold example changed from 2.84 to
4.69 ms matching, while cached matching changed from 2.13 to 2.08 ms. Exact-size
ring components avoid that preparation. These measurements establish the large
regression repair, not universal zero overhead or full-dataset throughput.

Binary publication now preserves shared typed fragments and matrices. For 16
distinct 374-atom placements sharing one source preparation, publication changed
from **0.933 to 0.0062 CPU s**, and peak Python allocation from **207.8 to
31.3 MB**. Explicit record export also shares source weight lists: 16 copies
became one. JSON remains an explicit expansion boundary. These allocation
measurements are distinct from whole-process RSS.

## Persistence and scope

The completed-result cache, immutable baseline/latest snapshot protocol and
shared-state transport are preserved. Binary publication excludes derived
lookup caches. The obsolete second finalization worker pool and duplicate cut
checkpoints are removed. Incompatible historical scientific schemas are
rejected explicitly; their producing revision remains necessary to inspect them.

[Measurements](measurements.json) retain the tested source/engine identities,
phase timings, memory scopes, cache statistics and limitations. The prior full
Golden/EnzymeMap results remain evidence for their recorded revisions. This
repair has not rerun those datasets.
