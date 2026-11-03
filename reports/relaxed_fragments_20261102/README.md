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

The final performance revision passed **833/833 repository tests** in 113.09
seconds with `python -m pytest -q --tb=short`. New checks cover independent
weighted placement oracles, shared-cache ownership and eviction, dense bond
neighborhood equivalence, and once-only internal fragment preparation. All ten
complete graph and decoded-record audits also match the initial integration.

The initial integration passed **824/824 repository tests** in 109.20 seconds
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

## Paired comparison with the previous implementation

The correct previous-implementation baseline is `2e99e989`, the parent of the
integration commit `64ee9bb5`. The initial measurements below compared two stages
of the new implementation; they did **not** compare against that previous release.

The paired comparison found avoidable regressions and removed them: repeated
global domain/predicate checks, unchanged parent propagation, copied fragment
records, duplicated cache preparations, and repeated internal decoder work.
Shared preparations remain inside the unchanged **64 MiB** retention budget;
last-owner eviction releases their index entries. Cache snapshots now use schema
`mappa.fragment-cache/8`. No mapping rule, hierarchy cap or tolerance changed.

Ten complete search graphs and ten complete decoded-record hashes match the
unoptimized `64ee9bb5` results exactly. This checks the optimizations against the
approved new algorithm; it does not claim identical scientific outputs between
the old and new algorithms.

The final comparison used 92 isolated worker processes, run serially, and 184
timed pipeline runs on the same Apple Silicon Mac. Each process performed a cold run and one
warm repeat. Revision order alternated, with five repetitions for the eight
smaller cases and three for the two large cases. No tests or builds ran alongside
these timings. All 92 processes completed successfully.

Both revisions used one seed, random seed 42, `iso_tol=1`, presence thresholds
0.2, ordinary/metal event thresholds 0.5/0.3, and the same cut and input bytes.
Branch cap was 2,000 for synthetic/Golden points and 100 for the EnzymeMap points.
The new fragment self tolerance remained 1. These are single-cut forward points,
not full sweeps. Two-ring targets have different external F/Cl attachments;
the attached-ring example is a six-atom path on two C6 rings joined by a bridge.

Median CPU seconds cover planning, preparation, complete matching, reduction,
decoding, record export, binary serialization and restoration. Imports and input
loading are excluded. Warm repeats reuse matching preparation and its cache;
reduction and decoding are performed again.

| Case | Previous cold | Current cold | Previous warm | Current warm | Families: previous → current |
|---|---:|---:|---:|---:|---:|
| C6 → C6 | 0.00464 | 0.00387 | 0.00100 | 0.00072 | 1 → 1 |
| C6 → two C6 rings | 0.00634 | 0.00465 | 0.00238 | 0.00114 | 12 → 2 |
| C24 → two C24 rings | 0.02717 | 0.00876 | 0.01955 | 0.00444 | 48 → 2 |
| P6 → attached C6 rings | 0.00713 | 0.00869 | 0.00262 | 0.00506 | 11 → 18 |
| Shared root, 32 arms | 0.08827 | 0.03274 | 0.04054 | 0.01965 | 1 → 32 |
| Golden 1001 | 0.00512 | 0.00428 | 0.00141 | 0.00092 | 1 → 1 |
| Golden 785 | 0.01122 | 0.01545 | 0.00518 | 0.01083 | 3 → 18 |
| Golden 1518 | 0.25378 | 1.26174 | 0.20254 | 1.24537 | 262 → 2,000 |
| EnzymeMap 10383, cut (131,132) | 2.92197 | 1.88792 | 0.75202 | 0.68608 | 100 → 100 |

Six of these nine points are faster; three are slower. All three slower cases
retain more families under the corrected target-set definition. Golden 1518
performs about eight times as many cache requests and hierarchy transitions;
its returned families are stalled terminal outputs, and the current run also
has cap evidence. Its 4.97-fold cold and 6.15-fold warm total-time increases
remain real absolute slowdowns. Equal caps do not imply equal work, and these results do not
support a universal no-regression claim.

For 10383, cold matching is about 1.24 CPU s and warm matching 0.08 s. Its cache
retains 30.28 MiB with zero evictions, versus repeated eviction before this
optimization. Median cold peak RSS is 174.38 MiB, compared with 174.73 MiB for
the previous release. Golden 1518's larger output raises cold peak RSS from
81.91 to 175.41 MiB. RSS is the complete process lifetime high-water mark,
including imports and inputs; cache retention is a separate measurement.

The tenth point, **37270**, is excluded from returned-output speed comparisons:
the previous version returned 100 families in 9.90 CPU s, whereas the current
version returns none after its explicit 101 > 100 fragment cap (0.254 CPU s).
That is neither recovery nor a decoding speedup. Its cap/output behavior is
unchanged by these optimizations.

[Paired measurements](paired-performance.json) contain every timed phase,
repetition, output count, cap summary, cache counter, RSS observation, input hash
and source/engine identity, without private operational paths. This study makes
no full-dataset accuracy or cluster-throughput claim.

## Initial integration diagnostics

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
