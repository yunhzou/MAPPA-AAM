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

The user subsequently approved lossless compression across equivalent target
placements. A representative continuation retains a shared global relationship
that moves its entire continuation jointly. Certification preserves the complete
target WBO matrix, active topology, source-node compatibility and every prior
lock. Fragment-local two-sided relaxed relationships remain separate. Compact
set-orbit exclusions avoid expanding all equivalent target sets before merging;
caps count representative continuations. Uncertified relaxed exchanges cannot
discard branches. This replaces the explicit target-set branch contract of the
historical comparison below.

On the final measured source, EnzymeMap 37270's first fragment retains all 288
previously established target sets through two representatives, each covering
144 sets. The 132 source atoms
and 131 accepted edges are unchanged. Each old set belongs to exactly one saved
orbit; cap 100 no longer truncates this first fragment. This check is not a
full-cut or reference-recovery result.

## Placement-compression validation

The final integrated source passed **880/880 tests in 110.21 seconds** on the
local Apple Silicon Mac with `python -m pytest -q --tb=short`. This includes
independent Python/native mapping oracles, correlated occupancy after later
fragments, original-WBO local/boundary event checks, external reference queries,
RAM/bundle worker pools, cache exchange and saved-result roundtrips.

The C6/H example retains all 768 mappings in one continuation instead of 64
explicit atom-set records. A multi-fragment ring/H example retains all 144
independently enumerated mappings in three continuations; a ring plus isolated
carbon retains all 24 mappings and their original bond events. Prior locks,
non-transitive relaxed relationships, custom compatibility, active cut topology,
unequal continuous WBOs and subthreshold contacts have separate checks.
Mutable graph inputs invalidate prepared native views and placement certification;
prepared reaction sessions keep their direct reuse path.

Placement relationships survive hierarchy transitions, catalogue reduction,
conversion, decoding, reference queries, checkpoints and replay. Catalogue schema
is `mappa.final_fragment_branches/v3`, search-graph record schema is `/3`, and
cache snapshots use `mappa.fragment-cache/11`; incompatible records are rejected.
The existing bounded cache and baseline/latest snapshot publication remain.
The obsolete expanded-orbit discovery and its image counter are removed.

## Current placement-compression performance

This comparison uses **`9cd89b0f`**, the immediate pre-compression implementation,
as its baseline. It is distinct from the historical `2e99e989` comparison below.
The [sanitized measurements](placement-performance.json) retain both complete
source/engine identities, input hashes, settings and all individual results.

The study ran 60 isolated processes serially on the same Apple Silicon Mac:
three repetitions of ten cases and two revisions. Each process performed one
cold run and one warm repeat, giving 120 timed pipeline runs. All processes
finished without execution errors. No builds or tests overlapped the timings.
Each point used one seed, random seed 42, discovery/self tolerance 1, presence
threshold 0.2 and event thresholds 0.5/0.3. The branch cap was 2,000 for
synthetic/Golden points and 100 for EnzymeMap points. These are single-cut
forward runs, not full sweeps. The two-ring targets retain their distinguishing
external attachments.

The table gives median total CPU seconds from preparation through matching,
reduction, decoding, export and binary roundtrip; cold runs also include
planning. Imports and input loading are excluded. Warm runs reuse matching
preparation and cache, then repeat the remaining stages.

| Case | Baseline cold | Placement cold | Baseline warm | Placement warm | Families: baseline → placement |
|---|---:|---:|---:|---:|---:|
| C6 → C6 | 0.00443 | 0.00383 | 0.00087 | 0.00077 | 1 → 1 |
| C6 → two C6 rings | 0.00438 | 0.00435 | 0.00105 | 0.00120 | 2 → 2 |
| C24 → two C24 rings | 0.00862 | 0.00834 | 0.00433 | 0.00451 | 2 → 2 |
| P6 → attached C6 rings | 0.00838 | 0.00531 | 0.00481 | 0.00162 | 18 → 4 |
| Shared root, 32 arms | 0.03249 | 0.01893 | 0.02024 | 0.01013 | 32 → 1 |
| Golden 1001 | 0.00459 | 0.00425 | 0.00098 | 0.00102 | 1 → 1 |
| Golden 785 | 0.01517 | 0.00901 | 0.01078 | 0.00471 | 18 → 3 |
| Golden 1518 | 1.26478 | 0.09664 | 1.24574 | 0.08905 | 2,000 → 94 |
| EnzymeMap 10383, cut (131,132) | 1.87397 | 3.39598 | 0.68567 | 2.06830 | 100 → 100 |
| EnzymeMap 37270, cut (65,66) | 0.25017 | 4.41074 | 0.07222 | 4.36885 | 0 → 100 |

Golden 1518 is 13.1 times faster cold and no longer hits its hierarchy cap.
Its uncut reference query still returns `not_recovered`, with no unknown outcome,
as at both earlier revisions; this speedup is not a new reference recovery.
The attached-ring and shared-root points improve by 1.58 and 1.72 times cold,
respectively. Smaller warm-run increases remain on three points; these are
0.036–0.188 ms in this study.

**10383 remains slower: 1.81 times cold and 3.02 times warm.** Both runs reach
cap 100, but the saved placement relationships retain wider mapping scope in
the new run. Equal representative counts do not mean equal retained mappings.
Its decoder prepares broader endpoint domains and boundary responses; it does
not enumerate complete group elements or bijections. Largest observed cold RSS
increases from 175.1 to 227.7 MiB. This is a real total-time increase, not a
universal no-regression result.

**37270 is not an equal-output speed comparison:** the baseline stops with zero
families, whereas placement compression returns 100 families and preserves later
cap evidence. The largest observed RSS across this study is 475.5 MiB on its
warm run. RSS includes the complete process lifetime and earlier repeats; it is
not a memory bound for other reactions.

The final support-sharing and decoder-preparation optimizations preserve the
complete 10383 decoded-record SHA exactly
(`e9d2f120db5a306dc4c5e941a30bdd233871d66fda09740264dbc584095ab5ac`).
Independent relation and event oracles validate compression separately from
timing. These point measurements establish neither full-dataset reference
recovery nor cluster throughput, and do not support a universal speedup claim.

## Historical validation before placement compression

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
retain more families. Golden 1518 exposes a change in branch representation:
under the same parent locks, the old fragment represents source H90 mapping to
any of target H29, H30 or H31 in one record. All three target hydrogens are free
and bonded only to C0 with WBO 1, so their exchanges preserve the original
endpoint and every lock. The current matcher emits three separate fragments and
child branches. These are three different labelled target atom sets, so the
separate branches follow the approved design. The user reaffirmed this rule
after inspecting the example. Global symmetry can share discovery work without
combining their branch identities. The old grouping across target sets explains
part of the count difference; that difference alone establishes neither a
matcher defect nor improved reference coverage.

The previous 262-family run has no cap records. The current run reaches the
2,000 combined-live-branch cap and performs about eight times as many cache
requests and hierarchy transitions. Both fill all 46 target atoms; their
terminal `stalled` status reflects the 104-atom source, not a failure to fill the
target. The current 4.97-fold cold and 6.15-fold warm total-time increases remain
real measured slowdowns. They do not by themselves justify changing branch
identity. The 2,000 families reuse 76 unique fragments; process memory also
includes the retained hierarchy, catalogue and exported records. The live-branch
cap is not a byte limit on all those structures.

Own-revision reference verification of this exact uncut point returns
`not_recovered` for both versions, with no timeout or unknown outcome. The
historical full cut sweep recovered Golden 1518 in both directions, but is a
different search scope. This diagnosis does not establish a new
reference-recovery loss on this point. These results
do not support a universal no-regression claim.

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

A focused diagnosis of 37270's first fragment, on the same cut and seed order,
establishes its complete placement count: **288 distinct target atom sets** for
one 132-atom source fragment (49 C, 11 N, 13 O, 59 H), within 661-atom endpoints.
Four independent target hydrogen pairs and two triples each contribute one
selected hydrogen, giving `2^4 × 3^2 = 144` images for a fixed heavy-atom set.
Every hydrogen group has identical complete original-WBO rows, so these choices
preserve all discovery constraints. Two heavy-atom sets differ only in which
oxygen receives source O248: target O169 has a double bond to C168, whereas
O170 has a single bond to that carbon. Both satisfy discovery tolerance 1,
giving `144 × 2 = 288` placements. These are alternative placements of one
fragment, not 288 fragments composing a reaction.

Explicit Python discovery exhausts at 288 unique images; native completion at
diagnostic cap 300 returns exactly the same images, source atoms and 131
accepted edges. Each witness was checked for injectivity, elements and accepted
bond constraints. The native first-fragment call took 0.145 CPU s in this single
diagnostic, with 206.56 MiB process-lifetime peak RSS including imports and raw
record export. This was not a complete-reaction run or a new reference-recovery
measurement. Production settings and source were unchanged. The original
101 count was the cap-plus-one stopping point, not the complete candidate count.

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
