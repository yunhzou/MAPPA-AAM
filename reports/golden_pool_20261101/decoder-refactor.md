# Unique-fragment decoder refactor

The initial v3 refactor and its measurements are recorded below. The subsequent
v4 change removes automatic whole-reaction reporting; see the final section.

The decoder now collects all unique saved reactant/product fragments before
event-response work. Local results are shared across branch occurrences only
when their mapping, ordered actions and cut context agree. Branches retain their
original action bindings, constraints and boundary relationships. Matching,
symmetry choices, hierarchy caps and event thresholds are unchanged.

The exporter preserves the v3 record values and ordering while sharing repeated
term dictionaries and lists. Binary publication therefore saves shared objects
once. Scheduling and bundle formats are unchanged; no per-fragment files or new
worker-memory allowances were introduced.

## Golden 1561 forward on Mac

The same saved catalogue was decoded before and after the refactor in separate
processes. The baseline is source `adfb626f7978ddaaa0b0d5ff4ae7e14f26edc592`.
Both complete expanded JSON records have SHA256:

`f415b641434ac0bf5e1a5211c43260731f0f75273bffd66de1e71b7bb1e9dfa3`

This comparison includes all families, fragment partitions, constraints, action
dependencies, provenance, internal/boundary terms, response tables, annotations,
statistics and representative descriptions. It is stronger than comparing
representative event counts alone.

| Measurement | Before | After |
|---|---:|---:|
| Complete record built, including catalogue load | 43.43 s | 43.11 s |
| Peak RSS through record hashing | 1.511 GiB | 0.947 GiB |
| Logical bond-term occurrences | 1,006,967 | 1,006,967 |
| Stored term dictionaries | 1,006,967 | 13,342 |
| Term dictionaries alone, excluding contents | 261.2 MiB | 3.46 MiB |
| Full binary publication under 2 GiB | Guard terminated | Completed |

The matched record-building measurements show essentially unchanged runtime,
with 37.3% lower peak RSS. These are single-process point measurements, not a
general speed guarantee. The subsequent publication replay completed in 46.67 s
including decoding, with **1.187 GiB peak RSS** and a **4,592,149-byte compressed
output**. Its parent sampled RSS every 20 ms under the unchanged 2 GiB guard.
That replay overlapped the test run and is not a paired timing comparison.

Evidence: [baseline record](decoder-baseline-digest.jsonl),
[refactored record](decoder-refactored-digest.jsonl),
[publication replay](decoder-refactored-publication.jsonl),
[original failure diagnosis](decoder-memory.md). The same
[profile tool](profile_decoder.py) performs both record hashing (`--digest-only`)
and publication measurement. Input provenance is in
[decoder-memory-source.json](decoder-memory-source.json).

## Correctness validation

The affected decoder, pool and benchmark tests passed. New tests check collection
before event work, shared local results across distinct branches, rebinding when
global action positions differ, cut annotations, required-edge constraints,
alternative/coupled partitions, and shared serialization. Every allowed choice
in the small branch examples is checked against an independent dense event oracle.

The full suite ran 742 tests: 741 passed; one subprocess test could not import
its test helper because `tests/` was absent from the inherited `PYTHONPATH`.
That test passed when rerun with the test directory on `PYTHONPATH`. No scientific
or execution failures remain. The suite includes the smoke checks. Pool tests
required local socket and process-inspection access outside the filesystem sandbox.

## DFW publication check

Job **21136551** completed decoding and durable publication from the saved
catalogue using source **5fd12602addfb597bcad422fc5abc132594dd7d4**. One compute
worker retained its **2 GiB guard**; the two-CPU, 4 GiB allocation also housed the
coordinator. Peak worker RSS was **0.931 GiB during computation** and **1.141 GiB
during publication**, versus the previous fresh-process retry's 1.632 and
2.021 GiB (terminated). The published compressed result is 4,592,181 bytes.

The successful task used 94.23 CPU seconds for computation and 1.76 CPU seconds
for publication; elapsed task time was 96.08 s and allocation time 1 min 45 s.
The old fresh-process computation used 89.85 CPU seconds on a different node.
Different-node measurements cannot isolate a code runtime change. The same-input
Mac measurements above showed essentially unchanged runtime.

The output and audit byte ranges were checksum-verified before recording the
recovery through the existing ledger API. The actual published DFW result was
also loaded on Mac and its complete expanded record hashed: it matches the
pre-refactor SHA256 above exactly. This checks the persisted result across the
Python 3.10 cluster / Python 3.12 Mac binary round trip; see
[published-record comparison](decoder-published-digest.json).

The failed ledger row was retained in
the retry history, and every other stage row was checked unchanged. All five
ledgers now report **3,702 completed directions in each of matching,
deduplication, decoding and verification**, with no remaining operational
failures. Matching and reference verification were not repeated. Coverage is
unchanged: 1,830 of 1,851 reactions recovered.

[Launch and source revision](refactor-decoder-launch.json) ·
[Published artifact, resource audit and completed counts](decoder-refactor-completion.json)

## Why the total decoding time barely changed

A stage-timed replay of the same saved catalogue confirms that cross-cut/seed
deduplication has already finished. The input is a `CatalogueResult` with 14,132
families and 14,132 distinct complete family keys: there are no exact duplicate
family records. Its 14,132 input paths group into 11,899 fragment-set branches.
The family key includes policy, witness, required edges and ordered actions;
different keys are retained even when branches share structural fragments.

The complete record took 43.30 s to build in this instrumented replay:

| Work | Elapsed time |
|---|---:|
| Unique-fragment collection | 0.21 s |
| Collection plus action annotations | 0.35 s |
| All family/partition table assembly | 8.87 s |
| Representative descriptions for every family | 33.15 s |

Rows are inclusive: collection is inside collection/annotations. Table assembly
includes fragment lookups, boundaries and dependency assembly. Other time covers
input loading and record wrapping. These timers add measurement overhead and are
for attribution, not a new paired performance claim.

The shared layer has 2,441 structural fragments, 4,472 distinct local decoding
contexts and 252,442 occurrences. All fragment lookups and internal decoding
together cost 1.24 s, included in the 8.87 s table-assembly time. Conditional
contexts remain separate when witnesses, ordered actions or cut annotations differ.

`decoded_events_record()` still calls `EventInterpretation.evaluate()` once for
every family after exporting the local tables. That call uses a separate
whole-reaction description path instead of composing the already prepared local
tables. Its 33.15 s includes 2.63 s of bond classification, 25.60 s of event-ID
construction/description, 3.39 s of dense count cross-checks and 0.84 s of
representative validation. Event-ID construction canonicalizes the whole
event-marked graph when its existing cache misses. This is a representative report
for each saved family; it does not enumerate that family's symmetry choices.

Thus about 77% of the recorded time remains in a reporting pass untouched by
fragment sharing. The previous decoder also already cached local bond responses;
sharing their surrounding records primarily removed storage duplication. Faster
end-to-end decoding must address the representative-reporting path while preserving
the required output contract. The completed refactor does not establish that all
branch-level reporting now composes shared fragments.

[Timing evidence](decoder-timing-breakdown.jsonl) was collected with
`profile_decoder.py --timings --digest-only`. The complete record still has the
same SHA256 as the pre-refactor baseline. No production algorithm was changed for
this diagnosis.

## Compose local results without automatic representative reports

On 2026-11-01 the user confirmed that normal decoding should derive event results
from the saved fragments and boundary relationships, without another global
event-pattern deduplication pass. Export now returns `mappa.decoded_events/v4`:
the automatic `representatives` section is removed. The catalogue, partitions,
local responses, dependencies, constraints and provenance are unchanged.
`contributions()` composes requested valid choices; `evaluate()` remains an
explicit optional whole-reaction report with a canonical ID. No new scheduling,
cache, matching or event-counting algorithm was introduced.

The complete new Golden 1561 record was compared with the actual published DFW
v3 result, changing only its schema version and removing `representatives`.
Every remaining value and ordering is identical. Both expanded JSON records,
after this normalization, have SHA256:

`6f67159f8a07143ba94a67f5534159e52d543a36b99f7be492cad3980802bad6`

The instrumented replay built the record in 10.60 s, versus 43.30 s previously,
and made zero calls to representative descriptions, canonical event IDs,
whole-reaction classification or dense whole-reaction count checks. It retained
all 14,132 families, 2,441 structural fragments, 4,472 local decoding contexts,
and 13,342 shared term dictionaries. This replay overlapped the test suite;
the stage timing is attribution evidence rather than a controlled paired trial.

A subsequent replay with no tests running completed decoding and binary
serialization in **11.34 s**, at **1.004 GiB peak RSS**, under the unchanged
2 GiB guard. The output compressed to **3,001,904 bytes**. The previous v3
publication replay took 46.67 s at 1.187 GiB; it overlapped tests, so that
approximately fourfold timing difference is a point observation, not a general
throughput guarantee. Existing completed Golden artifacts were not replaced.

The full suite ran 743 tests: 742 passed, and one updated CLI test initially
compared Python tuples with JSON lists. Its expected representation was corrected
and that test passed on rerun. Existing independent dense event oracles and
explicit canonical-ID tests remain. A new regression test forbids whole-reaction
reporting during normal decoding/export and verifies local choice counts.

Evidence: [normalized complete-record comparison](decoder-composition-digest.jsonl)
and [publication measurement](decoder-composition-publication.jsonl). Reproduce
the comparison with `profile_decoder.py --timings --compare-record` and the trusted
single-blob v3 baseline; the default command measures binary serialization.
