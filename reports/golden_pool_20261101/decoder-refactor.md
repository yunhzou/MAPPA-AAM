# Unique-fragment decoder refactor

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

The failed direction will be retried from its saved complete catalogue with one
worker guarded at 2 GiB. Its matching and reference-verification results remain
reused. Completion and measured cluster resources will be recorded here after
the point run.
