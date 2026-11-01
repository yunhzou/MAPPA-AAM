# Direct fragment completion — Mac validation

The implementation follows `design.md`: one witness during growth; direct
compressed completion from the accepted constraints; the existing hierarchy and
caps consume those choices. Candidate-growth replay, candidate deduplication and
symmetry-certificate retention are removed. Obsolete matching switches, unused
source-orbit work and merge-history output are removed. Completed-fragment cache
reuse and its approved snapshot recovery remain.

Local pools use matching predicates and locks. A component label alone cannot
split an interchangeable free-atom pool. Exact global target actions still
preserve molecular components and correlations. The immutable global action is
prepared once per native target; a subgroup solve is unnecessary when that action
already fixes every eligible pool class. Custom policies retain their compatibility
with all source atoms, including later fragments.

## Correctness

The direct-completion repair passed **753/753**, with no skips or failures, in
**55.22 s**. The subsequent support-check optimization is validated below.

```sh
PYTHONPATH=src:tests MAPPA_NATIVE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -m pytest -q
```

The critical tests are in `tests/test_fragment_matching_contract.py`. They check:

- C6 into C6: one candidate, exactly twelve rotations/reflections.
- C6 into two C6 rings: two candidates, twelve mappings per component.
- The same existence solver finds one witness or proves absence; resource limits
  never become a negative matching result.
- A tolerance pool containing **14!** valid assignments stays one candidate even
  when the exact target automorphism group is identity.
- Five interchangeable free hydrogens form the three valid pools rather than five
  unnecessary singleton choices.
- Internal pools, locks, custom policies, finalization, serialization and hierarchy
  caps preserve their scientific meaning.

Independent permutation checks also passed **400/400** small random queries,
including partial fragments, locks and custom policies: exact mapping coverage
and disjoint completed choices. The full suite includes the unchanged 4,608-mapping
decoding oracle and local worker-pool/cache/storage tests. Worker tests ran outside
the Mac sandbox because it blocks localhost listeners and process-tree inspection.

## Fragment performance

Mac M4 Pro; one compute process, one thread, no completed-fragment cache. Each
cell is the median of nine timed batches in a fresh process. Prepared calls reuse
endpoint preparation; public calls include the public fragment API preparation.
The archived baseline native engine has SHA-256
`13b2eb922245601cb0ded87e86a4c662b75ec84ad593ef2b976a8df3a1838e41`.

| Case | Prepared baseline → current | Public baseline → current |
|---|---:|---:|
| C6 → C6 | 48.97 → 16.45 µs | 114.44 → 48.43 µs |
| C6 → two C6 rings | 117.55 → 31.79 µs | 277.48 → 114.82 µs |
| 14! tolerance pool | 10 s guard → 25.34 µs | 10 s guard → 352.30 µs |

Ring candidate counts remain one and two. The old 14! run returned no candidate
before the external guard and reached **10.07 GiB** OS-reported peak RSS. Current
fragment runs peak below **62 MiB**, including imports. The timeout is a failed
completion measurement, not an absence result or an exact baseline runtime.

Fixed existence queries also pass the search-state and time gates. Earlier point
measurements of the retained existence kernel were 3.42 ms native for the saved
73-atom positive query, 1.58 ms for the independently proved 317-atom negative,
and 0.018 ms for the competing-destination negative. Both negative cases are
proved at the initial propagation, before branching on assignments.

## EnzymeMap failure cases

Identical saved inputs, cut order, seeds, tolerance and branch cap. Times below
cover the first cut start through final publication, excluding process startup.
RSS is sampled every 20 ms and includes endpoint preparation, matching,
finalization and publication. These are local point measurements, not cluster
throughput measurements.

| Scope | Previous → current time | Previous → current peak RSS |
|---|---:|---:|
| 03855 forward, cuts 0–262 | 79.51 → 60.02 s | 2.455 → 0.494 GiB |
| 37270 forward, cut 144 | 18.45 → 12.38 s | 3.212 → 1.893 GiB |
| 16175 reverse, cut 0 | unfinished at 180 s → 1.00 s | current 112.8 MiB |

All **23,446** returned full representative mappings across these diagnostic
scopes satisfy injectivity, element compatibility and every accepted edge.

Strict annotated-heavy-atom reference queries have no unresolved outcomes:

| Scope | Previous reference result | Current reference result |
|---|---|---|
| 03855, all 263 cuts | recovered in 84 cuts | recovered in 84 cuts |
| 37270, cut 144 | not recovered | not recovered |
| 16175, cut 0 | not recovered | not recovered |

The 03855 positive cuts overlap in 82 cases: cuts 55 and 224 gain recovery;
113 and 225 lose it. The combined reaction still recovers the annotated reference. Capped searches
can retain different paths after replacing history-dependent packing and ordering;
this report does not claim identical per-cut output or establish accuracy for the
entire EnzymeMap dataset. The two other scopes are individual diagnostic cuts,
not complete reaction sweeps. All matching and decoding scientific oracles pass.

The native completed-fragment snapshot format is updated; old snapshots are
explicitly incompatible. No legacy algorithm or error-triggered fallback remains.
This repair has been validated locally; it has not been deployed to a cluster.

Inputs and original failure measurements are alongside this report.
`summary.json` retains the original investigation and a `repair_validation`
section with final source/native identities, raw timing samples, reference hashes,
reference mappings and per-cut accuracy evidence. `memory_probe.py` in
`reports/native_transport_20261031/` is the existing replay/measurement entry point.

## Native support-check optimization, 2026-11-01

Compared with commit `85e1fcb`, the native constraint reducer now marks the allowed
destinations in a temporary lookup array and checks only actual bonded neighbors.
It reuses that scratch array within the reduction, then releases it. Matching
rules, search order, compressed completion, caps and persistent caches are
unchanged. Packed support tables and binary-search variants were tested and
removed because they added overhead to small queries.

The final full suite passed **763/763**, with no skips or failures, in **54.39 s**.
Ten added critical cases cover locks/occupied destinations, exact conditional
tolerance boundaries, and a positive/negative path pair with an independent
articulation-based absence proof. All fourteen performance probes retain identical
answers and search counters.

Fresh-process medians on the same Mac, one thread, eleven timed batches per probe:

| Prepared query | Baseline → current |
|---|---:|
| Saved 73-atom positive | 3.151 → 0.812 ms |
| Saved 317-atom negative | 1.473 → 0.483 ms |
| C5 into C66, impossible | 53.66 → 40.17 µs |
| Fourteen-atom path through three joined C6 rings, impossible | 155.15 → 148.51 µs |
| Small competing-destination positive | 14.85 → 16.10 µs |

C6 completion remains about **16 µs** for one ring and **32 µs** for two rings,
producing exactly one and two compressed matches. The small positive query costs
about 1.25 µs more; this is not a claim that every individual call becomes faster.
The large-query improvements repeat in both final timing rounds.

| Saved real scope | Native search baseline → current | Full cut scope baseline → current | Peak RSS baseline → current |
|---|---:|---:|---:|
| 03855 forward, cuts 0–262 | 28.87 → 26.28 s | 59.94 → 56.55 s | 508.8 → 502.3 MiB |
| 37270 forward, cut 144 | 8.53 → 8.38 s | 11.70 → 11.47 s | 1914.9 → 1901.0 MiB |
| 16175 reverse, cut 0 | 0.829 → 0.097 s | 0.970 → 0.240 s | 104.0 → 102.8 MiB |

Every field of all **265/265 finalized search graphs** is identical to the saved
baseline, including contexts, roots, states, transitions and stops. Native search
includes graph export; full scope runs from the first cut start to the last
publication. RSS is sampled every 20 ms. These are local point measurements,
not cluster throughput or a uniform runtime guarantee for hard absence queries.

`summary.json` → `repair_validation.propagation_experiment` retains final source,
header and engine identities, raw probe samples, exact per-cut comparisons and
phase timings. The earlier 400-query permutation check remains historical; the
full regression suite and these exact comparisons validate this further change.
