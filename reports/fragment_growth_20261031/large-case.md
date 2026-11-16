# Previous large EnzymeMap case

Historical experiment. The comparison modes and scripts were removed from the current checkout on 2026-11-01; reproduce these commands at their recorded Git revision.

This is the **pre-cache integration** experiment at revision `30f0bc9`.
Its diagnostic build commands require that revision. The later archived `replay_cut.py`
uses the installed engine and bounded fragment-result caching for both policies.
See [the cache integration report](cache-integration.md) for the subsequent results.

EnzymeMap **10383 forward, cut 334**, the archived 4 GiB failure: 374 atoms
per endpoint, cut `(131, 132)`, one original seed order, random seed 42,
matching tolerance 1.0, branch limit 100. The original trusted input is
large-case-input.pkl.gz (raw archive retained privately).

On the user's M4 Pro Mac, existence growth substantially improves the first
search, with exactly equal complete finalized search graphs. It regresses on
an identical repeated cut because it does not use the production fragment replay
cache. It is therefore **not yet a performance-safe replacement for reused
workers**.

| Measurement | Incremental | Existence |
|---|---:|---:|
| First search + symmetry finalization, CPU | 438.60 s | 35.54 s |
| First search + symmetry finalization, wall | 438.79 s | 35.56 s |
| First-run peak process RSS | 4.193 GiB | 0.586 GiB |
| Identical repeated cut + finalization, CPU | 0.442 s | 12.046 s |
| Peak process RSS through repeated run | 4.232 GiB | 0.617 GiB |
| Retained certificate entries | 48,978 | 4,549 |
| Retained certificate/key payload | 3.447 GiB | 0.321 GiB |

The first run is **12.34× faster with 86.0% less peak RAM**. The repeated cut is
**27.24× slower**. These are one cold and one repeated observation per policy,
not statistical estimates or a complete-reaction throughput measurement.

## Exact results

All four complete finalized graphs compare equal: contexts, roots, ordered
states, transitions and stops. Each contains **9,340 states, 9,339 transitions
and 100 terminal states**. The comparison includes the saved witnesses,
compressed symmetry records, fragment boundaries, accepted/deferred bonds,
branching and cap evidence. It does not expand atom bijections.

Both existence runs also compare equal to the complete graph saved by the earlier
Mac replay of this failed cut. The checked graphs are after symmetry finalization;
this experiment did not separately save the raw pre-finalization graphs or each
individual fragment-call return. Decoding and reference verification were not rerun.

## What accounts for the difference

There are 2,476 fragment calls in each cut. On the first incremental run,
1,019 calls reuse a saved fragment result. On the repeated incremental run,
**all 2,476 calls reuse saved results**. Existence performs every call again;
it cannot currently use the candidate-based replay protocol.

Existence makes 8,235 whole-fragment feasibility searches per cut. Its completion
stage rejects 21,774 intermediate compressed candidates using the finished
fragment constraints. Completion costs 24.79 s on the first run and 1.65 s on
the repeated run; the other approximately 10.4 s remains largely in witness
growth. Its certificate cache does benefit from the repeated run, but this is
different from reusing whole fragment results.

The profile's maximum of one growth candidate describes the growing witness.
It does **not** assert that final completion retains only one compressed candidate.
No branch cap, certificate-cache limit or scientific constraint was changed.

## Test boundary and reproduction

Source revision: `e9481be`. [archived replay script](https://github.com/yunhzou/MAPPA-AAM/blob/e9e97f255e8d0f5b1695a2ace858ecc6948e84a6/reports/fragment_growth_20261031/replay_cut.py) builds a separate
diagnostic engine. It selects the existing growth policy at the native hierarchy's
fragment call and adds counters. Branch admission, deduplication, hierarchy caps
and native transport are unchanged. Incremental keeps its production replay
cache; existence calls its existing growth implementation without candidate replay.
The installed native engine's hash was confirmed unchanged.

Each policy ran in a fresh, single-thread process and repeated the same cut in
that process. The two processes overlapped on the 14-core Mac; each was observed
using one CPU. Measurements include native search, transport and symmetry
finalization, excluding graph preparation, imports and result serialization.
RSS includes the process and its retained caches. The repeated-run RSS peak can
also include serialization of the first output.

With the project Python environment and `OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `MAPPA_NATIVE=1`, `PYTHONPATH=src`:

```sh
python reports/fragment_growth_20261031/replay_cut.py \
  --build out/fragment-growth-large-20261031/build

# Run once for each POLICY: incremental, existence.
python reports/fragment_growth_20261031/replay_cut.py \
  --engine out/fragment-growth-large-20261031/build/lib/mappa/_engine.cpython-312-darwin.so \
  --input reports/fragment_growth_20261031/large-case-input.pkl.gz \
  --output out/fragment-growth-large-20261031/10383-POLICY --policy POLICY

python reports/fragment_growth_20261031/replay_cut.py \
  --compare out/fragment-growth-large-20261031/10383-incremental \
            out/fragment-growth-large-20261031/10383-existence \
  --output out/fragment-growth-large-20261031/comparison.json
```

The native filename above is for the measured Mac/Python version. The build prints
its actual filename. [large-case-results.json](large-case-results.json) retains
measurements, exact equality assertions, build/source hashes and historical
comparison provenance. This test did not change or deploy production code.
