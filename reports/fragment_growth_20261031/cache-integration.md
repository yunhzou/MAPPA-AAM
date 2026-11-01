# Bounded result reuse for existence growth

Historical experiment. The comparison modes and scripts were removed from the current checkout on 2026-11-01; reproduce these commands at their recorded Git revision.

At this revision, both growth policies used the existing complete-fragment result cache in
native reused execution. This removes the experimental bypass measured in
[the earlier large-case report](large-case.md). The growth, compression, branch
admission and hierarchy-cap rules are unchanged.

The cache retains completed compressed answers and the conditions required for
reuse. It has a **64 MiB accounted retention budget**, evicts least-recently-used
entries, and returns oversized results without retaining a copy. Preserved bonds
and nested vector capacity are now included in the accounting. Each instance is
bound to one growth policy. Candidate-prefix checkpoints remain separate.

## Previous large raw case

EnzymeMap **10383 forward, cut 334**, with the original 374-atom endpoints, seed
order, random seed 42, tolerance 1.0 and branch cap 100. Measurements used the
installed implementation on the user's Apple M4 Pro Mac. Each policy ran in its
own single-thread process; matching and symmetry finalization are included.

| Measurement | Incremental | Existence with cache |
|---|---:|---:|
| First run, CPU | 306.085 s | 32.094 s |
| First-run peak process RAM | 4.194 GiB | 0.565 GiB |
| Repeated-run median CPU, separate processes | 0.312 s | 0.355 s |
| Repeated CPU, same process with balanced order | 0.361 s | 0.362 s |
| Retained complete-fragment entries | 1,457 | 1,457 |
| Accounted fragment-cache memory | 23.027 MiB | 23.027 MiB |
| Configured fragment-cache budget | 64 MiB | 64 MiB |

The measured first run is **9.54× faster with 86.5% less peak RAM**. There is one
cold sample per policy. The processes overlapped on separate CPUs; the full test
suite also ran during part of these first searches. The repeated medians use two
incremental repeats and nine existence repeats. The separate-process difference
motivated the paired timing check recorded with the measurements.

That check alternated both policies in one process for 15 pairs, with independent
prepared caches and exact graph comparisons after timing. The second operation
in a pair was consistently slower regardless of policy. Giving the two execution
orders equal weight across **all 15 samples** yields 0.360573 s for incremental
and 0.362101 s for existence, a **0.42% measured difference**. Both caches hit on
all warmed fragment requests. The former 12.05-second repeated-run regression
has been removed; the separate-process 0.312/0.355-second comparison alone did
not control execution environment and ordering sufficiently.

Across all ten existence runs, the cache retained exactly 1,457 entries and
24,145,284 accounted bytes. Every repeated fragment call was a cache hit. Total
process peak RAM settled at **0.601 GiB**, including retained library/allocator
memory and earlier output serialization. The cache itself did not grow.

All **26 raw/finalized graph comparisons** passed across the thirteen runs.
Every graph contains the same 9,340 states, 9,339 transitions and 100 terminals.
Equality includes complete ordered states, compressed fragment witnesses/pools,
accepted/deferred bonds, contexts and cap stops. It checks both the raw matching
graph and the graph after symmetry finalization, without enumerating bijections.

## Validation and use

**714 tests passed in 53.44 seconds.** Coverage includes forced cache eviction,
oversized-result rejection, changed cuts, locked mappings, island membership,
target occupancy and deferred boundaries. It also includes both policies through
reference/native/shared-seed execution, local multiworker RAM/bundle transport,
cross-cut deduplication and decoding. The large-case timing itself does not rerun
decoding or reference verification.

```python
from mappa import AAMSearchConfig, search_aam

config = AAMSearchConfig(growth_policy="existence")
aam = search_aam(problem, config, execution="reused_native")
```

Pass the same config to `pool.search(cases, config)`. Incremental remains the
default. The reference backend supports either growth policy without result
caching. No Slurm jobs or worker memory requests were changed.

The cache budget applies to retained completed-fragment entries; it is not a
total process-memory limit. In particular, the existing target symmetry-certificate
cache is separate. The measurements above establish this case's behavior, not a
universal memory bound for all reactions.

## Reproduction

Build the current native engine with the project Python environment. Set
`PYTHONPATH=src`, `MAPPA_NATIVE=1`, `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`
and `MKL_NUM_THREADS=1`. Use [archived replay script](https://github.com/yunhzou/MAPPA-AAM/blob/e9e97f255e8d0f5b1695a2ace858ecc6948e84a6/reports/fragment_growth_20261031/replay_cut.py):

```sh
python reports/fragment_growth_20261031/replay_cut.py \
  --input reports/fragment_growth_20261031/large-case-input.pkl.gz \
  --policy existence --repeats 10 \
  --output out/fragment-result-cache-20261031/existence

# Repeat in a fresh process with --policy incremental and its own output directory.
python reports/fragment_growth_20261031/replay_cut.py \
  --compare out/fragment-result-cache-20261031/incremental \
            out/fragment-result-cache-20261031/existence \
  --output out/fragment-result-cache-20261031/comparison.json

python reports/fragment_growth_20261031/replay_cut.py \
  --input reports/fragment_growth_20261031/large-case-input.pkl.gz \
  --paired 15 --output out/fragment-result-cache-20261031/paired.json
```

[cache-results.json](cache-results.json) contains timings, memory counters, exact
graph comparisons, source/native hashes and validation scope. The historical
uncached diagnostic build remains in Git; the maintained runner now uses the
production API directly.
