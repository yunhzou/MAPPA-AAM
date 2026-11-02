# Golden seed scaling — provisional local results

Snapshot: 2026-11-02T06:02:52+00:00. **104/128 requested direction tasks are complete; the trial grid is ongoing.** 52 complete cache-on/off pairs have identical catalogues, seed/cut scope, per-cut graph counts and cap evidence. Unfinished cells remain missing; no values are zero-filled or extrapolated.

![Matching CPU across seed counts](scaling.png)

The eight reactions were selected at recorded matching-cost quantiles using [selection.json](selection.json). Selection costs came from the earlier one-seed, cap-100 study; the new measurements use cap 2,000. This cost-stratified sample is not an estimate of benchmark-wide speed or accuracy.

## Fixed six-case comparison

Cases **1001, 1850, 1022, 785, 912 and 1518** have all 16 requested settings complete. Every row below sums the same six reactions and both directions. One trial was run per setting; there are no repeats or uncertainty estimates.

| Seed orderings | Cache on: matching CPU s | Cache off: matching CPU s | Off/on CPU ratio | On versus one seed | Recovered reactions, on/off |
|---:|---:|---:|---:|---:|---:|
| 1 | 10.820 | 26.353 | 2.44× | 1.00× | 6/6 · 6/6 |
| 2 | 24.328 | 58.959 | 2.42× | 2.25× | 6/6 · 6/6 |
| 5 | 63.521 | 144.041 | 2.27× | 5.87× | 6/6 · 6/6 |
| 10 | 131.502 | 290.334 | 2.21× | 12.15× | 6/6 · 6/6 |

The off/on ratio measures the completed-fragment cache effect for the same requested work. It is not a comparison with an older algorithm. Millisecond-scale points are sensitive to one-run timing noise.

## Timing boundaries

Matching uses only the measured matching stage. Preparation, saved-symmetry finalization, family reduction, reference verification and diagnostic fingerprinting are recorded separately below. **No event decoding was run.** CPU values sum worker execution; summed direction wall times in the JSON are not parallel reaction latency.

Cache-on totals for the same fixed six-case cohort:

| Seeds | Preparation CPU s | Matching CPU s | Finalization CPU s | Reduction CPU s | Verification CPU s | Audit CPU s (excluded) |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.028 | 10.820 | 10.338 | 18.055 | 3.985 | 22.285 |
| 2 | 0.028 | 24.328 | 20.009 | 35.198 | 4.436 | 36.006 |
| 5 | 0.026 | 63.521 | 49.877 | 101.292 | 8.367 | 98.880 |
| 10 | 0.026 | 131.502 | 93.818 | 197.128 | 2.981 | 187.496 |

## Scope and checks

- Hardware: Apple M4 Pro. Each case/direction/seed-count/cache-setting starts in a fresh worker. The fragment cache stays resident across that direction’s cuts and seeds. There is no cross-task cache import.
- Seeds: 1, 2, 5 and 10; random seed 42; branch cap 2,000; full single-edge cut sweep; graph and cut floors 0.2; isomorphism tolerance 1. Cache-off disables completed-fragment retention, not every internal prepared structure.
- A completed task means its requested seed/cut work finished. Cap evidence remains in the report; completion does not establish exhaustive mapping coverage.
- Reference recovery checks the exact heavy-atom relation including unmatched atoms, modulo endpoint chemical symmetry. These selected reactions and settings do not substitute for the separate 21-miss follow-up cohort.
- The first batch stopped normally after five completed reactions; the remaining three continue separately. No completed task was repeated or overwritten.
- Local memory protocol: the initial tasks use a 4 GiB guard. The running batch is being drained normally; only remaining large-case settings are planned at 16 GiB and 1,800 s per direction. Scientific settings and the 64 MiB fragment-cache budget remain unchanged. Actual guards are retained per task, and paired guard equality is recorded.
- The batches have the same scientific source identity and configuration but two fingerprint implementations. The second hashes compact catalogue records. Both runner hashes are retained below and in the JSON. Equality is checked within a version’s cache pairs; cross-version hash equality is not asserted.
- Audit CPU is excluded from matching and the other scientific stages. Whole-worker CPU and RSS include audit work. The older audit can increase those totals; whole-worker peaks cannot be assigned to matching or used as a clean cross-version memory comparison. RSS is sampled every 0.2 s, alongside the process-lifetime peak; a short unsampled phase is not evidence of zero memory use.

## Provenance

- Scientific source identity: e210faf401a2d33d679c8194a47137237478bbb7074532e4b827b14eeee9d1ca.
- Selection SHA256: 23efb3e68c513b4e345a1aaf65902fb2adba295f30336fc8f8e977bcef83ae68.
- scaling_first5: runner 15aa015b06ad77235075db050206ebea773a30a335fafd6c065699cf079aa7bc; catalogue-signature/v1; result snapshot d44e4af66de33d231aef87caaab0e5426b3742c8f479aa5e7cf8664f4d71dae1.
- scaling_remaining3: runner 4b46f2d7b2a55b8f902defac80b19d2a6f033628327dd31412a7a5995705e043; catalogue-signature/v2; result snapshot 4509666b9ce452e94ced3838580a376873e68dc6c8777a1cf23e9fd56376f87f.

[summary.json](summary.json) preserves all 128 requested cells, per-direction stage measurements, cache counters, completed scope, cap counts, catalogue signatures, reference outcomes and source-input hashes. [scaling.svg](scaling.svg) is the vector plot. This report will be updated from the completed trial grid.
