# MAPPA decoder optimization

All 140 saved searches decode completely with the same 300 reaction-local event classes.
All 110,477 previously completed per-family class-support records match.
The two previously interrupted runs now finish. The saved catalogue still has 120,052 families.

| CPU seconds | Before (same 138 cases) | After (same 138 cases) | After (all 140) |
|---|---:|---:|---:|
| Mean | 5.666 | 3.196 | 4.647 |
| Median | 0.345 | 0.276 | 0.288 |
| 95th percentile | 32.013 | 20.332 | 24.361 |
| Total | 781.942 | 441.008 | 650.544 |

New timing source: a fresh search-plus-decoding run; full stage timings are preserved in the accompanying end-to-end report.

Paired decoder CPU reduction: 43.6% (1.77x faster in aggregate).
Case 25 now takes 204.621 decoder CPU seconds; case 31 takes 4.914.
Their old attempts hit the 300-second **whole-process** watchdog; these are censored outcomes, not completed decoder times.

Four independent reaction processes, one numerical thread each; decoder CPU is not divided by four.
Both decoder measurements include lazy initialization and per-family journaling; search and catalogue construction are excluded from decoder stage times. The summary records whether the new decoding followed a fresh search or a saved checkpoint.
Thresholds (.5/.3), matching tolerance, saved families, event windows, and exact class identity are unchanged. No candidate cap or bijection/group-element enumeration was added.

Changes:
- Reuse exact raw-weight pair response tables and symbolic subset selectors.
- Express event totals as Boolean cardinality constraints instead of arithmetic search.
- Omit a redundant global injectivity constraint: each pool/group action is already a permutation.
- Share source edges across families; screen absent bonds in NumPy before scalar witness validation.

The 16,384-entry local-table cache is a memory policy, not a solution limit: evicted entries are recomputed exactly.
All original families and symmetry-query information remain available.

Reproduction:

```sh
python bench/benchmark_saved_decoding.py --archives /path/to/final_end_to_end_run --output /path/to/new_decode_run --workers 4
python bench/report_decoder_optimization.py --run /path/to/new_decode_run --output reports/decoder_optimization_20260916
```

Source hashes, per-reaction class IDs, timings, and validation flags are in `summary.json`; exact per-family results are in `family-journals.jsonl.gz`.
