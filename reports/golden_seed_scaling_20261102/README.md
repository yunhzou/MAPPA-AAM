# Golden seed scaling and recovery — local results

Final snapshot: 2026-11-02T09:19:02+00:00. **Seven reactions have complete balanced grids: 112 trials.** Five additional case 1799 trials also finished, giving **117 executed trials out of 128 planned cells**. All 117 recovered their references. The 11 optional cells were not started after the active trial drained: seven cases satisfy the requested five-to-ten-case scope, and memory repair took priority. All 58 cache-on/off pairs have identical catalogue and scope hashes, per-cut records, graph counts, cap evidence and verification outcomes. Missing cells are not imputed.

![CPU per seed relative to the one-seed baseline](scaling.png)

The eight reactions were selected at recorded matching-cost quantiles using [selection.json](selection.json). Selection costs came from the earlier one-seed, cap-100 study; these measurements use cap 2,000. This cost-stratified sample does not estimate benchmark-wide speed or accuracy.

## Fixed seven-case comparison

Cases **1001, 1850, 1022, 785, 912, 1518 and 709** have all 16 requested settings complete. Every row sums the same seven reactions and both directions. There is one trial per setting, with no repeats or uncertainty estimates.

| Seeds | Matching on CPU s | Matching off CPU s | Off/on | Through reduction on CPU s | Through reduction off CPU s | Off/on |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 35.022 | 81.170 | 2.32× | 149.255 | 194.869 | 1.31× |
| 2 | 79.929 | 183.974 | 2.30× | 322.000 | 424.598 | 1.32× |
| 5 | 198.384 | 443.270 | 2.23× | 782.807 | 1024.824 | 1.31× |
| 10 | 396.839 | 880.683 | 2.22× | 1531.983 | 2026.124 | 1.32× |

Cache-off/on CPU ratios are **2.22–2.32×** for matching and **1.31–1.32×** through saved-symmetry finalization and family reduction in this completed cohort. These ratios compare the same algorithm and requested work, with completed-fragment retention on or off. They exclude preparation, external reference verification and audit. **No event decoding was run.** Millisecond-scale points remain sensitive to one-run timing noise.

All seven reactions and all 14 directions recover their references at every completed seed/cache setting. Cap evidence remains present: completed work does not establish exhaustive mapping coverage.

## Work across seed counts

**This cohort does not show sublinear matching CPU overall.** With the cache enabled, ten seeds take 11.33 times the one-seed matching CPU. The cache still reduces matching CPU relative to disabling completed-fragment retention, but that is a separate comparison.

The figure shows **T(N)/(N × T(1))** for each reaction and cache setting. A value below one means less average CPU per seed than at one seed; the horizontal line at one marks linear scaling. Each reaction sums both directions. Curves share the same vertical scale, and cache-on and cache-off are each normalized to their own one-seed time.

| Seeds | Matching on CPU s/seed | Matching off CPU s/seed | On normalized | Off normalized | Through reduction on CPU s/seed |
|---:|---:|---:|---:|---:|---:|
| 1 | 35.022 | 81.170 | 1.000 | 1.000 | 149.255 |
| 2 | 39.964 | 91.987 | 1.141 | 1.133 | 161.000 |
| 5 | 39.677 | 88.654 | 1.133 | 1.092 | 156.561 |
| 10 | 39.684 | 88.068 | 1.133 | 1.085 | 153.198 |

Reaction sizes span **9–104 explicit atoms per endpoint**. Cut counts below include each direction’s no-cut task. This cost-stratified selection covers different sizes, but atom count alone does not explain search cost.

| Case | Explicit atoms R/P | Cut tasks R/P | On, 1 seed CPU s | On, 10 seeds CPU s | On, T(10)/(10 T(1)) | Off/on at 10 seeds |
|---:|:---:|:---:|---:|---:|---:|---:|
| 1001 | 9/9 | 8/8 | 0.0034 | 0.0094 | 0.277 | 1.04× |
| 1850 | 26/23 | 26/24 | 0.0203 | 0.0928 | 0.456 | 1.24× |
| 1022 | 26/26 | 28/28 | 0.1448 | 2.0569 | 1.420 | 1.19× |
| 785 | 40/40 | 41/41 | 0.3064 | 2.0633 | 0.673 | 1.30× |
| 912 | 44/44 | 42/44 | 3.2001 | 40.1999 | 1.256 | 1.43× |
| 1518 | 104/46 | 104/50 | 7.1449 | 87.0802 | 1.219 | 2.62× |
| 709 | 98/98 | 102/102 | 24.2024 | 265.3365 | 1.096 | 2.22× |

Three individual reactions have lower CPU per seed at ten seeds, including the two smallest timings. Four have higher CPU per seed. These single observations do not establish a universal scaling law. The increased number of seed searches, cache calls and retained families remains visible despite high cache-hit fractions.

| Seeds | Fragment calls | Full hits | Logical extensions | Reused extensions | Evictions | Saved states | Families |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 5,678,515 | 97.03% | 2,459,295 | 82.06% | 96,665 | 958,444 | 245,953 |
| 2 | 12,399,940 | 97.30% | 5,403,975 | 83.93% | 248,382 | 1,977,113 | 494,553 |
| 5 | 30,476,642 | 97.70% | 12,963,418 | 85.58% | 589,173 | 4,906,678 | 1,186,799 |
| 10 | 56,934,197 | 97.78% | 24,706,626 | 86.24% | 1,119,778 | 9,504,322 | 2,304,520 |

Counts sum the same fixed cohort. Saved states and families are summed across directions, not deduplicated across different reactions. High hit fractions do not remove the cost of cache calls, new seed orders, finalization or reduction. These measurements do not establish eviction as a bottleneck.

The raw result checksums, all 58 complete cache pairs, and all 112 balanced task scopes were rechecked for this analysis. Git objects for the matching/search code, native kernels and maintained runner are identical between the measured study commit `1134e99` and audited commit `06f7941`. The later decoder changes are outside this timing scope. No new matching run was needed.

## Timing boundaries

The cache-on stage totals below use the same seven cases. CPU values sum worker execution; summed direction wall times in the JSON are not parallel reaction latency. Verification and fingerprinting are external to the matching-through-reduction comparison.

| Seeds | Preparation CPU s | Matching CPU s | Finalization CPU s | Reduction CPU s | Verification CPU s | Audit CPU s (excluded) |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.046 | 35.022 | 45.320 | 68.913 | 8.216 | 51.121 |
| 2 | 0.041 | 79.929 | 92.227 | 149.845 | 11.332 | 96.971 |
| 5 | 0.039 | 198.384 | 208.072 | 376.351 | 173.868 | 240.880 |
| 10 | 0.043 | 396.839 | 395.438 | 739.707 | 523.181 | 465.310 |

## Scope and checks

- Hardware: Apple M4 Pro. Every case/direction/seed/cache setting starts in a fresh worker. The fragment cache stays resident across that direction’s cuts and seeds. There is no cross-task cache import.
- Seeds: 1, 2, 5 and 10; random seed 42; branch cap 2,000; full single-edge cut sweep; graph and cut floors 0.2; isomorphism tolerance 1. Cache-off disables completed-fragment retention, not every internal prepared structure.
- Reference verification checks the exact heavy-atom relation including unmatched atoms, modulo endpoint chemical symmetry. The 21-miss follow-up is reported separately below; neither study is a full 1,851-case rerun.
- The six named scaling result groups below supply the grid. Preflight and startup-only results are excluded. Recovery is a separate cohort. The first two groups drained normally before their remaining work continued in separate groups. No completed grid task was repeated or overwritten.
- Local resource protocol: 105 completed tasks used 4 GiB and 600 s per direction; 12 further tasks used 16 GiB and 1,800 s. No scaling task reached its guard. Case 709, five seeds, forward is the one cache pair with unequal guards: on used 4 GiB/600 s, off used 16 GiB/1,800 s. Their catalogue and work-scope evidence still match. The other 57 completed pairs have equal guards. These local limits do not alter cluster reservations or scientific settings; the enabled fragment-cache budget remains 64 MiB.
- All groups have the same scientific source identity, configuration and overlapping input hashes, with two fingerprint formats and three runner hashes. V2 hashes compact catalogue records; the later runner adds task-selection controls. Equality is checked within each cache pair’s fingerprint format; cross-format hash equality is not asserted.
- Audit CPU is excluded from matching and the other scientific stages. Whole-worker time and RSS include verification and audit, and the older audit can increase them. Whole-worker peaks are not matching-only memory measurements or clean cross-runner comparisons. RSS is sampled every 0.2 s alongside process-lifetime peaks; an unsampled short phase is not zero memory.

## Extra case 1799

All four one-seed trials and the two-seed cache-on forward trial finished and recovered the reference. Both one-seed cache pairs have identical catalogue and search-scope evidence. These five trials are excluded from the fixed seven-case aggregates; the remaining 11 optional trials were not started.

| Seeds | Cache | Direction | Matching CPU s | Through reduction CPU s | Reference |
|---:|---|---|---:|---:|---|
| 1 | on | forward | 63.259 | 224.173 | Recovered |
| 1 | on | reverse | 10.981 | 24.115 | Recovered |
| 1 | off | forward | 333.100 | 497.235 | Recovered |
| 1 | off | reverse | 25.423 | 38.792 | Recovered |
| 2 | on | forward | 131.841 | 442.723 | Recovered |

The two-seed forward trial overlapped a separate local decoder memory diagnostic from 09:02:18 to 09:04:34 UTC. Its timing is not an isolated-host comparison. The balanced seven-case study and case 1799’s one-seed trials preceded that overlap.

## Recovery of 21 previous misses

**Eight of 21 cases now recover a reference.** Twelve completed both directions without recovery; case 590 remains unresolved because its forward direction reached the memory guard. Across all 42 directions, 40 finished: 13 recovered and 27 did not. Two reached the memory guard; there were no unknown verdicts or timeouts. All retained outcomes use 10 seeds, cap 2,000 and the completed-fragment cache.

| Case | Forward | Reverse | Case outcome |
|---:|---|---|---|
| 7 | Not recovered | Not recovered | Not recovered |
| 19 | Not recovered | Not recovered | Not recovered |
| 20 | Recovered | Recovered | Recovered |
| 44 | Recovered | Recovered | Recovered |
| 590 | Memory limit | Not recovered | Unresolved |
| 867 | Recovered | Recovered | Recovered |
| 964 | Recovered | Recovered | Recovered |
| 986 | Not recovered | Not recovered | Not recovered |
| 1034 | Recovered | Recovered | Recovered |
| 1227 | Recovered | Not recovered | Recovered |
| 1228 | Not recovered | Not recovered | Not recovered |
| 1285 | Not recovered | Not recovered | Not recovered |
| 1314 | Not recovered | Not recovered | Not recovered |
| 1358 | Not recovered | Not recovered | Not recovered |
| 1377 | Not recovered | Not recovered | Not recovered |
| 1475 | Not recovered | Not recovered | Not recovered |
| 1553 | Not recovered | Not recovered | Not recovered |
| 1568 | Not recovered | Recovered | Recovered |
| 1606 | Memory limit | Recovered | Recovered |
| 1739 | Not recovered | Not recovered | Not recovered |
| 1740 | Not recovered | Not recovered | Not recovered |

Case 590 forward stopped during matching after 56/58 cuts, with sampled peak 8.040 GiB. Case 1606 forward stopped during reduction after 131/153 cuts, with sampled peak 8.003 GiB; its reverse direction recovered the reference. Partial scope and consumed stage measurements are retained in the summary; full raw outputs remain archived separately. Memory limits are not negative reference verdicts.

This follow-up used disjoint case groups, at most three simultaneous local one-worker jobs, 8 GiB and 1,800 s per direction. It initially overlapped decoder profiling; recovery timings are not part of the controlled scaling comparison. No event decoding or automatic retries were performed. Catalogue fingerprinting was skipped because this cohort has no cache-off control.

The earlier miss list comes from the frozen one-seed, cap 100 validation. This follow-up changes both seed count and cap, so the eight recoveries cannot be attributed to seeds alone and do not establish a new full-benchmark score. A completed capped search without recovery also does not prove that no valid mapping exists.

## Provenance

- Scientific source identity: e210faf401a2d33d679c8194a47137237478bbb7074532e4b827b14eeee9d1ca.
- Selection SHA256: 23efb3e68c513b4e345a1aaf65902fb2adba295f30336fc8f8e977bcef83ae68.
- scaling_first5: runner 15aa015b06ad77235075db050206ebea773a30a335fafd6c065699cf079aa7bc; catalogue-signature/v1; result snapshot d44e4af66de33d231aef87caaab0e5426b3742c8f479aa5e7cf8664f4d71dae1.
- scaling_remaining3: runner 4b46f2d7b2a55b8f902defac80b19d2a6f033628327dd31412a7a5995705e043; catalogue-signature/v2; result snapshot dfa415b87a4a983b4a64493b714cb2d1b49da9e5e83a348fe98f64be003140bf.
- scaling_709_seed5_on_reverse: runner af0c5d4c5e4c82f0f8b745201f56c3c1bafdd99a574b0b0fb805ca15fef49b19; catalogue-signature/v2; result snapshot 688e13812c991cb0dc9400d93e0e09fb70bbd7f030284ce739dcd1ebb357c9e8.
- scaling_709_seed5_off: runner af0c5d4c5e4c82f0f8b745201f56c3c1bafdd99a574b0b0fb805ca15fef49b19; catalogue-signature/v2; result snapshot 586476ea23df3e3abb121e13d09e2c0763b36809c4757d3833ffd8b91000d077.
- scaling_709_seed10: runner af0c5d4c5e4c82f0f8b745201f56c3c1bafdd99a574b0b0fb805ca15fef49b19; catalogue-signature/v2; result snapshot af63ad14e50a2be3d67c3427ca250b517b6b7c83ad8ba2bdd7073cf9777bab66.
- scaling_1799: runner af0c5d4c5e4c82f0f8b745201f56c3c1bafdd99a574b0b0fb805ca15fef49b19; catalogue-signature/v2; result snapshot 0b53a4b5778acaeb0b3c8d949593671a91ccb8ce69900a30b0017ecd77fd4276.

[summary.json](summary.json) preserves all 128 requested cells, per-direction stage measurements, cache counters, completed scope, caps, catalogue signatures, reference outcomes, actual guards and source-input hashes, plus all 42 recovery directions and their partial failure records. [scaling.svg](scaling.svg) is the vector plot. This planned local study is closed; unrun optional cells remain explicit.
