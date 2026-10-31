# Local Mac point test

Reaction **5735, forward, cut 111**, the exact archived cut that crossed the
2 GiB worker limit. Endpoints contain 387 atoms each. One seed, branch limit 100,
random seed 42, matching tolerance 1.0; all settings and seed order are unchanged.

Machine: Apple M4 Pro, 48 GiB RAM, Python 3.12. One CPU thread per process.
Historical transport: `733adfc`; corrected transport: `ec0bc91`. The old native
binding and Python graph adapter were rebuilt from Git in a temporary build
directory. The scientific C++ growth implementation is unchanged. No cluster
timings are used in this comparison.

Three fresh processes per build, run sequentially in old/new, new/old, old/new
order. Each process executes the same cut twice to distinguish the first search
from reuse of prepared matching caches. Medians:

| Measurement | Before | After |
|---|---:|---:|
| First search, including native export | 3.576 s | 2.491 s |
| First matching plus symmetry finalization | 3.909 s | 2.825 s |
| Same cut with prepared caches, matching plus finalization | 1.933 s | 0.869 s |
| First-run peak resident RAM | 2.444 GiB | 0.834 GiB |
| Saving the result | 0.247 s | 0.284 s |
| Cross-path grouping and family deduplication for this cut | 0.353 s | 0.499 s |
| Local event decoding | 2.138 s | 1.961 s |
| Exact annotated heavy-atom reference query | 0.912 s | 0.905 s |
| Sum of first matching, saving, deduplication, decoding and verification | 7.649 s | 6.477 s |

Matching plus finalization is 27.7% faster on the first run and 55.1% faster
with prepared caches. First-run peak RAM is 65.9% lower. The sum of measured
stages is 15.3% lower. This sum excludes interpreter startup and checkpoint
reload; it is not an end-to-end complete-reaction benchmark. Independent row
medians do not necessarily add up to the median of their per-trial sums.

**The candidate does not yet meet a strict no-stage-regression requirement.**
Saving costs about 0.037 s more, and deduplication about 0.146 s more. These
regressions remain open; the lower overall runtime does not erase them.

All **12 complete finalized search graphs** compare equal, including contexts,
24,758 states, 24,757 transitions, 100 terminal states, fragment choices,
correlations, joins and stop evidence. All **six catalogue and complete decoded
records** also compare equal; each retains 100 families. This comparison checks
the records, not merely representative mappings or event counts.

Reference verification is **not recovered** in both versions under exact
annotated heavy-atom indices without endpoint-symmetry relaxation. That is a
statement about this one cut's retained families, not the complete reaction
sweep. No scientific result difference was observed.

Raw timings, comparison assertions and machine metadata:
[mac-point.json](mac-point.json). The exact small input artifact is
`mac-point-input.pkl.gz`; it contains trusted internal records. The diagnostic
runner is [mac_point.py](mac_point.py), and the historical build preparation is
[build_mac_baseline.py](build_mac_baseline.py). The historical import selection
exists only in this test, not in the production API.
