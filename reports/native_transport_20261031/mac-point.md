# Local Mac point test

Reaction **5735, forward, cut 111**, the exact archived cut that crossed the
2 GiB worker limit. Endpoints contain 387 atoms each. One seed, branch limit 100,
random seed 42, matching tolerance 1.0; all settings and seed order are unchanged.

Machine: Apple M4 Pro, 48 GiB RAM, Python 3.12. One CPU thread per process.
Historical transport: `733adfc`; corrected transport: `e86fee4`, including whole-sequence
sharing and direct immutable-sequence pickle reduction. The old native
binding and Python graph adapter were rebuilt from Git in a temporary build
directory. The scientific C++ growth implementation is unchanged. No cluster
timings are used in this comparison.

Three fresh processes per build. The initial comparison alternated old/new,
new/old, old/new; after removing duplicate island containers the final build
was measured in three more fresh processes. Each process executes the same cut twice to distinguish the first search
from reuse of prepared matching caches. Medians:

| Measurement | Before | After |
|---|---:|---:|
| First search, including native export | 3.576 s | 2.481 s |
| First matching plus symmetry finalization | 3.909 s | 2.815 s |
| Same cut with prepared caches, matching plus finalization | 1.933 s | 0.847 s |
| First-run peak resident RAM | 2.444 GiB | 0.833 GiB |
| Saving the result | 0.247 s | 0.238 s |
| Cross-path grouping and family deduplication for this cut | 0.353 s | 0.485 s |
| Local event decoding | 2.138 s | 1.967 s |
| Exact annotated heavy-atom reference query | 0.912 s | 0.907 s |
| Sum of first matching, saving, deduplication, decoding and verification | 7.649 s | 6.417 s |

Matching plus finalization is 28.0% faster on the first run and 56.2% faster
with prepared caches. First-run peak RAM is 65.9% lower. The sum of measured
stages is 16.1% lower. This sum excludes interpreter startup and checkpoint
reload; it is not an end-to-end complete-reaction benchmark. Independent row
medians do not necessarily add up to the median of their per-trial sums.

**Accepted by the user on 2026-10-31 as the sole production version.**
The measured stage tradeoff remains recorded:
Deduplication costs about 0.132 s more in this sequence of operations. A separate
instrumented run measured 0.138 s in one generation-2 Python garbage collection
during the 0.486 s deduplication step. The cost is included in the timings;
the lower overall runtime does not erase it. The first candidate's saving
slowdown was removed by sharing identical whole island sequences (58 values,
instead of 24,758 containers) and serializing shared blocks directly.

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
runner is [mac_point.py](https://github.com/yunhzou/MAPPA-AAM/blob/f10999c75d7b232e8996cedf996f4d9b3a5e94ef/reports/native_transport_20261031/mac_point.py), which now runs only the approved version.
Historical build preparation and import selection are retained in Git at
`246aa3d`, under this report directory. Those switches and the temporary old
native build were removed after acceptance; the comparison data remain saved.

Full local validation: **654 tests passed in 45.26 s**. This includes native/Python
graph equality, independent compression/decoder oracles, RAM and bundle handoffs,
worker scheduling, recovery and the public CLI. The earlier sandboxed run could
not bind worker sockets or inspect child processes; validation was repeated with
those local operations permitted. One outdated CLI assertion was updated from
the retired refinement record to the current compressed-decoder contract.
