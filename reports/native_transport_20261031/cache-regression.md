# Cache-limit experiment rejected

The 512-entry certificate-cache experiment is not regression-free. It was never
installed in production. The temporary native-export patch was discarded;
the accepted native transport and scientific implementation at `3b43327` remain
unchanged.

EnzymeMap reaction 5735 forward, cut 147, run twice in one workspace on the Mac:

| Measurement | Accepted cache | 512-entry cache + shared deferred records |
|---|---:|---:|
| First native growth | 22.60 s | 22.38 s |
| Repeated native growth | 1.44 s | 21.54 s |
| Repeated cut through publication | 4.66 s | 23.92 s |
| Repeated cache hits | 17,171 | 512 |
| Repeated cache misses | 0 | 16,659 |

Both fresh and repeated complete finalized graphs compare equal: 61,756 states,
61,755 transitions, 100 terminals, and equal cut/seed provenance. The regression
is recomputation, not a changed search result. Native growth is measured before
Python export/sharing. Both processes used the same instrumented engine and ran
concurrently, one compute thread each. These diagnostic timings include probes;
they are not a whole-benchmark throughput estimate.

The repeated identical cut deliberately tests cache reuse; normal successful
benchmark cuts are not scheduled twice. Earlier runs of distinct cuts 111, 104,
147 showed similar growth times with and without the cap. That does not remove
the demonstrated regression in a supported reused session. No new full suite or
downstream decoder comparison was completed after this performance gate failed.

Raw events and exact graph comparisons: [cache-regression.json](cache-regression.json).
Reproduce at the measured revision using [the historical memory probe](https://github.com/yunhzou/MAPPA-AAM/blob/f10999c75d7b232e8996cedf996f4d9b3a5e94ef/reports/native_transport_20261031/memory_probe.py) with `--cuts 147 147`, first with the
instrumented accepted cache, then with `--memo-entries 512 --share-deferred`.

The user chose a resource increase instead: **4 GiB per benchmark worker**.
Native algorithms, cache policy, branch caps and saved scientific results did
not change. Raw operational records and allocation settings remain private.

Matching resumes from existing ledgers, explicitly retrying archived memory
failures and retaining completed cuts and attempt history. Old active pools
receive the existing drain signal before replacement. Downstream passes depend
on successful completion of every replacement matching bin.

Validation of the resource-only change passed 23 smoke
checks and 5 benchmark/resume/completed-cut preservation checks. Shell syntax,
CLI help and Git whitespace checks passed. The core source tree has no diff
from `3b43327`; the installed native engines were not rebuilt or replaced.
