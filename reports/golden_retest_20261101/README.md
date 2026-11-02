# Golden rerun: accuracy and time

Commit `096a010961f53944ddd4ddc433bd8d56f282a2f2` completed all **1,851 reactions /
3,702 directions**, including matching, reduction, decoding and reference
verification. There were no failed
tasks, retries, unknown verification results or cache errors.

## Reference recovery

| Measure | Previous Golden pool run | Current run |
|---|---:|---:|
| Recovered reactions | 1,830 / 1,851 (98.87%) | **1,830 / 1,851 (98.87%)** |
| Recovered directions | 3,617 / 3,702 | **3,618 / 3,702** |
| Unknown / failed verification | 0 | **0** |

The same 21 reactions remain unrecovered. No previously recovered reaction or
direction was lost. **Case 1351 reverse** is the sole changed outcome, moving from
not recovered to recovered; its forward direction already recovered the reference.
Recovery means the annotation exists among the retained alternatives, not that
an automatically selected first mapping is correct. A missed capped search does
not prove that no valid mapping exists.

## Time and memory

CPU seconds below are summed across workers and exclude result-publication CPU:

| Stage | Previous compute CPU seconds | Current compute CPU seconds | Change |
|---|---:|---:|---:|
| Matching | 8,378.07 | **5,136.07** | **−38.7%** |
| Reduction | 2,822.63 | 2,913.60 | +3.2% |
| Decoding | 7,996.45 | 8,044.71 | +0.6% |
| Reference verification | 1,306.12 | 1,334.67 | +2.2% |
| All stages | 20,503.27 | **17,429.04** | **−15.0%** |

The improvement is concentrated in matching; the other stages did not improve
in this point comparison. Both runs executed exactly **177,552 matching cuts**.
Matching's sampled worker RSS fell from **414.2 to 364.0 MiB**. The current
pipeline peaked at **990.75 MiB** sampled worker RSS under the unchanged **2 GiB**
guard. Sampling is not a proof that every transient allocation was observed.

The completed-fragment cache recorded a **91.83% hit rate**, no baseline
recovery loads, no errors, and at most **40.95 MiB** retained in a snapshot.

## Local point test of the previous failure

Golden **1561 forward** previously exceeded 2 GiB while publishing decoded
output. A fresh Mac test at report commit `26f0a5b` used the same scientific code
as this full rerun and retained the 2 GiB guard.

The exact historical catalogue now decodes and serializes successfully in
**11.66 seconds**, peaking at **1.005 GiB** sampled RSS. It contains the original
14,132 families across 192 cuts. This diagnostic measures publication using the
production serializer settings without retaining its output.

A separate fresh run started with empty results and caches on **one persistent
Mac worker**. All 192 cuts, reduction, decoding, durable publication and reference
verification completed in **49.09 seconds**, with no failures or retries. Peak
sampled worker RSS was **1.145 GiB**. The complete 45-heavy-atom reference was
recovered. Reading the saved result confirmed a v4 decoded record with 14,130
family tables and 13 shared responses.

| Fresh local stage | Compute elapsed seconds | Publication seconds |
|---|---:|---:|
| Matching, summed over 192 cuts | 14.70 | 6.34 |
| Reduction | 13.44 | 0.25 |
| Decoding | 10.76 | 0.70 |
| Verification | 1.80 | 0.01 |

The 49.09-second total also includes pool startup and stage coordination. The
isolated replay and persistent-worker pipeline have different resident memory
histories; their RSS values need not match. These Mac timings are not directly
comparable to the distributed-run timings above. [Point-test evidence](point-1561.json)
contains input checksums, commands, stage audits, verification and the complete
saved-catalogue replay log. No algorithm or configuration was changed.

## Comparable scope and limits

The comparison uses the previous [Golden pool run](../golden_pool_20261101/README.md),
whose matching commit was `6397d461156f86b877efeb46ea40e4fe6056fc18`, and its latest
v4 decoder measurements from `2f009516016efe3e1dbf97079908ac2f712e0327`.
The current rerun includes the intervening matcher repairs as well as the latest
support-check optimization; it does not isolate that final optimization alone.

Both runs use the same 1,851 archived inputs, original five batch memberships and
direction order, one seed, random seed 42, cap 100, full bidirectional cut sweeps,
matching tolerance 1, graph/cut floor 0.2, and event thresholds 0.5/0.3. Golden
verification preserves the complete heavy-atom relation, including unmatched
atoms, modulo chemically colored endpoint symmetries. Its per-family query limit
remains 10 seconds. All 5,553 unique prepared input/reference blobs were verified
against their saved checksums. The new run imported no old results or caches.

All matching and reduction finished before decoding began. Private execution
configuration and raw operational records are excluded from this public report.

There is no clean previous end-to-end wall-time baseline: the original run
included a decoder failure and later repair, and the latest decoder-only pass
reused one previously completed point result. Its CPU figure above explicitly
includes that point, covering all 3,702 outputs; its elapsed window
excludes the earlier point run. CPU totals
provide the stage comparison; elapsed time is not presented as an isolated
algorithm speedup. These are single full runs, not repeated timing trials.

[summary.json](summary.json) retains a scientific projection of the source
identity, input checksums, stage CPU and memory measurements and comparison.
[reference-results.json.gz](reference-results.json.gz) retains every case and
direction outcome, stage measurements and verification evidence checksums;
private artifact locations are omitted. Original records are retained privately.
The audit reads saved results and performs no further scientific computation.
