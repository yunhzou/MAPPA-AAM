# Golden rerun: accuracy and time

Commit `096a010961f53944ddd4ddc433bd8d56f282a2f2` completed all **1,851 reactions /
3,702 directions**, including matching, reduction, decoding and reference
verification. All ten DFW allocations exited successfully. There were no failed
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

The complete fresh run took **11 minutes 56 seconds**, from the first allocation
start to the last allocation finish. Matching and reduction occupied a
**5 minute 25 second** allocation window; decoding and verification occupied
**5 minutes 57 seconds**, with 34 seconds between the windows. These durations
include staggered allocation starts and pool startup, and exclude preparation,
compilation and the subsequent read-only audit.

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

Each batch used all 100 matching workers and reached 100 overlapping tasks.
Two pools could run concurrently; observed global task overlap peaked at 192
during matching and 200 during decoding. Average total worker CPU consumption
within each matching batch was 16.5–18.2 cores. Worker concurrency does not mean
all allocated CPUs were continuously busy. The completed-fragment cache recorded
a **91.83% hit rate**, no baseline-recovery loads, no errors, and at most
**40.95 MiB** retained in a worker snapshot.

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

The same four-node / 100-worker layout was retained: 28 allocated CPUs and
54 GiB per node (112 CPUs per pool), with a 2 GiB guard per compute worker. Five disjoint batches
ran with at most two pools concurrently. All matching and reduction finished
before decoding began. Jobs were `21159701` and `21159702`.

There is no clean previous end-to-end wall-time baseline: the original run
included a decoder failure and later repair, and the latest decoder-only pass
reused one previously completed point result. Its CPU figure above explicitly
includes that point, covering all 3,702 outputs; its 225-second allocation window
does not include the earlier point run. The old matching-only 373-second window
also differs in scope from the new matching-plus-reduction window. CPU totals
provide the stage comparison; elapsed time is not presented as an isolated
algorithm speedup. These are single full runs, not repeated timing trials.

[summary.json](summary.json) retains source/build identities, input hashes,
launch commands, Slurm accounting, per-stage/per-batch CPU and memory evidence,
cache statistics and the complete comparison. [reference-results.json.gz](reference-results.json.gz)
retains every reaction/direction outcome, compact per-case stage metrics and
references to the original verification evidence. The audit reads saved results;
it performs no further scientific computation.
