# Golden: tolerance 2, cap 200 and minimum plus 5

A fresh full Golden benchmark was submitted on 4 November 2026 using
[`52c223f`](https://github.com/yunhzou/MAPPA-AAM/commit/52c223f361b26c2cc8e24ac708817bea619ece1b).
It covers **1,851 reactions / 3,702 directions**, with no reused matching,
decoding or verification results. All **5,553 distinct immutable input/reference
blobs** passed checksum verification. Scheduling and operational evidence remain
in the private campaign directory.

| Setting | Value |
|---|---|
| Discovery / fragment self tolerance | 2 / 1 |
| Seeds / random seed | 1 / 42 |
| Branch cap | 200, with recursive subtree pruning |
| Search scope | Both directions; uncut plus full single-edge sweep |
| Graph / cut presence threshold | 0.2 / 0.2 |
| Event / metal event threshold | 0.5 / 0.3 |
| Decoder window | Independent-fragment model minimum through minimum + 5 |
| Input weights | Original prepared bond-weight matrices; no binary transform |
| Requested compute workers | 32 in one persistent pool |
| Requested CPUs, including execution overhead | 40 |
| Compute-worker memory guard | 2 GiB |
| Direction compute budget | 600 seconds; no automatic fresh retries |

Shared reaction caches are enabled. The same pool runs matching, reduction,
decoding and external reference verification as separate stages. All matching
and reduction must finish before decoding starts. Requested capacity does not
establish worker readiness or utilization.

The decoder uses shared local minima to find the smallest composition lower
bound `L`, resolves sweep contributions only within `L + 2` to find `M`, then
reports conditional responses through `M + 5`. This is the accepted
independent-fragment count approximation, not a global physical minimum.
Reference verification queries the retained mapping catalogue separately; it
does not measure reference membership in the reported event window.

Preflight validation: the decoder repair passed 984 Mac tests. The subsequent
benchmark-only policy handoff passed 63 campaign tests and three final targeted
checks. Linux preflight passed **73 tests in 14.75 seconds**, including actual
worker outputs with `event_window == minimum_events + 5`. The runner records
the explicit decoder policy in its input, manifest, plan and saved outputs.

## Regression assessment

No completed full run has identical tolerance-2/cap-200 settings. Compare
reaction-level reference gains and losses against the latest full
[tolerance-1/cap-2,000 benchmark](../golden_placement_20261103/README.md)
(1,831/1,851 recovered), and label the settings when comparing with the older
[tolerance-1/cap-100 run](../golden_retest_20261101/README.md)
(1,830/1,851). The selected multi-setting recovery union is not a full-run
baseline. A changed outcome cannot automatically be attributed to a code
regression because discovery tolerance and pruning also differ.

Report matching plus reduction separately from decoding, with publication and
external verification separate. Preserve failures, timeouts and unknown outcomes
in the full denominator. Compare timings on the same successfully completed
reaction scope and report both compute CPU and elapsed time.

The initial 110-worker allocation remained pending and performed no computation.
At the user's request, it was cancelled and replaced with a 32-worker pool using
a 2 GiB guard. The input index, scientific settings, source identity and compute
budgets are unchanged. The replacement started at 06:30 UTC on 4 November.
All 32 workers registered in 6.16 seconds. At 06:30:48 UTC, 69 of 3,702 directions
had completed matching, with 29 tasks running at that snapshot and no reported
failures. The live plan confirmed tolerance 2, cap 200 and event slack 5; the
progress record confirmed the 2 GiB guard. Reduction, decoding and verification
were not yet complete. No benchmark accuracy or performance conclusion is
available yet.
