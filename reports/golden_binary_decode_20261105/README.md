# Full Golden decoding with binary fragment symmetry

The final decoder at `deadf374a97de7b171e86eee85e07f975deb6d38` completed all
**1,851 reactions / 3,702 directions**, with **zero failed directions and zero
retries**, under the unchanged **2 GiB / 600-second** per-direction limits.
The full repository suite passed **1,127 tests in 175.00 seconds**; Linux
preflight passed **179 tests in 13.39 seconds**. The final 36-direction pilot
passed completely before the full restart.

This is decoding of the saved full Golden matching catalogue: discovery tolerance
2, one seed, branch cap 200 and the full single-edge sweep. Matching is not
repeated. The decoder retains the agreed independent ordinary-boundary model,
binary fragment symmetry, original WBO event thresholds and minimum-plus-five
window. Historical repair attempts follow the final measurements below.

## Final full-set timings

All 32 workers registered in 3.061 seconds; queue snapshots showed 32 distinct
workers simultaneously processing 32 distinct directions. The allocation
completed successfully in **305 seconds (5m05s)**; the pool measured **300.542
seconds**, including startup and publication. Conversion was reused and is
reported separately below.

| Decoding measure | Median | 90th percentile | 95th percentile | 99th percentile | Maximum |
|---|---:|---:|---:|---:|---:|
| One direction, wall seconds | 1.135 | 5.296 | 8.405 | 20.885 | **137.192** |
| Both directions summed per reaction, wall seconds | 2.410 | 11.169 | 16.754 | 36.724 | **147.157** |

The maximum direction is **1606 reverse**; the maximum two-direction sum is
reaction **1606**. Direction timings include catalogue loading, decoding and
record construction. The reaction sum is the total work duration of both
directions, which can execute concurrently. Publication is separate: median
0.039 seconds and maximum 3.001 seconds per direction.

Total decoding CPU time is **9,153.399 seconds (2.543 CPU hours)**, with another
239.553 CPU seconds for publication. Decoding CPU divided by the 32-worker pool
capacity over its measured runtime is **95.2%**. The allocation reserved 40 CPUs
including execution overhead; this percentage uses the 32 compute workers.
The largest sampled task RSS was **1.574 GiB**, below the 2 GiB guard.
[Per-direction measurements and full distributions](timings.json) retain every
direction, its attempt count, guard, compute/publication timing and sampled memory.

All 3,702 terminal ledger rows, task audits, input hashes, source/engine identities
and empty saved queues passed the maintained timing audit. A separate read-only
audit verified all 3,702 published outputs: their complete matching catalogues,
metadata and configurations match their inputs. Both directions are present for
every one of the 1,851 reactions, with exactly one attempt per direction.

The outputs report **3,701 proved model minima and one best-found model minimum**.
These refer to the selected fragment symmetries and agreed approximate boundary
model, not a physical global minimum over every discovery. All outputs retain
that distinction and the minimum-plus-five query window. External reference
verification was not rerun; these measurements do not provide a new
reference-recovery score. The read-only output audit took a separate 380-second
allocation on eight workers and is excluded from decoding timings.

Pilot comparisons do not establish a universal per-case speedup. On the 23
directions completed by the early preparation-only pilot, final-pilot CPU fell
from 568.161 to 464.135 seconds in total, while individual cases varied. The
immediately preceding endpoint pilot completed 35 directions; those directions
took 834.529 seconds there and 849.895 seconds in the final pilot. Its remaining
120-second timeout, 607 forward, completed in 45.652 seconds. Pilots used
separate allocations; the final full pass above is the throughput measurement.

## Input conversion and earlier attempts

An explicit benchmark input conversion removes the retired fragment
self-tolerance filter and preserves the saved binary generators, original WBOs,
discovery constraints, atom sets, placement actions and branch histories. The
original catalogues remain unchanged. Converted inputs do not constitute a
fresh benchmark of the current matcher. Conversion evidence retains both input
and output checksums and the original matching configuration.

Conversion is independently checked against the producing revision's complete
record export for Golden 625 forward. All fields match after the declared schema
and symmetry-policy changes. Full decoding and a save/load check completed in
0.051 seconds locally. Golden 569 forward exceeded a separate 240-second local
watchdog during full decoding; this is not a fragment-only timing.

The full cluster test retains event thresholds 0.5 / 0.3 and minimum plus 5,
using the current default five-second fragment optimization allowance. Each
direction has a 600-second compute watchdog and a 2 GiB memory guard, matching
the original full campaign's limits. Conversion, decoding compute and output
publication are measured separately. Timings will retain failed directions in
the denominator and identify the largest completed decoding time separately
from censored timeout measurements. No matching rerun or recurring monitor is
started by this experiment.

## Stopped at the observed performance issue

The user directed that an observed issue must stop broad validation immediately.
The cluster allocation was cancelled during input conversion, before any cluster
decoding started. The source matching campaign is unchanged. There is no full
Golden decoding result or measured full-set maximum from this experiment.

A bounded local profile of 569 forward spent 59.35 of 60 seconds in symmetry
searches for composed historical actions. The first slow contribution already
had lower bound **1** and an actual witness with **1** event, but the decoder
entered another minimum search instead of accepting that certificate. Its 793
variable bond terms each carried all 138 permutation positions from two saved
actions as their dependency scope. Repeated traversal of these scopes consumed
41.16 profiled seconds. These are two concrete sources of unnecessary work;
they do not establish that every slow direction has the same cause.

A process-local diagnostic accepted attained local bounds directly, without
changing production code. It avoided 273 redundant local searches but still
reached its 60-second watchdog. The next isolated branch had five final fragment
pieces with proved minima `0, 1, 0, 0, 0`, fixed boundary contribution 4, and one
sweep bond. Its lower/attained bounds were 5 and 7. Resolving that sweep still
constructed a search over 795 bond terms, nine saved actions and 154 permutation
positions. The branch program did not retain the local response's event-invariant
compression: its invariant subgroup is explicitly the identity. This explains
why removing only the first redundant search does not resolve the integration
bottleneck. These diagnostics do not establish that local minima can be added
without checking their shared choices or endpoint conditions.

At that stopping point, production scientific code was unchanged. The subsequent
repair and bounded local validation are recorded below. Fragment-level beta
measurements alone do not establish full branch-decoding performance.

Conversion and campaign validation passed 66 local tests; Linux preflight passed
44 tests. The reusable timing summary was checked against actual worker outputs
in one further integration test. It remains available for a later completed pass
as `python bench/performance.py decoding-timings CAMPAIGN OUTPUT.json`; it refuses
an undrained or incomplete pass. No full-set timing distribution is reported here.

## Prepared-response reuse repair

The repair accepts an attained local lower bound immediately. Branch assembly
now references prepared fragment cost projections instead of expanding them
back into raw bond terms. Weight-preserving placements retain internal costs;
they still move endpoints. Endpoint dependencies follow the ordered actions and
contain only positions that can determine their images. Local event-invariant
compression is reused only after protecting the possible images read by other
costs, including the sweep. Shared actions and injective occupancy remain joint.
No event model, allowed mapping relation, tolerance, time allowance or cap was
changed, and no alternative decoder was added.

All six formerly timed-out directions completed full local decoding and record
construction with a 60-second point-test watchdog. These are sequential Mac
measurements, including the minimum/window preparation and output record, but
excluding the separately measured explicit input conversion and save/load check.

| Direction | Conversion (s) | Decoding and record (s) | Model minimum |
|---|---:|---:|---:|
| 569 forward | 0.041 | 5.868 | 7 |
| 569 reverse | 0.104 | 2.313 | 7 |
| 1361 forward | 0.176 | 24.616 | 7 |
| 1361 reverse | 0.633 | 9.579 | 7 |
| 575 reverse | 0.566 | 2.910 | 11 |
| 837 forward | 0.039 | 4.941 | 6 |

Success is **6/6**, with no point-test timeout; the measured maximum is **24.616 s**
on 1361 forward. Each model minimum was proved within the retained decoded orbit.
All six matching catalogue records remained unchanged, and decoded record
save/load/save equality passed. Process peak RSS reached 784.3 MiB across the
sequential process, including periods retaining both original and restored
outputs; this is not a separate per-direction peak or an enforced memory guard.

The directly comparable pre-repair 569-forward test exceeded 240 seconds on the
same converted input; it now takes 5.868 seconds, over 40 times faster than that
timeout lower bound. The earlier six cluster timeouts used the former decoder
and symmetry policy, so their 600-second limits are not an isolated speedup
baseline for this repair. No full Golden maximum, recovery result or universal
no-regression claim follows from these six cases.

Critical tests require a ring's sweep search to contain only the two endpoint
costs, retain an attained contextual minimum without another search, and check
all small ordered fragment/placement choices against independent scalar event
counts. Queries retain endpoint/internal dependence, occupancy, unattainable
count gaps, differing source/target connectivity and saved-response restoration.

The final fixed-source repository suite passed **1,091 tests in 162.61 seconds**.
An earlier interrupted suite is excluded: editing a core file during that run
caused its execution-layer source/cache identity checks to reject worker tasks.
The affected test passed independently, and the complete suite above ran with
unchanged source.

## Full-set restart

The user authorized a fresh decoding pass with the repaired version. Linux
preflight passed **80 tests in 6.89 seconds**. The restart covers all 1,851
reactions / 3,702 directions, with 32 compute workers and the same 2 GiB
per-worker guard, 600-second direction allowance and minimum-plus-5 window.
The source and engine remain fixed throughout the run. Output goes into a
fresh campaign; the original matching catalogues and stopped attempt are kept.

The allocation started on 2026-11-05 at approximately 07:40 UTC. Explicit input
conversion precedes worker startup; reserved worker capacity is not decoding
utilization during that preparation. The full pass was stopped for the further
performance issue below. No full-set result or timing distribution is claimed.

## Further overhead found during the restart

All 3,702 catalogues converted successfully in 813.129 wall / 794.835 CPU
seconds, with 402.6 MiB peak RSS. Their checksums, complete two-direction scope
and declared settings passed inspection; none had zero families. Conversion is
separate from decoding. All 32 workers registered in 2.898 seconds, and saved
queue observations confirmed 32 distinct workers executing 32 distinct directions.

Case 862 exceeded four minutes per direction, versus 10.309 / 6.320 seconds in
the original run. The changed binary-symmetry policy prevents attributing that
comparison solely to implementation, but the slowdown warranted investigation.
New admission was stopped after 304 directions entered the pool. Accepted tasks
drain under their original limits; failures and partial results are retained.

A 45-second profile of 862 forward identified repeated work: 19.043 seconds
finding the minimum remaining atom pair, and 14.559 seconds checking placement
group membership. The first repeatedly scanned a shrinking set even when most
pair orbits were singletons. The second rebuilt group machinery even for the
identity placement. These are preparation/composition costs, not evidence that
an exhaustive bond-event table is required.

The follow-up visits potential pair representatives in sorted order once and
accepts the identity directly, while retaining membership checks for every
nonidentity permutation. Pair-orbit ordering and the represented choices are
unchanged. Independent small-group tests check complete pair orbits and reject
an isolated swap when only the correlated joint swap is allowed. Identical
80- and 160-atom singleton-pair results took 0.1133 / 1.7752 CPU seconds before,
and 0.00189 / 0.02047 seconds afterward (60–87 times faster). These component
measurements are not full-catalogue speedups.

The unchanged follow-up source passed 1,098 repository tests: 1,037 in the first
suite and 61 worker tests rerun with process/socket access after macOS sandbox
permission failures. The 108 focused symmetry/decoder tests also passed. Full
catalogue point validation and a completed full-set pass remain outstanding.

The stopped restart ultimately drained with 291 completed directions and thirteen
600-second compute timeouts; the other 3,398 directions were never admitted.
There were no memory-limit failures or retries. A separate 36-direction pilot
of the pair/identity fix completed 23 directions and reached its stricter
120-second watchdog in thirteen. Case 862 finished in 47.915 / 33.493 seconds;
the largest completed pilot direction was 1654 reverse at 73.394 seconds.
Pilot failures are censored at 120 seconds, not failures of the original
600-second allowance. The full pass stayed stopped.

A subsequent 45-second profile of 1264 forward spent 15.017 seconds rediscovering
independent symmetry factors. Ordered branch programs had discarded the binary
graphs already attached to their fragment responses. The follow-up carries
those graphs into their separate action frames, allowing the existing core to
reject impossible factor splits directly. Each allowed generator preserves this
graph; event terms still join dependent choices. This derived graph is prepared
only when needed and is not duplicated in persisted records. The correlated-arm
test checks structural invariance, dependence without stabilizer construction,
and a valid two-event query. The final unchanged source passed all **1,099 tests
in 157.84 seconds**. Point validation is still required.

Repeating the same pilot with fragment graph reuse alone again completed 23/36
directions; thirteen reached the 120-second diagnostic limit. A bounded trace
then isolated a 104-position global placement whose factorization still lacked
structural information. The next change also supplies a placement's internal
weight matrix when its generators preserve that matrix. Otherwise no structural
certificate is inferred; the original group and event solver remain intact.
Factorization metadata is now separate from point-selection metadata, preserving
the established search ordering. Both local and placement correlated-arm tests
check the reused certificate, with non-preserving choices retained unchanged.
The final fixed-source suite passed **1,100 tests in 162.29 seconds**; the 126
focused symmetry, query and decoder tests also passed.

The placement-certificate pilot still completed 23/36 directions, with thirteen
120-second timeouts and no memory failures. Its maximum completed direction was
862 forward at 90.652 seconds. Completed-scope CPU time was 578.234 seconds,
versus 568.161 seconds for the first pilot; this was not an overall speedup.

A low-overhead 120-second stack sample of 1264 forward then identified the
remaining repeated work. After local preparation, the branch window attempted
thousands of symmetry search states for programs containing only two sweep
costs. Each unresolved branch could consume another five seconds. The following
repair reuses the native joint endpoint projection and retains one ordered
permutation proof for the best reachable pair. Internal costs that vary remain
in the existing correlated event calculation; only a purely endpoint-dependent
cost uses this direct projection. Fixed changes of starting discovery occur once,
and successive actions retain their order. No complete permutation group or
mapping table is expanded. Independent small-graph enumeration validates counts
and witnesses for both-sided symmetry, changed discoveries, occupancy and ordered
placements; the ring test now forbids a full symmetry search for its sweep.

Revision `719225b` passed 177 focused Linux checks. Its identical pilot completed
35/36 directions, with only 607 forward reaching 120 seconds; no memory failures
or retries occurred. The largest successful direction was 1606 reverse at
105.259 seconds. Complete-scope pilot CPU was 834.529 seconds over the 35
successes, with publication measured separately. These totals have a different
successful denominator from the earlier 23-direction pilots.

A 60-second sample of 607 forward found two further repetitions: an already
prepared variable fragment minimum was solved again after weight-preserving
placements, and an endpoint-only component inside a mixed-cost program did not
use pair projection. The follow-up reuses that fragment's minimum, witness and
proof status directly. Pair projection also accepts an independent search
component with its current fixed choices; it requires a generator certificate
that the remaining group factors across ordered action frames. Other components
keep their choices, and cross-action constraints remain joint. Small complete
coset oracles and a variable weighted-ring example test these conditions.
