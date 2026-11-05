# Full Golden decoding with binary fragment symmetry

This study tests the repaired integrated decoder from `4d2e0d7` on the saved full Golden
matching catalogue: 1,851 reactions, 3,702 directions, discovery tolerance 2,
one seed, branch cap 200 and the full single-edge sweep. Matching is not repeated.

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
