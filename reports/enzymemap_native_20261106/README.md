# Fresh EnzymeMap benchmark of the accepted milestone

**Decoding is complete: 95,948 / 95,948 directions have successful saved outputs.**
The independent scope audit confirms all 47,974 reactions, both directions,
42 disjoint batches, unchanged input hashes and original ledgers. The final
three recoveries finished in **408.51 seconds wall time** using three concurrent
workers. Their measured peaks were **7.42–7.59 GiB**, below eight GiB even though
this final recovery reserved a 12-GiB guard. Original attempts, consumed budgets,
scientific settings and matching outputs remain intact. These source changes
share storage; they do not alter the mapping relation or event model.

**External catalogue verification is now conclusive for all 95,948 directions.**
Reaction-level reference recovery is **47,798 / 47,974 (99.6331%)**, accepting
recovery in either direction. The other **176 reactions** are not recovered;
there are **zero unresolved reactions, unknown directions or outstanding stage
failures**. Direction-level counts are 94,823 recovered and 1,125 not recovered.
These totals combine preserved conclusive outputs, registered failure recoveries
and separately saved certificates resolving the earlier unknown outputs.
The independent decoded-reference audit has **not run**: these are catalogue
membership results, not proof of membership in the saved decoded event window.

### Shared-fragment reference verification

Revision `5b29531` checks fragment membership once for each distinct local query
and reuses those answers across branch compositions. Impossible compositions
are rejected before solving. Independent choices compose directly; overlapping
placements retain their original joint order. Unannotated atoms still require
feasible full assignments. Every accepted composition is replayed through its
saved family, retaining a concrete injective mapping and its choices. Mapping,
reduction, decoding, reference policy and permitted relationships are unchanged.

All **17 remaining verification timeouts recovered** in **54.13 seconds of pool
wall time**, with 18 registered workers and a four-GiB guard. The longest task
took **48.77 seconds including catalogue loading**; maximum sampled peak memory
was **1.41 GiB**. This recovery explicitly granted another 120 seconds per selected
task and preserved all previous attempts and consumed budgets. The independent
scope audit checked all 17 outputs, empty queues and unchanged original ledgers.
An initial launch stopped during preflight, before worker startup, because the
preserved history required its earlier execution-error authorization; the
corrected launch retained that recorded authorization.

Two bounded point measurements separate input loading from verification:

| Case | Original result | Catalogue loading | Verification and certificate replay |
|---|---|---:|---:|
| 10393 reverse | 600-second timeout | 2.155 s | 1.515 s |
| 10383 reverse | 600-second timeout | 17.315 s | 0.00685 s |

All **26 previously unknown references** also recovered in separate queries of
their unchanged catalogues. Twenty-one resolved with one second per shared
component query. The final five used ten seconds per component, retaining
a 120-second case watchdog; their complete query times were 11.60–16.85 seconds.
The original unknown records remain intact alongside the new certificates.

Validation found no differences in **2,400 small explicit-permutation oracle
comparisons** or **16 real conclusive controls** (eight recovered and eight
not recovered). Local reference/Golden tests (77), campaign tests (47), smoke
checks (26) and Linux reference tests (34) passed. This is scoped equivalence
validation, not a fresh timing or rerun of all previously conclusive queries.
See [verification evidence](fragment-reference-validation.json) for individual
measurements, certificate hashes, scope and final coverage. Operational details
and full certificates remain private.

### Measured stage times

Saved audits cover every reaction and all three completed stages. Times below
sum successful worker computation across both directions and the full cut sweep,
with one seed, branch cap 100, discovery tolerance two and decoding slack five.

| Stage | Mean CPU seconds per reaction | Median CPU seconds per reaction |
|---|---:|---:|
| Matching | 3.02 | 1.25 |
| Reduction | 1.39 | 0.89 |
| Matching plus reduction | **4.41** | **2.15** |
| Decoding | **8.75** | **3.46** |
| Matching, reduction and decoding | **13.16** | **5.80** |

Including output publication and measured worker overhead, the mean is
**14.60 CPU seconds per reaction** (194.51 CPU hours across the dataset).
Historical failed attempts add approximately 27.93 sampled CPU hours, bringing
the observed cost including those attempts to about 16.69 CPU seconds per reaction.
Verification, diagnostic runs, scheduler waits and coordinator overhead are
excluded. These are summed CPU costs, not parallel campaign wall times.
The combined computation distribution has a 95th percentile of 41.86 seconds
and a 99th percentile of 128.58 seconds per reaction.

This measurement combines original successful tasks with successful recoveries;
it is not a fresh full-dataset timing of one final performance revision. See
[stage timing evidence](stage-timings.json) for separate computation, publication,
direction/reaction distributions and failed-attempt accounting. Original audit
references and source-ledger checksums remain in private evidence.

### Decoder storage correction and recovery

The 4-GiB recovery completed 471 of 545 failed decodes in 628.94 seconds,
with 74 memory-limit failures remaining (69 during computation and five during
publication), and no new timeouts or crashes. This brings successful decoding
to **95,874 / 95,948 directions**. The next recovery completed 71 of the 74,
bringing decoding to **95,945 / 95,948 directions**. The remaining three are
37269 forward/reverse and 37270 forward, which exceeded eight GiB. Verification
recovery waits for decoding to finish. Original results and resource failures
remain recorded.

A bounded profile reproduced original decoder peaks of 4.08 GiB for 6723 forward
and 4.13 GiB for 33459 forward; 23481 forward exceeded an 8-GiB diagnostic guard.
The latter retained 2.93 GiB in prepared response arrays and 817 MiB in expanded
pair-weight summaries near termination. The storage correction shares identical
immutable weights, uses shared uniform thresholds, and stores exact weight
multiplicities. Sweep endpoint summaries retain only incident counts and derive
nonincident counts from the fragment summary. No scientific settings or mapping
permissions change. The final storage revision also releases unused Python action-domain copies and
shares identical native domain rows without closing overlapping relationships.
All 1,345 local tests and 196 affected Linux tests passed. The complete saved
decoded record for 16297 forward is unchanged on the final native build. Final
real-case validation passed, including complete saved-record equality. A smaller
allocation started with all 24 registered workers simultaneously decoding the
remaining 74 directions; the larger unstarted request was cancelled first.
That exception pool used an 8-GiB guard because the measured largest corrected
case still exceeded 4 GiB. Matching and successful final-stage results are retained.
The three remaining catalogues load at approximately 4.5–5 GiB. In 37269 reverse,
6,889 retained weight arrays occupy 2,891,764,048 bytes; exact value sharing reduces
this to 1,504 arrays and 2,190,036,504 bytes. The additional storage correction
shares these equal matrices without replacing fragment objects, merging distinct
atom sets or anchors, or changing comparison keys. Existing keys also release
the duplicate matrix allocation. All 1,346 local tests and 297 affected Linux
tests passed on revision `058ec735`. The complete saved response for 16297 forward
still compares equal to its earlier successful output. The additional sharing
reduced 37269 reverse from 8.24 to 7.60 GiB and 37270 forward from 8.52 to 7.46 GiB
in instrumented full-catalogue checks. A 37269 forward diagnostic crashed before
this correction; its debugger repeat with Python memory checks completed normally.
That intermittent crash's cause remains unconfirmed. Its instrumented memory
measurement includes debugger and debug-allocator overhead and is not comparable
to an ordinary worker's peak. The final three-direction recovery uses a 12-GiB
guard and explicit additional 600-second allowances, preserving previous costs.

| Diagnostic direction | Before peak | After peak | Outcome after matrix correction |
|---|---:|---:|---|
| 6723 forward | 4.08 GiB | 2.94 GiB | Complete |
| 33459 forward | 4.13 GiB | 2.35 GiB | Two traced repeats complete; one earlier SIGSEGV remains unexplained |
| 23481 forward | Exceeded 8 GiB | 5.84 GiB | Complete |

These are instrumented saved-catalogue measurements, including loading and
publication, rather than new end-to-end benchmark timing. See the sanitized
[memory validation evidence](decoder-memory-validation.json). The crash is not
presented as fixed; diagnostic capture remains enabled during recovery.


**Earlier recovery launch:** The user authorized continuing diagnosis and
resource retries until decoding and verification close, while preserving the
core algorithm. The validated storage correction (`2e54202`) keeps shared numeric
matrices in binary decoded output; explicit JSON export is unchanged. An audit
compatibility correction (`ac990ab`) hashes the two storage forms identically.

The recovery runner now covers either final stage, preserves every prior attempt
and consumed budget, and records explicit added time and memory. It validates
original input hashes, scientific settings, catalogue/reference identity, drained
queues and exclusive ownership. Successive retries use the latest registered
failure; completed work is never repeated.

Validation: 260 affected storage/query/pool tests passed. The full suite passed
1,328 tests and exposed one audit serialization incompatibility; after correction,
all 11 targeted audit tests passed. Another 27 recovery tests and a successive
resource-grant test passed, including exhausted budgets, changed memory,
reference validation, source immutability and duplicate-owner rejection.

The first recovery tier selects exactly **545 failed decodes** and **21 failed
catalogue verifications** from the terminal finishing campaign. It uses 4 GiB
workers and an explicitly recorded additional 600-second compute allowance per
selected stage, retaining all original costs. Decoding registered 96 workers;
the first running-task snapshot showed all 96 active. Verification is scheduled
after that pool releases ownership. Matching, reduction and the minimum-plus-five
scientific settings are unchanged. These are active attempts, not recovered
results or a new accuracy claim. Operational evidence remains private.

**Status at 2026-11-07 02:10 UTC: the bulk pass has ended, but the benchmark
is not complete. All 95,948 directions matched and reduced successfully.
Decoding has 95,403 successful outputs and 545 resource failures. Catalogue
verification has 95,927 completed checks and 21 resource failures. The independent
decoded-reference audit has not run.**

| Stage | Successful tasks | Memory-limit failures | Compute timeouts | Full denominator |
|---|---:|---:|---:|---:|
| Matching | 95,948 | 0 | 0 | 95,948 |
| Reduction | 95,948 | 0 | 0 | 95,948 |
| Decoding | 95,403 | 508 | 37 | 95,948 |
| Catalogue verification | 95,927 | 4 | 17 | 95,948 |

All 42 finishing batches have terminal decoding and verification records and
empty queues, with no active benchmark owners. A completed verification task
can report recovery, absence or unknown; the table does not report accuracy.
Original failure records, results and consumed budgets are intact. At that observation, the user had
requested diagnosis before any extra timeout allowance; none had been granted.
See the [timeout investigation](#timeout-investigation).

**Earlier observation, 2026-11-07 00:43 UTC:** decoding had 1,441 successful
outputs and 73 memory-limit failures; verification was pending.

The 95 reduction memory failures were recovered without repeating matching:
92 completed at 4 GiB in a 290-second pool pass, and the remaining three at
8 GiB in a 46-second pass. Original attempts and consumed compute budgets were
preserved; no extra compute time was granted. The 98 original and subsequent
failed attempt records remain available alongside the successful outcomes.
These are stage completions, not benchmark accuracy results.

After all earlier pools and saved queues drained, revision `2a06a15` deployed
the validated storage correction with the native engine unchanged. Only the
22,320 untouched matching directions were rebalanced. A later resize moved
7,932 of those directions before their allocations started. The matching pass
then completed with zero failed directions. The two measured source identities
remain separately recorded; completed scientific outputs were retained.

The final-stage handoff independently verified the complete, disjoint input
scope and successful matching/reduction for every direction. It balances 42
batches by saved catalogue bytes, copies only the small ledger records and
retains the original immutable result references. All 191,896 imported successful
stage rows and original input rows are unchanged; source ledgers remain intact.
The ordinary 95,853 directions use 2 GiB workers. The 92 and three previously
large directions use 4 and 8 GiB respectively. Every fresh decoding and external
verification stage retains the original 600-second allowance and decoding uses
minimum + 5. No algorithm or scientific setting changed during this adjustment.

The bulk matching resize reached 298 registered workers. The decoding observation
records 267 registered workers across ten active allocations, with 264 running
tasks in their latest saved snapshots and 308 allocated CPUs including overhead.
The 73 decoder failures are unresolved memory limits, not mapping absence. Sixteen
inspected audits exceeded their 2 or 4 GiB guards during computation; those records
do not identify the dominant allocation or establish an algorithm defect. Larger
retries must retain the original attempt history and remaining compute budget.
The global matching/reduction barrier was verified before any decoding job was
submitted.
See [scope, recovery and handoff evidence](throughput-adjustment.json). Operational
resource requests, ownership checks, hashes and scheduler records remain private.

**Earlier observation, 2026-11-06 22:45 UTC: bulk matching resumed at 2 GiB.**
The user authorized larger-memory retries for memory-limit failures. The first
retry tier will use 4 GiB, preserving completed results, original failure audits
and consumed compute budgets. Those retries have not been submitted yet.
All 33 unstarted matching batches were released, and continuations were submitted
for exactly the three drained batches (8, 28 and 29). The first resumed pool
registered 28 workers; its startup snapshot preceded task admission, so it does
not establish concurrent compute utilization. Other pools remain queued.
Decoding stays held until all matching and reduction results succeed. No
active-pool source changed or retry budget was reset. The stage counts below
are the last published counts at this observation.

| Stage | Successful directions | Failed directions | Full denominator |
|---|---:|---:|---:|
| Matching | 7,326 | 0 | 95,948 |
| Reduction | 6,131 | 5 | 95,948 |
| Decoding | 0 | 0 | 95,948 |
| Catalogue verification | 0 | 0 | 95,948 |
| Independent decoded audit | 0 | 0 | 95,948 |

Four failures occurred during reduction computation: cases **18749 and 18750,
both directions**. Each direction had already matched all **528 cuts** on
524-atom endpoints. **27769 forward** completed reduction computation but exceeded
the guard during publication; its 433-atom endpoints had matched all **466 cuts**.
Sampled peaks were **2.013–2.046 GiB**, reached after approximately 12–14 seconds.
These are resource failures, not missing mappings. The audits identify the
stage and phase; they do not yet identify the dominant allocation or establish
an algorithm regression.

The first four pools registered **88 compute workers**, with **80 simultaneous
running tasks** at the recorded startup snapshot. A fifth pool subsequently
started. Initial submitted concurrency targeted 344 compute workers; this was
requested capacity, not observed utilization. Per-array concurrency limits were
then relaxed to use available slots. A proposed resize was abandoned when its
safety check detected an active owner; no cancellation or input rewrite occurred.
After the failures, pending allocations were held and active coordinators were
asked to drain. The global barrier prevented every decoding allocation from
starting. No background monitor was created.

Evidence: [launch and settings](launch.json), [failure metrics](initial-failures.json)
and [verified paused scope](paused-observation.json). Private continuation evidence
records ownership checks, input hashes, pre-resume ledger hashes, empty saved
queues and scheduler actions.

The scope is **47,974 reactions / 95,948 directions**. Revision `c437ed6`
records the accepted milestone after indexed unique-fragment reuse. This run
starts with empty scientific outputs and caches. Only original endpoint and
reference input records are reused, with checksums verified during preparation.
The previous EnzymeMap campaign remains unchanged.

| Setting | Value |
|---|---|
| Discovery tolerance | 2 |
| Seeds / random seed | 1 / 42 |
| Hierarchical branch cap | 100 |
| Search scope | Both directions, uncut plus full single-edge sweep |
| Growth and cut presence thresholds | 0.2 |
| Ordinary / metal event thresholds | 0.5 / 0.3 |
| Decoding window | Best attained model count + 5 |
| Unique-fragment improvement allowance | 5 seconds |
| Worker memory guard | 2 GiB |
| Direction compute allowance | 600 seconds per stage |
| Fragment retention | 64 MiB per prepared direction; shared reaction snapshots |

Matching and reduction must finish successfully across the complete dataset
before any decoding starts. Persistent pools execute the scientific tasks;
individual cuts and seeds do not submit scheduler jobs. Source and native engine
identity are checked before each allocation and at worker registration.

The current decoder retains binary fragment symmetry and original-weight events,
with the approved independent ordinary-boundary model and joint sweep correction.
Bounded optimization is not a claim of an exhaustive physical minimum.

Two separate reference measurements are planned. Catalogue verification queries
the retained matching relation. The independent decoded audit queries the saved
starting orbits and retained event window, replaying each recovered injective
witness and checking its model and physical event counts. Both use exact annotated
heavy-atom indices without endpoint-symmetry relaxation; hydrogen assignments must
admit a valid completion. Unknowns, scientific absence and execution failures
will be reported separately. Reaction recovery means either direction recovers.

Task archives retain matching, reduction, decoding, publication and external
verification CPU/wall times and memory samples. The independent decoded audit
records its own resource use. Parallel elapsed time will be kept separate from
summed compute CPU time.

The Linux native build passed **155 selected cache, execution and campaign tests**
and **26 scientific smoke tests**. Smoke collection excluded two unrelated viewer
and publication modules whose optional dependencies are absent from the compute
environment.
The previously accepted local release passed 1,318 tests; this launch makes no
scientific implementation changes. Operational launch commands, input locations
and scheduler records stay in the ignored private evidence directory.

The older full EnzymeMap run used tolerance 1 and different matching, decoding
and verification implementations. Comparisons will identify those differences;
its catalogue recovery does not establish current decoded recovery, and allocation
time does not measure isolated algorithm speed.

## Reduction memory investigation

The first 14 archived memory failures comprise nine during reduction computation
and five during publication. Matching completed successfully for these directions.
Three saved directions were reproduced in fresh isolated processes using the
unchanged scientific revision, a 4 GiB guard and a 180-second diagnostic limit.
No matching was repeated and no campaign result or retry budget was changed.

| Direction | Unique fragments | Families | Cold process peak, GiB | Reduction + inspection + publication, s |
|---|---:|---:|---:|---:|
| 18749 forward | 3,041 | 4,065 | 2.35 | 14.53 |
| 27769 forward | 3,374 | 8,045 | 1.99 | 13.53 |
| 10383 forward | 7,688 | 29,780 | 2.03 | 20.65 |

All three completed under the 4 GiB diagnostic guard. These instrumented timings
are diagnostic evidence, not a replacement for benchmark stage timings.
The cold reproductions show that retained matching-worker state is not required
to explain the failures in 18749 and 10383. Case 27769 sits close to the original
guard; worker history and sampling can affect which side of it a run reaches.

There is confirmed avoidable storage overhead. In 18749, retained fragment-side
bond matrices occupy 717,461,344 bytes, and their cached comparison keys contain
another 717,461,344 bytes of matrix copies. Equivalent side preparations across
different fragment pairs also repeat 100,751,648 bytes of matrix storage. These
comparison keys are separate from the bounded growth-reuse cache. Reduction
keeps one fragment per atom-set pair and performs no symmetry/bijection search;
these measurements do not indicate a combinatorial matching explosion.
The matrix-copy representation predates the current cache milestone.

10383 uses less matrix storage but its memory peaks during serialization.
Its 29,780 family witnesses contain **11,137,720 separate assignment-pair tuples**
representing only **14,833 distinct pair values**. The pair tuples alone occupy
623,712,320 bytes; their enclosing mapping tuples occupy another 90,292,960 bytes.
These are repeated records such as the same `1→a` assignment in many family
witnesses, not an enumeration of 11 million complete mappings. The serializer's
own allocated memo/buffer storage measures 536,871,112 bytes (approximately
512 MiB). The diagnostic's identity-counting sets temporarily increase RSS, so
the table above uses the first cold run, without those sets, for process peaks.

The appropriate correction is storage sharing and compact comparison keys,
preserving each retained fragment and every branch's mapping values. This does
not require changing growth, symmetry choices, branch caps or decoding rules.
No scientific correctness defect has been demonstrated by this memory
investigation; that is not a new accuracy validation. The bulk implementation
remains frozen, and the three diagnostic successes are not campaign retries.

Sanitized measurements: [memory investigation](memory-investigation.json).

## Validated storage correction

Revision `170dc2a` removes the retained matrix copies from comparison keys and
shares equal immutable assignment pairs while building the catalogue. Keys
retain exact byte equality, hashing and ordering. Fragment identity, retained
witnesses, constraints, ordered actions and provenance remain unchanged; matching,
branch pruning, decoding and verification policies are unchanged.

The correction passed **143 focused tests** and the **full 1,321-test suite**.
Old/new reduction was then compared on the same saved matching inputs for the
three directions above, in fresh serial processes on one Linux node. All three
produced identical scientific-content digests, including complete retained
fragment data, original matrices, constraints, mappings, action order, history,
provenance and branch composition. All six save/reload checks preserved those
digests. This checks retained scientific content, not merely branch counts; it
does not claim a new full-dataset reference-recovery measurement.

| Direction | Previous peak, GiB | Corrected peak, GiB | Reduction CPU, old → new, s | Publication CPU, old → new, s |
|---|---:|---:|---:|---:|
| 18749 forward | 2.359 | 1.440 | 9.331 → 9.315 | 3.502 → 3.118 |
| 27769 forward | 1.993 | 1.037 | 8.614 → 8.834 | 3.241 → 2.630 |
| 10383 forward | 2.031 | 0.660 | 14.711 → 17.637 | 4.342 → 2.014 |

Peak RSS covers reduction, publication, content comparison and reloading, sampled
every 0.1 seconds. Each process had a 4 GiB diagnostic guard and a 240-second
watchdog; all corrected peaks were below 2 GiB. This does not establish that every
direction will fit below 2 GiB in a persistent worker. Reduction CPU increased
from 32.656 to 35.786 seconds across these cases, while publication decreased
from 11.085 to 7.762 seconds. Combined CPU was essentially unchanged:
**43.742 → 43.548 seconds**. These are single paired measurements, not a throughput
or speedup estimate.

The isolated candidate imported exactly the two changed modules from committed
Git objects; all other modules and the native engine came from the frozen
baseline. The recorded filesystem identity therefore identifies the baseline,
while the candidate commit and individual module hashes identify the tested
correction. The active benchmark checkout was not changed beneath running pools,
and no campaign results, errors or budgets were overwritten. The correction was
subsequently deployed at a verified drained execution boundary, as recorded
above; these diagnostic successes are not retries of the campaign's failed
reductions.

Sanitized validation: [storage comparison](storage-comparison.json).

## Timeout investigation

The 37 decoding and 17 catalogue-verification timeouts exhausted their original
600-second compute allowances. Sampled CPU time averaged **599.66 seconds** for
decoding and **598.30 seconds** for verification. Every timeout occurred during
computation. These observations rule out scheduling waits as their explanation;
they do not establish mapping absence.

Nine representative failures were profiled from their original archived task
inputs, using the unchanged scientific revision and original memory guards.
Diagnostics requested interruption after 60 seconds and used an external
90-second watchdog. Eight produced interrupted profiles; one required the hard
watchdog. None completed the scientific stage or replaced a benchmark result.
Profiling overhead means these partial measurements are attribution evidence,
not replacement stage timings or speedup estimates.

| Direction | Measured bottleneck in the partial diagnostic |
|---|---|
| 40793 reverse | One symmetry-restriction call occupied 59.18 s despite the five-second fragment allowance. |
| 24793 forward | Two calls to the same symmetry-restriction helper occupied 57.79 s. |
| 16297 forward | Repeated stack snapshots remained inside native symmetry restriction; the external watchdog stopped the process at 90 s. |
| 15959 reverse | Fragment pair-weight lower bounds took 18.28 s; contextual weight-bag bounds took another 16.12 s in the initial 60 s. |
| 13414, 18744 and 18749 forward | Initial profiles also concentrated on fragment bounds and ordered-response preparation. |
| 10393 reverse, verification | Compiling 19 family queries took 54.70 s, including 41.25 s in finite-constraint propagation. |
| 10383 reverse, verification | Loading and reindexing took 42.35 s; the first query then spent 23.93 s preparing symmetry factors. |

The native restriction operation computes which symmetry operations remain when
specified atoms must stay fixed. It has no internal deadline check. The caller's
cooperative five-second deadline cannot stop it before it returns, although the
execution layer still enforces the overall stage limit. Capturing the first
restriction input for 40793 revealed 83 generators on 269 positions, with all
145 movable positions fixed. Its answer must be the identity group. Rebuilding
a general group chain for that query is avoidable work, not a requirement to
enumerate possible atom mappings. The other two captures stop at their first
restriction call and do not identify every later expensive input.

The Python event-bound paths repeatedly classify and compare atom pairs across
distinct fragments and contextual pieces. Unique-fragment preparation remains
cached by atom-set identity; these profiles do not show that cache being bypassed.
They identify expensive early preparation, not a complete attribution of each
600-second attempt. Catalogue restoration also repeats assignment and action
interning, which is especially visible in the 26,589-family verification case.

For verification, the one-second family-query allowance is checked **after**
constraint compilation. Expensive compilation can overrun it before the solver
is reached. The observed propagation work includes millions of Python domain
lookups and choice hashes. Timeout remains unknown; no constraint should be
weakened to turn it into recovery.

The indicated corrections are compact native symmetry restriction, cooperative
checks inside its long operations, and faster preparation/reuse of the same
event bounds and verification constraints. Their scientific equivalence and
performance still require validation before deployment. No algorithm changes,
extra retry time or failed-stage reruns were made during this investigation.
The separate memory failures still require recovery.

Sanitized evidence: [profiles and terminal stage counts](timeout-investigation.json).

### Isolated symmetry-restriction diagnosis

The captured restriction for **40793 reverse** confirms a performance regression
in the native replacement introduced by `040f817`. The preceding implementation
used SymPy's `pointwise_stabilizer(..., incremental=False)`. On the same captured
input on the local Mac, that historical operation completed in **65.95 ms**;
the current native operation did not finish within a **15-second process
watchdog**. These are isolated operations, not full historical/current decoding
runs.

The native helper fixes atoms one at a time. After each effective restriction,
it reconstructs a complete generator-reduction chain, including unrelated
independent symmetry components. Each insertion into that chain also rebuilds
intermediate representatives. The helper does not recognize the immediate
identity answer when all movable positions are fixed. The main event search
separates independent components, but this helper again processes the full
event-preserving symmetry group.

An isolated calculation using the **existing native independent-factor and
restriction functions** first separated those components, restricted each one,
and collected their generators. It returned the exact identity for the captured
case in **1.35 ms median across seven repeats**, including factorization and
result construction. This experimental composition is not integrated into
production. No SymPy fallback or change in allowed mappings was introduced.

The same calculation passed **483 checks against independent exhaustive small
groups**, including a paired exchange that must move together and a case where
combining generators permits a restriction that no individual generator supplies.
For 64 independent pair swaps with all but the last pair fixed, the current
operation took **2.323 s**, while separate restrictions took **0.695 ms** and
retained exactly the one free swap. This also demonstrates avoidable cost when
the answer is not the identity.

These measurements identify redundant native group construction as a concrete
cause. They do not establish the performance of a production repair, its effect
on all 37 decoding failures, or a full-decoder speedup. The scientific code,
existing benchmark outputs and retry allowances remain unchanged.

### Native restriction repair

The production helper now restricts disjoint generator-support components in
compact local frames. All positions moved by one generator remain joined, so
coupled exchanges cannot become independent swaps. Fixing every movable position
returns the identity immediately. The remaining computation uses the same exact
pointwise restriction; no allowed operations are added or removed.

Long native restriction, factorization and group-chain operations now receive
the existing remaining search allowance and check it internally. Interrupted
optimization retains an actual incumbent with incomplete status; interrupted
queries remain unknown. No partial group is returned as an exact result. The
external process watchdog remains responsible for hard resource limits.

The captured 40793 reverse operation takes **0.0650 ms median over nine fresh
groups** on the Mac, and returns the expected identity. This is a helper timing,
not a full-decoder speedup. Compared with the measured historical SymPy helper
(65.95 ms), it is approximately 1,014 times faster on this particular input; the
unrepaired native helper exceeded its 15-second watchdog. Keeping one free swap
among 64 independent pairs takes **0.0329 ms median over nine repeats**, versus
**2.323 s** before, with exactly the same remaining swap.

Both native extensions were rebuilt. The full local suite passed **1,325 tests
in 126.03 s**, including independent exhaustive small-group comparisons, coupled
and composed restrictions, large independent swaps, deadline interruption,
fragment witnesses, queries and complete stage pipelines. Both native extensions
also rebuilt on Linux, where **247 focused tests passed in 21.86 s**. The repaired
scientific source is `ea4d9e6`.

Three isolated full-catalogue checks used the original two-GiB guards, original
decoding settings including `minimum + 5`, and no profiler. Their original
600-second campaign attempts had timed out. The new checks progressed past the
previous symmetry restriction, but **all three then reached the memory guard**:

| Direction | Wall time to memory stop | Sampled peak GiB | Complete decode |
|---|---:|---:|---|
| 16297 forward | 86.10 s | 2.093 | No |
| 40793 reverse | 73.03 s | 2.000 | No |
| 24793 forward | 85.99 s | 2.106 | No |

These times include validation startup and catalogue loading; they are neither
CPU times nor completed-stage timings. Memory is sampled, so a process can exceed
the guard between observations. Stack samples show fragment event-bound
preparation, local weight-bag comparisons and response/action construction rather
than the previous persistent restriction call. This does not yet establish which
objects dominate memory, a full-decoder speedup, or recovery of all 37 timeouts.

The specific helper regression is repaired and tested. Full-case recovery remains
blocked by memory in these checks, with other preparation costs still visible.
Original campaign errors, outputs and consumed budgets are unchanged. No fresh
retry time was granted, and these validation measurements do not count as
benchmark recovery.

### Four-GiB decoding diagnosis

The user authorized testing four GiB rather than treating two GiB as a scientific
limit. The same repaired revision (`ea4d9e6`), saved catalogues, fragment budget
and `minimum + 5` settings were used in fresh isolated workers. Lightweight
measurements sampled prepared arrays and cache sizes. These runs do not change
the campaign ledger or replenish consumed retry budgets.

| Direction | Outcome under 4 GiB | Peak GiB | Decode and record construction, wall / CPU s |
|---|---|---:|---:|
| 16297 forward | Completed | 2.609 | 91.536 / 91.375 |
| 24793 forward | Completed | 2.549 | 92.045 / 91.935 |
| 40793 reverse | Memory guard | 4.082 sampled | Incomplete |

The two completed responses retain model minimum four and window maximum nine.
16297 has one incompletely optimized fragment; this status remains explicit,
and does not prevent the model minimum from being established by other retained
choices. These checks do not establish external reference recovery. Timings
include catalogue loading and instrumentation, and exclude durable publication.

Loading the three catalogues uses 0.53–0.92 GiB. They contain 1,533–2,145 unique
fragments. The prepared fragment arrays and contextual responses grow with that
retained work. The bounded native domain cache remains small (under two MiB in
the recorded samples); it does not account for the multi-GiB footprint.

There is also confirmed avoidable duplication in output construction. An isolated
`catalogue.to_record()` on 40793 reverse, without decoding, converts 1,899 shared
fragment-side matrices containing 57,936,100 numeric cells into Python lists.
The arrays occupy **463,488,800 bytes (442 MiB)**. The additional lists, pointers
and individual floating-point objects occupy **1,869,959,400 bytes (1.742 GiB)**,
while the arrays remain live. This is repeated storage of bond-weight matrices,
not additional candidate mappings. Total RSS rises from **0.922 to 3.063 GiB**
in this output-only measurement, a **2.141 GiB** increase; the list allocation
accounts for most, but not all, of that increase. Output construction takes
3.917 seconds.

A separate diagnostic tests retaining numeric arrays in binary output while
checking every matrix against the normal list representation. Production modules
and public JSON behavior remain unchanged. Its first attempt exited with signal
11 before reaching output construction. The subsequent debugger run completed
decoding in **90.178 s**, publication in **1.883 s**, and a save/reload check of
every matrix and the saved event window. The worker peak was **2.564 GiB**
(**2.731 GiB** including debugger overhead). Output construction increased RSS
by only 30.1 MiB in that run; source and event-preparation memory remained live.
The response has model minimum four and maximum nine, with one explicitly
incompletely optimized fragment. No external reference audit was performed.

A normal repeat also completed: **90.516 s** for decoding and record
construction, **1.876 s** for publication, and **94.807 s wall / 94.802 s CPU**
including digest and reload checks. Its worker peak was **2.564 GiB**. Both
completed runs produced the same matrix digest, model minimum/window and
fragment/table counts. The binary published result occupies 10,511,020 bytes.

The initial pre-export crash did not recur in either repeat; its cause remains
unidentified and is not treated as fixed. The duplicate storage is established
independently of that process failure. Four GiB is sufficient for these tested
cases with shared numeric-array output, but the production serializer is still
unchanged. Integrate this storage correction through the binary transport owner
while preserving explicit JSON export behavior, rather than raising this case
above four GiB to accommodate duplicate lists. This is a storage correction,
not a change to matching, allowed symmetries or event counting. These checks do
not establish a whole-dataset memory bound or new reference-recovery accuracy.
