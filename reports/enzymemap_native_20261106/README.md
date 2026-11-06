# Fresh EnzymeMap benchmark of the accepted milestone

**Status at 2026-11-06 22:45 UTC: bulk matching resumed at 2 GiB.**
The user authorized larger-memory retries for memory-limit failures. The first
retry tier will use 4 GiB, preserving completed results, original failure audits
and consumed compute budgets. Those retries have not been submitted yet.
All 33 unstarted matching batches were released, and continuations were submitted
for exactly the three drained batches (8, 28 and 29). The first resumed pool
registered 28 workers; its startup snapshot preceded task admission, so it does
not establish concurrent compute utilization. Other pools remain queued.
Decoding stays held until all matching and reduction results succeed. No
scientific source changed or retry budget was reset. The stage counts below
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
