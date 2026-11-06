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
