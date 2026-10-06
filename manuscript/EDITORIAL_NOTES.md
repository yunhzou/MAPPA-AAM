# Author information remaining before submission

The scientific text has been rewritten around the final algorithm and supported results. The following information still requires author input; none has been invented:

1. Corresponding-author details, individual contributions, and competing-interest declarations. The existing author order and supplied affiliations have been preserved. Placeholder declarations are omitted from the compiled paper.
2. Provenance and selection of the 140 coordinate reaction pairs, the electronic-structure method/software used to obtain their cached WBO matrices, and redistribution permissions. The paper calls this a development collection without annotated mappings, not an independent accuracy benchmark.
3. Permanent data/archive identifiers if required for submission. The code URL points to the existing repository on `main`.

Golden and coordinate configurations are reported separately. Golden tests the fragment-search stage; the published coordinate protocol also tests full event-window decoding. Competition is excluded from both published protocols. Resource-limited outcomes are explicitly distinguished from verified misses; excluded or unsearched continuations are not credited.

The paper's direct signed WBO-change score uses 0.5 for ordinary pairs and 0.3 for metal-involving pairs, including H. Superseded scoring definitions, abandoned reconstruction variants, and experiment diaries are excluded from the manuscript and compact bundle.


## Molecular figure conventions

Reserve molecule colors for fragment membership. Keep original black single/double bond strokes. Mark bond cleavage/weakening with a vivid red × directly on the affected bond, and formation/strengthening with paired bright green inward arrows. Keep the marks compact so fragment colors and black bond lines remain visible. Match fragment colors across mapped endpoints; different candidate partitions may assign colors differently and must say so. The selected discussion example is Golden case 9, with the verified heavy-atom reference labeled. Ground truth refers to the dataset annotation, not a validated mechanism. Count events from the evaluation matrices, including H, while distinguishing heavy-only depictions and unannotated H identities.

TS motivation concerns querying jointly allowed reaction-core atom assignments and building mapping-dependent endpoint interpolations. Do not illustrate this as arbitrary rigid translation/rotation of a fragment. Keep these as downstream justifications, not validated TS results, and distinguish illustrative symmetry queries from a demonstrated core shuffle.

Describe the default cut sweep as source-edge masking that changes growth routes and relaxes preservation constraints. It counteracts greedy fragment commitment without forcing a final broken bond or a permanent boundary. One seed still means one ordering for every uncut/single-edge subrun; final events always use original matrices.

A single search-edge mask can redirect later growth and produce multiple bond events. Never equate the number of sweep cuts with the number of final events, or imply that multi-event alternatives require enumerating event combinations. Search completeness is a separate issue.

Default MAPPA uses one seed ordering and cut sweep, followed by separate decoding. The experimental competition module remains optional and off by default; exclude it from the paper, pseudocode, and film. Use baseline-only coordinate evidence from `reports/published_baseline_20260915/`. Golden already measures this pipeline.

Keep Algorithm1 at the level of branches, fragments, cut conditions, conditioned growth, and separate decoding. Cross-reference Figure1 panels; keep isomorphism and automorphism internals in their own Methods subsections. Search caps preserve truncation status.

Credit the cut-sweep strategy as introduced in this work and integral to default MAPPA. SLAP sweep is our extension to the prior algorithm. Distinguish the original SLAP literature score (86.9% on 1,758 reactions), our strict default binary score (85.95% on 1,851), the paired uncut binary/weighted bidirectional union (89.74%), and that union with our sweep (97.03%). Attribute the sweep benefit to the paired comparison, not the cross-protocol literature difference.


## Hardware disclosure for publication

Report only the CPU model when describing compute hardware, including future
benchmark updates. Do not publish site names, machine identifiers, login or
network addresses, account names, absolute storage paths, scheduler commands,
queue or allocation details, system memory, or job-concurrency settings.
Keep algorithm parameters, measured timings, completion status, and the
distinction between timing cohorts. Author affiliations are separate from
compute-hardware disclosure. Check the PDF and every submission-bundle member
before publication; raw execution manifests are not submission artifacts.

## Reproducibility in the paper

Keep the paper focused on the algorithm, evaluation definitions, essential settings, results, and limitations. Put exact versions, hashes, configuration fields, commands, run bookkeeping, and detailed audit records in the code repository, linked through `manuscript/REPRODUCIBILITY.md`. Do not reproduce an execution manifest in prose. Retain settings that change interpretation (tolerances, seeds, caps, sweep/direction policy), timing scope, denominators, and completeness limits. Apply the CPU-only disclosure policy to the guide and submission bundle as well.

## Branch-cap labeling

Keep cap 100 labeled as the default. Label cap 2,000 as the larger-budget comparison. Do not relabel the earlier coordinate results as cap 100. Count an empty capped search separately from incomplete decoding; a completed empty catalogue does not provide a mapping or a minimum. Report the paired 140-case evidence and link its code report.
