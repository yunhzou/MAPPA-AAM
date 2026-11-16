# Why the eleven ten-seed Golden references are absent

This report audits the final saved search families at matching tolerance 1, ten seed orderings per cut, cap 100, both directions, and the uncut plus independent single-edge sweep. Recovery remains **1,840/1,851**, with eleven definitive nonrecoveries. No diagnostic witness is added to that benchmark score.

## What the saved families establish

| Cases (zero-based) | Observation | Consequence |
|---|---|---|
| 7, 1358, 1377, 1475 | Every terminal has exactly one more heavy-atom pair than the reference: 30 vs 29, 37 vs 36, 23 vs 22, and 19 vs 18, respectively. | A target permutation cannot delete a mapped pair. The retained maximal matching families cannot reproduce the reference's unmatched-atom decisions. |
| 19, 590, 871, 986, 1228, 1553 | No terminal has a compatible reference atom-pair orbit signature, even when all permitted target actions within a cut are pooled into an overapproximation. | These references cannot be hidden among unexpanded permutations of the retained families. Different atom choices or fragment placements are required. |
| 1285 | Some representatives pass the necessary independent-orbit test, but the full joint correspondence is absent. | Independent orbit compatibility does not establish a compatible correlated mapping. Carbon/oxygen attachments matter jointly. |

All scheduled raw cuts were present and inspected in both directions. Cardinality and source-choice invariants apply to entire families. The target-orbit test deliberately pools actions across different paths and ignores constraints, enlarging the permitted action set: rejection is conclusive, acceptance is only a necessary condition. Exact family verification, already completed for the benchmark, supplies the final negative verdict where that relaxation passes. There was no atom-bijection or permutation-group enumeration.

All eleven have cap stops in at least one direction. That observation cannot establish that the cap removed the reference. Case 986's ten-order reverse search is uncapped. In eight additional reference-blind checks of cases 986 and 1285, increasing the cap to 1,000 or reducing matching tolerance to 0.5 did not recover either reference; five directional checks were uncapped. Those bounded sensitivity checks are recorded separately in `probes/` and do not establish failure under every possible seed schedule.

## Causal traces of the final growth implementation

The new traces reproduce the previously documented mechanisms using the frozen final benchmark engine, not merely the historical summaries. Ordinary Python and native graphs agree exactly, and both ordinary trajectories are uncapped at cap 2,000.

**Case 986: unchanged pyridine ring versus reference label exchange.** The original RDF swaps adjacent pyridine labels 15 and 16. The ordinary mapping preserves the ring, with three changed heavy-atom pairs; the reference has seven, including four caused by that swap. The traced accepted aromatic edge has weight 1.5 but its reference image has weight zero, so the reference cannot survive in that fragment at tolerance 1. A diagnostic reversal of the two labels is equivalent to the ordinary result under the existing chemical certificate. The actual benchmark reference is unchanged. Reference-constrained node compatibility allows a different fragment decomposition and represents all sixteen reference pairs.

**Case 1285: correlated oxygen/carbon attachments.** The original oxalic-acid reference attaches each retained carbonyl oxygen to the other retained carbon. Some saved candidates choose both reference oxygen atoms, so this is not simply choosing a different oxygen donor. Ordinary growth accepts a C=O connection whose reference image is a nonbond (weights 2 versus 0). Ordinary witnesses have five or six changed heavy-atom pairs; the reference has nine. Reference-constrained matching represents all fourteen pairs. A separate ordinary ten-order run with two C–O attachments and one O–H bond omitted also recovers the reference. These three edges were selected from the reference for diagnosis; the successful cut set is absent from the independent single-edge schedule. This is a feasibility demonstration, not blind recovery, a minimal-cut proof, or a recommendation to enumerate arbitrary triple cuts.

The constrained probes only restrict local node compatibility; they do not preload a completed reference mapping. Their success establishes representability with different choices, not that the ordinary blind search should have found the same path. They locate failures in selected greedy trajectories and do not prove that every possible seed ordering must fail. For the other five correspondence misses (19, 590, 871, 1228, 1553), the audit proves absence from the retained families but does not isolate a unique causal contribution of cap truncation versus unsampled growth choices.

## Reference and comparator checks

Independent reconstruction from the original RDF agrees with the prepared reference in all eleven cases, including full endpoint stereochemical certificates. There is no detected atom-index conversion error. `original-reference-audit.json` retains the mapped reaction strings, RDF hash, and reference-induced changes. The ring-label exchange and oxygen-attachment exchange are present in the supplied RDF; they were not introduced by the canonical input ordering.

Only case 871 is recovered by any comparator configuration examined: RXNMapper and SLAP (including the sweep). None of the other ten is recovered by any of the seven released configurations or by the bidirectional SLAP sweep under the same strict relation criterion. Shared misses warrant inspection of annotations and search objectives; they are not proof that the reference labels are wrong. In particular, heavy-atom bond-change counts alone are not a chemical truth criterion. All eleven remain counted as failures to recover the supplied reference.

## Implications

The four cardinality mismatches require explicit decisions about unmatched atoms if exact partial relations are to be represented. The traced examples require retaining a branch that declines an otherwise successful extension, or selectively reopening coupled boundaries. More symmetry enumeration within an existing fragment cannot turn a bond-breaking reassignment into a fragment-preserving automorphism. More seeds, higher caps, and single-edge cuts increase coverage but do not make the fragment partition search exhaustive. Competition was not enabled in the Golden configuration, so this audit does not measure whether it would recover these eleven.

## Scope and execution

The complete saved-family audit used one worker per task, with a five-minute command watchdog. The successful audits took 178.99 CPU seconds in total and peaked at 221.84 MiB. Eight small sensitivity probes used at most four concurrent workers. Two short causal traces were completed locally with a 60-second subprocess watchdog each. No full benchmark rerun was performed.

Initial diagnostic output-writing errors are recorded in the Slurm accounting and logs; these were analysis-script serialization errors, not AAM failures. Original traces and search archives remain under the dedicated HPC validation project. Compact audit outputs, probe outcomes, current causal traces, independent reference checks, scripts, source identities, and accounting are included here. The summary is validated against all 22 directional audits and the unchanged final list of eleven misses before updating the paper.
