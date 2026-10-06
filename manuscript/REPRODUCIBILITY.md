# Reproducing the paper

Use this guide for implementation settings and archived protocol details. The paper reports the algorithm, essential evaluation settings, and results; this repository retains the full technical record.

## Start here

- [Public Python API and simple notebook](../docs/AAM_SIMPLE.ipynb): matching, unique event candidates, and conditional symmetry queries.
- [Source-to-result index](evidence/paper_sources.json): the exact evidence snapshots and report paths used by the paper. Follow these records rather than choosing an older experiment from the reports directory.
- [Golden campaign and coordinate pipeline](../reports/current_validation_20260912/): recorded configurations, component identities, and final-stage results.
- [Golden reference evaluator](../bench/golden_evaluation.py), [checkpoint evaluator](../bench/golden_checkpoint_evaluation.py), and [released-mapper adapters](../bench/golden_competitors.py).
- [No-sweep controls](../reports/golden_unswept_20260913/), [directional recovery](../reports/golden_direction_size_20260914/), [timing analysis](../reports/golden_timing_detail_20260914/), and [published coordinate minima](../reports/published_baseline_20260915/).

These paths are relative to the repository. In the compact manuscript bundle, use the same paths in the linked code repository to access benchmark scripts and reports. The bundle contains the evidence needed to regenerate the paper's tables and figures, not every raw search archive. The coordinate WBO-generation method and redistribution permissions remain unspecified; the current evaluation does not establish regeneration of those inputs from molecular structures.

## Rebuild tables, figures, and the PDF

From the repository root, install the dependencies listed in `manuscript/requirements-build.txt` and Tectonic, then run:

```sh
MANUSCRIPT_PYTHON=python3 bash manuscript/scripts/build.sh
```

This reads the checked-in evidence and runs no mapping searches. It checks numerical claims, figure chemistry, source hashes, references, and publication disclosure rules. For a fresh benchmark, use a new output directory and the source/input pins and configuration recorded for that campaign; current library defaults are not a substitute for those settings. Keep reference verification separate from blind search and preserve incomplete outcomes.

The default branch cap is 100. The original Golden tables use cap 100; the original coordinate tables and Figure 4 use cap 2,000 and are labeled as the larger-budget comparison. The paired coordinate cap study reports both settings separately. The completed same-CPU Golden cap ablation is reported separately in the appendix and in `reports/golden_controlled_20260915/final/`.

## Benchmark source layout

The current entry points and evaluators are indexed in `bench/README.md`.
Historical campaign drivers are preserved byte-for-byte under
`bench/archive/campaigns/`; original paths and SHA256 hashes are recorded in
`docs/publication_cleanup_20260916.json`. Measurements and witnesses remain in
`reports/`. For historical reproduction use the recorded source/input pins;
current library defaults are not a substitute. Superseded visual artifacts
remain available at the Git revision recorded in the cleanup manifest.

## Search and comparison protocols

The active-edge floor is $b=0.2$ and the matching tolerance is $\tau_{\mathrm{iso}}=1.0$. Default MAPPA uses one seed ordering per cut and the sweep. Both published evaluations exclude the experimental competition extension. Additional seeds and the no-sweep ablation are explicitly labeled. Seed streams are deterministic with root seed 42. Sweeps retain the uncut search and individual source-edge deletions. Bidirectional Golden output is a union of returned relations; branches from opposite directions are not spliced. Growth and branch caps bound compressed candidates rather than explicit mapping counts.

Golden uses dataset revision `793475e`. Its heavy-atom evaluator compares correspondence relations modulo endpoint chemical symmetry, retains unmatched atoms, and rejects altered endpoint chemistry. Hydrogen atoms participate in search but their individual identities are not reference labels. All 1,851 records remain in the denominator. All reported MAPPA configurations were rerun against the final algorithm source, commit `3a9ef9a`. Algorithm configurations and component hashes accompany the per-case results. Golden evaluates the fragment-search stage, while the coordinate protocol also evaluates event decoding.

The no-sweep MAPPA control retains only the original uncut checkpoint at one seed in each direction. Its 3,702 saved graphs are finalized and evaluated using the same reference-family criterion, yielding 1,489 recoveries and 362 completed non-recoveries. For SLAP, the control retains ordinal-zero calls from both directions and both bond modes; incomplete sweep records are checked at that ordinal separately. Case 1723 has no committed reverse uncut output and remains unresolved. These are paired controls of the final campaign, not substituted historical configurations. Re-evaluation time is excluded from search-cost claims.

Search and evaluation are bounded by time and resource limits. Interrupted phases can resume completed cuts through bounded continuations. Completed raw cuts can also be finalized and verified individually, without new search. A positive result requires a verified witness; it can establish recovery even when another direction or the combined archive is incomplete. A negative result requires completed verification of the saved output; other cases remain unknown.

The paired timing comparison uses the Apple M4 Pro CPU. Python 3.12.14, RDKit 2026.3.6, pynauty 2.8.8.1, SymPy 1.14.0, and Z3 5.1.0.0 are fixed across the evaluation. Source-hash checks and search/reference parity checks verify consistent execution. Numerical-library calculations use one CPU thread. CPU averages use 1,821 common reactions with both one- and two-order searches completed on that CPU and complete comparator output. Search-call timing includes checkpoint I/O and excludes imports, preparation, reference evaluation, and interrupted work; SLAP's instrumented workflow excludes I/O. Higher-order searches were not all measured on the comparison CPU and are excluded from this timing comparison.

The SLAP comparator uses SLAPMapper 1.0.0, revision `ea248fd`, with both directions and binary and weighted modes. Each sweep variant uses one original ordering and rebuilds initialization for the uncut input and each single-edge deletion. Returned labels are restored to the original endpoint structures before evaluation. Its 1,796 verified Golden recoveries leave 52 completed misses and 3 unresolved records; all remain in the denominator. The sweep perturbs SLAP's input graph and assignment objective, whereas MAPPA relaxes fragment-preservation constraints.

On the coordinate collection, the exact original XYZ-derived SLAP adjacency matrices are restored from the saved graph records and verified against every source-edge deletion. Atom order and atomic numbers are retained; mapping labels are reset. Fresh native and bidirectional single-edge SLAP runs use binary graphs and the original initialization. The XYZ-to-graph preparation itself is not rerun. MAPPA uses the cached WBO matrices as supplied. The available benchmark metadata do not specify their electronic-structure method or charge/spin preparation; these matrices therefore define the computational inputs, while regeneration from molecular structures is not established by this evaluation. Each returned SLAP heavy-atom assignment receives a symbolic hydrogen refinement under the comparator's original bond-change objective, and its witness is rescored on the original WBO matrices using the signed-event definition in Methods. Tied hydrogen optima can give different signed-event classes; the reported comparator set comprises the returned witnesses, not the entire hydrogen label relation. This compares the stated endpoint workflows rather than isolated solvers on identical graph representations.

## Default comparator reproduction

The default comparator predictions were generated with RXNMapper 0.4.3, LocalMapper 0.1.5, Chython 2.18 with chython-rxnmap 2.0, Indigo 1.46.0, SLAPMapper 1.0.0 at revision `ea248fd`, and RDT 4.0.0 at revision `0d663225`. Chython's released ONNX model is identified by its saved hash; it is not represented as a reproduction of a historical GraphormerMapper checkpoint. RXNMapper and LocalMapper use their released default models on CPU. SLAP is evaluated over all returned explicit candidates, without ranking, with its default heavy-atom symmetry handling, with both binary and weighted settings retained in the archive. Table 2 selects the binary default, which has higher recovery; its bidirectional SLAP rows combine both directions and both modes. Indigo uses its native automapping interface; RDT uses the public `RDT.map` API.

The complete archive retains mapping strings, explicit candidate order, version/model hashes, and per-record evaluations. All 12,957 method/reaction records were rescored with the current evaluator and agree with the archived outcomes. Dataset and label-free input hashes match the final MAPPA campaign. Agent-field components returned by a mapper are retained on the input side, preserving their original atom labels. Changed endpoint chemistry and duplicate heavy-atom labels invalidate a prediction. All 1,851 records remain in each denominator.

RDT's final records include execution retries after a worker failure; every final call returns an output. Invalid outputs are not repaired. These default implementations were timed separately from the Apple M4 Pro cohort, so their archived call times are not included in the paired comparison in the main Golden timing figure. Table 2 reports their mean CPU times with an explicit marker for unrecorded hardware; full timing and failure records accompany the reproducibility archive.

## Published coordinate output

The published pipeline is sweep search followed by separate decoding. `reports/published_baseline_20260915/` contains the baseline-only evidence, provenance, and decoding checks. The optional experimental module is not called by this pipeline.

## Event windows and output equivalence

A signed event class contains the nonzero atom-pair signs from the signed-event definition in Methods, modulo source permutations preserving the full event-response coloring relative to the opposite endpoint. This coloring retains every threshold comparison needed to transport a pattern without changing its identity. It is separate from the intrinsic fragment matching symmetry in the fragment-response definition in Methods. The exported witness records both the atom mapping and the formed/strengthened and broken/weakened pair lists.

The 140-case comparison uses fixed per-case windows from the evaluation protocol. Each limit is one event above the larger of the native-SLAP and SLAP-sweep minima under the common signed-event definition. All saved baseline families are exhausted in these windows without a class-count cap; higher-event patterns are outside the completeness claim.

For the saponification example in Figure 1, the direct assignment $(1,2,3,4,5,6)\mapsto(a,b,c,d,e,f)$ loses O4--C5 and gains C5--O6. Exchanging the images of O3 and O4 additionally weakens C2--O3 and strengthens C2--O4: two versus four heavy-atom events. Both assignments preserve the internal matching constraints of fragment ${1,2,3,4}$ at $\tau_{\mathrm{iso}}=1$. Formal bond orders define this illustration; the displayed event counts include only heavy atoms. Exchanging the acetate oxygens also changes the localized resonance representation, so the event signatures do not establish distinct mechanisms. Both product acetate oxygens have zero attached hydrogens. The figure builder verifies these mappings and event lists. For the separate alcohol-site tree, it verifies every attached hydrogen and every bond inside each local matched fragment, the distinct target-site environments, and injectivity along all six displayed paths.

## Audit of Golden nonrecoveries

We inspected every saved cut in both directions for the eleven ten-order misses, preserving the original reference and search outputs. Independent reconstruction from the original RDF agrees with the prepared reference in all eleven cases. Mapping cardinality is invariant within each saved family. For six of the seven remaining cases, necessary source/target orbit-pair counts exclude the reference even after pooling permitted target actions into an overapproximation; case 1285 requires the full joint-correspondence check. These tests concern the retained output, not every possible search trajectory.

All eleven cases have cap stops in at least one direction, which does not establish cap truncation as their cause. Separate ten-order sensitivity checks on cases 986 and 1285 at cap 1,000 or matching tolerance 0.5 remain nonrecoveries; several of those directional runs are uncapped. In the causal traces, ordinary Python and native graphs agree and are uncapped at cap 2,000. Reference-constrained matching represents all 16 and 14 annotated pairs, respectively. The successful case-1285 cut set omits two C--O attachments and one O--H bond; it is reference-selected, and neither its success nor the constrained witnesses contributes to the blind benchmark score. The audit, controls, traces, and source identities are archived with the benchmark report.

## Output multiplicity and the discussion example

A read-only audit of all 1,851 archived records per default comparator found no multiple-record outputs from RXNMapper, LocalMapper, Indigo, Chython, or RDT. Binary SLAP returned multiple records in 522 reactions and weighted SLAP in 417; these are explicit record counts, not counts after expanding label symmetries or grouping events. Figure 2 uses three selected witnesses from a current-source, one-seed MAPPA cut sweep for Golden case 9, with branch cap 2,000 and $\tau_{\mathrm{iso}}=1$. Reference labels were removed before search. The reference witness matches all 15 original heavy-atom pairs, independently checked against the source RDF; the two alternatives change oxygen origins and, in the third witness, carbon correspondence. All 33 explicit atoms are included in event checks at threshold 0.5; H atoms are implicit only in the drawing. The five, six, and seven heavy-atom changes become nine, eight, and nine total events. For each fixed heavy mapping, the sum of absolute differences in attached-H counts bounds the number of H-involving changes from below; each displayed witness attains that bound. These selected witnesses are not an exhaustive decoded set. Matrices, mappings, fragment partitions, and source hashes are in `golden_case9.json`; original-reference and H-bound checks are recorded in `golden_case9_verification.json`. A separate six-atom example in `docs/AAM_SIMPLE.ipynb` illustrates the public API interactively.

### Coordinate search cap control

The archived paired search control at caps 100 and 2,000 returns complete baseline mappings for 138 and 140 reactions, respectively. Cases 123 and 125 have no complete mapping at cap 100. The published event-pattern comparison uses baseline searches at cap 2,000. The old combined search-and-extension decoding tables in `reports/coordinate_cap_consistency_20260915/` are historical and are not used as baseline-only pattern or decoder-time measurements.

### Final Golden cap ablation

All 1,851 reactions were processed under twelve configurations on one AMD EPYC 9J14 CPU thread per search. Caps 100/2,000 each use 1/2/3/10 seed orders with the sweep and one uncut order; SLAP uncut and swept baselines combine binary/weighted modes and both directions. Exact settings and source pins remain in `reports/golden_controlled_20260915/`. The completed results, per-case records and direction tables are in its `final/` directory.

The same 1,403 completed-search cases are used for every paired CPU mean, median and 95th percentile. Recorded search CPU over all 1,851 selected attempts is reported separately, including interrupted calls. Superseded checkpoint-conflict attempts and the initial pilot are retained operationally but are excluded from these timing estimates; these are not total campaign expenditure. Reference verification and process imports are outside the primary search timer. The new timings are not pooled with the main Mac cohort.

At caps 100/2,000, confirmed swept recovery is 1,834/1,834 for one seed, 1,834/1,834 for two, 1,837/1,837 for three, and 1,840/1,835 for ten. Unknown counts at cap 2,000 are 2/4/6/12, respectively; cap 100 has none. Five ten-seed references become unresolved at the larger cap (cases 833, 850, 1568, 1691 and 1786), with no swept gain or loss among cases resolved at both caps. Uncut MAPPA gains 19 confirmed recoveries at cap 2,000 and leaves five cases unresolved. SLAP sweep has two additional unresolved previously recovered cases (1665 and 1806); the new count of 1,794 is kept separate from the main sweep count of 1,796.

## Final end-to-end timing

`reports/final_end_to_end_optimized_20260916/` records the fresh search-through-decoding pass used in the paper. All 140 coordinate inputs complete at one seed and cap 2,000, reproducing all 300 archived event classes and all saved family counts. The script `bench/final_end_to_end.py` uses the public decoder and published event windows, with four single-thread processes, no class-count cap, and five-minute per-reaction watchdogs. `bench/summarize_final_timing.py` reconstructs stage distributions from per-case records.

The earlier run in `reports/final_end_to_end_20260916/` had two decoder watchdogs and is retained as historical evidence. `reports/decoder_optimization_20260916/` compares identical completed cases and verifies all 110,477 previously completed family-support records. Exact timed source changes and hashes are supplied with the new report.
