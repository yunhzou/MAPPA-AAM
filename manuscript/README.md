# MAPPA preprint

The manuscript presents the final algorithm: weighted continuous fragment growth, correlated automorphism actions, conditional saturated-placement branching and cut sweep, unordered final fragment pairs, and symbolic signed-event decoding.

- `manuscript.pdf`: compiled paper.
- `preprint.tex`: LaTeX entry point, retaining the supplied author list and preprint style.
- `includes/paper.tex`: methods, evaluation, and discussion.
- `includes/supplement.tex`: evaluation scope, completeness argument, and supplementary results.
- [REPRODUCIBILITY.md](REPRODUCIBILITY.md): full technical settings and links to code, input/source pins, and per-case evidence.
- `figs/fig1_algorithm.{pdf,svg,png}`: the branching algorithm overview: RDKit vector depictions of alcohol sites and methyl acetate saponification illustrate hydrogen-preserving local matches, three conditional assignment levels, automorphism compression, cut sweep, and event decoding.
- `evidence/competitors.json`: audited Golden results for RXNMapper, LocalMapper, Chython, Indigo, RDT, and default SLAP configurations; Table 2 combines their returned-output reference recovery and mean CPU call times with all evaluated MAPPA seed/sweep configurations, without assigning a ranking to SLAP alternatives.
- `evidence/molecule_example.json`: checked SMILES, mappings, and heavy-atom event lists for the constructed illustration.
- `figs/fig2_golden.{pdf,svg,png}`: Golden reference-family recovery and direction-specific mean CPU cost.
- `figs/fig3_coordinate.{pdf,svg,png}`: MAPPA and prior-method minimum-event comparisons, minimum-count distributions, and measured stage costs on 140 cases.
- `figs/fig4_alternatives.{pdf,svg,png}`: Golden case 9, with the recovered reference and two alternatives. Red × marks and green inward arrows identify bond changes; pastel colors identify actual saved fragments. The Methods demonstration motivates multiple event candidates; the Discussion describes potential downstream uses.
- `evidence/golden_case9.json` and `golden_case9_verification.json`: full witnesses, source fragment groups, source hashes, original-RDF reference verification, and H-event lower bounds for the selected heavy mappings.
- `evidence/multicandidate_example.json` and `evidence/output_multiplicity.json`: current-API example and archived comparator-output audit for the discussion.
- `manuscript_bundle.zip`: compact paper, source, figure, and evidence bundle.

The manuscript does not include abandoned variants, debugging history, or superseded scoring comparisons. Older research artifacts remain in Git for provenance but are not included in the paper or its compact bundle.

## Results and scope

`evidence/paper_sources.json` records the source commit, snapshot hashes and full report location. The figure and table builder requires complete fresh-campaign evidence; pending outcomes cannot silently enter a final table.

## Build

Install Python 3.12, Tectonic, and the packages in `requirements-build.txt` in an isolated environment. Then:

```sh
MANUSCRIPT_PYTHON=/path/to/python bash scripts/build.sh
```

The build regenerates figures and tables from checked-in evidence, compiles the paper, verifies numerical claims and references, and creates the compact bundle. It performs no mapping searches. The first Tectonic build may download TeX packages. A saved `preprint.bbl` accompanies the source bundle.

For visual review, render the PDF with Poppler:

```sh
pdftoppm -scale-to 1400 -png manuscript.pdf build/page
```

`EDITORIAL_NOTES.md` records the remaining author-supplied metadata and dataset provenance needed for submission. These are not fabricated in the paper. This is an arXiv-oriented manuscript, not a submitted arXiv record.

Paired no-sweep controls are included in the main Golden results table. Evidence: `reports/golden_unswept_20260913/`; MAPPA 1,489/1,851, SLAP 1,661/1,851 (one unresolved). No new mapping searches were run. The molecular figure now checks alcohol-site matches including attached H and shows three conditional assignment levels.

Runtime detail is in Table 3 (paired mean/median/95th-percentile CPU by endpoint-size orientation, paired with recovery), Figure 3, and the unified comparison in Table 2. `timing_comparison.json` retains per-case measurements and matched-cohort identities. Archived default-call CPU times are marked with a superscript in Table 2 because their CPU model is not recorded; no cross-CPU speed ranking is implied.

The 140-reaction coordinate collection compares bond events under unverified mappings; it is not a mapping-accuracy benchmark. The cap control returns complete mappings for 138 reactions at cap 100 and all 140 at cap 2,000; cases 123 and 125 are the two search-completion losses. Table 4 reports the baseline-only bond-event comparison at cap 2,000, Table 5 reports archived stage costs, and Table A4 reports fresh full-pipeline timing. Against SLAP sweep and native SLAP, lower/equal/higher event counts are 4/136/0 and 15/125/0, respectively. Final event-pattern evidence is in `reports/published_baseline_20260915/`; older competition-enabled decoding reports are excluded.

Figure 2 uses Golden case 9: the recovered reference and two alternative mappings. Compact red × marks for losses and paired green inward arrows for gains sit on the affected bonds in Figures 1 and 2, retaining fragment colors and the original black bond lines. Golden labels concern heavy atoms; totals explicitly include unannotated H assignments. The three selected witnesses are not an exhaustive decoded set. Evidence: `golden_case9.json` and `golden_case9_verification.json`.

## 3D research preview

[Grow. Branch. Decode.](animations/mappa_research_preview/index.html) is a 22-second baseline-only 3D film of Golden case 15. Conditional fragment placements grow into a tree and the two shown branches decode to two distinct event classes. The complete baseline catalogue has 44 unordered branches, 96 families and nine classes through six events; two classes are shown. [MP4](animations/mappa_research_preview/mappa-grow-branch-decode.mp4) · [GIF](animations/mappa_research_preview/mappa-grow-branch-decode.gif).

The completed same-CPU Golden cap ablation is in the appendix (`evidence/golden_cap_ablation.json`). It includes cap 100 and cap 2,000 accuracy, unknown outcomes, and search CPU for all four swept seed settings, both uncut controls, and SLAP baselines. All paired timings use the same 1,403 reactions; all-attempt recorded cost is shown separately. Full results are in `reports/golden_controlled_20260915/final/`. The 3D film uses the verified baseline Golden example in `reports/published_baseline_20260915/film/`.

## Published pipeline

The published method uses one-seed cut sweep followed by separate decoding. Competition is excluded. Golden results are unchanged. Baseline-only coordinate evidence and decoding checks are in [the publication report](../reports/published_baseline_20260915/). At cap 2,000, MAPPA recovers 162/168 SLAP-sweep patterns and 151/160 native-SLAP patterns within the fixed windows. At equal minima, the corresponding counts are 158/164 and 138/140. The output has 300 event patterns, including 166 at MAPPA's own minima.

Final end-to-end timing is recorded in `evidence/end_to_end_timing.json` and Table A4: all 140 reactions complete within the five-minute watchdog, preserving all 300 archived event classes and all 120,052 saved families. Mean/median CPU seconds are 1.257/0.300 for search, 4.647/0.288 for decoding, and 6.073/0.812 end to end. These are per-reaction CPU times, not divided by the four workers. Golden timings retain their explicitly labeled search-only scope. The final paper results, evidence, and bundle are ready; no further run is pending for these reported claims.

## Authors

Yunheng Zou; Olalla Nieto Faza; Shifa Hussain; **Varinia Bernales**<sup>†</sup>;
**Alán Aspuru-Guzik**<sup>†</sup>. † Principal investigators. Shifa Hussain is
assigned to the Department of Biology, University of Toronto Mississauga.
[Affiliation provenance](evidence/author-affiliations.md) records the supplied
assignments and official institutional sources. `main` is the publication branch.

## Unified Golden comparison

Table 2 presents all Golden recovery results and mean CPU times in one continuous
table, without separate default-mapper or timing panels. MAPPA 1/2/3/10 seeds,
MAPPA no sweep, and bidirectional SLAP use the paired 1,403-reaction AMD EPYC 9J14
timing cohort. Superscript `a` marks archived default-call times whose CPU model
was not recorded; these cannot support a cross-hardware speed ranking. The
default SLAP row uses binary, which had higher recovery than the weighted default.
The bidirectional SLAP baseline still combines both directions and both modes,
as stated in the caption. Single-output methods are labeled “One bijection.”
Detailed archived call counts and wall-time distributions remain in
`evidence/timing_comparison.json`; the duplicate appendix timing table is removed.
