# MAPPA preprint

The manuscript describes one-witness fragment growth, direct compressed completion under global target symmetry, conditional fragment placement and cut sweep, family reduction, and separate decoding from a complete inventory of unique fragments. Whole-reaction event-window queries are distinct from ordinary local decoding.

- `manuscript.pdf`: compiled paper.
- `preprint.tex`: LaTeX entry point, retaining the supplied author list and preprint style.
- `includes/paper.tex`: the scientific stages, molecular examples, completed current-method validation, and discussion.
- `includes/supplement.tex`: solver and representation details, evaluation protocols, and historical comparisons.
- [REPRODUCIBILITY.md](REPRODUCIBILITY.md): full technical settings and links to code, input/source pins, and per-case evidence.
- `figs/fig1_algorithm.{pdf,svg,png}`: the integrated molecular overview: RDKit vector depictions connect one-witness growth, conditional completed placements, aggregation across cuts, and shared fragment/boundary event responses. Reference verification belongs to the evaluation and is outside this MAPPA overview.
- `evidence/current_golden.json` and `current_golden_cases.json`: publication-safe summary and per-case outcomes for the completed Golden validation.
- `evidence/competitors.json`: archived Golden results for RXNMapper, LocalMapper, Chython, Indigo, RDT, and default SLAP configurations. The historical comparison retains its original seed/sweep settings and timing cohorts.
- `evidence/molecule_example.json`: checked SMILES, mappings, and heavy-atom event lists for the constructed illustration.
- `figs/fig6_current_golden.{pdf,svg,png}`: completed current-method Golden reference recovery and mean total CPU cost for matching, family reduction, local decoding, and verification. Both directions are included; CPU cost is not elapsed latency. The builder checks the aggregate values against all per-case records.
- `figs/fig2_golden.{pdf,svg,png}`: supplementary archived Golden reference-family recovery and direction-specific mean CPU cost.
- `figs/fig3_coordinate.{pdf,svg,png}`: supplementary archived MAPPA and prior-method minimum-event comparisons, minimum-count distributions, and measured stage costs on 140 cases.
- `figs/fig4_alternatives.{pdf,svg,png}`: Golden case 9, with the recovered reference and two alternatives. Red × marks and green inward arrows identify bond changes; pastel colors identify actual saved fragments. The Results illustration explains why retaining alternative correspondences matters; the Discussion describes potential downstream uses.
- `evidence/golden_case9.json` and `golden_case9_verification.json`: full witnesses, source fragment groups, source hashes, original-RDF reference verification, and H-event lower bounds for the selected heavy mappings.
- `evidence/multicandidate_example.json` and `evidence/output_multiplicity.json`: saved API example and archived comparator-output audit for the discussion.
- `manuscript_bundle.zip`: compact paper, source, figure, and evidence bundle.

The main paper explains the method through an integrated molecular overview, then presents the completed current-method validation. Historical comparator, cap, seed, and coordinate studies remain in the supplement with their original scientific scope. They do not supply a speed comparison for the current implementation.

## Results and scope

The completed Golden validation at source `096a010961f53944ddd4ddc433bd8d56f282a2f2` processed all **1,851 reactions / 3,702 directions** through matching, family reduction, local decoding and separate reference verification. It recovered **1,830 reactions (98.87%)** and **3,618 directions**, with no failed tasks or unknown verification outcomes. Recovery means that the annotation belongs to a retained mapping family, not that a first-ranked mapping is correct. These results apply to the recorded version; they do not certify a full rerun of every subsequent revision.

The September seed, cap, comparator and coordinate studies remain archived evaluations. Their default Golden result of 1,834/1,851 is not substituted for the separate current-design validation. Their timing cohorts and event-window results are not presented as measurements of the current implementation.

The expanded EnzymeMap benchmark is still running. No final EnzymeMap recovery score or completion time is included in the paper. Completing the Golden validation does not establish completion of all new benchmark work.

`evidence/paper_sources.json` records the source versions, snapshot hashes and report locations. Pending outcomes cannot enter a completed-results table. [REPRODUCIBILITY.md](REPRODUCIBILITY.md) explains each study's protocol and timing scope.

## Build

Install Python 3.12, Tectonic, and the packages in `requirements-build.txt` in an isolated environment. From this `manuscript` directory:

```sh
MANUSCRIPT_PYTHON=/path/to/python bash scripts/build.sh
```

The build regenerates figures and tables from checked-in evidence, compiles the paper, verifies numerical claims and references, and creates the compact bundle. It performs no mapping searches. The first Tectonic build may download TeX packages. A saved `preprint.bbl` accompanies the source bundle.

For visual review, render the PDF with Poppler:

```sh
pdftoppm -scale-to 1400 -png manuscript.pdf build/page
```

`EDITORIAL_NOTES.md` records the author-supplied metadata and dataset provenance still needed for submission. These are not fabricated in the paper. This is an arXiv-oriented manuscript, not a submitted arXiv record. The build checks the PDF and compact bundle for execution-site disclosures; only CPU models are reported as hardware information.

## Archived studies and examples

The archived Golden comparison retains one-, two-, three- and ten-seed settings, uncut controls, and released comparators. Its cap study uses 1,403 common completed reactions on an AMD EPYC 9J14; its orientation study uses 1,821 common completed reactions on an Apple M4 Pro. Default-call times with an unrecorded CPU model are marked separately. No cross-CPU or current-implementation speed ranking follows from these tables.

The 140-reaction coordinate collection compares bond events under unverified mappings. Its archived cap control returns complete mappings for 138 reactions at cap 100 and all 140 at cap 2,000. At cap 2,000, explicit event-window queries retained 300 event patterns across 120,052 saved families, including 166 patterns at MAPPA's own minima. These are historical event-query results, not a mapping-accuracy benchmark or an exhaustive structural search. Their evidence and timing scopes remain in [the coordinate report](../reports/published_baseline_20260915/) and `evidence/end_to_end_timing.json`.

The Golden case 9 illustration uses three selected archived witnesses: the recovered heavy-atom reference and two alternatives. Event totals include H, whose individual identities are not reference labels. The witnesses are neither an exhaustive decoded set nor validated mechanisms. Their mappings, constraints and checks remain in `golden_case9.json` and `golden_case9_verification.json`.

## 3D research preview

[Grow. Branch. Decode.](animations/mappa_research_preview/index.html) is an archived 22-second 3D film of Golden case 15. It shows two event classes from a saved catalogue containing 44 unordered branches, 96 families and nine classes through six events. It illustrates conditional placements and event alternatives, not a trace of the current growth implementation. [MP4](animations/mappa_research_preview/mappa-grow-branch-decode.mp4) · [GIF](animations/mappa_research_preview/mappa-grow-branch-decode.gif).

## Authors

Yunheng Zou; Olalla Nieto Faza; Shifa Hussain; **Varinia Bernales**<sup>†</sup>;
**Alán Aspuru-Guzik**<sup>†</sup>. † Principal investigators. Shifa Hussain is
assigned to the Department of Biology, University of Toronto Mississauga.
[Affiliation provenance](evidence/author-affiliations.md) records the supplied
assignments and official institutional sources. `main` is the publication branch.
