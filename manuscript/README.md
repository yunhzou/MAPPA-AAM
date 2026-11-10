# MAPPA manuscript

The Methods section opens with the complete workflow: grow a fragment, discover
its placements, continue conditional branches, reduce repeated fragments and
compositions, then decode bond events. Explicit Matching Stage and Decoding Stage subsections introduce the two
scientific stages. Reduction is the final subphase of matching.
Constraint-discovery philosophy is integrated into growth, sweep and boundary
interpretation. Detailed matching equations and full pseudocode are in the
supplement; parallel execution has its own final Methods subsection.
Reference verification is external evaluation, outside the MAPPA pipeline.
Growth matching, saturation matching and hierarchical branching have separate
subsections. A zero/one/many transition table distinguishes the existence decision
from completed placement construction. Decoding likewise separates local symmetry,
internal events, connectivity-mismatch handling, ordinary boundaries, sweep
correction and window selection. Each operation is paired with its inputs,
decision and molecular example.
The Introduction develops AAM's purpose and applications, ambiguity, taxonomy
and prior work, the computational challenge of retaining alternatives, and
MAPPA's contribution and principal findings. Figure 2 supports the ambiguity
discussion before the output taxonomy. Alternative-generating solvers and
single-answer predictors are compared in prose and a concise table; a later
selector connects their objectives. MAPPA and SLAP belong to the first category.
Methods introduces base and target as search roles: either endpoint can be the
base, and reverse search exchanges the roles. The notation and molecular diagrams
follow this convention.

## Current results

| Study | Catalogue recovery | Saved decoded min+5 recovery | Decoded absent |
|---|---:|---:|---:|
| Golden, one seed | 1,834 / 1,851 | 1,833 / 1,851 | 18 |
| Golden, ten seeds | 1,843 / 1,851 | 1,842 / 1,851 | 9 |
| EnzymeMap, one seed | 47,798 / 47,974 | 47,796 / 47,974 | 178 |

All requested mapping and decoding outputs are available. Decoded verification
uses exact annotated heavy-atom indices and replayed complete witnesses.
Golden catalogue verification additionally allows endpoint chemical symmetry;
EnzymeMap uses exact indices for both metrics. Catalogue evaluation has zero unresolved reactions in all three studies. Every reaction also has a definitive decoded-window outcome after the targeted verification closure. All positive checks require a full replayed mapping and an independent event count. Matching, decoding and their measured CPU times were unchanged.

Mean compute CPU seconds per reaction for matching (including reduction)
and decoding are respectively **1.31 / 1.88**, **11.64 / 9.28**, and **4.41 / 8.75**
for these three studies. Both directions and every cut are included. Publication,
external evaluation and failed-attempt accounting are separate. These are saved
stage measurements with their recorded versions and recoveries, not a new
uninterrupted run or a hardware-normalized speed ranking.

The event-offset figure subtracts the saved detected minimum from the event count
of a verified GT-consistent full mapping. It uses the lowest already recorded
verified count across directions, with no new optimization. Saved anchored
completions supply the 196 references absent from the original output. Every
reaction contributes once, using the same all-atom event model on both sides.
Each dataset uses one color; unoccupied negative ticks are omitted. Independent
coverage is unchanged.

## Comparator and minimum-choice results

The main Results compare MAPPA, SLAP, LocalMapper, RXNMapper, Chython, RDT and
Indigo by reference coverage. A single prediction counts as a singleton candidate
set. Current MAPPA scores, paired SLAP protocols and re-evaluated archived default
outputs retain their provenance; only the paired study supports speed ratios.

The separate paired study gives MAPPA catalogue coverage of **1,839/1,851**,
expanded SLAP coverage of **1,794/1,851**, and standard bidirectional binary SLAP
coverage of **1,620/1,851**. Default forward-only SLAP gives **85.95%** inclusion,
consistent with the original paper's **86.9%** on a different Golden version.
The original score counts recovery among alternatives, not top-1 selection.
These separate runs are not pooled with the seed-study table above.

On 1,834 reactions completed by all paired protocols, mean search CPU is
**0.111 s** for standard SLAP, **1.202 s** for MAPPA and **12.475 s** for the
expanded SLAP sweep. The broader standard-SLAP/MAPPA pair has 1,849 reactions
and means **0.487/1.259 s**. Speed claims retain their protocol and denominator.
The expanded wrapper changes weight mode and search scope; it is not a
sweep-only ablation. All comparator verdicts passed an independent audit.

Minimum-event branch counting covers all **1,851 Golden reactions**. Distinct
fragment compositions have a median of **31**, mean **51.06**, 95th percentile
**120** and maximum **2,250**. Internal shuffles count once. The calculation uses
saved choices without further search and takes **36.7 seconds** with eight
processes. Branch compositions can overlap in their represented mappings, so
these counts are not a random-selection accuracy measure.

The paper reports approximately linear matching cost over the measured one-to-ten
seed range and omits the cache-control comparison. No new experiments were run.

## Sources and build

The title, full affiliations, abstract and institutional logos retain their
separate opening page. Scientific figures do not appear beneath the abstract.
The molecular overview appears on the first main-text page, in the Introduction.
Other figures are placed at their first discussion where space permits, or on
the adjacent page; the final PDF is checked against those citation locations.

- `preprint.tex`: entry point; author list and supplied affiliations are retained.
- `includes/paper.tex`: main Methods, Evaluation, Results and Discussion.
- `includes/supplement.tex`: technical details and clearly separated archived studies.
- `evidence/current_benchmarks.json`: sanitized coverage and stage timing evidence.
- `evidence/reference_offsets.json.gz`: verified mapping counts, original detected minima, offsets and certificate hashes for every reaction.
- `evidence/branch_offsets.json.gz`: superseded branch-minimum analysis, retained as historical evidence and the original-baseline source.
- `evidence/paper_sources.json`: evidence hashes and source-to-claim scopes.
- [REPRODUCIBILITY.md](REPRODUCIBILITY.md): settings, source pins and evaluation definitions.
- `figs/fig1_algorithm.*`: molecular overview. The lower panels show both endpoint atom sets for each fragment, then compare both sides of the same fragment under two oxygen assignments. Boundary changes are visible on the full endpoints. Verification remains outside the AAM pipeline.
- `figs/fig6_current_benchmarks.*`: full-output catalogue coverage and matching/decoding CPU costs.
- `figs/fig7_reference_offsets.*`: verified GT-consistent mapping counts relative to the saved detected minimum, with one color per dataset.
- `figs/fig12_reference_percentages.*`: appendix companion with each dataset normalized to 100%, using the same bins and logarithmic axis.
- `figs/fig11_saturation.*`: one carbon ring, two distinct target ring placements and two branches; internal rotations remain compressed.
- `figs/fig8_design_examples.*`: molecular drawings of the sweep mask, cyclopropanol/propan-2-ol assignments and acetonitrile/ethylamine boundary events. These are heavy-atom illustrations with hydrogen constraints omitted.
- `figs/fig9_paired_slap.*`: paired comparator coverage and CPU costs.
- `figs/fig10_minimum_choices.*`: recorded branch compositions at the detected minimum.
- `figs/fig4_alternatives.*`: accepted EnzymeMap case 13175. Carbon 1 enters the three-carbon molecule in the reference and CO₂ in the alternative. The drawings show four versus two heavy-atom events; exhaustive hydrogen completion of each fixed heavy-atom mapping gives six total events for both. This illustrates ambiguity, not a strict all-atom advantage. All nine heavy atoms are shown, with consistent atom-origin colors and only carbons numbered. Portable mappings and count checks are in `evidence/enzymemap_case13175.json`.
- `figs/fig2_golden.*` and `figs/fig3_coordinate.*`: archived supplementary studies.

Install the dependencies in `requirements-build.txt`, Tectonic and Poppler. From
the repository root:

```sh
MANUSCRIPT_PYTHON=/path/to/python bash manuscript/scripts/build.sh
```

The build regenerates figures and tables, compiles `manuscript.pdf`, validates
numeric evidence and references, checks publication privacy, and writes
`manuscript_bundle.zip`. It performs no matching or decoding. Visual inspection
is recorded separately against the final PDF checksum.

The decoder's ordinary-boundary model is approximate; its selected starting
orbits and bounded improvement restrict coverage. Model-window membership is
not physical event-bin accuracy. Default-comparator coverage is included in the main Results with its recorded
evaluation scope. Historical MAPPA, cap, direction, coordinate and timing results
retain their original protocols in the supplement; they do not establish
current-version speed ratios. The coordinate collection lacks
annotated mappings and complete WBO-generation provenance.

The archived [3D research preview](animations/mappa_research_preview/index.html)
illustrates its saved historical catalogue, not an execution trace of the current
matcher. [EDITORIAL_NOTES.md](EDITORIAL_NOTES.md) records remaining author metadata
and dataset-provenance decisions; those have not been invented.

## Authors

Yunheng Zou; Olalla Nieto Faza; Shifa Hussain; **Varinia Bernales**<sup>†</sup>;
**Alán Aspuru-Guzik**<sup>†</sup>. † Principal investigators. Shifa Hussain is
assigned to the Department of Biology, University of Toronto Mississauga.
[Affiliation provenance](evidence/author-affiliations.md) records the supplied
assignments and official institutional sources. `main` is the publication branch.
