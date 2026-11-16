# MAPPA manuscript

This directory is the single Overleaf-compatible source project. Edit `main.tex`
and its `includes/` files here. `mappa_overleaf.zip` is a generated upload copy;
no separate editable export is maintained.

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

| Study | Reference recovery | Not recovered |
|---|---:|---:|
| Golden, one seed | 1,833 / 1,851 (99.028%) | 18 |
| Golden, ten seeds | 1,842 / 1,851 (99.514%) | 9 |
| EnzymeMap, one seed | 47,796 / 47,974 (99.629%) | 178 |

Coverage refers to verified mappings within the decoded minimum + 5 window.
The check preserves annotated heavy-atom identities and requires a complete,
replayed witness and an independent event count. Every reaction has a definitive
outcome. The paper reports this single coverage measure; the underlying benchmark
reports retain the full audit evidence. Matching, decoding and their measured
CPU times are unchanged.

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

Expanded SLAP covers **1,794/1,851** reactions and standard bidirectional binary
SLAP covers **1,620/1,851**. Default forward-only SLAP gives **85.95%** inclusion,
consistent with the original paper's **86.9%** on a different Golden version.
The original score counts recovery among alternatives, not top-1 selection.
The paired study supplies the computational-cost comparison; Figure 9 reports
timing only. MAPPA coverage comes from the completed decoded-output audits above.

On 1,834 reactions completed by all paired protocols, mean search CPU is
**0.111 s** for standard SLAP, **1.202 s** for MAPPA and **12.475 s** for the
expanded SLAP sweep. The broader standard-SLAP/MAPPA pair has 1,849 reactions
and means **0.487/1.259 s**. Speed claims retain their protocol and denominator.
The expanded wrapper changes weight mode and search scope; it is not a
sweep-only ablation. All comparator verdicts passed an independent audit.

Reaction-core counting covers all **1,851 Golden reactions**. The saved
minimum-event mappings give **1,252 reactions with one observed core (67.64%)**,
a median of **one**, mean **2.97**, 95th percentile **four** and maximum **522**.
Cores identify signed bond changes, including hydrogen events; equivalent atom
labels and different fragment partitions with the same pattern count once.
The complete analysis took **15.66 seconds** on the local Mac with eight processes,
including loading and replaying the saved choices. No matching, decoding or
additional symmetry search was performed. The earlier branch-composition count
is retained only as historical evidence.

The paper reports approximately linear matching cost over the measured one-to-ten
seed range and omits the cache-control comparison. No new experiments were run.

## Sources and build

For Overleaf, upload `mappa_overleaf.zip` as a new project, select
`main.tex` as the main document and use **XeLaTeX**. The ZIP includes the
source, supplement, bibliography, figures, class and custom font; no Python
processing is needed. The maintained bundle command also creates this archive.

The title, full affiliations, abstract and institutional logos retain their
separate opening page. Scientific figures do not appear beneath the abstract.
The molecular overview appears on the first main-text page, in the Introduction.
Other figures are placed at their first discussion where space permits, or on
the adjacent page; the final PDF is checked against those citation locations.

- `main.tex`: sole entry point; author list and supplied affiliations are retained.
- `includes/paper.tex`: main Methods, Evaluation, Results and Discussion.
- `includes/supplement.tex`: current algorithm details, benchmark settings, measurement definitions and the percentage event-offset figure.
- `evidence/current_benchmarks.json`: sanitized coverage and stage timing evidence.
- `evidence/reference_offsets.json.gz`: verified mapping counts, original detected minima, offsets and certificate hashes for every reaction.
- `evidence/branch_offsets.json.gz`: superseded branch-minimum analysis, retained as historical evidence and the original-baseline source.
- `evidence/paper_sources.json`: evidence hashes and source-to-claim scopes.
- [REPRODUCIBILITY.md](REPRODUCIBILITY.md): settings, source pins and evaluation definitions.
- `figs/fig1_algorithm.*`: molecular overview. The lower panels show both endpoint atom sets for each fragment, then compare both sides of the same fragment under two oxygen assignments. Boundary changes are visible on the full endpoints. Verification remains outside the AAM pipeline.
- `figs/fig6_current_benchmarks.*`: decoded-output reference coverage and matching/decoding CPU costs.
- `figs/fig7_reference_offsets.*`: verified GT-consistent mapping counts relative to the saved detected minimum, with one color per dataset.
- `figs/fig12_reference_percentages.*`: appendix companion with each dataset normalized to 100%, using the same bins and logarithmic axis.
- `figs/fig11_saturation.*`: one carbon ring, two distinct target ring placements and two branches; internal rotations remain compressed.
- `figs/fig8_design_examples.*`: molecular drawings of the sweep mask, cyclopropanol/propan-2-ol assignments and acetonitrile/ethylamine boundary events. These are heavy-atom illustrations with hydrogen constraints omitted.
- `figs/fig9_paired_slap.*`: paired comparator CPU costs.
- `figs/fig10_minimum_choices.*`: observed reaction cores among saved minimum-event mappings.
- `figs/fig4_alternatives.*`: accepted EnzymeMap case 13175. Carbon 1 enters the three-carbon molecule in the reference and CO₂ in the alternative. The drawings show four versus two heavy-atom events; exhaustive hydrogen completion of each fixed heavy-atom mapping gives six total events for both. This illustrates ambiguity, not a strict all-atom advantage. All nine heavy atoms are shown, with consistent atom-origin colors and only carbons numbered. Portable mappings and count checks are in `evidence/enzymemap_case13175.json`.

Install the dependencies in `requirements-build.txt`, Tectonic and Poppler. From
the repository root:

```sh
MANUSCRIPT_PYTHON=/path/to/python bash manuscript/scripts/build.sh
```

The build regenerates figures and tables, compiles `manuscript.pdf`, validates
numeric evidence and references, checks publication privacy, and writes
`mappa_overleaf.zip` directly from the same source project. It performs no matching or decoding. Visual inspection
is recorded separately against the final PDF checksum.

The decoder's ordinary-boundary model is approximate; its selected starting
orbits and bounded improvement restrict coverage. Model-window membership is
not physical event-bin accuracy. Default-comparator coverage is included in the main Results with its recorded
evaluation scope. Historical MAPPA, cap, direction, coordinate and timing studies are excluded
from the appendix and Overleaf export. Their evidence remains in repository
records and Git. The appendix retains only support for the current paper.

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
