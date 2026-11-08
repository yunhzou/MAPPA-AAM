# MAPPA manuscript

The main manuscript presents the current scientific design: one-witness growth,
compressed completion, hierarchical branching and cut sweep, atom-set reduction,
shared binary-symmetry event responses, branch composition and event windows.
Portable execution and fragment reuse have their own Methods subsection.
Reference verification is external evaluation, outside the MAPPA pipeline.
The introduction and molecular illustrations are largely retained.

## Current results

| Study | Catalogue recovery | Saved decoded min+5 recovery | Decoded absent / unknown |
|---|---:|---:|---:|
| Golden, one seed | 1,834 / 1,851 | 1,832 / 1,851 | 16 / 3 |
| Golden, ten seeds | 1,843 / 1,851 | 1,841 / 1,851 | 8 / 2 |
| EnzymeMap, one seed | 47,798 / 47,974 | 47,796 / 47,974 | 178 / 0 |

All requested mapping and decoding outputs are available. Decoded verification
uses exact annotated heavy-atom indices and replayed complete witnesses.
Golden catalogue verification additionally allows endpoint chemical symmetry;
EnzymeMap uses exact indices for both metrics. Unknowns are not absences.

Mean compute CPU seconds per reaction for mapping (matching plus reduction)
and decoding are respectively **1.31 / 1.88**, **11.64 / 9.28**, and **4.41 / 8.75**
for these three studies. Both directions and every cut are included. Publication,
external evaluation and failed-attempt accounting are separate. These are saved
stage measurements with their recorded versions and recoveries, not a new
uninterrupted run or a hardware-normalized speed ranking.

The branch-event histogram is explicitly **preliminary**: 45,737 exact offsets,
3,891 pending offset checks, 194 absent references and three unresolved
memberships. A GT-containing branch's unconditioned minimum is compared with
the overall detected count using the same saved all-atom model. The completed
subset overrepresents zero and cannot establish the final distribution.

## Sources and build

- `preprint.tex`: entry point; author list and supplied affiliations are retained.
- `includes/paper.tex`: main Methods, Evaluation, Results and Discussion.
- `includes/supplement.tex`: technical details and clearly separated archived studies.
- `evidence/current_benchmarks.json`: sanitized coverage, timing and cache-control evidence.
- `evidence/branch_offsets.json.gz`: complete reaction accounting for the preliminary figure.
- `evidence/paper_sources.json`: evidence hashes and source-to-claim scopes.
- [REPRODUCIBILITY.md](REPRODUCIBILITY.md): settings, source pins and evaluation definitions.
- `figs/fig1_algorithm.*`: molecular overview, with no verifier in the AAM pipeline.
- `figs/fig6_current_benchmarks.*`: current decoded coverage and separate stage CPU costs.
- `figs/fig7_reference_offsets.*`: preliminary branch-event distribution.
- `figs/fig4_alternatives.*`: selected archived Golden case 9 witnesses, labeled as illustrations.
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
not physical event-bin accuracy. September comparator, cap, direction and
coordinate results retain their original protocols in the supplement. They do
not establish current-version superiority. The coordinate collection lacks
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
