# Reproducing the manuscript

The paper separates scientific design, measured implementation versions and
external evaluation. Rebuilding it reads saved evidence and performs no mapping.

## Build and inspect

From the repository root:

```sh
MANUSCRIPT_PYTHON=/path/to/python bash manuscript/scripts/build.sh
pdftoppm -scale-to 1400 -png manuscript/manuscript.pdf manuscript/build/page
```

Python dependencies are in `requirements-build.txt`; Tectonic and Poppler are
also required. The build checks evidence hashes, numeric totals, citations,
figure inclusion and publication privacy. The Overleaf export contains the same `main.tex`, shared source files and
figure PDFs used for the local build. Evidence and figure builders remain in
this repository. The visual-review record identifies the final PDF checksum.

## Scientific settings

The paper calls the endpoint on which fragments grow the base and the other
endpoint the target. Each direction assigns these roles independently; reverse
search exchanges them. Existing saved-data keys and atom-index frames are
unchanged.

Current full studies use original prepared formal bond orders with explicit H,
discovery tolerance 2, graph/cut threshold 0.2, branch cap 100, random seed 42,
both directions and uncut plus every eligible single-edge cut. Golden compares
one and ten seed orderings per cut; EnzymeMap uses one. Bond events use original
bond orders and empirical change thresholds: 0.5 for organic bonds (the default
for non-metal pairs), and 0.3 for bonds involving at least one metal atom. An event
requires a change strictly greater than the applicable threshold. The decoding
window extends through the best attained model count plus five, with a five-second
improvement allowance per unique fragment. A bounded result is not automatically
a proved minimum.

[design.md](../design.md) is the current scientific contract.
[The Python API](../docs/PYTHON_API.md) and
[execution guide](../docs/shared-execution.md) describe the separate stage APIs.
The manuscript introduces the scientific modules; those guides retain interface
and persistence details.

## Golden evidence

The full source reports are in
[the Golden campaign](../reports/golden_native_seeds_20261106/README.md).
Matching was measured at `27b18abad691648cf2a87d2e5ae43db117d80b66`.
Reduction and decoding use `7b248506a4d541cf2c50540e41ce6204bc87bf80`, including
the accepted one-record-per-atom-set reducer. Every saved matching result was
reused; snapshot conversion is separately measured and is not search time.
These full matching passes predate indexed containment reuse. Their cache-hit
ratios must not be attributed to the later cache implementation.

Reference recovery within the decoded minimum + 5 window is **1,833/1,851**
with one seed and **1,842/1,851** with ten. Every reaction has a definitive
outcome. Verification preserves exact annotated heavy-atom indices and permits
completion of unannotated hydrogens. The final targeted closure uses
external-verifier revision `dc92f469269f0e1dfcf0db9f76abd7ed19136903`, with
unchanged saved mapping and decoding outputs. Each positive has an injective,
replayed witness and an independently counted event total inside the saved window.
The benchmark report preserves the complete external-audit history.

## EnzymeMap evidence

[The completed EnzymeMap report](../reports/enzymemap_native_20261106/README.md)
covers 47,974 prepared reactions and 95,948 directions. Fresh matching began
from the indexed-reuse milestone `c437ed6`; only original immutable inputs were
reused from earlier campaigns. Subsequent native and storage corrections
preserved the declared mapping relation and event model. Final scientific
source/engine identity is
`ba841d410beb9356c1bf1d86c9c34467a42e06d9e82a66cfd45c4aa65fa92e6b`.
The final external decoded audit uses `91eea28`.

Reference recovery within the saved minimum + 5 window is **47,796/47,974**,
with 178 absent and no unresolved reaction-level outcomes. All 94,820 positive
direction certificates passed replay, injectivity and independent counting.
Negative outcomes require conclusive exclusion: 133 were established by direct
decoded queries and 45 by exclusion from the containing matching relation.
Bounded unknown direction queries remain recorded separately and are never
relabeled as direct negative answers. The source report retains this evidence.

The accepted independent-fragment model differs from the physical whole-reaction
count for 2,801 positive directional witnesses. Decoded coverage refers to the
model window and is not a claim of physical event-bin accuracy. EnzymeMap's
curated mappings are reference annotations, not independent mechanism proofs.
Original input identities and complete scope are preserved in the report's
source evidence; no new preparation or dataset filtering was performed for this
manuscript rewrite.

## Timing boundaries

`evidence/current_benchmarks.json` contains successful compute CPU totals for
fragment discovery, reduction and decoding. Published matching time sums discovery
and reduction; reduction is a subphase, not a separate scientific stage. Means
divide the sum over both directions by the reaction denominator. Stage input loading and result construction are
included. Publication, external verification, diagnostics and
historical failed attempts are separate. Medians of stages are not additive.
CPU is accumulated computation, not parallel elapsed latency.

Golden reuses matching and corrected postprocessing from separate passes.
EnzymeMap combines original successful tasks and successful recoveries. Neither
is a new uninterrupted full run of one final performance revision. CPU models
were not retained in the exported full-campaign summaries. There is no
hardware-normalized comparison with the archived baselines. EnzymeMap's
successful computation averages 13.16 CPU seconds per reaction, or 14.60
including publication and measured worker overhead. Earlier failed attempts
add approximately 27.93 sampled CPU hours; the independent decoded audit has
separate accounting in its source report.

## Reference-mapping event offsets

[The offset report](../reports/reference_event_offsets_20261107/README.md) and
`evidence/reference_offsets.json.gz` record the counts used in Figure 7:

`n = verified GT-consistent mapping event count - saved detected minimum`.

Every reaction contributes one concrete count difference. Both terms use the
saved all-atom decoder model. For recovered references, select the lower event
count among the final saved forward/reverse verification witnesses, with forward
winning ties. The numerator is the count of that actual mapping, including its
hydrogen completion. It is not an optimized count over all reference completions.
The denominator is the original detected minimum across directions, unchanged.
No new matching, decoding, membership or minimum-event search is run for this plot.

The Golden counts come from `decoded_verification_closure_directions.json.gz`.
The EnzymeMap counts were extracted from the final audited direction records;
the snapshot retains their source identity and each selected certificate,
decoded-output and reference hash. The export contains counts and identifiers,
not molecular inputs. `build_benchmark_figures.update_reference_offsets` performs
the deterministic selection and subtraction from these saved records.

Saved anchored completions supply the 18 Golden and 178 EnzymeMap references
absent from the original decoded output. `evidence/anchored_offsets.json` retains
the original detected minima and both audited completion counts. These use one
seed per direction, no sweep, tolerance two and cap 100. All 392 witnesses
preserve their supplied pairs and pass independent counting; physical and model
counts agree. Eleven EnzymeMap references required a documented index translation.
These completions do not increase independent recovery.

The merged histograms are Golden `{0:1630, 1:90, 2:88, 3:24, 4:11, 5:6, 7:1,
12:1}` and EnzymeMap `{-2:9, 0:44086, 1:1384, 2:2190, 3:122, 4:171, 5:9, 6:3}`.
They contain all 1,851 and 47,974 reactions. Each dataset uses one color, with no
separate completion legend. Unoccupied negative ticks are omitted. The nine
negative values are actual anchored completions below the original detected
minimum; no lower bounds are plotted.

The appendix percentage view uses the same bins and plotting function, scaling
each count by `100 / dataset size`. The denominators are 1,851 and 47,974, so
each dataset sums to 100% before displayed rounding. The logarithmic axis and
dataset colors are shared with the main figure.

The prior branch-minimum analysis remains in `evidence/branch_offsets.json.gz`
for provenance and original-baseline checks. It does not supply the plotted
numerator. Replacing its values and lower-bound substitutions moves 126 Golden
and 2,196 EnzymeMap reactions to higher offsets. Coverage, saved decoding minima
and the separate minimum-choice analysis are unchanged. Case 268 in Golden now
has offset four: verified count 16 minus detected minimum 12.

## Evidence retained outside the paper

Default-comparator coverage appears in the main Results: the returned output is
a candidate set, with single-answer predictions treated as singleton sets. Saved
outputs were re-evaluated on all 1,851 inputs, without changing any per-case verdict.
Historical MAPPA, direction, seed/cap, 140-coordinate and comparator timing studies
are excluded from the appendix and Overleaf export. Their original records remain
in the repository and Git with their settings, source hashes and timing cohorts.
`evidence/paper_sources.json` indexes their unchanged snapshots. The coordinate
collection has no verified atom mappings and incomplete WBO-generation provenance;
its results compare event counts rather than mapping accuracy. Selected case 9
witnesses illustrate alternative atom origins, not an exhaustive current decoded
set or validated mechanism.

All publication data omit operational settings, addresses, account names,
resource inventories and private storage locations. Only recorded CPU model
names are used when identifying hardware.

## Paired SLAP comparison and minimum-event choices

The [paired comparison](../reports/golden_slap_paired_20261108/README.md)
uses identical prepared Golden inputs, one CPU model and common chemical
reference equivalence. MAPPA and standard SLAP have 1,849 common successful
reactions; requiring the expanded SLAP sweep as well leaves 1,834. Keep both
timing denominators explicit. The expanded wrapper changes cuts, directions and
weight modes relative to the standard forward-only call; its gain is not a
sweep-only ablation. `evidence/paired_slap.json` retains the scientific summary
and hashes of the full source report and independent verification audit.

The [reaction-core analysis](../reports/golden_minimum_branches_20261110/README.md)
uses the original one-seed Golden decoded outputs and their unchanged detected
minima. It combines saved qualifying witnesses across both directions, then
groups their signed physical bond-change patterns with MAPPA's existing
`signed_events/v3` identity. Hydrogen events are included, equivalent atom labels
count once, and fragment partitioning does not define a new core. All 1,851
reactions are represented. No matching, decoding or additional symmetry search
is performed. Evidence is in `evidence/minimum_reaction_centers.json` and
`evidence/minimum_reaction_centers_reactions.json.gz`. The older mapping-class
and branch-composition tallies remain historical evidence; neither supplies the
current figure. These observed patterns do not count all alternatives hidden
inside unqueried symmetry responses.

The constructed design figures are checked by direct pairwise bond comparison
in `build_molecule_figure.py`. They are graph illustrations without hydrogen
constraints, not new benchmark measurements or full chemistry simulations.
