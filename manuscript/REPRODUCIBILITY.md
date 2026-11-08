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
figure inclusion and publication privacy. The compact bundle includes the
source, checked evidence and figure builders. A visual-review record is included
only when its checksum matches the final PDF.

## Scientific settings

Current full studies use original prepared formal bond orders with explicit H,
discovery tolerance 2, graph/cut threshold 0.2, branch cap 100, random seed 42,
both directions and uncut plus every eligible single-edge cut. Golden compares
one and ten seed orderings per cut; EnzymeMap uses one. Bond events use original
weights and strict thresholds 0.5, or 0.3 for metal-involving pairs. The decoding
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

Catalogue recovery is 1,834/1,851 with one seed and 1,843/1,851 with ten. Its
endpoint chemical-equivalence policy preserves charge, isotope, H count and
stereochemistry. Exact-index decoded recovery is separately 1,832 with three
unknowns and 1,841 with two unknowns. One-seed targeted decoded verification uses
`b6d8c086cf92583d099509c7c7078daaf3cfc0af`; the ten-seed audit uses
`1591a8c90d810ba42b27273f19de3ea1278dd4e5`. Each positive has an injective,
replayed witness in the saved response and window. One ten-seed directional
catalogue timeout remains in its attempt history; the opposite direction
resolves that reaction's membership.

The three-reaction indexed-cache control uses cases 895, 1261 and 1841, both
directions, one and ten seeds, and an Apple M4 Pro CPU. It compares three cold
retention-on repetitions with one zero-retention control. All 48 measured
direction runs agree on catalogue fingerprints, scope and cap outcomes. The
control disables immutable sharing as well as existence reuse; it is not a
lookup-only speedup or a full-dataset scaling result.

## EnzymeMap evidence

[The completed EnzymeMap report](../reports/enzymemap_native_20261106/README.md)
covers 47,974 prepared reactions and 95,948 directions. Fresh matching began
from the indexed-reuse milestone `c437ed6`; only original immutable inputs were
reused from earlier campaigns. Subsequent native and storage corrections
preserved the declared mapping relation and event model. Final scientific
source/engine identity is
`ba841d410beb9356c1bf1d86c9c34467a42e06d9e82a66cfd45c4aa65fa92e6b`.
The final external decoded audit uses `91eea28`.

Exact-index catalogue recovery is 47,798 reactions. Saved decoded min+5 recovery
is 47,796, with 178 absent and no unresolved reaction-level outcomes after
combining direct queries with conclusive catalogue absence. The 250 bounded
unknown direction queries remain recorded as unknown. All 94,820 positive
direction certificates passed replay, injectivity and independent counting.
The two additional catalogue references are excluded by the saved event window.

The accepted independent-fragment model differs from the physical whole-reaction
count for 2,801 positive directional witnesses. Decoded coverage refers to the
model window and is not a claim of physical event-bin accuracy. EnzymeMap's
curated mappings are reference annotations, not independent mechanism proofs.
Original input identities and complete scope are preserved in the report's
source evidence; no new preparation or dataset filtering was performed for this
manuscript rewrite.

## Timing boundaries

`evidence/current_benchmarks.json` contains successful compute CPU totals for
matching, reduction and decoding. Means divide the sum over both directions by
the reaction denominator. Stage input loading and result construction are
included. Publication, external catalogue/decoded audits, diagnostics and
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

## Branch-event offsets

[The offset report](../reports/reference_event_offsets_20261107/README.md) and
`evidence/branch_offsets.json.gz` preserve the accepted preliminary snapshot.
For each reaction, find the smallest unconditioned model minimum among retained
branch contexts admitting the exact reference inside the saved window, and
subtract the overall detected count across both directions. The GT witness's
own event count is only an upper bound on that branch minimum.

Both terms use the same original all-atom decoder model. No heavy-only count or
new global physical optimization is substituted. Analysis revision `946fc6a`
uses saved response queries. Fifteen affected tests include a reference witness
with two events in a zero-minimum branch and a lower branch that excludes GT.
Unknown answers retain intervals outside exact bins. Zero-offset certificates
resolve first, so the partial histogram overrepresents zero. Its 3,891 pending
checks must finish before reporting the full distribution or percentages.

## Archived analyses

The September comparator, direction, seed/cap and 140-coordinate studies remain
in the supplement with their original settings, source hashes and timing cohorts.
`evidence/paper_sources.json` indexes their unchanged snapshots. The coordinate
collection has no verified atom mappings and incomplete WBO-generation provenance;
its results compare event counts rather than mapping accuracy. Selected case 9
witnesses illustrate alternative atom origins, not an exhaustive current decoded
set or validated mechanism.

All publication data omit operational settings, addresses, account names,
resource inventories and private storage locations. Only recorded CPU model
names are used when identifying hardware.
