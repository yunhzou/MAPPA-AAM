# Numerical manuscript evidence

`current_benchmarks.json` records the current Golden and EnzymeMap coverage,
stage CPU measurements. `branch_offsets.json.gz`
records the completed branch-event analysis, with exact offsets and bounded intervals.
`paper_sources.json` pins every bundled evidence snapshot by hash and labels its
scope. `paired_slap.json` records the current paired comparator study;
`minimum_choices.json` and its compressed per-reaction file record the minimum
choice count and incomplete outcomes. `design_examples.json` contains analytical
pair-count checks for the constructed graph illustrations. Other comparator, cap,
direction and coordinate snapshots describe
archived September studies, not current-version accuracy or performance.

The figure build reads these immutable snapshots and runs no matching or decoding.
Replace evidence only from measured, audited results; never fill incomplete bins
with guessed values. `PAPER_EVIDENCE.md` summarizes the claim boundaries.
