# Full Golden decoding with binary fragment symmetry

This study tests the integrated decoder from `4759d25` on the saved full Golden
matching catalogue: 1,851 reactions, 3,702 directions, discovery tolerance 2,
one seed, branch cap 200 and the full single-edge sweep. Matching is not repeated.

An explicit benchmark input conversion removes the retired fragment
self-tolerance filter and preserves the saved binary generators, original WBOs,
discovery constraints, atom sets, placement actions and branch histories. The
original catalogues remain unchanged. Converted inputs do not constitute a
fresh benchmark of the current matcher. Conversion evidence retains both input
and output checksums and the original matching configuration.

Conversion is independently checked against the producing revision's complete
record export for Golden 625 forward. All fields match after the declared schema
and symmetry-policy changes. Full decoding and a save/load check completed in
0.051 seconds locally. Golden 569 forward exceeded a separate 240-second local
watchdog during full decoding; this is not a fragment-only timing.

The full cluster test retains event thresholds 0.5 / 0.3 and minimum plus 5,
using the current default five-second fragment optimization allowance. Each
direction has a 600-second compute watchdog and a 2 GiB memory guard, matching
the original full campaign's limits. Conversion, decoding compute and output
publication are measured separately. Timings will retain failed directions in
the denominator and identify the largest completed decoding time separately
from censored timeout measurements. No matching rerun or recurring monitor is
started by this experiment.

Full-run results are pending. The fragment-level beta measurements do not
establish full branch-decoding performance.
