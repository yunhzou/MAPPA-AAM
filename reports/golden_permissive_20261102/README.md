# Golden permissive matching comparison

Historical study at revision `4ba1939`: the four-case pilot completed and the
full comparison was launched locally. This report preserves that observation;
it does not establish a currently running job or a completed full comparison.
The [current decoder validation](../golden_binary_decode_20261105/README.md)
uses a later algorithm and a separate saved-catalogue scope.

The user requested two full Golden runs, each with one seed, branch cap 2,000,
random seed 42, both directions and the full single-edge cut sweep:

| Arm | Matching matrices | Matching tolerance | Event matrices |
|---|---|---:|---|
| Binary | 1 where the original WBO is at least 0.2; 0 otherwise | 1 | Original WBOs |
| Weighted | Original WBOs | 10 | Original WBOs |

The scope is 1,851 reactions and 3,702 directions per arm. Atom indexing,
elements, reference annotations and event thresholds remain fixed. The ordinary
event threshold is 0.5 and the metal threshold is 0.3. Binary matrices affect
matching and its symmetry compression; they must not replace the event weights.
Saved matching constraints remain attached to their own input matrices.

The comparison will report reaction and direction reference recovery, paired
gains and losses, capped searches, unknown verdicts and resource failures.
Reference verification uses the existing Golden heavy-atom policy with endpoint
chemical equivalence. CPU time is recorded separately for matching, symmetry
finalization, family reduction, decoding and verification. Reported decoder
output is the retained compressed event response, not an exhaustive enumeration
or a proven minimum.

Both arms finish matching and reduction before decoding starts. A persistent
local work pool keeps each direction's completed-fragment cache across cuts.
Catalogues, decoded results and attempt evidence use aggregated bundles.
Terminal failures are retained; resuming does not grant new retry budgets.

The [earlier full benchmark](../golden_retest_20261101/README.md) supplies a
secondary regression comparison: 1,830/1,851 reactions recovered at tolerance 1,
one seed and cap 100. Its different cap and implementation timing must remain
explicit. The two new arms share cap 2,000; no third control run is requested.

Work, interpretation and manuscript text belong to the existing task-tracker
ticket 62. The current tolerance-1 benchmark remains the primary recorded score.

## Validated implementation and pilot

Scientific/runner revision: `4ba1939`. The Mac source/engine identity is
`7ccbb4b94beaefce637f6264a5c8941db143ae42747cb31cfd7d1199186f8ff3`;
the comparison runner SHA256 is
`02e386b82d90e30bab2bc5f83ad929c0dff4c3763f8e06856732bed07187c5ad`.
The explicit `event_problem` decoder input preserves matching constraints and
event matrices separately. Distinct weight sets export both identities in v6;
ordinary decoding retains v5.

Validation: 104 affected decoder/analysis tests, 34 campaign tests and 27 smoke
checks passed. After a final worker-lifetime adjustment, the affected two-arm
pipeline test passed again. Tests include an independent dense event oracle,
correlated choices, retained matching connectivity, persistence of both weight
sets, global stage barriers and preservation of failed-stage budgets on resume.

The pilot used cases 0, 1228, 1553 and 1739: eight directions per arm, with all
48 stage records successful. Both arms recovered cases 0, 1553 and 1739;
case 1228 remained unrecovered at one seed. Elapsed time was 144.53 seconds,
peak overlapping attempts was four, and maximum sampled worker RSS was
1.294 GiB. This selected pilot does not estimate full-dataset coverage or timing.
[Pilot measurements](pilot.json) retain the separate stage timings and settings.

The fresh full run uses four persistent workers on an Apple M4 Pro, each with an
8 GiB memory guard and a 1,800-second budget per direction/stage. The coordinator
has a 24-hour admission allowance. There are no automatic retries. Each
direction starts with its own 64 MiB completed-fragment cache, reused across
that direction's cuts; results and caches are not imported between the arms.
