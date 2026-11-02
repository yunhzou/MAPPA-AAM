# Golden permissive matching comparison

Status: implementation and validation before launch. No full comparison result is available yet.

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
