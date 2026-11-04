# Decoder repair validation

The current repair implements the approved independent-fragment event model:
shared conditional fragment responses, fixed ordinary boundary events, one
sweep correction, and a minimum-plus-slack window (default six). Original
mapping constraints remain attached. The reported minimum and bins belong to
this approximation; explicit physical event evaluation remains separate.

Boundary lower bounds exclude compositions before local response preparation.
Only a fragment or the two sweep endpoint fragments can require a further
minimum query. Constant endpoint responses settle the sweep without that query.
There is no eager attainable-count pass, complete-bijection enumeration or
whole-branch correlation solver. Alternative cut histories remain separate.
Final catalogue provenance now associates each cut with its fragment partition;
this records existing search results without changing matching or pruning.

Local validation reused the saved case 1654 catalogues after the accepted
recursive pruning fix: one seed, tolerance 1, branch cap 100, full cut sweep.
Matching was not rerun. Both directions ran sequentially in one Mac process:

| Direction | Saved families | Retained compositions | Model minimum | Window | CPU seconds | Wall seconds |
|---|---:|---:|---:|---|---:|---:|
| Forward | 2,289 | 78 | 3 | 3–9 | 5.25 | 5.27 |
| Reverse | 2,488 | 1,181 | 3 | 3–9 | 11.17 | 11.18 |

Timing includes decoder preparation and record construction, excluding catalogue
loading and file publication. Peak process RSS across this sequential pair was
1.48 GiB. These are single point measurements, not a full Golden rerun or an
accuracy/performance comparison against the previous exact-boundary decoder.
No cluster benchmark was launched for this repair.

Validation: **979 tests passed on the Mac in 143.82 seconds**. Independent
oracles cover local joint constraints, both self rearrangements, action order,
fixed ordinary boundaries, strict WBO/metal thresholds, original sweep counting,
alternative cut histories, model minima, window exclusion before preparation,
shared-cache reuse and exact physical evaluation. Process-pool publication,
resume ownership and consumed-budget checks also pass. No search, hierarchy
pruning, native matcher or scheduling algorithm was changed.

## Withdrawn experiment and rollback

The user requested a complete rollback on 2026-11-04. The decoder rewrite in
`8641033`, including its follow-up optimizations, is withdrawn. Code, tests,
APIs, exports, native source and design documentation are restored to `0349a8c`,
the revision immediately before that rewrite. Its separate recursive branch-cap
repair remains intact. The new local count solver, approximate event bins and
minimum-plus-slack decoder API are removed; no alternate decoder mode remains.

The three-case cluster pilot used saved forward catalogues without repeating
matching. Cases 00001 and 01001 completed in 7.795 and 0.646 CPU seconds,
excluding publication. Case 01654, with 59,119 families, exceeded the explicit
300-second limit during computation; sampled peak memory was 1.13 GiB.
A separate bounded profile located substantial work in the added local minimum
queries. These measurements describe the withdrawn implementation only.

The prepared full Golden study (discovery tolerance 3, branch cap 200, one seed,
both directions and complete cut sweep) was never launched. Existing pilot
artifacts remain historical evidence and must not be resumed under the restored
decoder. No benchmark is started as part of this rollback. Operational evidence
and preserved uncommitted follow-up work remain private.

Rollback validation: the restored source, tests and documentation match
`0349a8c` exactly. The native module was rebuilt, and the restored API was checked
to exclude the withdrawn event-window argument and native helper. The full
local regression suite passes: **956 tests in 136.02 seconds**. Unrelated working
reports retain their original checksums. Rollback commit `b36baa2` is pushed to
main and deployed to the existing cluster checkout through Git. Its native
module was rebuilt there; **89 affected Linux tests pass in 8.25 seconds**.
No matching or decoding benchmark was launched.
