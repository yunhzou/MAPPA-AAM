# Decoder repair validation

The decoder now follows the approved sequence explicitly: combine shared local
fragment minima with fixed boundary contributions to find the smallest lower
bound `L`; resolve the sweep only within `L + 2` to find the model minimum `M`;
then retain conditional responses through `M + n` (default `n=6`). The preceding
repair, `a3888ba`, ordered work by boundary bounds but did not impose the initial
`L + 2` window.

Cheap bounds exclude compositions before unnecessary local preparation. Numeric
fragment minima are shared. A bounded sweep query returns an attained count or
no result; its upper limit is never assumed to be a witness. The final reporting
pass reuses exact bounds already available or conservative local bounds, without
additional minimum queries. Its conditional count predicate remains attached to
the shared responses. Constant sweep responses still avoid a joint query.

The independent-fragment event model, fixed ordinary boundaries, original local
relations and separate physical evaluation are unchanged. There is no eager
attainable-count pass, complete-bijection enumeration, whole-branch correlation
solver or fallback outside the initial window. Alternative cut histories remain
separate. Matching and hierarchy pruning were not changed.

Local validation reused the saved case 1654 catalogues after the accepted
recursive pruning fix: one seed, tolerance 1, branch cap 100, full cut sweep.
Matching was not rerun. A fresh before/after measurement used one Mac process
per version, with forward then reverse running sequentially:

| Direction | Saved families | Retained compositions | Model minimum | Window | Before CPU seconds | After CPU seconds | After wall seconds |
|---|---:|---:|---:|---|---:|---:|---:|
| Forward | 2,289 | 78 | 3 | 3–9 | 5.11 | 5.24 | 5.27 |
| Reverse | 2,488 | 1,181 | 3 | 3–9 | 11.08 | 11.31 | 11.32 |

Timing includes decoder preparation and record construction, excluding catalogue
loading and file publication. Peak process RSS across the sequential pair was
1.47 GiB before and 1.39 GiB after. Both versions retained the same composition
counts and model minimum/window. These single measurements show similar runtime,
with a small observed increase; they do not establish a speedup or absence of
regression over the full dataset. No full Golden or cluster benchmark was run.

Validation: **984 tests passed on the Mac in 145.86 seconds**. New molecular
regressions distinguish the initial `L + 2` window from the larger reporting
window, retain its inclusive upper endpoint, require attainable counts, and
prevent minimum queries on reporting-only fragments. Existing independent
oracles cover local joint constraints, both self rearrangements, action order,
WBO/metal thresholds, sweep counting, cut histories, cache reuse and physical
evaluation. Process-pool publication, resume ownership and consumed-budget checks
also pass. A trial that eagerly minimized reporting-only fragments exceeded the
60-second point-test limit; that work was removed before the measurements above
and is now prohibited by a regression test.

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
