# Withdrawn decoder experiment

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
