# Branch-cap and unresolved-verdict audit

The live frontier cap is enforced per seed ordering and cut. The evidence does not show that cap being exceeded. The implementation retains the histories of all runs, including capped/dead states, combines the cut graphs, and then finalizes symmetry and evaluates the combined result. Consequently the live frontier cap is not a bound on the resident archive or verifier memory.

For Golden case 590, reverse, ten orderings and 56 cuts produced 52,738 terminals and 533,903 capped stops with a recorded maximum of 100 retained terminals per context. These terminals are raw search output, not final unordered-fragment deduplication counts. For case 1358, forward, 85 cuts produced 29,739 terminals and 473,410 capped stops. The native scheduler bounds its admitted frontier at 100, but commits proposed children to the trace before combined-frontier admission. Cap-rejected history therefore remains in the trace.

Source locations in the producing checkout:

- `native/src/fragment_scheduler.h:166`: frontier admission and cap checks; child records are committed before admission at lines 219–224.
- `src/rxn_core/aam.py:111`: independent seed-order graphs are combined per cut.
- `src/rxn_core/aam.py:242`: all cut graphs retained/restored before a combined graph is built at line 298. Both original graph containers and the combined containers are resident during subsequent finalization.
- `src/rxn_core/search_graph.py:423`: combination copies state/transition/stop wrappers and retains all histories.
- `src/rxn_core/search_symmetry.py:66`: finalization retains a Python cache per conditioned transition, independent of the bounded native coloring cache.
- `bench/golden_evaluation.py:225`: whole-graph ranking and representative/certificate collections; subsequent path verification retains deduplication keys and typed fragment caches. This is not enumeration of all automorphism-related bijections.

The original combined phases exceeded the external 6,144 MiB worker limit. Saved telemetry does not isolate an exact allocation or stack frame at termination, so the precise dominant allocation is not established. However the existing per-cut verifier completed all four directions using 200–304 MiB peak memory. This directly demonstrates that the reference verdict did not require the failed whole-graph working set.

## Incorrect unresolved aggregate

`reports/current_validation_20260912/drivers/score_checkpoints.py:76` assigns `unknown` whenever there is no positive witness, even if every expected cut was checked and returned `not_recovered`. This was intentionally a positive-only recovery adapter, but its result is overly conservative once full cut coverage is independently certified.

The audit recomputes the expected cuts from the frozen reaction inputs and cut floor 0.2, verifies exactly the expected cut indices, all negative outcomes, successful per-cut worker exits, driver hashes, and the exact-graph finalization proof hash. All 284 expected cut checks are present:

| Case | Direction | Expected / checked | Peak MiB |
|---|---|---:|---:|
| 590 | Forward | 58 / 58 | 303.7 |
| 590 | Reverse | 56 / 56 | 278.5 |
| 1358 | Forward | 85 / 85 | 240.4 |
| 1358 | Reverse | 85 / 85 | 199.6 |

All cut verdicts are negative. Since combining independent cuts is a disjoint union and finalization is cut-local, both reaction-level reference verdicts are **not recovered**, not unknown. This certifies absence from saved capped families, not all possible atom mappings. The failed monolithic phases remain failed execution records.

Corrected ten-seed Golden counts: **1,840 recovered, 11 not recovered, 0 unknown, denominator 1,851**. Recovery percentage remains 99.41%. The paper and the corrected seed_comparison.json in this report use these counts. The original campaign summary remains unchanged as historical evidence.

## Implementation direction

Preserve the cap and search semantics. Verify finalized cuts sequentially and aggregate their verdicts. Allow a negative certificate only if the manifest matches, every scheduled cut is present and checked, and every cut returned a definitive negative verdict. Missing cuts or solver timeouts must remain unknown. Keep detailed historical traces on disk rather than requiring them all in memory for a coverage query; retain complete ancestry for every reported family. Do not infer a global search guarantee or expand automorphism groups into explicit bijections.

Run `audit.py` with the existing `work/aam-event-env/bin/python`. This investigation performs no new search or decoding. It audits saved results and does not re-read the remote raw archives. `audit.json` records the source evidence hashes; guard checks reject missing/duplicate cuts and unknown outcomes.

## Implemented fix and validation

`search_aam_checkpoints` persists raw cuts without assembling the global history. It returns checkpoint paths, aggregate search metrics, and capped status. It preserves the existing search semantics and full historical records on disk. Its working set includes the current cut per worker; one unusually large cut can still require substantial memory. The original `search_aam` API continues to return a materialized graph for callers that need it.

`bench/golden_checkpoint_evaluation.py` validates the input/configuration manifest and expected cut indices, independently finalizes each cut with a fresh symmetry workspace, checks reference membership, and releases the cut. Positive witnesses retain their cut identity; their rankings are explicitly local to that cut. Definitive absence requires every scheduled cut to be checked with a negative verdict. The Golden validation runner and both fallback adapters now use this path and accept complete negative certificates. Coordinate search behavior is unchanged.

The tests cover exact graph equality with the existing search on a real molecule, serial and parallel checkpoint parity, release of previous cut graphs, resume without new matching or global combination, manifest/context errors, missing/duplicate cuts, unknown verdicts, and all 284 recorded negative cut verdicts for the two ten-seed cases.

Fresh saved-archive verification using the new code, one seed, one CPU, both directions of cases 590 and 1358: 27.69 seconds combined verifier wall time, 25.83 CPU seconds, peak 162.67 MiB. A 3 GiB / 300-second external guard was applied per process. All four verdicts remained not_recovered. See fresh-seed1-case*.json. No full benchmark searches were rerun; the ten-seed correction uses the independently audited original per-cut records.

For new Golden campaigns, deploy the updated src and bench modules together with the updated runner scripts. Existing frozen campaign archives and original execution records are not rewritten. The standalone saved-archive integration script accepts the existing campaign directory as an argument.

Raw operational records and full artifacts are retained privately.
