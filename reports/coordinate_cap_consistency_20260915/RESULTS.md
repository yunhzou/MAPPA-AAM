# Branch-cap comparison on 140 coordinate reactions

Complete paired evaluation.

Both configurations use one seed, the cut sweep, fragment competition, and the same event windows. Only the branch cap changes: 100 or 2,000. The cap applies to initial growth and competition. Search and event decoding remain separate.

| Result | Cap 100 | Cap 2,000 |
|---|---:|---:|
| Complete saved-family decoding windows | 140 | 140 |
| Reactions with a certified minimum in the saved full families | 138 | 140 |
| Reactions without a complete mapping | 2 | 0 |

Equal minimum counts on 138 reactions. Higher cap-100 minima: none; lower: none. Cap-100 searches without complete mappings: 123, 125. Incomplete paired decoding: none.

Cap-2,000 minimum patterns lost at cap 100: case 123: 1; case 125: 1. All-window pattern losses: case 123: 1; case 125: 1. Additional cap-100 window patterns: none.

Fresh cap-2,000 cases disagreeing with saved minimum counts or omitting a saved minimum pattern: none. Full-window differences from the saved cap-2,000 output: none.

## Comparison with the saved SLAP results

Mappings in this collection are unverified. These results compare event counts and patterns, not mapping accuracy. Minima refer to returned candidates, not global optima.

| MAPPA cap | Comparator | Fewer events | Equal | More | No full mapping | Unresolved | Comparator minimum patterns covered in window |
|---|---|---:|---:|---:|---:|---:|---:|
| 100 | native_slap | 15 | 123 | 0 | 2 | 0 | 153/160 |
| 100 | slap_sweep | 4 | 134 | 0 | 2 | 0 | 164/168 |
| 2000 | native_slap | 15 | 125 | 0 | 0 | 0 | 155/160 |
| 2000 | slap_sweep | 4 | 136 | 0 | 0 | 0 | 166/168 |

## Timing

CPU: AMD EPYC 9J14 96-Core Processor. Same completed-mapping cohort, one pinned CPU thread, both caps serial in shuffled order; worker CPU includes imports and persistence. Decoder includes all passes. Search/competition internal stage CPU excludes imports.

| Stage | Cap 100 mean CPU s/reaction | Cap 2,000 mean CPU s/reaction |
|---|---:|---:|
| Sweep search (n=138) | 1.209 | 2.158 |
| Fragment competition (n=138) | 1.389 | 2.189 |
| Decoding worker (n=138) | 12.216 | 19.461 |
| Total worker CPU (n=138) | 15.778 | 24.793 |

All processed worker CPU (including empty results and all recorded decoder passes): cap 100: 2180.294 seconds over 140 reactions; cap 2000: 3427.310 seconds over 140 reactions.

The worker total includes process imports, so it is not the sum of import-excluded search/competition and decoder rows. Timing uses identical reactions with completed mappings at both caps; failed searches and unfinished decoding are excluded from this paired latency summary. Raw stage records preserve their computational cost. SLAP times from a different CPU are not used to construct a speed ratio.

## Interpretation

Case 123 reaches the cap before obtaining a complete parent mapping at cap 100. Competition has no complete parent to revise, so it cannot repair this search. At cap 2,000 the initial search returns 73 complete terminals and the final pipeline recovers the saved three-event pattern. This is a search-budget failure, not an unresolved decoder or proof that the reaction has no mapping.

Decoding certificates cover the retained complete families within the common event window. Neither branch cap establishes exhaustive search. A lower cap may also change which competition candidates receive its fixed budget, so the decoded sets are compared directly rather than assumed to be nested.

Case 125 is a second confirmed empty search at cap 100; cap 2,000 returns its saved four-event solution. Both cases must count as search failures under a uniform cap-100 protocol.

Source and input hashes are frozen in manifest.json. The experimental driver is run.py; analyze.py interprets complete empty searches separately from unresolved decoding. Raw mappings and private execution details are excluded from this report.
