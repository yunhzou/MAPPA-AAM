# Prior-method baselines in the paper

LocalMapper is labeled as a historical accuracy SOTA (2024), supported by Chen et al., Nature Communications, DOI 10.1038/s41467-024-46364-y. Its reported 89.8% used 1,758 retained Golden reactions and manual calibration; this is explicitly separate from our strict 1,605/1,851 (86.71%) re-evaluation. SLAP sweep is labeled the strongest evaluated prior-method configuration for returned-alternative coverage, not a published SLAP SOTA score; our bidirectional/cut wrapper is disclosed.

Figure 2 includes LocalMapper and SLAP prior baselines. Table 3 marks both roles. Table 4 includes mean-CPU ratios relative to the same-direction SLAP baseline on the same 1,821 Mac reactions. Table A3 pairs strict recovery with archived default call costs and marks LocalMapper; RXNMapper is the lowest observed mean call-time comparator. Those cluster timings are not used for a cross-host MAPPA speed ranking.

Figure 3 includes MAPPA, native SLAP, and SLAP-sweep minimum-event distributions on all 140 coordinate reactions. Tables 5–6 distinguish the prior mapping baseline from full-family decoding; no equivalent prior full-family decoder timing is available. These mappings remain unverified, and the comparison is not mapping accuracy.

All numerical values derive from existing saved evidence. No searches, decoding runs, or mapper calls were rerun. The literature source was checked on 2026-09-14. The new script assertions verify the strongest measured single-output baseline and all timing ratios.
