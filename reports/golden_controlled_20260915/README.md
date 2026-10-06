# Controlled Golden branch-cap ablation

Status: **complete**. All 1,851 reactions were processed for twelve configurations on the same CPU model, with one numerical-library thread per search. The [final comparison](final/COMPARISON.md) and manuscript Table A4 report cap 100 and cap 2,000 for MAPPA, with paired SLAP controls.

With the sweep, cap 2,000 confirms the same recovery at one, two and three seeds: 1,834, 1,834 and 1,837. Among cases resolved at both caps, no swept setting gains or loses a reference recovery. Ten seeds at cap 2,000 confirm 1,835 recoveries and leave 12 cases unresolved; five were recovered at cap 100. This is an unresolved verification difference, not a demonstrated loss of representable mappings. Without the sweep, cap 2,000 adds 19 confirmed recoveries.

The [raw numerical records](final/controlled-results.json.gz), [directional table](final/directional.csv), [paired audit](final/paired-cap-audit.json), and [runtime proof](final/runtime-proof.json) accompany the summary. Host identifiers are removed from the public numerical records. [Provenance](final/provenance.json) distinguishes original downloaded hashes from hashes of the published files.

## Comparison scope

All recovery percentages use all 1,851 reactions, retaining unresolved cases in the denominator. MAPPA uses both directions and is measured before fragment competition. SLAP unions binary and weighted modes in both directions. The common completed-search cohort contains the same 1,403 reactions for all twelve configurations; compare CPU means within that cohort. Search CPU includes graph preparation, mapping and output persistence, excluding initial imports and reference verification. This timing scope differs from the earlier Mac study, so those timings are not pooled.

The selected-attempt CPU column averages recorded expenditure over all 1,851 reactions, including interrupted calls. It is not completed-search latency or total campaign expenditure: 16 superseded checkpoint-conflict attempts and the initial pilot are excluded. The legacy CSV/JSON field names `all_attempted` refer to this selected-attempt measure.

The SLAP sweep confirms 1,794 references under this study's limits, versus 1,796 in the main evaluation. Cases 1665 and 1806 change from recovered to unknown; no other verdict changes. The cap-100 MAPPA counts reproduce the main evaluation. Resource-limited or inconclusive verification is reported as unknown, never as a verified exclusion.

## Reproduction

The frozen source is commit `1ba2d3706a214df2d394e8f2b941ecbfc76dbafa` plus the recorded uncut checkpoint-verifier fix. The existing manifest, payload and runners preserve the submitted protocol. Root seed 42, iso tolerance 1 and graph floor 0.2 are shared; MAPPA uses 1/2/3/10 seed orders per swept cut and one seed without the sweep. Search and reference verification have independent five-minute watchdogs and a 6 GiB worker RSS limit. See the runners for the full settings and adapt execution paths before running.
