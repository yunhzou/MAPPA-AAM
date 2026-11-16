> Correction (2026-09-13): the ten-seed Golden verdict is 1,840 recovered, 11 not recovered, and zero unresolved. Complete per-cut verification resolves cases 590 and 1358. See [the audit and corrected summary](../checkpoint_memory_fix_20260913/README.md). Original campaign JSON and execution records below remain historical evidence; the driver scripts now use bounded cut-by-cut verification.

# Final-source validation, September 2026

The algorithm source is commit 3a9ef9a. The search/evaluation source is unchanged from eb1a7d1; the final commit adds exact adjacent nested-subgroup absorption during final-family decoding. The source hashes and version proof record this boundary.

Golden includes all 1,851 original records, dataset revision 793475e. Search uses both directions, floor 0.2, iso tolerance 1, cap 100, root seed 42, uncut plus single-edge sweeps, and 1/2/3/10 seed orderings. Competition is disabled. Reference labels enter only verification. Recovered / not recovered / unknown counts are 1834/17/0, 1834/17/0, 1837/14/0, and 1840/9/2. The pinned SLAP sweep returns 1796/52/3. These are reference-family recoveries, not top-1 accuracy. The remaining ten-seed unknowns are cases 590 and 1358.

The coordinate evaluation uses all 140 original WBO inputs, forward one-seed search at cap 2,000, and the 128-completion competition policy. The final two-seed coordinate pipeline was not measured. Complete decoding covers all 236,653 saved flat families within fixed event windows, producing 336 reaction-specific classes. Coverage of returned SLAP-sweep minimum-event classes is 166/168, with every compared class recovered in 139/140 reactions. Native SLAP coverage is 155/160. This measures comparator-output coverage, not chemical accuracy. XYZ-derived comparator adjacency is restored and verified; XYZ preparation itself is not rerun. SLAP hydrogen label families are not exhaustively decoded.

Final decoding takes 25.70 recorded CPU-minutes and 11.60 elapsed minutes on three Mac workers, including continued passes and certificate journals, excluding preceding search and competition. The 1,821-case same-Mac timing cohort gives mean search/workflow CPU seconds of 2.518 (one seed), 4.953 (two seeds), and 4.819 (SLAP sweep). Reference verification is excluded. AAM search includes checkpoint IO; SLAP graph/mapping/export instrumentation excludes IO. Three- and ten-seed mixed-host times are not used for speedup claims.

The remaining campaign used distributed workers. All 6,706 work items were processed; the initial direction-level statuses include 72 memory limits and 9 timeouts. These are not counts of unresolved reactions because an opposite direction can verify recovery. Separate positive-only checks of saved cuts recovered ten of twelve remaining unknown reactions, without new search. Only verified positive membership is credited from partial archives; absence is never inferred from an interrupted search. Each subprocess retains a five-minute watchdog; continuations use at most 6 GiB sampled RSS.

The main release-to-summary duration was 42.96 minutes. Resource accounting
records about 47.40 CPU-hours, including multiple configurations, checks and
retries; this is not one configuration's successful search CPU. Completed task
results and raw cut checkpoints were preserved during scheduling corrections.
Raw operational records are retained privately.

Compressed records retain configurations, archive hashes, explicit verification witnesses, decoder patterns, initial resource limits, and subsequent cut proofs. All 140 coordinate audits and full-window decoder checks passed. Linux checks comprise 52 decoder tests, eight native tests, two Mac/Linux search/reference parity cases, and end-to-end AAM/SLAP/audit checks. The original Mac archives and large HPC raw checkpoints remain in their campaign directories, identified in the evidence; they are not duplicated in Git. Historical driver sources retain their original scientific meaning. Use the portable launcher and private site configuration for new runs; raw operational records remain private.

The manuscript contains the final algorithm and measured results. Operational details and development audits belong in this report. The coordinate data's original WBO provenance and redistribution terms still require author confirmation before publishing the underlying dataset.
