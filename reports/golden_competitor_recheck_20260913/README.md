# Golden default-comparator recheck

All 12,957 saved method/reaction records from the seven completed default configurations were checked using the current strict heavy-atom relation evaluator. Mapping implementations and models were not rerun. Every per-case evaluation agrees with the archived result. The labeled audit and label-free input hashes match the final MAPPA/SLAP campaign exactly. Seven evaluator regression tests pass.

`competitors.json` supplies the paper table; compressed per-case files, `audit.json`, original version/model provenance, and the recheck drivers retain the evidence. Original predictions and failed attempts remain in `reports/golden_competitors_20260908/saved_outputs.tar.gz`, identified by SHA256. Signature evaluation preserves unmatched atoms and endpoint chemistry and handles agent-field spectators on the input side. Duplicate heavy labels are invalid.

The table distinguishes the first explicit mapping from recovery in any explicit returned mapping. It is not exhaustive SLAP symmetry-family decoding. MAPPA family recovery and the broader bidirectional SLAP sweep use separate search protocols. Original cluster mapping timings are retained in JSON but are not pooled with same-Mac MAPPA timing. RDT retries already completed in the original run; no final mapping failures remain for RDT.

Other archived tests are not conflated with these results. The coordinate-to-SMILES trials were cancelled and withdrawn; no scores from them are included. SAMMNet has no completed reproduction in this archive and receives no score. Superseded MAPPA experiments remain research history; the paper describes the final algorithm and current measured configurations.

The two Python drivers are preserved exactly as executed for audit provenance. They reference the original local staging directories; reproducing them requires adapting those paths and extracting the source archive. They are not installed command-line entry points. `audit.json` records the executed evaluator and driver hashes.
