# 80-atom MAPPA case: branch cap 10, decoding threshold 1.0

The requested configuration finished in **1.49 seconds**, but retained **zero complete mappings**. All 91 cut checkpoints were examined independently: there were no terminal mappings and 275 recorded cap stops. Consequently the decoder had zero families to process. This is a resource-capped search outcome, not a successful mapping or a proof that no mapping exists.

| Measurement | Result |
|---|---:|
| Search, including cut checkpoint writes | 0.818 s |
| Catalogue processing | 0.629 s |
| Decoder call (empty catalogue) | 0.00172 s |
| Through decoding, including setup and catalogue saving | 1.476 s |
| Entire instrumented run, including audits | 1.490 s |
| Peak process resident memory | 57.53 MiB |
| Cut checkpoints | 91 |
| Recorded cap stops | 275 |
| Retained families / output candidates | 0 / 0 |

Search used three seeds, branch cap 10, matching tolerance 1.0, and graph/cut floors 0.15. The signed-event threshold was 1.0 for ordinary and metal-containing pairs in both search metadata and the decoder. The decoder used no event-count cutoff and a 120-second per-family budget. A single cluster CPU and the same 5 GiB virtual-address limit were used. Scheduler elapsed time was seven seconds, which includes shell/interpreter startup and job overhead; queue time is excluded from the measured computation.

A cap hit stops the affected branch; it does not guarantee that ten complete mappings survive. Lowering the cap therefore does not merely return fewer candidates from the cap-100 result. With no retained mappings, threshold 1.0 did not affect the output here. A matching tolerance and a decoding threshold also govern different operations: allowed fragment differences versus reported bond-order changes.

For context, the earlier cap-100 run retained 27,295 complete families and decoded 23,956 event classes. Its search took 139.38 seconds and its serial decoding took 771.21 seconds. Those earlier decoder thresholds were 0.5 for ordinary pairs and 0.3 for metal-containing pairs, so this is not a controlled timing comparison isolating branch cap alone.

## Reproducibility

- Cluster job: 7572965, completed with exit code 0.
- [Exact configuration and input hashes](result/configuration.json)
- [Timing and outcome](result/result.json)
- [Search counters](result/search.json)
- [Independent checkpoint audit and hashes](checkpoint-audit.json)
- [Runner](run_decode.py), [submission](run.sbatch)
- The frozen benchmark source is the parent report's `source.tar.gz` and `source_sha256.json`; this experiment does not change its engine.

The raw decoder marks an empty retained catalogue complete. The audit explicitly records that this does **not** mean a complete molecular mapping was found.
