# 80-atom MAPPA result: branch cap 100, decoding threshold 1.0

The 80-atom PMC11894599_10 case completed in **357.34 seconds (5 min 57 s)** on one cluster CPU, including checkpointing, catalogue preparation, output saving and witness audits. All **27,295 retained families** were decoded without an event-count cutoff. The run emitted **342 signed-event classes**, and every emitted mapping passed independent WBO rescoring and producing-family certificate replay.

| Phase | Measured time |
|---|---:|
| Search, including 91 cut-checkpoint writes | 132.06 s |
| Catalogue processing | 62.09 s |
| Catalogue serialization and other preparation before decoding | 68.28 s |
| Full decoding, including per-family report saving | 93.55 s |
| Through decoding | 355.98 s |
| Independent witness audit | 0.81 s |
| Entire run, including remaining output writes | 357.34 s |

Peak process resident memory was **1,783.84 MiB (1.742 GiB)**. The Linux virtual-address-space limit remained **5 GiB**; resident memory and virtual address space are different measurements. There was no memory-limit failure.

## Mapping result and scope

All 27,295 family reports were complete, with zero incomplete paths. The minimum was **6 signed bond-order changes** within the retained cap-100 families. Three event classes tied at that minimum; each contained four WBO decreases and two WBO increases meeting the threshold of 1.0. This is not a global minimum over branches omitted by the search cap, nor a stereochemical correctness verdict.

A threshold of 1.0 ignores smaller bond-order differences. Therefore this minimum of 6 is not directly comparable to the previous minimum of 24 measured at ordinary/metal thresholds of 0.5/0.3; it does not establish that the underlying chemical mapping improved.

## Controlled comparison with cap 10

| Setting | Cap 10 | Cap 100 |
|---|---:|---:|
| Matching tolerance | 1.0 | 1.0 |
| Ordinary and metal decoding thresholds | 1.0 / 1.0 | 1.0 / 1.0 |
| Completed retained families | 0 | 27,295 |
| Output event classes | 0 | 342 |
| Search | 0.82 s | 132.06 s |
| Decoder call | 0.00172 s (empty input) | 93.55 s |
| Entire instrumented run | 1.49 s | 357.34 s |

Only the branch cap changed in the runner. Both runs used the same unchanged inputs and frozen engine, three seeds, all 91 cut checkpoints, graph/cut floors 0.15, one CPU, the same memory limit, and the same unrestricted decoder with a 120-second per-family budget. Cap 10 stopped all complete paths; its short time is not a successful faster mapping.

## Evidence

- [Measured result](result.json) and [exact configuration/input hashes](configuration.json).
- [Cluster audit](independent-audit.json): all 342 witnesses rescored and checked against constructive family certificates.
- [Local validation](local-validation.json): repeated WBO rescoring of all 342 witnesses; verified all 27,295 complete family reports, all 91 cut checkpoints, input hashes and archive checksum.
- [Full saved evidence](evidence.tar.gz) and [checksum](evidence.sha256).
- [Runner](run_decode.py), [submission](run.sbatch), [scheduler completion](scheduler-final.txt), cluster job **7573145**.
- The exact frozen engine is the parent report's `source.tar.gz` and `source_sha256.json`.
