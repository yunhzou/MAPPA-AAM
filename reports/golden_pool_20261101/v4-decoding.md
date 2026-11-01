# Completed Golden v4 decoding

All **1,851 reactions / 3,702 directions** have current `mappa.decoded_events/v4`
outputs. Job array **21137977** decoded 3,701 saved catalogues and reused the
validated 01561-forward v4 result. Every direction has exactly one successful
decoding attempt in the final set. There were no failures or memory retries.

The pass used source `2f009516016efe3e1dbf97079908ac2f712e0327` and the existing
campaign runner. New ledgers retain the original matching, deduplication and
verification rows exactly. Those stages were not recomputed. Old decoded
artifacts remain intact in `golden_pool_20261101`; the current set lives in
`golden_decode_v4_20261101` under the DFW MAPPA root.

Five disjoint 100-worker pools used the existing DFW `cpu_short` allocation shape:
four nodes, 112 allocated CPUs, and 54 GiB per node. Each compute worker retained
its **2 GiB guard**. Two allocations could overlap under the node quota. Each
pool used all 100 workers and reached 100 overlapping tasks; these concurrency
measurements do not imply continuous full CPU utilization.

| Measurement | Result |
|---|---:|
| First allocation start to last allocation finish, including queue gaps | 3 min 45 s |
| Individual allocation elapsed time | 44–79 s |
| New decoding CPU time | 7,974.99 s |
| New publication CPU time | 272.60 s |
| Peak worker RSS | 1,102,872,576 bytes / 1.027 GiB |
| Complete decoded outputs | 3,702 |
| Compressed result bytes | 815,213,754 |
| Result bundle files, including reused point result | 501 |

Average worker CPU consumption over each pool's task span was 25.7–47.7 cores.
This report measures the completed pass; it makes no new scaling or comparative
speedup claim. The runtime includes waiting between allocations, and excludes
the earlier point test, development, preparation and subsequent integrity audit.

The final audit checked the exact denominator and disjoint reaction ownership,
identical preserved stage rows, decoding-only task operations, original catalogue
references, event thresholds, successful single attempts, and unchanged guards.
It verified SHA256 for every published result byte range without loading the
dataset into RAM. The combined manifest has SHA256:

`746bf4a86074d6e604a4233d931ff4e566cd1052fd06cb3ce1cf3aaed965e391`

The manifest references existing DFW bundles; retain those files and the reused
point-result bundle. All five Slurm allocations exited successfully. Existing
reference coverage remains **1,830 / 1,851 reactions (98.87%)**, because the
verification results and underlying catalogues are unchanged.

[Launch](v4-decoding-launch.json) · [audit and measurements](v4-decoding-completion.json)
· [combined decoded manifest](v4-decoded.json) · [decoder validation](decoder-refactor.md)
