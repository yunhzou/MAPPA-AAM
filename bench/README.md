# Benchmarks

This directory contains the current Golden search, anchored-fragment conversion,
and saved-result analysis pipeline. Both reaction directions and seed prefixes
1–10 use the same public MAPPA implementation.

| Module | Responsibility |
|---|---|
| `prepare_golden_benchmark.py` | Pinned Golden dataset preparation |
| `golden_evaluation.py` | Reference mapping and strict scoring |
| `publication_analysis.py` | Representative ranking and coverage summaries |
| `publication_timing.py` | Search stage CPU and I/O measurements |
| `rebenchmark/campaign.py`, `search.py`, `shared.py` | Search planning and execution |
| `rebenchmark/analysis.py`, `saved_analysis.py` | Saved search evaluation and input indexing |
| `rebenchmark/decoding.py`, `decoder_campaign.py` | Public anchored conversion on saved AAM graphs |
| `rebenchmark/decoder_inventory.py` | Exact task accounting, including unresolved outputs |
| `rebenchmark/check_conversion_outputs.py` | Independent saved conversion checks |
| `rebenchmark/decoder_metrics.py`, `metrics_campaign.py` | Per-case decoded metrics and scheduling |
| `rebenchmark/report.py`, `timings.py` | Search/reference and timing aggregation |
| `rebenchmark/evidence.py`, `relocate_index.py` | Verified archives and explicit input mount relocation |
| `rebenchmark/run_conversion.sbatch` | Shared pilot/coordinator/worker launcher |
| `rebenchmark.sbatch` | Golden search drain launcher |
| `run_typed_ts_array.sbatch` | Explicit-reaction TS batch launcher |

Use `PYTHONPATH=src:bench`. Launches require explicit plans, source hashes and
new output directories. Saved results and paper provenance remain in `reports/`;
reproduce historical campaigns from their recorded Git revisions.

Do not rerun matching to inspect saved outputs. Preserve every expected task
key, direct anchors, fragment atom sets, incomplete flags, resource accounting
and publication hashes. Conversion minima describe retained converted branches.
Missing outputs and limits are unresolved, never negative results.

No retired solver campaign, experiment runner or benchmark-specific viewer is
kept here. Viewers use the shared MAPPA tools. Independent test fixtures live
under `tests/`.

Cluster launchers live here. Their partition, QoS, account and resource settings
must match the destination cluster; moving a launcher does not change its settings.
