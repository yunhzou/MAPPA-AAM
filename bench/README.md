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

Fresh full Golden runs use `rebenchmark.shared prepare/coordinator/worker` with
`run_conversion.sbatch` (seventh argument `rebenchmark.shared`). Preparation
verifies the original input manifest and packs the immutable search, finalization
and conversion task specifications into bundles. Seed batches and cuts share the
same executor; finalized graphs are persisted before conversion. A pilot selects
cuts from this exact index. Production imports verified completed pilot outputs
and charges their elapsed directional allowance, without repeating them.

Every task records CPU and wall time, sampled process-tree peak RSS separately
for computation and publication, the sampling interval, and its phase on failure.
Successful tasks also report publication CPU and the compute process's lifetime
high-water RSS (not a per-task peak). Conversion additionally logs loading, setup,
decoding, independent event-count checks and record construction. Node allocation
RSS/CPU comes from Slurm accounting and is kept separate from per-process data.

Fresh-run workers stop claiming tasks one directional allowance plus three
minutes before their Slurm allocation ends. Already claimed work keeps its
original deadline and publication guard. The retirement receipt permits a new
allocation to use the same queue without retrying accepted or active work.
Pilots from an earlier operational revision require an explicit accepted commit
list and recorded source diff; task specifications and input index must match.
