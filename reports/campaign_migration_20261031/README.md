# Campaign migration validation

Source revision: `56c4b6c`. The producing and destination checkouts were
synchronized through Git. Source writers stopped before immutable inputs and
saved outputs were transferred with checksum verification. Existing scientific
results were preserved.

Five disjoint batches partition reaction IDs modulo five: **47,974 reactions /
95,948 directions**, retaining **3,682 completed matching directions**. Bundle
byte ranges, checksums, denominator and scientific settings remained unchanged.
Storage-root relocation resolves the original references without changing data.

Validation passed **50 native-engine, scheduler and campaign tests** and **23
smoke tests**. Two real reactions completed matching, deduplication, decoding and
verification in both directions.

Decoding required successful completion of every matching batch. Incomplete or
drained matching work could not release the downstream pass. This report records
migration validation, not completion of the full benchmark. Raw operational
records and storage locations are retained privately; new campaigns use the
[portable launcher](../../bench/README.md).
