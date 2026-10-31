# DFW campaign migration

Source revision: `56c4b6c`. Local canonical checkout and the DFW checkout were
synchronized through Git bundles. NRT array 195448 was stopped before transfer.
Immutable inputs and stopped outputs were copied with checksum-verified rsync;
no private keys were copied. Temporary forwarding agents were terminated.

The five DFW bins partition reaction IDs modulo 5: 95,948 distinct directions
across 47,974 reactions, with 3,682 completed matching directions preserved.
Original bundle byte ranges and checksums remain unchanged; storage-root mapping
resolves NRT references on DFW. The original denominator and scientific settings
are preserved. There is no overlap between bins.

DFW validation on CPU nodes:
- Job 21087919: 50 native-engine, scheduler and campaign tests passed.
- Job 21087969: 23 smoke tests passed; two real reactions completed matching,
  deduplication, decoding and verification in both directions.

Matching array: 21088167. Five bins each request two nodes, 84 CPUs per node,
80 compute workers per node and 180 GiB per node: 800 compute workers and
840 allocated CPUs total. Partition cpu_short, account nemotron_omni_vision.
Downstream array: 21088172, dependent on successful completion of the entire
matching array. Incomplete/time-drained matching jobs do not release decoding.

The user authorizes scaling up to 2,000 CPUs when eligible capacity permits.
Current DFW cpu_short limits the user to ten nodes; cpu_dss rejects this account.
Requested workers are not evidence of actual concurrency; consult Slurm and the
per-bin progress checkpoints for the current running count.

DFW root:
`/lustre/fsw/portfolios/nemotron/projects/nemotron_omni_vision/users/yunhengz/mappa`
Campaign: `enzymemap_pool_20261030/dfw160/{0..4}/result` beneath that root.
