# Fresh EnzymeMap benchmark

`enzymemap.py` runs matching, global per-reaction/direction deduplication,
original-space decoding, then exact heavy-atom reference queries as separate
passes. It calls the public MAPPA pool; no matching or decoding algorithm lives
here. A later pass never starts when the preceding pass has failed directions.
Search caps and solver unknowns remain scientific result metadata.

The input manifest covers 47,974 reactions in both directions. Only prepared
endpoint/reference records are read from the existing immutable input bundles.
Saved legacy task lists are used to locate endpoint records; their operations
are never executed. Old scientific outputs are not reused.

The current run uses one-witness growth, one seed, branch cap 100, tolerance 1.0,
random seed 42 and the full cut sweep. Verification asks for exact annotated heavy-atom indices, with
hydrogens unconstrained and no endpoint-symmetry relaxation. It does not claim
stereochemical validation. The manifest records this policy explicitly.

`prepare_batches.py fresh OUTPUT INDEX... --workers COUNT...` makes a completely
new campaign. It preserves immutable endpoint/reference inputs, explicitly uses
the current algorithm, and imports no matching results, ledgers or cache snapshots.
It keeps both directions together and balances prepared-input byte sizes per
planned worker. This is a workload estimate, not a measured CPU cost. Larger
inputs run first. The preparation rejects duplicate or missing directions.

Each independently scheduled allocation starts one persistent work pool. Workers
prefer their current reaction and retain bounded caches; idle workers can take
other available cuts. Shared snapshots remain enabled. Each batch owns its output
directory and ledger. Submit matching with `--stop-after matching`, then submit
the remaining passes with `--resume` and `afterok` on **all** matching jobs.
A drained allocation exits 75 and does not release downstream work. Resume an
unfinished batch explicitly; no task invokes Slurm.

`prepare_batches.py resume INDEX PREVIOUS OUTPUT` is the separate operation for
repartitioning a stopped campaign while preserving its results. It is not used
to initialize a fresh benchmark.

At most 128 reaction/directions per batch are active. Completed histories and attempt
records are appended to bundles before their tasks are removed from the dispatch
queue. Task polling reads states without repeatedly unpickling results. This
keeps the active queue bounded while preserving audit evidence. There are no
per-cut or per-seed files. The coordinator owns local SQLite databases and
publishes ledger/queue snapshots and progress every 30 seconds.

The DFW script stops admitting work after 2 hours 40 minutes and drains accepted
work within its 3-hour-45-minute allocation. Individual case budgets are one
hour. Slurm termination can still leave work incomplete. Reuse the
same output with `--resume` in a later explicit allocation; completed directions
in the last published ledger are reused. Work since that checkpoint or incomplete
directions can repeat. Failed directions are retained and stop transition to the
next pass; there are no automatic retries or automatic Slurm resubmissions.

`plan.json` fixes the complete denominator and input/configuration identity.
`progress.json` reports stage counts, ready workers and active work.
`ledger.sqlite` indexes outputs and archived task evidence. These files and the
bundle directory are the handoff; no dataset-wide graph is collected in RAM.

DFW uses `enzymemap-dfw.sbatch`, account `nemotron_omni_vision`, partition
`cpu_short`, and 2 GiB per worker. Set `MAPPA_BATCH_ROOT` and the validated
`MAPPA_EXPECTED_CODE_IDENTITY`; the script rejects a changed source/native build.
Override Slurm nodes, CPUs and memory together with `MAPPA_WORKERS_PER_NODE`
to fit available resources. The default is 100 workers across two nodes;
`MAPPA_ACTIVE_DIRECTIONS` controls the bounded admission window (default 128).
Requested workers are not evidence that they have started. Check Slurm and
`progress.json` for actual readiness, task overlap and worker CPU consumption.

The relocated campaign keeps aggregated bundles intact and sets
`MAPPA_BUNDLE_ROOTS` for original input references. Inputs and results are
verified by their original checksums when read. Old campaign outputs remain
untouched. Current DFW `cpu_short` limits the user to ten nodes across jobs.

## Golden with the same work pool

`golden.py` verifies and aggregates all 1,851 archived Golden inputs, prepares
endpoint-symmetry generators for verification, and partitions both directions
by reaction ID. It uses the same campaign runner and public MAPPA stages.
Golden's complete heavy-atom relation, including unmatched atoms and chemical
endpoint symmetries, is an explicit verification policy in its input manifest.
`golden-dfw.sbatch` requests five disjoint 100-worker pools at 2 GiB per worker.
Submit matching with `--stop-after matching`, then resume the ledgers with
`--resume` under an `afterok` dependency on the entire matching array.
Set `MAPPA_GOLDEN_ROOT` to the prepared DFW directory and, when relocating inputs,
`MAPPA_GOLDEN_INPUT_ROOT` to the original prepared directory.
