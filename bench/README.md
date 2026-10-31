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

The first run uses one seed, branch cap 100, tolerance 1.0, random seed 42 and the
full cut sweep. Verification asks for exact annotated heavy-atom indices, with
hydrogens unconstrained and no endpoint-symmetry relaxation. It does not claim
stereochemical validation. The manifest records this policy explicitly.

The campaign uses a Slurm array of four independent two-node allocations, each
with 100 compute workers and 104 allocated CPUs. Each allocation starts one
persistent pool. Slurm may start batches independently as capacity is available;
The current wave requests 400 workers. The cpu-short per-user limit is four
submitted jobs; reaching 800 requires increasing each allocation to 200 workers.

`prepare_batches.py` partitions reaction IDs modulo the requested batch count (four for this wave) and copies only the
corresponding ledger rows from the stopped original campaign. Saved bundle
references remain valid; keep the original bundle directory. Each batch has its
own output directory and ledger, with no concurrent writers across batches.
Submit matching with `--stop-after matching`, then submit the remaining passes
with `afterok` on the entire matching array. A drained allocation exits 75 so it
cannot release the dependency before completing matching. Resume unfinished
batches explicitly before releasing downstream work.

At most 128 reaction/directions per batch are active. Completed histories and attempt
records are appended to bundles before their tasks are removed from the dispatch
queue. Task polling reads states without repeatedly unpickling results. This
keeps the active queue bounded while preserving audit evidence. There are no
per-cut or per-seed files. The coordinator owns local SQLite databases and
publishes ledger/queue snapshots and progress every 30 seconds.

The run stops admitting work after 3.5 hours and drains accepted work. Individual
case budgets are one hour, so an unusually long in-flight case can outlast the
remaining allocation. Slurm termination then leaves it incomplete. Reuse the
same output with `--resume` in a later explicit allocation; completed directions
in the last published ledger are reused. Work since that checkpoint or incomplete
directions can repeat. Failed directions are retained and stop transition to the
next pass; there are no automatic retries or automatic Slurm resubmissions.

`plan.json` fixes the complete denominator and input/configuration identity.
`progress.json` reports stage counts, ready workers and active work.
`ledger.sqlite` indexes outputs and archived task evidence. These files and the
bundle directory are the handoff; no dataset-wide graph is collected in RAM.

DFW uses `enzymemap-dfw.sbatch`: five independently scheduled two-node bins,
160 compute workers each (800 total), within cpu_short's ten-node limit.
The account is `nemotron_omni_vision`. The relocated campaign keeps aggregated
bundles intact and sets `MAPPA_BUNDLE_ROOTS` for original NRT references. Inputs
and results are verified by their original checksums when read. Matching runs
with `--stop-after matching`; downstream submissions depend on the entire array
completing successfully.
