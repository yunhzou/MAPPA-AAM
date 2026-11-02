# Benchmark campaigns

`enzymemap.py` runs matching, per-reaction/direction deduplication, original-space
decoding and reference verification as separate passes through the public MAPPA
work pool. It contains no matching or decoding algorithm. Failed directions block
the next pass; caps and solver unknowns remain scientific result metadata.

The full EnzymeMap input contains 47,974 reactions in both directions. Its
manifest fixes the scientific settings: one-witness growth, one seed, branch cap
100, tolerance 1.0, random seed 42 and the full cut sweep. Verification asks for
exact annotated heavy-atom indices with hydrogens unconstrained, without
endpoint-symmetry relaxation or a claim of stereochemical validation.

## Prepare fresh inputs

`prepare_batches.py fresh OUTPUT INDEX... --workers COUNT...` creates disjoint
batches from existing prepared input manifests. Only immutable endpoint/reference
records are reused. Saved task lists locate inputs; their operations are never
executed. No results, attempts or cache snapshots enter the new campaign.

Both directions stay together. Batches are balanced by prepared-input bytes per
planned worker, with larger inputs first; this estimates workload rather than CPU
time. Each count is the **total workers in one batch's allocation**, not workers
per node. Duplicate or missing directions are rejected.

For example, four batches planned for 16 workers each:

```bash
export PYTHONPATH="$MAPPA_SOURCE/src:$MAPPA_SOURCE/bench"
"$MAPPA_PYTHON" bench/prepare_batches.py fresh "$MAPPA_BATCH_ROOT" \
  /path/to/prepared/index.json --workers 16 16 16 16
```

Run preparation from the checkout after loading the private environment below.
The output directory must not already exist. `prepare_batches.py resume INDEX
PREVIOUS OUTPUT` instead repartitions a stopped campaign while preserving its
results; it does not initialize a fresh benchmark.

## Private configuration and one portable launcher

`enzymemap.sbatch` serves both EnzymeMap and Golden. Site accounts, partitions,
hostnames, mount paths, node names and operational logs must stay outside Git.
Keep them in a private file such as `~/.config/mappa/slurm.env`, never in a
site-specific copy of the launcher. Review generated reports before publishing:
manifests and logs can contain input paths and allocation details.
The repository privacy check also inspects compressed evidence and notebook
text. Keep any site-specific denylist private:

```bash
python manuscript/scripts/publication_privacy.py --repository . \
  --private-values /path/to/private-identifiers.json
```

The optional denylist is a JSON array of exact private identifiers. The check
reports file names and categories without echoing private values. It checks the
tracked working tree, including staged new files; it does not clean Git history.

The private environment must export:

```bash
export MAPPA_SOURCE=/path/to/MAPPA-AAM
export MAPPA_PYTHON=/path/to/venv/bin/python
export MAPPA_BATCH_ROOT=/path/to/new-campaign
export MAPPA_WORKERS_PER_NODE=8
export MAPPA_EXPECTED_CODE_IDENTITY=validated-source-and-engine-identity
```

Set the site's `SBATCH_ACCOUNT`, `SBATCH_PARTITION` and any required `SBATCH_QOS`
in that same private file. If immutable input bundles have moved, also export
`MAPPA_BUNDLE_ROOTS` there as a JSON object mapping old roots to their new roots.
Relocation preserves offsets and checksums; it does not rewrite scientific data.

Build and validate the native engine before submission. Obtain `code_identity()`
from `mappa.execution.identity` using the chosen interpreter and checkout, then
pin that value in the private environment. The launcher rejects a changed source
or native engine. It never builds code inside an allocation. Keep that checkout
unchanged while its workers or queued continuations use it.

## Run matching and reduction before decoding

Load the private environment and submit two passes. This illustrative example uses two nodes,
eight compute workers per node, a 2 GiB worker guard and explicit admission/runtime
settings. CPU and memory requests include room for supervisors and the coordinator;
size that headroom for the machine and workload. The worker guard is not the
allocation's total memory request.

```bash
source ~/.config/mappa/slurm.env
cd "$MAPPA_SOURCE"
mkdir -p "$MAPPA_BATCH_ROOT/logs"

matching=$(sbatch --parsable --nodes=2 --cpus-per-task=10 --mem=20G \
  --time=02:00:00 --array=0-3%1 \
  --output="$MAPPA_BATCH_ROOT/logs/%A_%a.log" \
  bench/enzymemap.sbatch --memory-gib 2 --window 32 --run-seconds 1800 \
  --stop-after deduplication)

sbatch --nodes=2 --cpus-per-task=10 --mem=20G \
  --time=02:00:00 --array=0-3%1 --dependency="afterok:$matching" \
  --output="$MAPPA_BATCH_ROOT/logs/%A_%a.log" \
  bench/enzymemap.sbatch --memory-gib 2 --window 32 --run-seconds 1800 --resume
```

Allocation size, array scope/concurrency and runner resource flags are explicit
submission choices. The launcher forwards runner arguments unchanged. Each array
element owns `MAPPA_BATCH_ROOT/<array-index>/result` and starts one persistent pool;
no cut or seed submits a Slurm job. If matching is split across several arrays,
every finishing array must depend on **all** matching/reduction arrays. Thus no
decoding starts until the entire campaign has completed both preceding stages.

Workers prefer their current reaction, retain bounded caches and share reaction
snapshots; idle workers may take other ready cuts. `--window` bounds active
reaction/directions per allocation. Requested workers are not proof of execution:
check `progress.json`, task overlap and worker CPU usage for actual utilization.

`--run-seconds` or the advance Slurm signal stops admission and drains accepted
work. Allow time for draining within the allocation; each case has a one-hour
compute budget. A drained allocation exits **75**, so `afterok` keeps downstream
work blocked. Resume the same unfinished output explicitly with `--resume`, then
update downstream dependencies to the successful continuation jobs. Never run two
allocations against the same batch output at once. There are no automatic retries
or automatic Slurm resubmissions. Failures also block the next pass.

The coordinator publishes ledger/queue snapshots and progress every 30 seconds.
Completed directions in the published ledger are reused on resume; uncheckpointed
work or incomplete directions can repeat after termination. `plan.json` fixes the
denominator and input/configuration identity, `progress.json` reports stage counts
and workers, and `ledger.sqlite` indexes outputs and archived attempts. Aggregated
bundles preserve histories without per-cut or per-seed files. These outputs form
the handoff; no dataset-wide graph is collected in RAM.

## Golden and cache pilots

`golden.py INPUTS MANIFEST OUTPUT --batches COUNT` checks and aggregates the original
1,851 Golden inputs and prepares both directions. Its manifest selects complete
heavy-atom reference verification, including unmatched atoms and chemical endpoint
symmetries. Point `MAPPA_BATCH_ROOT` at that output and use the same launcher and
two-pass barrier above, with the matching array range. No separate Golden wrapper
or runtime policy override is needed.

`cache_pilot.py` compares worker-only and shared reaction caches on a small input
set, checks exact finalized graph equality, then runs reduction, decoding and
verification. Run it once inside an existing allocation using the same private
environment and single-thread BLAS settings; it creates its own persistent pool:

```bash
export PYTHONPATH="$MAPPA_SOURCE/src:$MAPPA_SOURCE/bench"
export MAPPA_NATIVE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
"$MAPPA_PYTHON" bench/cache_pilot.py --index /path/to/prepared/index.json \
  --output /path/to/new-pilot --count 10 --workers 8
```

The pilot uses a 4 GiB guard per compute worker and a one-hour case budget; reserve
additional allocation memory for its supervisor/coordinator processes. Repeated
`--index` arguments combine manifests. For a normal computer, add `--local` to
either Python runner; `--workers` then means total local workers. Build and test
before launching a pilot, using the checks in the repository's main README.
