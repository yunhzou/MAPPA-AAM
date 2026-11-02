# Benchmark campaigns

`enzymemap.py` runs matching, per-reaction/direction deduplication, original-space
decoding and reference verification as separate passes through the public MAPPA
work pool. It contains no matching or decoding algorithm. Matching or reduction
failures block later passes. Once decoding is terminal, reference verification
continues from the catalogues in the same pool, even if decoding failed; those
failures remain recorded and prevent overall success. An explicit stop after
decoding still stops there. Caps and solver unknowns remain scientific metadata.

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

## Verify completed catalogues separately

Reference verification reads the reduced catalogue, not decoded events. After
matching and reduction finish, `--resume --verify-only` runs this external
evaluation even if decoding failed. Use the original index and output, unchanged
MAPPA source/engine identity and the same resource flags, for example inside an existing
allocation:

```bash
"$MAPPA_PYTHON" bench/enzymemap.py "$BATCH/index.json" "$BATCH/result" \
  --resume --verify-only --workers "$MAPPA_WORKERS_PER_NODE" \
  --memory-gib 2 --window 32 --run-seconds 1800
```

The portable launcher forwards `--verify-only` in the same way. The selected
batch must have no active or queued continuation owning its output. Other
disjoint batches may continue normally. Verification refuses a published queue
that still contains unarchived tasks or attempts, so their budget evidence cannot
be overwritten. Every selected direction must already have a completed catalogue;
input/configuration and saved MAPPA source/engine identity are checked before
submission. That fingerprint does not cover the benchmark verifier; keep its
implementation unchanged. The existing verifier and its per-family timeout are unchanged.

Completed verification records are reused. Failed verification records are also
left untouched: this option grants no retry or fresh budget to them. Decode
failures, consumed budgets and attempt audits remain recorded. A successful
verification-only pass reports `pass_finished: verification`; it does not claim
that decoding or the full benchmark finished. No decoding is retried.

## Recover decoding after a repair

After validating a decoder repair, retry a completed campaign's archived decoding
memory failures into a **separate output directory**:

```bash
python bench/enzymemap.py /path/to/original-campaign /path/to/decoding-recovery \
  --recover-decoding --local --workers 4 --memory-gib 2 --window 16
```

The original campaign root must contain the manifest and batch directories made
by `prepare_batches.py fresh`. All original stages must be terminal and their
queues drained. Recovery validates input/configuration identities, saved
catalogue references and original decoder thresholds. Its plan records the
original campaign identities and the current decoder identity. It restores only
archived decoding memory failures, with their original task IDs, attempts,
memory guards and remaining time budgets. It never reruns matching, reduction
or reference verification, grants no fresh time and leaves the original scientific
artifacts unchanged. Catalogue payloads are loaded only by compute workers.
Original batch ownership stays locked for the recovery pool's lifetime, preventing
concurrent source continuations or a second recovery of the same tasks.
A single `decoding-recovery.json` registration beside the original manifest binds
the campaign to this recovery output. A second output or fresh restart is rejected
so the original pre-retry budget cannot be reused.

Within an existing Slurm allocation, omit `--local` and set `--workers` per node.
This campaign-wide command uses one recovery pool; the batch-array launcher is
for the ordinary per-batch passes. Use the same worker guard as the original run.
After a normal drain, repeat the recovery command with `--resume`. Resume requires
the same recovery code/input plan and a drained saved queue; unarchived attempts
require explicit diagnosis. Terminal recovery rows, including failures, are
never retried automatically. Success reports `pass_finished: decoding`, without
rewriting the original benchmark's status or accuracy evidence.

If a later validated repair is needed, pass the **latest terminal recovery** as
the input and choose another separate output:

```bash
python bench/enzymemap.py /path/to/decoding-recovery /path/to/next-decoding-recovery \
  --recover-decoding --local --workers 4 --memory-gib 2 --window 16
```

This selects only its remaining memory failures, retaining every prior attempt
and the latest consumed budget. Successful directions stay in their existing
outputs. All predecessor scopes must be terminal, drained and free of active
owners. Their registered plans, ledgers, original inputs and scientific settings
are checked, and ownership remains locked through the new run. Each recovery
registers at most one successor; copied or older sources cannot restart a chain
or restore earlier budgets. A predecessor with a successor cannot resume. A new
decoder identity requires a new successor output; resuming that output still
requires its unchanged decoder identity and plan.

## Golden and cache pilots

For a bounded local study of selected Golden cases, use the prepared Golden
indexes and the current verifier directly:

```bash
export PYTHONPATH=src:bench
python bench/golden.py diagnose --index /path/to/prepared/index.json \
  --cases 1001 1850 1022 785 912 1518 709 1799 \
  --seeds 1 2 5 10 --cache on off --branch-limit 2000 \
  --direction-seconds 300 --run-seconds 7200 --memory-gib 4 \
  --output /path/to/new-local-diagnostic
```

Case IDs are explicit, zero-based selections; this command does not select the
whole dataset implicitly. Repeat `--index` when the chosen cases span several
prepared indexes. Both directions run by default; `--directions forward` or
`--directions reverse` selects one for an explicit follow-up in a new output
directory. Each selected case/direction must be present exactly once.
The example cases span matching-cost quantiles of the saved current Golden run;
they are a diagnostic sample, not a new benchmark denominator. For a study of
previous misses, supply those case IDs with `--seeds 10 --cache on` and the same
cap. A relocated input bundle can use the existing `MAPPA_BUNDLE_ROOTS` setting.
Pass `--cpu-model` with a verified hardware model when recording timing results.
The report includes platform information, core code identity and the diagnostic
runner's SHA-256; each task rejects a runner changed after planning.

Each reaction/direction and seed/cache setting starts a fresh process. The same
native `SearchSession` then retains its prepared state across the full canonical
cut sweep and all requested seeds. Seed-order prefixes, random settings and
input checksums remain fixed. `--cache off` sets completed-fragment retention to
zero; it preserves the native algorithm, dependency tracking and other caches.
The enabled control uses the existing 64 MiB completed-fragment cache. These
controls do not compare different search algorithms or shared-file cache I/O.

The runner performs matching, symmetry finalization, family reduction and the
existing Golden reference-membership query; it does not decode events. CPU and
wall time are reported separately for each stage, preparation and parity audits.
Task totals also include progress/result publication, while process startup is
reported separately. The existing local worker records sampled process-tree RSS
and the fresh process's lifetime peak RSS. Reports include cache calls/hits,
retention statistics, weighted cap/stop counts and completed cut/seed scope.
Cache controls compare the saved scientific family records and branch
associations, exact search scope, and per-cut graph/stop counts. This is catalogue
parity, not a claim of byte-identical search-history storage.
Two-cache studies automatically retain catalogue fingerprints. For separate
single-cache runs that will be compared later, pass `--audit-catalogue` to each
run. Single-cache studies otherwise skip that fingerprint calculation.

One `results.json` aggregates all requested settings; one atomically replaced
`progress.json` reports the current direction and stage. No cut or seed files
are created. Direction deadlines and the memory guard are enforced by the
existing worker supervisor. The run budget, SIGINT or SIGTERM stops further
admission; active work drains within its direction deadline and publication
allowance. Timeout, memory failure, unknown verification and nonrecovery remain
distinct. There are no automatic retries. Incomplete runs exit with status 75
and preserve completed records; use a new output directory for an explicit
follow-up. These times describe sequential local work, not distributed latency.

`golden.py INPUTS MANIFEST OUTPUT --batches COUNT` checks and aggregates the original
1,851 Golden inputs and prepares both directions. Its manifest selects complete
heavy-atom reference verification, including unmatched atoms and chemical endpoint
symmetries. Point `MAPPA_BATCH_ROOT` at that output and use the same launcher and
two-pass barrier above, with the matching array range. No separate Golden wrapper
or runtime policy override is needed.

### Minimum-event analysis of saved results

`analytics.py` queries completed decoder outputs without rerunning matching or
decoding. Supply an explicit JSON manifest with a `directions` list. Each row has
`dataset`, `case`, `reversed`, `decoded` and `reference`; the last two are saved
blob references (`path`, `offset`, `size`, `sha256`). Optional `source_provenance`
records the measured implementation and cohort. Select the latest successful
decoded output after recovery, retaining each direction exactly once.

```bash
PYTHONPATH=src:bench python bench/analytics.py /path/to/manifest.json \
  /path/to/new-analysis --workers 1 --seconds 30 --task-seconds 120 --memory-gib 2
```

The statistic is the annotated heavy-atom reference's bond-event count minus
the certified minimum over that direction's saved compressed mapping relation.
Both counts use heavy–heavy original WBO changes and the saved event thresholds;
all original matching constraints, including hydrogen constraints, remain in
force. Partial mappings follow the decoder's mapped/unmapped boundary rule.
Reports retain witness coverage so differing mapped atom sets remain visible.
This is not a minimum over all chemically possible mappings or both directions.

Saved local responses first provide conservative bounds. Exact feasibility
queries certify any remaining improvement, preserving joint symmetry choices.
Alternative fragment partitions are never added together. `--seconds` is a soft
per-direction analysis budget; input restoration is measured
separately. `--task-seconds` and the memory guard bound the complete worker task.
Unknown minima retain their bounds and never enter a numerical histogram.

One aggregate `results.json` records whole-input checksums, problem identities, resource
measurements and positive, zero or negative differences for certified directions.
`terminal_complete` means all tasks ended; `successful_analysis` means all tasks
published an analysis; `minimum_scope_complete` means every requested minimum
was certified. Unknown analysis and worker resource failures remain distinct.
Use a new output directory for each explicit run. Begin with a bounded pilot
before selecting a complete benchmark cohort.

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
