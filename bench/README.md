# Benchmark campaigns

## Reference-anchored completion

`anchored_completion.py` fixes the supplied reference pairs and completes the
remaining assignments through the existing MAPPA API. It uses one seed in each
direction, no cut sweep, tolerance 2, cap 100 and a minimum-plus-five event window.
This is reference-conditioned analysis, separate from benchmark recovery and
minimum queries over previously retained branches.

```bash
export PYTHONPATH=src:bench OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python bench/anchored_completion.py /path/to/manifest.json /path/to/new-output \
  --slurm --workers 8 --memory-gib 4 --task-seconds 180
```

The manifest's `directions` contain `dataset`, `case`, `reversed`, `prior_status`,
and checksummed `input` and `reference` bundle references. Input endpoints are
already oriented for that direction; reference pairs retain the original frame.
Include both directions exactly once per reaction. Omit `--slurm` locally.
The existing pool guards each direction and aggregates outputs. Results retain
the completed mapping, anchor checks, replay evidence, model and physical event
counts, minimum proof status, any added heavy-atom pairs, and separate matching,
full decoding and witness-audit timings. Empty or resource-limited results remain
explicit; they are not interpreted as impossible reference completions.

## Current Golden comparison with SLAP

`compare_slap.py` pairs current MAPPA (one seed, cap 100, tolerance 2) with pinned
upstream SLAP (`ea248fd`), binary plus weighted modes. Both consume the same
prepared explicit-hydrogen endpoint graphs, both directions and every single-edge
cut. SLAP retains its upstream isolated element-balancing atoms. No reference
labels enter matching. MAPPA retains its reaction fragment cache across cuts.

```bash
export PYTHONPATH=src:bench OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python bench/compare_slap.py /path/to/golden/index.json /path/to/new-comparison \
  --slurm --workers 16 --memory-gib 4 --seconds 600
```

Run in one allocation on one CPU model. Omit `--slurm` for local execution.
Use `--cases 0 7 19` for a pilot; omitting it selects all 1,851 reactions.
Install the pinned SLAP package and RDKit in the benchmark environment first.
Imports are outside scientific timings; graph preparation, matching and reduction
form search cost. Input loading, publication, MAPPA decoding (minimum + 5), and
external reference verification are measured separately. Both methods' coverage
uses complete heavy-atom correspondence modulo chemical endpoint symmetry.
MAPPA's additional decoded audit uses exact annotated indices and is separate.

The global matching barrier precedes decoding. One persistent pool alternates
method admission by reaction, records CPU model, source hashes, software versions,
budgets and per-direction timings, and publishes through aggregate bundles and
the existing campaign ledger. SLAP failures remain explicit incomplete comparator
outcomes; they are neither retried nor interpreted as mapping absence. Any MAPPA
failure stops admission and drains active work for investigation. `--resume` requires a drained owner and
identical scientific settings; completed and failed attempts are preserved without
automatic retries. `terminal_complete` means every required task has a recorded
outcome; `complete` additionally requires every task to succeed. Public evidence
must exclude the private operational paths.

`slap_audit.py` independently checks completed SLAP label records and reference
verdicts using a direct colored-graph VF2 oracle, rather than the comparison's
nauty certificates. It also times default binary, uncut `SlapAAM.get_maps` on
the identical prepared graphs for every reaction and direction. Graph preparation
and output extraction are timed; file loading and SMILES conversion are outside
this scope. This standard baseline is separate from the expanded two-mode sweep
above. It saves actual returned mappings and positive correspondence witnesses. Neither reference assignments nor atom-map hints enter the mapper.

```bash
python bench/slap_audit.py /path/to/completed-comparison /path/to/new-audit \
  --workers 16
```

Run inside one allocation with the same CPU model, pinned SLAP installation and
one library thread. The audit uses a 4 GiB guard and 600-second task watchdog;
resource failures stay explicit. It requires a drained source campaign, preserves
its results, and records its source hashes. Report standard forward-only and
bidirectional calls separately; neither is the full sweep protocol.

## Compare revisions locally

`performance.py` runs identical single-cut cases serially in separately built
source trees. Each process starts cold, then repeats with its own prepared
session and cache. The manifest is `{"cases": [{"name": "case-name", "input":
{"path": "/path/to/inputs.bundle", "offset": 0, "size": 123, "sha256": "..."},
"config": {...}, "cut": [[0, 1]]}]}`. Each verified input blob contains a plain
problem record accepted by `problem_from_record`, not a saved search object.
Use common configuration fields only; omit `self_tolerance`. An empty cut is `[]`.

```bash
python bench/performance.py --baseline /path/to/baseline --candidate /path/to/candidate \
  --cases /path/to/cases.json --output /path/to/new-comparison.json \
  --repetitions 3 --warm-repeats 1 --timeout 120
```

Use the same interpreter and separately built engines for both revisions.
`--case NAME` selects a case and can be repeated. Optional
`--baseline-revision` and `--candidate-revision` record commit provenance for
exported trees; actual source/engine hashes and loaded module paths are checked
independently. Worker environments select only that tree's `src` and use one
BLAS/OpenMP thread. Revision order alternates between repetitions.

The atomic JSON report retains raw timings and timeout/error output. CPU and
wall time cover planning, session preparation, complete cut matching (including
the revision's finalization), catalogue reduction, decoding, record export,
binary serialization and restoration. Input loading and imports are excluded.
Export materializes local response tables. RSS is the process lifetime peak,
including imports and inputs; it does not reset for warm repeats. Caps, terminal
counts, output sizes and cache counters remain explicit. A capped search with no
outputs is not a successful decoding benchmark, and equal settings alone do not
establish equal scientific work between revisions.

For diagnosis, `--profile-directory /path/to/profiles` writes one cProfile file
per worker. Profiled timings include profiling overhead and must be reported
separately. `--audit-output` hashes the complete decoded record outside the timed
stages to check exact output preservation. Run this audit separately from memory
measurements because hashing expands the record temporarily. Keep input paths
and raw reports in a private area when necessary.

For a postprocessing comparison, each case instead supplies a `matching` blob
containing a `MatchingResult`, a `reference` blob, `reversed`, the original
`config`, and `decode_config`. Both revisions read the same immutable matching
outputs; matching is not repeated. Loading, reduction, catalogue verification,
decoding, decoded-window reference queries and publication preparation are timed
separately. Decoding includes lazy window construction. The decoded reference
check uses exact annotated atom indices, independently of the catalogue's
endpoint-equivalence policy.

These comparisons stop on task errors, changed branch compositions, lost
reference recovery, or reduction/decoding CPU exceeding both 0.1 seconds and
1.25 times the paired baseline. Report smaller slowdowns too; this stop threshold
is a diagnostic guard, not proof that every case is faster. Saved-matching mode
checks decoded record transport equality and branch-composition hashes; it does
not claim complete decoded-output equality across revisions. The process timeout
is a wall-clock watchdog, not a memory allowance.

`enzymemap.py` runs matching, per-reaction/direction deduplication, original-space
decoding and reference verification as separate passes through the public MAPPA
work pool. It contains no matching or decoding algorithm. Matching or reduction
failures block later passes. Once decoding is terminal, reference verification
continues from the catalogues in the same pool, even if decoding failed; those
failures remain recorded and prevent overall success. An explicit stop after
decoding still stops there. Caps and solver unknowns remain scientific metadata.

`--seed-batch N` selects the number of original seed orders executed together
for each cut (default one). Use the full seed count to reuse the worker's cache
directly across seeds and avoid an extra cut-combination task. This changes
scheduling granularity, not seed selection or branch caps. Nondefault batching
is saved in the run plan and must match on resume. Task audits include per-seed
CPU/wall time and cache counters alongside total cut-task and publication costs.

The full EnzymeMap input contains 47,974 reactions in both directions. Its
manifest fixes the scientific settings: one-witness growth, one seed, branch cap
100, tolerance 1.0, random seed 42 and the full cut sweep. Verification asks for
exact annotated heavy-atom indices with hydrogens unconstrained, without
endpoint-symmetry relaxation or a claim of stereochemical validation.

`--reduce-only` accepts a fresh index whose directions contain `case`, `reversed`,
`matching` and original `catalogue` blob references, plus `config` and a
`reduction_source` provenance object. It runs three separate measured passes in
one pool: republish unchanged matching records into current binary transport,
reduce through `WorkPool.deduplicate`, and compare every catalogue's complete
scientific values and provenance. No matching, decoding or reference queries run.
Conversion and equality checking are excluded from reduction time. Outputs use
a new directory; original ledgers and bundles are read-only. Any failed task
stops further admission. Resume preserves completed passes, settings, source and
runner identity and requires a drained queue.

`--postprocess-only` uses the same saved-matching index, with each direction's
original reference blob and verification policy. It runs transport conversion,
reduction, decoding and external reference verification as separate measured
passes in one pool. Every reduction finishes before decoding starts. Original
catalogue blobs are neither required nor loaded: this mode validates a changed
reducer against references, rather than asserting equality with an old schema.
It records fresh budgets and source identity in separate output, preserves the
original campaign, and stops admission on any failed task. Saved matching times
and new postprocessing times must be reported separately; they do not establish
an uninterrupted end-to-end wall time.

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

The input index can declare `decode_config`, for example
`{"threshold": 0.5, "metal_threshold": 0.3, "event_slack": 5}`. Both a full run
and a decoding-only run use that policy. Preparation and execution record it in
the manifest and plan; recovery and output audits preserve it. Changing the
policy requires a fresh run rather than resuming an existing result directory.

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
work. Allow time for draining within the allocation; `--case-seconds` sets the
compute budget per fresh direction (default 3,600 seconds). A drained allocation exits **75**, so `afterok` keeps downstream
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

## Measure decoding again from saved catalogues

Use a fresh decoding-only campaign to measure a new decoder on every saved
catalogue, including directions whose previous decoding succeeded. Matching,
reduction and reference verification do not run again.

```bash
python bench/prepare_batches.py decoding /path/to/original-campaign /path/to/new-decoding \
  --workers 16 16 --memory-failures-workers 8
python bench/enzymemap.py /path/to/new-decoding/0/index.json /path/to/new-decoding/0/result \
  --decode-only --local --workers 16 --memory-gib 4 --case-seconds 1800
```

Preparation requires complete matching and reduction, drained source queues and
consistent scope/configuration. It records the original campaign identities and
catalogue references without importing old decoding results. The optional memory
failure worker list creates separate batches for directions with previous
decoding memory failures; it does not choose or increase memory automatically.
Set each allocation's memory explicitly. The portable Slurm launcher accepts
`--decode-only` with the same arguments.

The new plan records decoder identity, event thresholds and fresh time/memory
allowances. `--resume --decode-only` continues that unchanged measurement and
preserves its completed outputs and errors. Original results and consumed budgets
remain in the source campaign. Compare operation CPU on identical successful
directions, report publication separately, and keep failed or missing results in
the full denominator.

After both campaigns are terminal and drained, compare every saved decoded
record outside the timed stage:

```bash
python bench/performance.py audit-decoding --baseline-campaign /path/to/original-campaign \
  --candidate-campaign /path/to/new-decoding --output /path/to/output-audit \
  --workers 4 --memory-gib 8 --task-seconds 1800
```

Add `--slurm` inside an existing allocation. The audit streams complete record
and metadata hashes, releases each old artifact before loading its replacement,
and reports missing or failed outputs separately. `--resume` requires a cleanly
drained audit with unchanged source files, implementation and budgets.

## Recover final stages after a repair

Recover archived decoding or verification resource failures into a separate
output, preserving the original campaign and successful results:

```bash
python bench/enzymemap.py /path/to/original-campaign /path/to/decoding-recovery \
  --recover-stage decoding --local --workers 4 --memory-gib 4 --window 16 \
  --retry-seconds 600
```

Use `--recover-stage verification` for failed catalogue verification tasks.
The original campaign root must contain its manifest and batch directories.
All stages must be terminal and queues drained. Recovery checks input hashes,
scientific settings, catalogue/reference identities and task operations before
writing its plan. Compute workers alone load catalogue payloads.

The optional `--retry-seconds` explicitly grants that additional compute allowance
from the already consumed time, and permits the specified memory guard. Every
old attempt and resource grant remains in the task audit. Without this option,
recovery requires the original guard and uses only its remaining time. Core
scientific settings are always preserved.

After diagnosing and correcting an execution error, `--retry-error-reason TEXT`
explicitly selects that error for recovery and records the reason in the plan.
It preserves the failed attempt and remaining budget; it grants no extra time.
Without this option, recovery selects only memory and time limits.

Original batch ownership remains locked for the pool's lifetime. A stage-specific
registration binds each source to one recovery output, preventing duplicate
ownership or restoring older budgets. In an existing Slurm allocation omit
`--local`; `--workers` counts workers per node. Resume a normally drained allocation
with the same command and `--resume`. A nonempty saved queue needs explicit recovery,
not a fresh restart. Terminal failures are never retried automatically.

For another diagnosed resource retry, use the latest terminal recovery as input
and a new output directory. The runner validates its lineage, every prior attempt,
consumed time and resource grants, and selects only remaining failures. Successful
outputs stay in place. Changed code requires this new successor; resume requires
unchanged source identity and plan. Completion reports the selected stage pass,
not overall benchmark completion. Combine the original and successor ledgers for
the final scope and accuracy audit.

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
input checksums remain fixed. `--cache off` sets unique-fragment retention to
zero; it preserves the native algorithm and other prepared state.
The enabled control uses the 64 MiB unique-fragment pool. These
controls do not compare different search algorithms or shared-file cache I/O.

The runner performs matching, family reduction and the
existing Golden reference-membership query; it does not decode events. CPU and
wall time are reported separately for each stage, preparation and parity audits.
Task totals also include progress/result publication, while process startup is
reported separately. The existing local worker records sampled process-tree RSS
and the fresh process's lifetime peak RSS. Reports include existence queries,
fragment hits and witness checks,
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

### Paired Golden matching study

`golden.py compare` runs exactly two arms: binary matching weights at tolerance 1,
and original matching weights at tolerance 10. Binary weights are one where the
original WBO is at least 0.2 and zero otherwise, with a zero diagonal. Both arms
use one seed, cap 2,000, random seed 42, both directions and every single-edge cut.
Both decode events from the **original WBOs** at thresholds 0.5 and 0.3. Binary
matching endpoints and their constraints remain stored separately from the event
endpoints; replacing the catalogue's problem would change its scientific meaning.

```bash
PYTHONPATH=src:bench python bench/golden.py compare \
  --index /path/to/full-golden/index.json --output /path/to/new-paired-study \
  --workers 4 --memory-gib 8 --stage-seconds 1800 --run-seconds 86400
```

Omitting `--cases` requires all 1,851 cases exactly once in each direction: 3,702
directions per arm. For a separate pilot, add `--cases 0 1228 1739 1553`. The
memory guard applies to each worker; leave room for the coordinator and other
applications. The default admission window is twice the worker count.

One persistent local pool processes reaction/directions. Each matching task owns
its cache and immediately adds each completed cut’s compressed families to the
direction catalogue. Individual cut graphs are released rather than published.
All matching/reduction tasks in both arms become terminal before decoding starts.
The next pass decodes completed catalogues; the final pass independently verifies
the original heavy-atom references against matching catalogues. Failed matching
blocks only that direction's later stages. A decoding failure does not suppress
reference verification. Memory limits, timeouts, unknown reference queries and
completed nonrecoveries remain distinct.

Typed catalogues and decoded results go into aggregated bundles. Workers return
small descriptors with their blob references and measurements, so the coordinator
never loads molecular catalogues. Large-artifact publication is timed explicitly
inside the supervised stage budget. Matching, reduction, decoding,
verification and publication CPU/wall are separate; process-tree memory samples
and archived attempts retain resource failures. No whole-reaction mapping
enumeration or per-cut/per-seed files are introduced.

`plan.json` freezes both arms, input hashes, source and runner identities, guards
and per-stage budgets. `progress.json` and the ledger checkpoint every 30 seconds.
`summary.json` contains paired outcomes and the archived 1,830/1,851 baseline,
explicitly labelled as a different cap-100 study. `terminal_complete` means every
stage row is terminal, including blocked/error rows; `complete` additionally
requires every stage to succeed. Neither means every reference was recovered.

The run budget or a termination signal stops admission and drains accepted work.
Resume that unchanged output with `--resume`; completed and failed stage rows are
both retained, with no automatic retries or fresh failed-task budgets. Resume
requires a durably recorded drain and an empty archived queue. An unexpected
coordinator loss requires an explicit audit instead of silently repeating work.

### Minimum-event analysis of saved results

To plot where a GT-containing branch starts, use `analytics.py --positions`.
Both counts use the saved decoder model: the smallest minimum of retained
GT-containing contexts minus the overall detected minimum over both directions.
Ground truth identifies branch membership; the branch's event query is then
unconditioned on the reference. Thus a GT witness's event count is only an upper
bound, not automatically the branch minimum. Use exact annotated indices and the
saved decoded window for membership, consistent with the decoded audit.

Reuse existing verified witnesses at the overall minimum to establish offset
zero immediately. For other reactions, inspect saved context bounds and query
only branches that can improve the current value. The analysis changes no
matching/decoding output. Bounded unknowns retain offset intervals outside exact
numeric bins; missing references are reported separately.

The manifest has `source_files`, pending `directions` and already certified
`completed` rows. Pending directions carry decoded/reference blob references,
`overall_minimum`, optional `saved_minimum` and verified `certificates` identifying
the family, partition, cut and witness event count from the same decoded output.

```bash
PYTHONPATH=src:bench python bench/analytics.py /path/to/positions.json \
  /path/to/output --positions --slurm --workers 8 --seconds 30 \
  --task-seconds 120 --memory-gib 8
```

Render numeric output with `analytics.plot_positions(results_path, figure_stem)`
for PNG, PDF and SVG: integer offsets on x, logarithmic reaction counts on y.

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
set, checks exact completed graph equality, then runs reduction, decoding and
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

## Audit saved decoded recovery

Each reference first supplies a cheap lower bound on its forced model events.
Contexts above the saved window are rejected before preparing their symmetry
query. Unspecified atoms and sweep corrections remain conservative; a surviving
reference still needs the complete replayed witness. Audit rows record filter
time, examined/rejected contexts and symmetry-query counts.

```bash
PYTHONPATH=src:bench python bench/golden.py audit-decoding /path/to/batch \
  /path/to/separate-audit --workers 4 --seconds 10 --task-seconds 60 --memory-gib 2
```

The batch contains `index.json` and `result/ledger.sqlite`. Add `--cases 0 19`
for a pilot. This external audit loads existing successful decoded records; it
does not rematch, optimize starting correspondences or retry scientific stages.
Every recovery retains a full mapping and checked symmetry/placement choices in
one aggregate evidence bundle. It checks exact annotated heavy-atom indices,
including unmatched atoms; hydrogens must admit a valid injective completion.
This differs from endpoint-symmetry-normalized catalogue recovery.
The shared implementation is `decoded_reference.py`. It also accepts an explicit
`mappa.decoded_reference_inputs/v1` JSON manifest with `directions`,
`source_files` (absolute filenames and SHA256 digests) and `unmatched_heavy`.
Each direction contains `name`, `case`, `reversed`, `input`, `reference`,
`decoded` and optional `verification` blob descriptors (use null without a hint).
Resolve successful recovery outputs exactly once when preparing the manifest;
include the original and recovery ledgers in its source hashes. For EnzymeMap,
set `unmatched_heavy` to false: unannotated pairs remain unspecified, while the
complete witness must still satisfy saved decoder constraints and injectivity.

```bash
PYTHONPATH=src:bench python bench/decoded_reference.py /path/to/manifest.json \
  /path/to/separate-audit --slurm --workers 16 --memory-gib 4
```

`--slurm` uses one persistent pool inside an existing allocation. Workers archive
their attempts and publish replayable certificates in aggregated bundles. No
source stage is modified. Resource failures remain unknown in this independent
audit; an invalid certificate stops the pass for investigation.
The report distinguishes the saved model's event count from physical bond changes,
retains resource failures as unknown, and checks that source records are unchanged.
An invalid certificate stops admission and drains active audit tasks.
Saved catalogue certificates can prioritize a branch or supply one candidate;
they cannot establish decoded recovery on their own. A candidate is rewritten
using the saved starts only when the required binary self operation is valid,
then fully replayed and counted against the saved window. This checks one concrete
witness and cannot establish absence. Missing or failed catalogue checks provide no hint
and do not block this independent decoded audit. Their original records and
budgets remain unchanged. An explicitly budgeted additional pass can use
`--unknown-from /path/to/previous-audit/results.json` and a new output directory.
It checks input identity and admits only unresolved directions, retaining the
previous report hash. Report both passes and their combined resource use.
Use `--seconds 0` in the shared CLI to check existing candidate certificates
without further membership search. Candidates that cannot be validated remain
unknown; they are never counted as absent.
The audit shares `decoded_constraints.py` with minimum-choice analysis. It encodes
the saved ordered choices directly, constrains the annotated pairs, and tests
the retained event interval. Fixed reference pairs simplify their bond-event
terms before solving. This changes no saved start, symmetry, event count or
window. A returned mapping is replayed and counted independently.
The shared CLI exposes `--max-conflicts` separately from `--seconds`; these bound
the external constraint query, not matching or decoding. An exhausted solver
allowance remains unknown. Preserve earlier attempts when granting more time.
Catalogue hints with different endpoint-equivalent labels are skipped before
replay; the exact-reference query still runs.

## Replay fragment event responses

`python -m bench.fragment_events INPUT.json.gz OUTPUT.json.gz` measures saved
Fragment response preparation, bounded count queries and shared serialization.
The input schema and measured original/xTB fixtures are in
[the integration report](../reports/fragment_optimization_20261104/README.md).
`--reference /path/to/symmetry-event-lab` interleaves preparation with the accepted
beta; this optional benchmark dependency is outside production MAPPA.

## Choices at the detected minimum

```bash
PYTHONPATH=src:bench python bench/minimum_choices.py /path/to/saved/batch \
  /path/to/choice-counts --slurm --workers 16 --seconds 60 --memory-gib 4
```

This external analysis uses both saved decoded directions at their lowest
detected model event count. It counts distinct heavy-atom mapping classes using
the same chemical endpoint equivalence as benchmark verification. Hydrogen
assignments need one valid completion and do not multiply the number of choices.
Repeated branches, cuts and directions are combined. No scientific stage is
rerun and no ground-truth information selects the minimum or its alternatives.

Every enumerated heavy mapping is checked against the saved ordered actions and
event response. Exhaustion certifies the count; a time or work limit reports only
a lower bound. For a complete count K, uniform-choice reference success is 1/K
when the reference class is present, otherwise zero. This is expected random
selection accuracy, not a learned ranking or a guarantee of physical minimum
events. Use `--cases` for a pilot; original source hashes and evidence remain in
the separate output directory.

`--remaining-from /path/to/previous/results.json` preserves every exact count
and processes only incomplete or unprocessed reactions from the same immutable
scope. Each pass records its own counting allowance and the prior report hash;
scientific-stage results and their budgets remain unchanged. `--plot` exports
sanitized per-reaction evidence, summary statistics and the choice histogram.
Incomplete counts are excluded from exact-count means and contribute bounds,
rather than point estimates, to full-scope random-selection accuracy.
