# Golden: ten seeds and cache scaling

This fresh study tests all **1,851 reactions / 3,702 directions** with ten
seeds, discovery tolerance 2 and branch cap 100. It starts with empty matching
results and caches, reusing only checksum-verified original endpoint/reference
inputs. The comparison is the [one-seed full run](../golden_tol2_cap100_20261105/README.md),
which recovered **1,834/1,851 reactions**. Source revision:
[`b51661d`](https://github.com/yunhzou/MAPPA-AAM/commit/b51661d766fc4549aece850a3e5f62e7c1ae230d).

**Status:** stopped at the user's request to repair observed overhead before
further benchmarking. All 3,702 directions completed matching and reduction.
Decoding produced 3,690 successful outputs and 12 memory-limit failures at the
unchanged 2 GiB guard. Reference verification published 3,701 of 3,702 checks;
the last check was interrupted. No retries were granted. The campaign is not
fully complete, and resource failures are not mapping absences. All 32 workers
registered in 6.05 seconds. The 384 cache-control runs and 192 paired scientific
checks are complete; their conclusions remain separate from this stopped run.

## Full-run protocol

| Setting | Value |
|---|---|
| Seeds / random seed | 10 / 42; stable prefixes of the original seed plan |
| Discovery tolerance / branch cap | 2 / 100 |
| Scope | Both directions; uncut plus full single-edge sweep |
| Input and event weights | Original prepared bond matrices |
| Graph/cut presence thresholds | 0.2 / 0.2 |
| Ordinary/metal event thresholds | 0.5 / 0.3 |
| Decoding window | Model minimum through minimum + 5 |
| Unique-fragment optimization allowance | Five seconds |
| Scheduling | Ten original seeds per cut task; cuts remain independently schedulable |
| Compute workers | 32 in one persistent pool; 40 allocated CPUs including overhead |
| Memory / compute allowance | 2 GiB per worker / 600 seconds per direction stage |
| Completed-fragment cache | 64 MiB per prepared worker session, with shared reaction snapshots |
| Expected work | 177,552 cut tasks; 1,775,520 seed searches |

Matching, reduction and decoding remain separate passes, with a global barrier
before decoding. Reference verification is external evaluation of the retained
catalogue. It does not establish reference membership inside the decoded event
window. Resource failures and unknowns remain in the full denominator.

The new logging retains each seed's original index, matcher CPU/wall time and
cache call/hit/dependency-check deltas. Small measurements accompany archived
task evidence, so analysis need not reload complete search graphs. Worker totals
also retain preparation, combination, publication, cache transport, sampled
memory and attempt intervals. Files are aggregated by worker, never per seed or
cut. Operational configuration stays private.

The execution change groups ten seeds on one worker for direct cache reuse and
avoids ten dispatches plus a cut-combination task. It does not change growth,
mapping permissions, pruning, seed orders or decoding. Grouped and separately
executed seeds were checked for exact graph equality. Three saved Golden point
cases also preserved their graphs against the preceding uninstrumented version;
raw timing repetitions are in [protocol.json](protocol.json).

## Controlled cache study

Eight reactions were selected before observing outcomes at explicit-atom size
quantiles 10%, 25%, 40%, 55%, 70%, 80%, 90% and 97%. The selection ranks maximum
endpoint atom count, then total atom count and case ID. Cases are **190, 1261,
895, 1137, 681, 1841, 1425 and 1779**, spanning maximum endpoint sizes **24–111
explicit atoms**. The full selection and input hashes are retained in the protocol
and private launch evidence.

Each case runs both directions at **1, 2, 5 and 10 seeds**, with the
completed-fragment cache enabled and disabled: **128 controlled direction runs**.
Every arm starts a fresh process and prepared session. Enabled retention is
64 MiB; disabled retention is zero. Other native preparation remains identical.
These controls measure the completed-fragment cache, not shared-file cache I/O.

The maintained `golden.py diagnose` runner measures matching, family reduction,
reference verification and audit costs separately. It compares complete catalogue
signatures, seed/cut scope and cap evidence between cache states, stopping on a
scientific mismatch. These controls do not repeat event decoding. The 2 GiB guard
and 600-second direction allowance remain explicit.

## Completed cache-control results

All **384 direction runs** completed. All **192 cache-on/off comparisons** agree
on complete catalogue signatures, seed/cut scope and cap evidence; those records
also agree across all three repetitions. Every directional reference was recovered.

| Seeds | Cache-on matching CPU (s) | Cache-off matching CPU (s) | Cache-on reduction CPU (s) |
|---:|---:|---:|---:|
| 1 | 3.258 | 5.169 | 1.250 |
| 2 | 6.315 | 10.308 | 2.633 |
| 5 | 15.679 | 26.904 | 6.538 |
| 10 | 31.958 | 55.665 | 13.746 |

Values are medians of totals across the 16 directions, measured three times.
At ten seeds, caching gives a **1.76× matching CPU speedup** (observed range
1.74–1.76×). One-to-ten-seed matching CPU increases **9.60×** with caching
(range 9.56–9.85×), versus **10.69×** without it (10.50–10.89×).
Including family reduction, the enabled ratio is **9.92×** (9.88–10.21×).
This shows useful absolute cache savings, but modest sublinear seed scaling.
The one-seed baseline already benefits from reuse across cuts.

All 66,376 ten-seed cache evictions in the first repetition came from case 1841;
the other seven cases fit within 64 MiB. A separate three-repetition capacity test completed all 144 arms under a 2 GiB
process guard, with every catalogue, scope and cap record equal to the original
64 MiB result. It changed no settings in the full benchmark. Each cold process constructs the
existing `SearchSession`, replaces its empty `FragmentRepair` with the requested
`cache_bytes`, then runs the same complete seed/cut plan and catalogue reduction.
The original compiled matcher is unchanged; no alternative algorithm is used.

| Cache limit (MiB) | Eight-case matching CPU (s) | Case 1841 CPU (s), both directions | Evictions per repetition |
|---:|---:|---:|---:|
| 64 | 31.892 | 13.023 | 66,376 |
| 256 | 31.205 | 12.402 | 2,070 |
| 1,024 | 31.266 | 12.371 | 0 |

Median matching CPU improves only about **2.2% at 256 MiB**. Eliminating every
eviction at 1 GiB gives no further aggregate improvement. Case 1841 improves
about 5%, but the larger cache adds only 546 hits across its 1.12 million
requests. Its maximum retained cache grows from 64 to 265.2 MiB. Most evicted
entries were not subsequently reused in these searches.

The original ten-seed sample's 92.5% request hit rate corresponds to **44.9% of
logical growth-extension steps reused**. Hits vary greatly in work avoided;
a high hit count does not imply that the same fraction of CPU is avoided.
[Capacity evidence](cache-capacity.json.gz) and its
[summary](cache-capacity-summary.json) retain every repetition and result check.

![Aggregate cache scaling](cache-scaling.svg)

![Scaling by selected case](cache-scaling-cases.svg)

Lines show medians and shaded bands the observed three-repetition range, not
confidence intervals. [Raw control evidence](cache-controls.json.gz) retains
phase timings, cache counts, memory samples, caps and scientific signatures.
The [summary](cache-summary.json) and PDF/PNG figure versions are included.
These controls exclude event decoding and shared-file cache transport.

The full-run allocation was reduced before starting to avoid an estimated
two-hour queue. Its replacement uses one 32-worker pool with identical scientific
settings and resource limits. The original owner was confirmed unstarted and
cancelled before replacement; no work was duplicated.

The two reactions with the largest full-run matching CPU totals, **1654 and
1606**, also completed all 12 targeted capacity arms with identical exact
per-cut catalogue, scope and cap signatures. These are separate exploratory
point tests, one measurement per arm, not part of the prespecified sample.
They retain only one cut's comparison catalogue at a time while the native
fragment cache persists across cuts.

| Case, both directions | 64 MiB CPU (s) | 256 MiB CPU (s) | 1 GiB CPU (s) |
|---|---:|---:|---:|
| 1654 | 101.978 | 101.157 | 97.435 |
| 1606 | 94.149 | 92.208 | 87.985 |

The 1 GiB limit reduces observed matching CPU by **4.5% and 6.5%**, respectively.
Some evictions remain at 1 GiB. This does not support a large capacity-driven
speedup even on the two slowest matching reactions. It does not establish that
64 MiB is optimal for every possible reaction. [Raw evidence](cache-capacity-heavy.json.gz)
and [summary](cache-capacity-heavy-summary.json) preserve timings and memory.

## Completed full-run matching and reduction

All **177,552 cuts / 1,775,520 seed searches** and all 3,702 directional reductions
are independently audited complete, without retries or duplicate task ownership.
Peak observed task overlap was 32. Matching used **15,482.407 CPU-seconds**;
reduction used **14,952.539 CPU-seconds**. Publication consumed a further
4,816.026 and 1,412.124 CPU-seconds, respectively, recorded separately.

Against the completed one-seed run, matching CPU is **7.00×**, reduction **13.30×**,
and their combined CPU **9.13×**, for ten times as many seed searches. This is an
observed campaign comparison, not an isolated cache effect: worker count and seed
batching also differ. Use the controlled cache study for causal cache comparisons.

Shared-cache transport recorded **26.449 worker-seconds loading** and **2,804.822
worker-seconds saving**, with 272,840,403 bytes loaded and 36,340,548,247 bytes
written across the run. These summed worker times are not elapsed campaign time.
There were no cache-load/save errors or baseline-recovery loads. Increasing
capacity can add snapshot cost, which the isolated capacity timings exclude.

## Reduction repair

The [baseline profile](reduction-profile.json.gz) separates saved-graph loading
from catalogue construction. In case 1654 forward, loading used 15.04 of 26.98
unprofiled CPU-seconds. Instrumented profiles locate repeated relation
construction checks during loading, plus repeated fragment hashing and action
preparation during catalogue ingestion. Profiling adds overhead; its absolute
times are not used in the comparison below.

Revision [`cb2a400`](https://github.com/yunhzou/MAPPA-AAM/commit/cb2a4005a33749f4c9f9e1bf433f37e1449a40a2)
restores trusted binary scientific fields without rerunning constructor checks,
reuses immutable fragment hashes and prepares each referenced graph transition
once per ingestion. Required-edge records are constructed only for new families.
Public construction still validates input. Fragment permissions, action ordering,
deduplication identity and every provenance entry are unchanged.

| Case, forward | Before CPU (s) | After median CPU (s) | Observed speedup |
|---|---:|---:|---:|
| 190 | 0.234 | 0.119 | 1.97× |
| 1841 | 10.598 | 8.293 | 1.28× |
| 1561 | 38.913 | 30.371 | 1.28× |
| 1654 | 26.977 | 8.899 | 3.03× |

All **12 complete-record comparisons** pass: four directions, three repetitions
each, retaining exactly the original branches, fragment constraints, action
order and provenance. Case 1561 contains 94,717 input paths, 91,763 retained
families and 6,631 unique fragments. This repair removes redundant preparation;
processing its large retained output still has a cost.

The baseline has one unprofiled observation per case; the repair reports three
repetitions. Both use the same CPU model in separate allocations. Searches were
not repeated. The unchanged saved graphs were republished to exercise the new
binary transport; that conversion and complete-record auditing are outside the
timed reduction. Existing old snapshots still invoke their original constructor
restoration until republished. Raw process memory includes conversion and audit,
not just reduction. Every validation retained the 2 GiB guard and completed.
These point timings do not establish a full-campaign speedup.

[Raw repair evidence](reduction-repair.json.gz) and
[summary](reduction-repair-summary.json) record the measurements. Validation also
passed 169 affected checks and 26 smoke checks; these are reported test scopes,
not a claim of exhaustive chemistry validation. The benchmark remains stopped;
this reduction repair does not resolve its 12 decode memory failures.

## Paper analysis

At the stop, saved verification establishes **1,843/1,851 reactions (99.57%)**
recovered in at least one direction, versus 1,834 with one seed: nine gains and
no losses. The interrupted check is 1351 reverse, whose forward direction already
recovered the reference, so it cannot change this reaction-level result. This
does not resolve the 12 decoding failures or establish recovery inside the event
window. The eight unrecovered reactions are 7, 19, 590, 986, 1285, 1358, 1377 and
1475. Directional checks remain 3,676 recovered, 25 not recovered and one unfinished.

The [stopped summary](stopped-summary.json) and
[direction records](stopped-directions.json.gz) retain phase CPU, publication,
memory and original seed measurements. Failed compute observations remain separate
from complete timings. The 44 distinct worker identities include replacements
after memory failures; they do not mean 44 concurrent workers. Pool capacity was
32 throughout.

Report reaction-level recovery, gains/losses from one seed, all failures/unknowns,
and the full denominator. Compare matching plus reduction separately from
decoding and external verification; retain CPU totals and per-reaction timing
distributions, including maxima. Report actual worker readiness and overlap,
not submitted capacity as utilization.

Plot CPU versus seed count for the size-stratified controls, with cache-on and
cache-off curves and a linear tenfold-work reference. Record cache hits, misses,
retained bytes, dependency checks and shared-snapshot I/O for the full run.
Compare only completed common scopes and show censored arms explicitly. Cache
hit rates alone do not establish a speedup, and scheduling/batching changes
prevent attributing the full-run wall-time difference solely to caching.

Validation before launch: **50 search/campaign tests**, **95
transport/cache/Golden tests**, **26 smoke tests** and **12 Linux seed tests**
passed. Cluster hosts, accounts, raw scheduler identifiers and paths are excluded
from public evidence.

## Reproducibility and interpretation

The full campaign archives every task and attempt, original seed/cut plan,
source/engine identity, stage handoff and reference answer. Per-seed measurements
are stored in the same aggregated task evidence as the matching results. The
public report contains sanitized summaries; original operational records remain
in the private campaign. A compute limit never becomes a negative mapping result.

For the cache control, each reaction/direction/count/cache arm starts a fresh
process, executes its entire cut sweep, reduces the resulting families and checks
reference membership. The reported matching CPU excludes preparation, family
reduction, signature auditing and reference verification; their separate times
remain in the raw records. The enabled and disabled arms share the same seeds,
original inputs, cap and compiled matcher. Their sole cache difference is
completed-fragment retention. Repeat measurements use the same fixed eight-case
sample; no case is added or removed based on speed or recovery.

The full run uses shared reaction snapshots and a persistent pool, whereas the
controls use a private in-memory cache within each fresh direction process.
Report these as distinct experiments. Full-run cache-hit counts alone do not
quantify how much CPU the snapshots saved.
