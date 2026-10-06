# Reuse validation and conditioned symmetry across competition branches

Repeated exact calculations are now shared across repairs and during event decoding. On the expensive case 11, measured competition-plus-decoding CPU fell from **57.223 to 32.186 seconds** for prefix competition and from **150.330 to 89.211 seconds** for maximal continuation. Search graphs, witnesses, symmetry generators and decoded classes matched the saved outputs exactly.

This is calculation reuse, not a new equivalence relation between mapping families. The large maximal-continuation result still contains 31,948 families and 31,909 different representative mappings. Those families cannot all be discarded as repeated copies. Event-class equality alone is insufficient for merging partial search states: their allowed future extensions and symmetry actions can differ.

## What changed

- One reaction-local conditioned-symmetry workspace is shared across competition repairs. Exact locked roles and candidate partitions remain part of the calculation. Its native cache has a 64 MiB budget.
- Representative validation propagates preservation-constraint unions and path counts through the search DAG. All histories ending at a state share its witness. If that witness satisfies the union, every incoming history is certified with one check. If it fails, the original pathwise validator determines precisely which histories fail.
- Fragment constraint masks and constraint arrays are cached across repairs and event-extraction calls. Each cache is limited to 4,096 entries. Event thresholds and output grouping are unchanged; decoding remains separate from search.

No bijection enumeration, event-based search pruning, new takeover policy, or changes to the completion budget were added. Existing exact state/completion deduplication remains in place. Distinct histories are retained for symmetry queries and audit. The optimization also benefits the default local competition policy; prefix and maximal continuation remain experimental ablations.

## Measured performance

Times below are process CPU seconds on the same Mac, measured by the same driver. Search includes loading the baseline and competition, excluding repair persistence. Decode includes catalogue construction, symmetry reconstruction, persistence, validation and event extraction. The initial saved sweep search is reused, so these are not end-to-end fresh reaction timings.

| Case 11 policy | Competition before → after | Decode before → after | Combined before → after | Reduction |
|---|---:|---:|---:|---:|
| Prefix | 26.604 → 14.100 | 30.619 → 18.086 | 57.223 → 32.186 | 43.8% |
| Maximal continuation | 66.145 → 45.862 | 84.185 → 43.349 | 150.330 → 89.211 | 40.7% |

These are single measured runs against saved baselines, not repeated statistical timing estimates. Small cases show smaller and occasionally noisy changes; `summary.json` lists every measured pair. Do not extrapolate these percentages to all reactions or the default policy. Four workers ran concurrently, each with one native numerical thread, a 300-second watchdog and 3,000 MiB memory guard. All 18 stages passed. Maximum recorded worker RSS was 2,456.6 MiB in the large case-11 search.

The prefix catalogue retains 12,336 families and three event classes; maximal continuation retains 31,948 families and one event class within the frozen case-11 window of at most five events. Every family was fully decoded within that window. This does not claim enumeration of every event count or every globally possible mapping.

## Correctness evidence

Eleven search comparisons cover local competition on cases 6, 11, 59, 64 and 101; prefix competition on those same five cases; and maximal continuation on case 11. Every accepted repair graph has the same canonical JSON SHA-256 as before, including its symmetry records. Existing counters, accepted offers, repair configurations and pending-work counts agree. Seven paired complete decodes have identical full event records, including witnesses/actions, catalogue sizes and class IDs.

`verification.json` contains all per-repair hashes. `results.json.gz` contains the before/after search and decode records. `execution.json` records watchdog outcomes and memory. `production-check.json` additionally checks the integrated default policy on case 64 against its saved local-mode control; diagnostic mode labels are omitted when comparing offers.

The production source passes 66 focused tests covering representative validation, native conditioned reuse, event extraction, final-family grouping and family scoring. The new tests check a shared DAG join, mixed valid/invalid histories, cuts, deferred edges, tolerance boundaries, wrong elements and duplicate target assignments. Initial broad test collection lacked the existing `rebuild_intrinsic` helper path; adding its benchmark directory resolved all five import failures.

## Reproduction

Use a separate source checkout at **807d36f**, install repository dependencies plus pybind11, setuptools, psutil and z3-solver, and invoke the scripts from this report directory. The adjacent prefix report supplies the isolated experimental engine; the adjacent unrestricted report supplies the five saved input archives. No working source is modified.

```sh
export MAPPA_EXPERIMENT_WORK=/tmp/mappa-reuse-before
python prepare.py --repo /path/to/source-at-807d36f --reference
python bench_reuse.py
export MAPPA_REFERENCE_WORK="$MAPPA_EXPERIMENT_WORK"
export MAPPA_EXPERIMENT_WORK=/tmp/mappa-reuse-after
python prepare.py --repo /path/to/source-at-807d36f
python bench_reuse.py
python verify_outputs.py
```

Choose fresh work directories. The optimization patch changes only Python code after the native engine is built. Experimental counters use nested diagnostic dictionaries; the integrated production API exposes additional numeric counters with `symmetry_` and `validation_` prefixes. The production validator is formatted differently but has the same behavior. The published historical timings use the original saved baseline runs; a fresh reproduction will have its own timings.

The separate [sweep/prefix ablation](../sweep_prefix_ablation_20260915/README.md) explains why prefix competition has not replaced cut sweep. Its timings predate this optimization and must not be presented as timings of the new version.
