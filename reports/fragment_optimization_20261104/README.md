# Integrated fragment symmetry responses

Production integration completed on 2026-11-05. The shared core now supplies
binary symmetry construction, bounded starting-correspondence improvement,
minimum witnesses and target-event queries. `decode_events` retains unique
branches referencing these fragment responses. The obsolete self-tolerance
filter and full attainable-count-table implementation are removed.

The matching catalogue keeps its full discovery relation. Decoding deliberately
retains the improved starting correspondence plus both-sided binary symmetry.
Other discoveries on the same atom sets may lie outside this orbit. Explicit
catalogue proofs remain inspectable; target-event queries report membership in
the selected decoded relation. Ordinary boundaries use the accepted independent
approximation, while the single sweep retains its joint endpoint/internal cost.
Historical overlapping actions keep their order and shared choices.

A five-second default bounds starting-map improvement. Best-found witnesses are
separate from proved minima. A compatible composition of proved local minima
provides the initial lower-bound-plus-two ceiling; otherwise the decoder uses an
attained upper bound and retains unresolved possibilities. Serialized v12 output
restores selected starts, responses and window decisions without rerunning
optimization or symmetry preparation. Runtime timing counters are not part of
its scientific record. RAM and aggregated bundle storage use the same API.

## Integrated validation

Full local release suite: **1,074 passed in 164.39 seconds**. The tests include
native/reference matching, branch caps, shared caches, direct and persistent-pool
execution, RAM/bundle equivalence, resume budgets, serialization and scientific
decoder oracles. The local Slurm transport tests simulate the allocation entry;
no actual cluster allocation was launched. Worker tests ran outside the sandbox
to permit local sockets and process-tree memory inspection.

The final focused optimizer/response/threshold checks passed **135 tests in
2.24 seconds** before the full release run. Earlier integration attempts exposed
retired API assertions, nondeterministic timing fields in scientific records and
missing new configuration fields in campaign fixtures; those were corrected
before the final passing run.

The critical independent tests cover all 384 mappings of weighted C(CH2)4,
including exactly 96 mappings at two events and gaps at odd counts; endpoint
conditions cannot assign two arms to the same position. Tiny unequal-connectivity
oracles check both-sided symmetry, anchors and original-WBO costs. Ordered action
tests reject incompatible local minima. Save/load tests forbid fresh optimization
and check sharing, selected coverage and window exclusions. Explicit limits return
unknown rather than absence. Existing native/reference matching, cap, cache,
transport, resume and work-pool tests remain part of validation.

## Saved-fragment timing check

Twelve saved Golden fragment placements, each with original and independently
calculated GFN2-xTB weights, give **24 responses and 72 target queries**. Every
minimum count and query outcome agrees with the accepted standalone beta. These
extracted fixtures retain atom sets, witnesses, anchors and local weights, but
lack their historical growth edges; this measures binary-orbit responses, not
full discovery improvement or a complete benchmark rerun.

Single-threaded local measurements, five repetitions:

| Operation | Minimum | Median | Maximum |
|---|---:|---:|---:|
| Response preparation | 1.502 ms | 6.790 ms | 57.095 ms |
| Target query | 0.0019 ms | 0.0718 ms | 0.1266 ms |

An interleaved comparison with beta commit `96c1f11` gives a median preparation
CPU ratio of **1.026**, with a median absolute overhead of **0.122 ms**. The
largest measured overhead is **0.509 ms**. This is not a claim of zero regression
on every input or of full-campaign throughput. MAPPA reuses structural generators
already stored in Fragment; beta constructs them during preparation.

The 69-atom unequal-connectivity fixture represents approximately **4.26e20
structural operation pairs**, not that many distinct mappings. Its event response
is constant; no complete permutations are expanded. Three shared query views,
including local matrices, occupy 5.4 KB compressed for original weights and
7.8 KB for xTB weights. These sizes are not a whole-reaction RAM bound. The
weighted arm test supplies variable counts that these real fixtures lack.

Inputs and measurements: [response-inputs.json.gz](response-inputs.json.gz),
[response-results.json.gz](response-results.json.gz). Replay without matching:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python -m bench.fragment_events \
  reports/fragment_optimization_20261104/response-inputs.json.gz \
  /tmp/fragment-response-measurements.json.gz \
  --reference ../symmetry-event-lab
```

Omit `--reference` to replay MAPPA alone. The optional comparison requires the
accepted standalone beta checkout; it is never a runtime algorithm alternative.
No cluster run or full Golden/EnzymeMap recovery rerun was performed for this
integration. Historical measurements below retain their original scope.

---

# Earlier bounded fragment optimization study

The core operation `optimize_fragment(fragment, time_limit=5.0)` improves an
actual Fragment's internal bond-event witness. It reuses the ported beta
symmetry optimizer and the existing discovery-constraint compiler. It does not
modify the saved relation, split candidates, or replace the production decoder.

## Algorithm and guarantee

1. Retain the saved valid mapping. Obtain a lower bound by pairing original bond
   weights within element-pair types, including zeros and subthreshold contacts,
   while ignoring atom consistency.
2. Optimize the saved mapping using its two-sided binary symmetry.
3. If needed, search the original growth constraints for a better starting map
   and apply the same symmetry optimizer. Share one deadline across all work.
4. Stop on bound attainment or budget exhaustion and return a valid witness,
   actual event count, lower bound, completion status and work counters.

The discovery search prunes using partial raw event counts. It can skip a
starting map whose raw count is worse but whose symmetry could improve it, so
exhausting that search alone does **not** prove a full-relation minimum. A minimum
is certified only by attaining a valid lower bound or completing symmetry that
provably covers every discovery. When growth uses every target binary edge,
target symmetry supplies that coverage certificate and avoids rediscovery.

The five-second allowance is per call and includes event preparation, discovery
search and symmetry search; constructing the input Fragment is outside it. Core
checks are cooperative. An opaque native operation can overrun until it returns;
strict preemption belongs to an external execution watchdog. No per-fragment
process creation was added. These local results do not establish a total memory
bound, an exhaustive event range, or reaction-level decoder accuracy.

## Measurements

Local macOS/ARM64, one numerical thread. The comparison below uses the unchanged
standalone beta revision `b4e5f52e05414bafe043e5f05bb8fd9376ff369a`.
Times include preparation and solution, exclude imports and Fragment creation,
and are medians of five interleaved runs after one warmup. The port consumes
the binary generators already stored in its actual Fragment; the standalone
beta constructs them during preparation. These are operation timings, not a
complete matching/decoding benchmark.

| Example | Atoms | Beta | Ported symmetry | Bounded search | Events, all three |
|---|---:|---:|---:|---:|---:|
| Four CH2 arms | 13 | 1.397 ms | 1.521 ms | 1.632 ms | 0 |
| Four CH2 arms, nonzero minimum | 13 | 1.643 ms | 1.828 ms | 1.920 ms | 1 |
| Weighted C6 ring | 6 | 0.387 ms | 0.464 ms | 0.499 ms | 2 |
| Weighted three-level tree | 40 | 17.814 ms | 18.933 ms | 19.544 ms | 6 |
| Three-level tree with one target edge removed | 40 | 464.261 ms | 423.074 ms | 423.763 ms | 6 |

The added bounded-search layer costs 0.035–0.689 ms over the ported symmetry
operation on these cases. This is a small point comparison, not a claim of no
runtime regression on every fragment.

Starting-map tests, each with a five-second allowance:

| Example | Best in saved symmetry | New witness | Lower bound | Elapsed | Status |
|---|---:|---:|---:|---:|---|
| Triangle with tail, incomplete growth edges | 2 | 0 | 0 | 0.574 ms | Minimum certified |
| Path into weighted double triangle | 3 | 2 | 2 | 0.547 ms | Minimum certified |
| C80 path with distinct additional chords on each side | 40 | 40 | 0 | 5.005 s | Best found; minimum unproved |

The C80 stress case shares the 79 path edges. Each side has 20 different extra
chords, selected without replacement with NumPy seed 42. It tests bounded work
when the lower bound is loose; the optimizer did not improve this witness.
The complete synthetic inputs, source hashes, individual timings and returned
witness proofs are in [measurements.json.gz](measurements.json.gz).

## Validation

- Full repository suite: **1,050 passed in 140.72 s**. This was run outside the
  filesystem sandbox because existing process-memory tests need OS process-tree
  inspection. The sandboxed attempt failed on that OS permission.
- Final focused optimizer suite: **47 passed in 1.01 s**, including one
  additional test after the full run for original event weights, metal thresholds
  and noncontiguous atom indices. Production code subsequently changed only in
  docstrings.
- Smoke suite: **26 passed**, 1,024 deselected.
- The independent tiny-case tests verify every returned mapping and event count,
  compare bounds against complete pair-assignment oracles, and check that a
  claimed minimum agrees with the complete mapping relation. They cover anchors,
  strict event thresholds, unattainable bounds and interrupted preparation.

## Reproduce

Use the project Python environment from the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python -m pytest -q -p no:cacheprovider tests/test_fragment_optimization.py
```

Replay the saved inputs directly, without rerunning matching or reconstructing
the examples by hand:

```python
import gzip
import json
from mappa.mappa_core.core.matching.representation.records import Fragment
from mappa.mappa_core.core.symmetry.optimization import optimize_fragment

with gzip.open("reports/fragment_optimization_20261104/measurements.json.gz", "rt") as stream:
    data = json.load(stream)
for row in data["comparison"] + data["starting_map_cases"]:
    fragment = Fragment.from_record(row["fragment"])
    answer = optimize_fragment(fragment, time_limit=data["budget_seconds"])
    print(row["name"], answer.event_count, answer.lower_bound,
          answer.complete, answer.elapsed_seconds)
```

Time-limited witnesses and work counts may vary by hardware and load. Only
certified minima carry an optimality guarantee.
