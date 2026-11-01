# Golden persistent-pool run, 2026-11-01

The user authorized a fresh high-parallelism Golden run with **2 GiB per worker**,
replacing the earlier 4 GiB resource choice for this run. Scientific limits remain
one seed, full bidirectional cut sweeps, branch cap 100, matching tolerance 1,
graph/cut floor 0.2 and random seed 42. Existence growth is explicitly selected.

All **1,851 original Golden records / 3,702 directions** are included. Each
prepared input was checked against the original per-record SHA256 manifest
(`652285008cea1d10a21ae80e476e0f6c0ff17df831f546dee426d1fd0306442b`).
Endpoints are copied unchanged from the trusted archived JSON inputs. No old
matching implementation or scientific result is executed or reused. Input and
reference records occupy one aggregate bundle, with five small batch indexes.

Golden verification retains the original complete heavy-atom relation, including
unmatched atoms, modulo chemically colored endpoint automorphisms. Generators
are prepared from the archived element/charge/isotope/H/stereo colors and bond
colors, then lifted to original atom indices. Reference data enters verification
only. This differs from EnzymeMap's exact-index verification policy; their
coverage percentages must not be conflated.

The existing campaign runner owns all stage orchestration. `bench/golden.py`
only prepares verified inputs and provides the Golden reference query using
MAPPA's existing correlated-family solver. The runner now accepts direct prepared
input references and a manifest-declared verification operation/policy. Core
matching, hierarchy caps, cache logic, deduplication and decoding are unchanged.
Eight affected tests passed in 10.41 seconds, including independent exhaustive
small colored-graph permutations and complete pipeline checks of unmatched-atom
semantics. The Slurm script passed shell syntax checking.

DFW `cpu_short` permits ten nodes per user, with 96 CPUs per node. The initial
two-node / 100-worker batch requests could not start: most spare CPUs belonged to
nodes whose memory was already reserved. A revised batch spreads 100 workers over
four nodes, requesting 28 CPUs and 54 GiB per node, with 25 compute workers per
node. Each worker remains guarded at 2 GiB; each node has 4 GiB above the sum of
worker limits. Five disjoint batches request 500 workers in total, but the
ten-node quota permits at most two revised batches (200 workers) simultaneously,
subject to other jobs. The allocation limit is one hour; admission stops after
50 minutes and drains active work. Requests are not proof of running parallelism.

Case IDs modulo five establish disjoint ownership. Each allocation starts one
persistent pool with one `srun`; tasks never invoke Slurm. Shared reaction caches
retain the permanent baseline and latest versions. Matching is submitted with
`--stop-after matching`; a dependent array resumes the same ledgers only after
the entire matching array succeeds, then performs global per-direction
deduplication, decoding and verification. Resource failures remain explicit and
block the next pass; no automatic memory increase or silent retry is enabled.

Launch and live-progress evidence are recorded alongside this report. Aggregated
inputs/results live under the DFW MAPPA root in `golden_pool_20261101/`.

The original Golden case 0 completed all four stages locally in both directions
in 15.05 seconds with two workers at 2 GiB each. Both reference queries recovered
the annotation; [local-preflight.json](local-preflight.json) records the checks.
The [input manifest](inputs.json) retains every input/reference hash.
