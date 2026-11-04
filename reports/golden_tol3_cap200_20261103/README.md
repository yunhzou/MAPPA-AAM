# Independent-fragment decoding and Golden tolerance-3 study

## Scientific scope

The approved decoder uses the independent-fragment approximation for both
selection and reported event bins. Each unique fragment retains its jointly
constrained internal event spectrum. Ordinary source and target boundaries
contribute fixed breaking and forming events. A single sweep correspondence
keeps its endpoint-dependent correction. Original WBOs supply event counts.

The model minimum is obtained from additive local bounds and a sweep contribution
between zero and two. The default window retains model counts through minimum
plus six, inclusively. Count gaps are preserved. No complete-branch minimum
solver or bijection enumeration is part of ordinary decoding.

Reported bins are approximate. Independent local choices need not admit a joint
mapping witness, and the model minimum is not a proved physical minimum.
The complete compressed catalogue remains available; explicit mapping requests
still enforce its original joint constraints and injectivity. Output schema v10
records these distinctions. Exact reference verification remains a separate
catalogue query.

The implementation reuses local numeric spectra, native endpoint-domain
bookkeeping, NumPy response evaluation and integer-bitset count composition.
Alternative cut histories retain their own partition association. Sweep
cancellation only removes events actually present in the fixed boundary ledger.

## Validation

Critical checks cover correlated internal count gaps, sweep/internal tradeoffs,
fixed boundary approximation, original versus transformed weights, strict event
thresholds, separate cut histories, inclusive windows, shared preparation and
unchanged joint mapping constraints. The frozen implementation passes the full
suite: **991 tests in 161.12 seconds**. This includes local worker/persistence
tests, native/Python matching comparisons and separate exact event analytics.
The latter explicitly recovers a physical minimum excluded by the approximate
decoder window, proving that the two contracts remain separate.

Local spectra retain one widest computed numeric response per fragment and
endpoint query. Narrower windows require no new solver calls. Solvers and models
are released after each local projection; equal boundary ledgers share storage
during window selection. Cluster timings remain to be measured.

## Planned full Golden run

Fresh matching of all 1,851 reactions in both directions (3,702 directions):
discovery tolerance 3, relaxed-self tolerance 1, branch cap 200, one seed with
random seed 42, and the complete cut sweep. Event thresholds remain 0.5 normally
and 0.3 for metals; decoding uses minimum plus six under the model above.

Matching and catalogue reduction must finish across the entire study before
decoding begins. Persistent worker pools execute disjoint input batches with
shared reaction caches. No earlier matching results or caches enter this run.
Operational configuration and scheduler evidence remain in ignored private
storage. The run has not yet been launched.

This changes both discovery settings and the event interpretation model, so
historical event-bin values are not accuracy or timing evidence for this run.
