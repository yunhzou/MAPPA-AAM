# Explicit execution cleanup

The user clarified that the prohibition on code fallbacks does not remove the
permanent first cache snapshot. Cache recovery remains unchanged.

Native growth now rejects unsupported tracing, node policies and orbit metadata
instead of selecting Python. The old capability-based dispatch predicate is
removed. Diagnostic trajectory replay explicitly invokes Python incremental
tracing and still verifies the saved fragment decisions.

Direct search, work pools and fragment APIs now default to the accepted Golden
settings: existence growth and reusable native reaction search. The explicitly
requested incremental and reference comparison modes remain available. There is
no error-driven selection of those modes. Growth, compression, hierarchy caps,
deduplication and decoding algorithms are unchanged.

The checkpoint loader no longer reinterprets the current `EventDecodeConfig` as
historical rematching configuration. A bundle round trip now preserves its exact
type and rejects obsolete rematching fields. Existing Golden bundles that were
already written with the wrong metadata class are preserved unchanged.

## Validation

- Rebuilt the native engine on the Mac.
- Full local suite: 725 passed; 23 tests needed explicit backend/policy selection
  after the default and dispatch changes. The 19 unchanged test IDs passed on
  rerun. The four ambiguity cases became eight cases testing both policies;
  those and the default-worker checks passed in a 24-test affected-module run.
  All 753 current cases are covered across these runs, including the independent
  symmetry oracles, exact Python/native comparisons, cache recovery, viewer,
  storage and worker-pool tests.
- Loaded the saved Golden `01561/forward` v4 bundle without recomputation.
  The complete decoded JSON record retains SHA-256
  `6f67159f8a07143ba94a67f5534159e52d543a36b99f7be492cad3980802bad6`.
- No benchmark rerun was needed for this cleanup.
- A distributed validation rebuilt the engine from commit `254b862` and completed
  successfully: 29 smoke checks and 52 affected-module tests passed. The latter
  cover native dispatch/parity, checkpoint types, default execution, branch caps
  and shared cache recovery. Build and tests took 67 seconds in one allocation.

Raw operational records and full artifacts are retained privately.
