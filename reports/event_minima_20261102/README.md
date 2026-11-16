# Event-count minimum pilot

The local feasibility pilot certified all four selected directions. This is not
the complete Golden or EnzymeMap distribution. The full analysis remains pending.

| Saved input | Direction | Reference events | Certified minimum | Difference |
|---|---|---:|---:|---:|
| Historical Golden smoke, case 0 | Reverse | 7 | 5 | +2 |
| Historical Golden smoke, case 0 | Forward | 7 | 5 | +2 |
| Golden point, case 1561 | Forward | 3 | 3 | 0 |
| EnzymeMap case 37270 | Forward | 2 | 2 | 0 |

Both counts include heavy–heavy original-WBO changes under the saved thresholds.
All original matching constraints, including hydrogen constraints, remain in
force. These four minimizing witnesses have the same mapped heavy-atom coverage
as their references. Positive and negative differences are permitted; unknown
minima do not enter a numerical bin.

The analysis reads each completed decoder output once. Its saved local responses
give conservative family bounds; alternative fragment partitions are never added
together. A valid witness attaining the global bound certifies the minimum. Only
unresolved families need the existing exact constraint queries. No matching,
decoding or fragment rematching is repeated. The minimum concerns the retained
direction catalogue, not every possible chemical mapping.

All four pilot directions were certified by saved responses without solver
queries. In the final implementation's EnzymeMap 37270 retest, all 8,309 families
were certified in 0.115 CPU seconds of analysis. Reading the saved input took
0.313 CPU seconds; total worker computation was 0.492 CPU seconds, with a sampled
peak RSS of 0.437 GiB under the unchanged 2 GiB guard. This is one local observation.
The three Golden observations used the same bounds with earlier preparation
overhead; they have not been retimed under the final runner.

The [machine-readable evidence](pilot.json) retains each development measurement,
input checksum, source identity, budget and result. The initial catalogue-only
prototype certified 1/4 directions within 30 seconds. Saved-response reuse certified
all 4; removing redundant fingerprinting and whole-catalogue restoration removed
preparation overhead. The superseded prototypes are not runtime modes.

All 12 focused tests passed, including independent dense count checks for ordered
overlapping choices, hydrogen constraints, partial occupancy, alternative
partitions, positive invariance, unresolved competing families and negative bins.
The supported command and manifest are documented in [bench/README.md](../../bench/README.md).

Next is a bounded random pilot on current saved outputs, followed by the planned
3,702 Golden and 95,948 EnzymeMap direction analysis if its resource use and
certification coverage are acceptable. Reaction-level aggregation requires an
explicit check of directional coordinate and counting equivalence.
