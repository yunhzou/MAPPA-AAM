# Saved minimum-event branch compositions in Golden

All **1,851 reactions** have at least one recorded fragment composition at the
saved detected minimum. Both directions are included and normalized to the same
endpoint indexing. No additional matching, optimization or symmetry search was
performed. The count pass took **36.665 seconds** with eight processes.

| Statistic | Branch compositions |
|---|---:|
| Median | 31 |
| Mean | 51.06 |
| 95th percentile | 120 |
| Maximum | 2,250 |
| Reactions with one composition | 8 |

![Saved minimum-event branch compositions](minimum-branches.png)

A branch is identified by its sorted collection of base/target fragment atom-set
pairs. Repeated histories and internal atom permutations count once. Different
fragment partitions remain different branches, including when they can realize
the same complete atom mapping. Consequently these are composition counts, not
chemically inequivalent mapping counts or random-selector accuracy.

## Calculation

Use the original one-seed Golden decoded outputs (tolerance 2, branch cap 100,
both directions and cut sweep). Keep their saved global detected minima fixed.
The older archives retain local optimized witnesses and event responses but not
the attained count of every branch. Replay those saved local choices in their
recorded action order, leave global placement choices unchanged, and directly
count internal events plus the saved boundary contribution and sweep correction.
Keep a composition only if that actual count equals the saved detected minimum.
This counts recorded choices; it does not establish every further symmetry tie.

The pass examined 244,386 retained contexts. No replay produced a count below
the saved detected minimum. No reaction is omitted from the figure. The original
input index, plan and result ledger hashes match the previous Golden analysis.
The sanitized per-reaction evidence includes the original decoded-output hashes.

## API and validation

Analysis revision: `8cf8392`. Decoder query revision: `8614e42`.

New decoder outputs retain the branch witnesses already computed during window
selection, using shared action references. `decoded.saved_branches()` queries
recorded compositions at the detected minimum; `event_count=5` selects another
count. `decoded.saved_minimum(...)` retrieves a branch witness and bound.
Neither lookup starts a search. Older archives explicitly lack the new saved
branch-count lookup; the read-only analysis above recounts their preserved local
choices rather than inventing branch minima.

Validation: 23 decoder-window tests, three branch-count tests and 26 smoke tests
passed. Tests cover duplicate compositions, internal shuffles, reverse indexing,
save/load and sweep corrections; solver calls are forbidden during lookup tests.
One hundred saved duplicate branch records returned one composition in 2.92 ms
per lookup on the local Mac. The four-reaction data check took 1.50 seconds.

Reproduce from saved data with `bench/minimum_branches.py`; its `--plot` command
uses the aggregate JSON. Operational configuration is private. The raw source
archive remains on its original host; no benchmark is rerun.
