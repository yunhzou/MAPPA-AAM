# Competition prefixes along one greedy growth sequence

**Retaining intermediate takeovers helps.** On the five selected failure cases, the prefix policy recovers the same 12/14 SLAP-sweep minimum-event classes as local competition, and retains the additional four-event minimum class in case 11. Keeping only maximal extensions previously recovered 10/14. Prefixes also add distinct event classes in case 64, although the two missing SLAP-sweep minimum classes remain missing.

This is an isolated experiment. The production default remains local competition.

## Exact proposal policy

Start with a current-policy local takeover, including the first eaten atom `a`. Save that placement. Continue native greedy growth once, and retain its successive prefixes `a`, `ab`, `abc`, through saturation. Only atoms acquired from other fragments count toward the eating sequence; growth through atoms of the original B owner does not create an extra takeover size. The initial local placement stays fixed during continuation. As in the existing local policy, that initial placement can relocate B and need not include every atom of the old B owner.

These are actual intermediate candidate states, not subsets cut out of a final witness. The native observer retains compressed candidate states without changing the growth order. Earlier states are finalized lazily. Every prefix goes through the same takeover planning, dependent-fragment release, deduplicated completion and validation as local competition. Different compressed placements can exist at the same prefix size. No arbitrary overlap subsets or atom-bijection permutations are enumerated.

The initial local branch is retained before expansion. Unlike the previous maximal-only experiment, a valid earlier prefix can also survive if subsequent growth reaches a branch cap. The observer extension adds support for a mapped seed; ordinary growth remains unchanged.

## Results at the same 128-completion budget

All cases use the same saved one-seed, cap-2,000 forward cut-sweep archives as the [previous experiment](../unrestricted_competition_20260915/README.md). All variants retain the original no-competition result. Cases were selected for known missing native or swept SLAP alternatives; this is not a new held-out accuracy estimate, and their mappings are not verified ground truth.

| Case | SLAP-sweep minimum classes | No competition | Local | Maximal extension only | Prefixes |
|---|---:|---:|---:|---:|---:|
| 6 | 4 | 3 | 4 | 3 | 4 |
| 11 | 2 | 1 | 2 | 1 | 2 |
| 59 | 1 | 1 | 1 | 1 | 1 |
| 64 | 5 | 2 | 3 | 3 | 3 |
| 101 | 2 | 1 | 2 | 2 | 2 |
| **Total** | **14** | **8** | **12** | **10** | **12** |

Prefixes recover four of the six SLAP-sweep classes missing without competition. They do not recover the two remaining classes in case 64.

Against native SLAP without the wrapper, local and prefix competition each recover 7/12 classes, but not the same seven: prefixes gain one in case 64 and lose the local recovery in case 59. Their union recovers 8/12. Native and swept SLAP targets overlap and have different per-case event minima; do not add these totals.

### What the intermediate stopping points add

- **Case 11:** a four-event minimum class absent from local competition is represented by saved prefixes after eating **3, 4 or 5 atoms**. This is the same additional minimum class previously found by anchored full growth. It is absent from the evaluated SLAP minimum classes. Its two broken pairs are `(0,49)` and `(50,51)`, and its two formed pairs are `(0,51)` and `(49,50)` (zero-based source indices).
- **Case 64:** three classes absent from all previously tested policies appear: two six-event classes at a **five-atom** prefix, and one seven-event class at prefixes of **six through eleven atoms**. The seven-event class is an additional native-SLAP minimum class; the swept-SLAP minimum is four events, so this does not close its remaining gap.

All these provenance examples are intermediate, unsaturated prefixes. They directly demonstrate that the chosen stopping point changes the available final event patterns. Full mappings, symmetry actions, exact family constraints and responsible proposals are in `novel-candidates.json.gz`. Each witness passed an independent exact-family membership query and canonical event-ID reconstruction.

## Runtime and decoding

CPU-seconds across the five cases, including input loading, catalogue construction and event decoding but excluding the original sweep search:

| Policy | Competition | Catalogue and decoding | Total |
|---|---:|---:|---:|
| Local, previous run | 3.410 | 11.392 | **14.802** |
| Maximal extension only, previous run | 71.487 | 95.009 | **166.496** |
| Prefixes | 28.099 | 38.465 | **66.564** |

These are same-Mac measurements with one native thread per worker. The batch used two workers. Prefixes cost about 4.5 times the local total and 40% of maximal-only extension here. Case 11 dominates: 57.2 CPU-seconds and 12,336 retained compressed families. The other four cases together cost 9.34 CPU-seconds. This is family growth and completion cost, not enumeration of every represented bijection. Startup/import time and search checkpoint writes are outside these timers.

All five final catalogues were fully decoded within the same pre-existing event windows: at most five events except case 64, which uses eight. Event thresholds are 0.5 and 0.3 for metal bonds; matching uses `iso_tolerance=1` and graph floor 0.2. There is no decoded-pattern count cap. All retained final families passed representative validation. Invalid proposed repair graphs were rejected before catalogue inclusion.

## Budget sensitivity

The unchanged 128-completion budget can interrupt a trajectory or prevent later local proposals from being visited. It is not a promise to complete every recorded prefix. Prefix sizes observed and finalized (excluding the initial one-atom proposal) were:

| Case | Observed sizes across sessions | Finalized sizes |
|---|---:|---:|
| 6 | 107 | 107 |
| 11 | 104 | 96 |
| 59 | 71 | 59 |
| 64 | 28 | 23 |
| 101 | 108 | 104 |

Repeated sessions/cache reuse and alternative placements are recorded separately. These counts describe prefix-size visits, not unique branches or completed mappings. Finalizing a size also does not guarantee that every placement at that size was completed before the budget was reached.

Because cases 59 and 64 were inexpensive, they were additionally tested with **512 completion calls**; the AAM branch cap remained **2,000**:

| Case | Event classes, budget 128 → 512 | SLAP-sweep recovery | Native-SLAP recovery | Total CPU-seconds at 512 |
|---|---:|---:|---:|---:|
| 59 | 9 → 10 | 1/1 → 1/1 | 0/1 → 0/1 | 1.490 |
| 64 | 21 → 30 | 3/5 → 3/5 | 2/6 → 2/6 | 1.056 |

The larger case-64 run has eight classes absent from local competition, versus three at budget 128. It still does not recover the remaining SLAP-sweep minima. Both additional decodes completed. These are targeted budget checks, not a five-case 512-budget comparison.

Parent limit eight, competition depth two and queue limit 512 are unchanged. Each main run allows up to 128 local proposal calls, 128 continuation sessions and 128 completion calls; all three call budgets become 512 in the targeted checks. Each worker has a five-minute watchdog and 3,000 MiB memory guard. No final run hit either resource guard. All seven searches reached their completion-call budget, so no global completeness claim is made.

## Checks and reproduction

Fourteen focused tests passed, including the existing fragment-choice tests and a mapped-seed chain test that checks every prefix. A separate local-mode control exactly reproduced production counters, accepted offers, pending work and all repair graphs on case 64. All accepted takeover records preserve their initial placement and agree with their recorded eaten-atom sequence; observed prefix sizes are contiguous.

`prefix.patch` applies only to an isolated copy. Its exact base and resulting source hashes are in `manifest.json`. The baseline archives and SLAP reference sets are reused from the adjacent `unrestricted_competition_20260915` report. With repository dependencies plus `pybind11`, `setuptools`, `psutil`, `pytest` and `z3-solver` installed:

```sh
export MAPPA_EXPERIMENT_WORK=/tmp/mappa-prefix-test
python prepare.py --repo /path/to/coordinate_alignment
PYTHONPATH="$MAPPA_EXPERIMENT_WORK/engine/src" python -m pytest tests/test_fragment_choices.py test_prefix.py -q
python check_control.py
python verify_candidates.py
python run.py
python run_budget.py
```

`prepare.py` builds both the regular extensions and the native growth engine in the copied source. It never edits production code. Choose a fresh work directory. Raw search offers, prefix audits, decoded patterns and per-family completion certificates are in `results.json.gz`; per-case counts and timings are in `summary.csv` and `summary.json`. Candidate proof records cover both budgets, so the same event class can have more than one saved witness record.

## Subsequent timing control

The total-time advantage above compared budget-limited searches visiting different
proposals. It is not an intrinsic speed advantage of prefixes. A matched-trajectory
control, completing every prefix of the same five trajectories, found prefixes
slower in every case after warmup. See the [six-arm sweep ablation and paired
control](../sweep_prefix_ablation_20260915/README.md). These experiments use the
pre-optimization source pinned in their manifests.
