# Can a competing fragment keep growing?

Allowing the challenger to absorb more atoms adds some distinct event patterns, but it does not replace the current local policy successfully on these cases. Keep the local takeover branch if further growth is explored. Production defaults are unchanged; the broader policies are an isolated experiment.

## What was tested

Cases **6, 11, 59, 64 and 101** from the 140 coordinate reactions were selected because the no-competition result missed minimum-event alternatives returned by native SLAP and/or SLAP with the edge-sweep wrapper. This is a failure-case development experiment, not a new accuracy benchmark. These coordinate mappings are unverified; matching a SLAP event class does not establish chemical ground truth.

All variants start from the same saved **one-seed, branch-cap 2,000, cut-sweep** result and retain that baseline. Four competition policies are compared:

- **Local (current):** grow on B plus one neighboring atom, allow B to relocate, then complete the remaining mapping.
- **Anchored full growth:** preserve B's current placement and let it grow on the whole source graph. Other fragments can be consumed.
- **Free full growth:** restart from B's boundary seed on the whole graph, allowing its placement to change.
- **Local then extend:** obtain a local takeover, fix that placement, then grow on the whole graph before completing the mapping. This experiment saves the extended proposal; it does not also save every intermediate local proposal.

Full growth has no added atom-count limit. It remains greedy and subject to matching constraints and the existing search budgets. Accepted proposals absorbed up to **72 atoms**, so this was not another one-atom restriction.

## Recovery of SLAP minimum-event alternatives

Entries are the number of distinct minimum-event classes recovered, after event-pattern deduplication. A class is a signed bond-event pattern modulo the endpoint symmetries used by the evaluator, not a raw atom bijection.

### SLAP with the edge-sweep wrapper

| Case | Available SLAP classes | No competition | Local | Anchored full | Free full | Local then extend |
|---|---:|---:|---:|---:|---:|---:|
| 6 | 4 | 3 | 4 | 3 | 3 | 3 |
| 11 | 2 | 1 | 2 | 1 | 1 | 1 |
| 59 | 1 | 1 | 1 | 1 | 1 | 1 |
| 64 | 5 | 2 | 3 | 2 | 2 | 3 |
| 101 | 2 | 1 | 2 | 2 | 1 | 2 |
| **Total** | **14** | **8** | **12** | **9** | **8** | **10** |

Of the **six classes missing without competition**, local takeover recovers four, anchored full growth one, free full growth none, and local-then-extend two. Combining all tested policies still recovers 12/14: the two remaining SLAP-sweep classes in case 64 were not recovered.

### Native SLAP, without the sweep wrapper

| Case | Available SLAP classes | No competition | Local | Anchored full | Free full | Local then extend |
|---|---:|---:|---:|---:|---:|---:|
| 6 | 2 | 2 | 2 | 2 | 2 | 2 |
| 11 | 1 | 0 | 1 | 0 | 0 | 0 |
| 59 | 1 | 0 | 1 | 0 | 0 | 0 |
| 64 | 6 | 0 | 1 | 0 | 0 | 0 |
| 101 | 2 | 1 | 2 | 2 | 1 | 2 |
| **Total** | **12** | **3** | **7** | **4** | **3** | **4** |

Native and sweep class sets overlap; these totals must not be added. Every broader variant retains the original baseline. Its lower score than local means it failed to reproduce alternatives that local competition added.

## New alternatives outside the SLAP targets

Anchored full growth adds **one new four-event minimum class in case 11**, absent from the current local result and the evaluated SLAP minimum classes. Its witness breaks pairs `(0,49)` and `(50,51)` and forms `(0,51)` and `(49,50)` (zero-based source indices). The responsible depth-two proposal, `takeover16`, grows B from atom 65 into a 21-atom placement, absorbing 20 atoms from other fragments. This is a distinct allowed mapping/event alternative, not a verified reaction pathway.

Local-then-extend adds **ten higher-event classes** relative to local: one five-event class in case 6, and nine classes with six to eight events in case 64. Free full growth adds none. Combining the saved local result with all broader results therefore adds **11 event classes**, including the one new minimum class. This union is an analysis of separate runs; it is not a jointly budgeted combined search.

Every new representative was independently queried against its exact saved mapping family and its canonical event identifier recomputed. All 11 checks passed. The witnesses, family constraints, symmetry actions and provenance are in `novel-candidates.json.gz`; `verify_candidates.py` repeats these checks.

## Cost

CPU-seconds summed over these five cases, on the same Mac with one native thread per worker:

| Policy | Competition | Catalogue + event decoding | Total |
|---|---:|---:|---:|
| No competition | 0.000 | 9.096 | 9.096 |
| Local | 3.410 | 11.392 | **14.802** |
| Anchored full | 2.032 | 9.949 | **11.981** |
| Free full | 2.502 | 10.155 | **12.657** |
| Local then extend | 71.487 | 95.009 | **166.496** |

These times exclude generating the original sweep archives, process startup and imports. Competition includes loading its input archive; decoding includes loading archives, deduplicating final branches, rebuilding shared symmetries, writing the catalogue and extracting event classes. Thus these are the added post-search costs, not complete AAM runtimes. Search output checkpoint writes are outside the search timer. CPU time is not elapsed wall time.

Local-then-extend is about **11.2 times** the local total here. Case 11 accounts for most of that cost: 66.1 search CPU-seconds and 84.2 decoding CPU-seconds, with 31,948 retained families. This grows the number of compressed families; it does not enumerate all atom bijections. Other full-growth variants are inexpensive on this subset, but recover fewer SLAP targets.

The batch used two workers, one numerical thread each, a five-minute per-stage watchdog and a 3,000 MiB per-worker memory guard. The largest observed worker used about 2.51 GiB. All final jobs finished without a resource stop. Per-case CPU, elapsed time, memory, pending work and budget counters are preserved in the result archive.

## Interpretation and limits

Larger takeovers change the constraints imposed on subsequent completion. A plausible explanation for the lost alternatives is that greedy extension commits to a larger preserved fragment before the completion stage can explore the intermediate placement. The results support retaining the local branch and treating further growth as an additional source of proposals. They do not prove that all unrestricted policies are inferior, or that the remaining alternatives are unreachable.

Each search uses the current parent ranking, at most eight parents, competition depth two, queue capacity 512, and budgets of 128 proposal calls and 128 completion calls. Local-then-extend additionally permits 128 full-growth expansion calls; it therefore has an extra matcher-call budget. Completed calls still use branch cap 2,000. Pending work and capped/rejected calls are recorded. There is no claim of exhaustive search over all takeover sequences or intermediate growth stopping points.

Matching uses `iso_tolerance=1` and graph floor 0.2. Event extraction uses bond thresholds 0.5 and 0.3 for metal bonds. All retained families were fully decoded **within the pre-existing event windows**: at most five events for cases 6, 11, 59 and 101, and at most eight for case 64. There was no event-pattern count cap. All 25 catalogue decodes completed with no invalid final families. This does not mean patterns above those windows were enumerated.

The experimental `local` policy was checked against the original implementation on case 64: counters, offers, pending queue and every repair graph agree exactly. An early extended-case-64 anchor assertion was corrected before its successful rerun; the audit note is retained in `manifest.json`.

## Reproduction

The production algorithm is not modified. `competition.patch` applies to the source commit recorded in `manifest.json`, and `prepare.py` creates an isolated engine. The five saved baseline inputs are included, so no expensive initial sweep is necessary. SLAP references enter comparison only, never proposal generation.

Using a Python environment with the repository dependencies, `z3-solver`, `psutil`, `pybind11`, and setuptools installed:

```sh
export MAPPA_EXPERIMENT_WORK=/tmp/mappa-full-growth
python prepare.py --repo /path/to/coordinate_alignment
python check_control.py
python verify_candidates.py
python run.py
```

`prepare.py` builds the native extensions in the isolated copy; `--skip-build` is available if compatible compiled extensions were already copied from the source tree. Choose a fresh work directory. `run.py` saves all search archives, complete decoded pattern dictionaries, per-family completion certificates and execution records there. To run one stage: `python run.py worker 64 anchored search` (direct worker invocation does not apply the parent watchdog; normal batch execution does).

Recorded evidence:

- `summary.csv` and `summary.json`: per-case counts and measured CPU times.
- `results.json.gz`: all final search records and offers, decoded patterns, completion certificates and execution records.
- `novel-candidates.json.gz`: all 11 additional witnesses and their exact families.
- `references.json`: fixed native and swept SLAP comparison sets.
- `manifest.json`: source revision, source/input hashes, settings and audit notes.
- `inputs/`: the five baseline compressed AAM checkpoints.
