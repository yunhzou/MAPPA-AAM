# Does takeover-prefix competition replace cut sweep?

**Not in its current post-search form.** In a fixed 25-reaction pilot, prefixes without sweep retain higher minimum event counts in cases **35, 64 and 69**. With sweep, both local and prefix competition recover the saved MAPPA minimum classes on all 25 cases. The production configuration is unchanged.

## Factorial pilot

The pilot includes the five development cases (6, 11, 59, 64, 101) and 20 further cases selected by a fixed SHA-256 ordering of case IDs, without consulting outcomes. The selection and inputs were frozen before running the pilot. This is still a development-set experiment on unverified coordinate mappings, not a mapping-accuracy benchmark or independent validation set.

The same input, forward direction, one seed, random seed 42, branch cap 2,000, matching tolerance 1 and event thresholds 0.5/0.3 were used for all six configurations. Competition has 128 completion calls, eight selected parents, depth two and queue limit 512. Prefix mode additionally retains intermediate stopping points along each started greedy takeover trajectory; its completion budget can interrupt that sequence.

| Cut sweep | Competition | Reactions matching saved MAPPA minimum count /25 | Saved MAPPA minimum classes /35 | SLAP-sweep minimum classes /36 | Search + competition CPU s |
|---|---|---:|---:|---:|---:|
| Off | None | 20 | 21 | 22 | 0.500 |
| Off | Local | 22 | 30 | 30 | 18.538 |
| Off | Prefix | 22 | 30 | 30 | 121.558 |
| On | None | 25 | 30 | 30 | 21.486 |
| On | Local | 25 | 35 | 34 | 72.457 |
| On | Prefix | 25 | 35 | 34 | 161.902 |

The saved MAPPA target is the previously decoded default sweep-plus-local result, not ground truth. These columns count recovery among all returned candidates within the frozen event window; they do not claim the target sets contain every valid alternative. The two missing SLAP-sweep classes are in case 64. Pooling local and prefix outputs does not improve either target-coverage total in this pilot, although their broader event sets can differ.

Without sweep, prefix competition's minimum is **2 versus 1** in case 35, **6 versus 4** in case 64, and **4 versus 3** in case 69. Thus it recovers much of the sweep result, but cannot yet replace it without regression.

## Seed-order control: cuts have an effect beyond extra starts

The sweep changes both the source constraints and the per-cut seed order. To separate those effects in case 64, all **21 exact seed orders** recorded by its sweep were replayed on the uncut graph. Prefix competition was then applied to that union.

All 42 returned full terminal states had one fragment. Their completely decoded event window still has minimum **six**, versus **four** with sweep. There was no second fragment for takeover competition, and no competition proposals were made. This is a concrete limitation of post-search competition: it revises existing fragment boundaries; it does not split a saturated one-fragment result to create its first challenger.

An initial-growth stopping policy could address a different search space, but that is not what this ablation implements. No claim is made that prefixes at every stage could never replace sweep.

## Correcting the earlier prefix timing interpretation

The earlier 66.6 versus 166.5 CPU-second result compared **budget-limited searches that visited different proposals**. It did not show that completing every prefix on the same trajectories was cheaper than completing only their maximal endpoints. In particular, the maximal-only case-11 run retained 31,948 families, while the truncated prefix run retained 12,336.

A separate paired control fixes the first started trajectory from each of the five development cases: the same parent, local placement, seed and contact in both arms. Further generations are disabled, and every compressed placement at every observed prefix is completed. Every maximal-result family is verified to occur in the prefix result by its full canonical family key; decoded maximal event sets are subsets too.

| Case | Maximal-only families | All-prefix families | Maximal completion calls | Prefix completion calls |
|---|---:|---:|---:|---:|
| 6 | 10 | 15 | 1 | 3 |
| 11 | 6 | 21 | 1 | 5 |
| 59 | 1 | 17 | 1 | 28 |
| 64 | 3 | 10 | 3 | 52 |
| 101 | 1 | 5 | 1 | 8 |
| **Total** | **21** | **68** | **7** | **96** |

After warming both arms, three repeats with alternating execution order showed higher median competition CPU and event-extraction CPU for prefixes in every case. Sum of per-case medians: maximal competition **0.193 s**, prefix **0.336 s**; maximal event extraction **0.035 s**, prefix **0.071 s**. These small matched-trajectory timings exclude archive loading and catalogue construction and are not whole-reaction performance estimates. The evidence is consistent with the expected additional work; the earlier speed reversal came from different budget allocation.

## Completeness, resources and timing scope

All **50 case/sweep pipelines**, representing **150 configurations**, completed their saved-family event-window decoding. Shared exact families were decoded once across the three competition arms, retaining each arm's family membership. Empty searches, unresolved decoding and patterns outside the chosen windows are distinct conditions; all 25 reactions had full mappings in every arm. Windows are the same frozen per-case windows used by the prior coordinate evaluation. They are not unrestricted enumeration of every possible event count.

The table reports fresh same-Mac **search + competition CPU**, including successful-stage archive persistence and excluding imports. Decoding is not included in those columns because the three arms share decoding work. Joint catalogue/decoding CPU was 125.013 s without sweep and 255.006 s with sweep. Successful measured stage CPU across the experiment was 732.489 s. CPU-seconds must not be read as elapsed wall time or added across columns that share a baseline search.

Two case pipelines ran concurrently, with one numerical thread per worker, a five-minute stage watchdog and a 3,000 MiB worker memory guard. Largest recorded worker RSS was about 971 MiB. No final stage hit a resource guard. Some initial competition jobs failed at output persistence due to a missing directory in the new driver; the directory creation was corrected and those stages rerun. Completed searches were reused. Failed attempts are excluded from the successful-stage table and were not algorithmic search failures.

## Reproduction and evidence

`manifest.json` fixes the pilot cases, input identities and protocol. `results.json.gz` contains every final search, competition and decode record, including arm memberships and per-family certificates. `summary.json` contains the per-case comparison. `seed-control/64` and `paired-control/` contain the two diagnostic controls and their saved catalogues; paired timing records preserve all repeats.

The experiment uses the isolated engine from the adjacent `prefix_competition_20260915` report. The additional `paired.patch` changes only a separate diagnostic module, forcing the fixed starting trial. It does not change the ordinary experiment or production policy. Use a source checkout at commit `807d36f` (the pinned pre-optimization implementation). Install the repository dependencies plus pybind11, setuptools, psutil and z3-solver, then:

```sh
export MAPPA_EXPERIMENT_WORK=/tmp/mappa-sweep-prefix
python prepare.py --repo /path/to/coordinate_alignment
python run.py
python summarize.py
python seed_control.py 64
python paired_control.py
python paired_timing.py
```

The regular batch enforces per-stage resource guards. The two diagnostic scripts are small direct runs and do not install a separate hard process watchdog. All controls use compressed matching states; no arbitrary overlap subsets or full-bijection expansion was introduced.
