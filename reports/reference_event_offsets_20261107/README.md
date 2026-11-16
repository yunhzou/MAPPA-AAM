# Verified reference-mapping event offsets

The current plot compares concrete saved counts:

`n = verified GT-consistent mapping event count - original detected minimum`.

[Figure](event-offsets.png) · [PDF](event-offsets.pdf) ·
[Per-reaction evidence](../../manuscript/evidence/reference_offsets.json.gz) ·
[Histogram](histogram.csv)

Use the lowest recorded verified count across both directions, with forward
winning ties. Both terms use the same all-atom decoder model. The 196 references
absent from the saved output use their previously audited anchored completions.
Every one of the 49,825 reactions contributes a concrete difference. This update
runs no matching, decoding or additional optimization. Independent coverage and
the original detected minima are unchanged.

The previous figure instead minimized GT-containing branches and inserted lower
bounds for unfinished minimum queries. That is a different quantity: a branch
can attain its minimum with a mapping that does not agree with GT. The corrected
figure moves 126 Golden and 2,196 EnzymeMap reactions to higher offsets. Golden
case 268 now lies at +4, from its verified count of 16 and detected minimum of 12.
The archived branch-minimum measurements below are not the plotted metric.

## Superseded branch-minimum analysis

[Archived branch-minimum data](event-offsets.json.gz) retain the original analysis.

The superseded analysis processed all 49,825 reactions, including the final
Golden membership closure.

| Dataset | Exact offsets | Bounded offsets | Reference absent | Total |
|---|---:|---:|---:|---:|
| Golden | 1,823 | 10 | 18 | 1,851 |
| EnzymeMap | 47,686 | 110 | 178 | 47,974 |

All reference memberships are resolved. The 120 bounded offsets are completed
queries whose 30-second allowances did not establish a single event count.
They are not pending work and are not relabeled as exact offsets.

Golden exact counts at offsets 0 through 5 are 1,722, 60, 28, 10, 1 and 2.
EnzymeMap exact counts at offsets 0 through 4 are 45,888, 960, 785, 10 and 43.
For each bin, the lower count is its established count; the upper count adds
each interval containing that bin. Upper counts across bins are not jointly
attainable and must not be summed as a histogram. Three EnzymeMap intervals
extend below zero, but no exact measured offset is negative. This original
bin-count analysis remains in the evidence; Figure 7 now uses verified mapping
counts as documented below.

The earlier partial snapshot and histogram remain as historical evidence. They
must not be used for the current manuscript figure.

Final analysis source SHA-256:
`1f5133f0f44b0bb503207c2e1a969a8bf85a9233d955cad3c5d74b237c618953`.
The numerical snapshot records the producing code identity and runner hash.
The [Golden membership closure](../golden_native_seeds_20261106/decoded-verification-closure.json)
settles case 852 at zero and excludes cases 1568 and 1606. No matching or decoding
was rerun for this publication update.

## Anchored hydrogen-completion trial

EnzymeMap case 13175 was rerun locally with all nine reference heavy-atom pairs
fixed as search anchors. Six hydrogen assignments were left unspecified. The
existing native API ran matching, reduction and decoding in both directions:
one seed, branch cap 100, discovery tolerance two, full single-edge sweep, and
minimum-plus-five decoding. Matching used the saved indexed elements and bond
matrices; coordinates were zero placeholders and do not enter these operations.

Both directions returned complete injective 15-atom mappings preserving every
anchor. Each has six directly counted bond events: four between heavy atoms and
two involving hydrogen. The hydrogen contribution is one O–H bond broken and
one O–H bond formed. Model and direct physical counts agree for these witnesses.
An independent check of all 720 hydrogen permutations, with the heavy mapping
fixed, confirms the six-event minimum; four assignments attain it.

Matching including reduction consumed 0.033438 CPU seconds across both directions.
The originally reported 0.001093 CPU seconds measured decoder construction only;
lazy minimum/window preparation was outside that timer, so it is not a full
decoding measurement. Interpreter startup and imports are excluded.
The independent 720-assignment check took 0.018340 CPU seconds and is
not part of MAPPA runtime. This is one small-case feasibility trial, not a timing
estimate for the full datasets.

[Trial evidence](anchored-completion-trial.json) records both mappings, event
lists, source identity, settings and stage timings. This completes a supplied
reference; it does not measure independent reference recovery. The manuscript
and benchmark coverage are unchanged.

### One seed per direction without a cut sweep

Repeating case 13175 with the same anchors and settings but `sweep_cuts=False`
produced one retained family per direction. Each direction ran only the uncut
condition and returned a complete 15-atom mapping preserving all nine heavy-atom
anchors. Both still attained six events: four heavy-atom events and two hydrogen
events, matching the independent exhaustive minimum and the earlier sweep trial.

| CPU time, both directions summed | Seconds |
|---|---:|
| Matching, including reduction | 0.005720 |
| Decoding, including minimum/window preparation | 0.008468 |

Witness retrieval took a further 0.001531 seconds of elapsed time in total.
These are single local measurements, excluding imports and startup. The sweep
was unnecessary for this case; one case does not establish that every anchored
completion reaches its minimum without it. The existing evidence file retains
the no-sweep configuration, mappings, direct event checks and separate timings.

## Anchored completion of all missing references

The user selected the 196 reactions whose reference was absent from the saved
decoded results: 18 Golden and 178 EnzymeMap. The 120 bounded branch-offset cases
were excluded. All selected reactions were run with supplied reference pairs as
anchors, one seed per direction, no sweep, tolerance 2, cap 100 and a
minimum-plus-five decoding window. The scientific MAPPA implementation was
unchanged. The portable runner is `bench/anchored_completion.py`, revision
`8a773c6`.

| Dataset | Reactions completed | Directions audited | Mean matching CPU s/reaction | Mean decoding CPU s/reaction |
|---|---:|---:|---:|---:|
| Golden | 18/18 | 36/36 | 0.0325 | 0.0724 |
| EnzymeMap | 178/178 | 356/356 | 0.0305 | 0.0632 |

Every supplied anchor was preserved. Every direction returned an injective,
element-preserving correspondence with maximum element-compatible atom usage.
All 178 EnzymeMap reactions map both endpoints completely. Eleven Golden
reactions have unequal endpoint element inventories and retain unmatched atoms;
none leaves an additional element-compatible assignment unfilled. Completion
adds previously unspecified heavy-atom pairs in Golden cases 7, 1358, 1377 and
1475; these additions are recorded explicitly.

Eleven EnzymeMap references use atom numbering different from the prepared
input: 6271, 6294, 43409, 43410, 43411, 43412, 43454, 43460, 43470, 43477 and
46384. Their saved index-agreement flag was false. Before anchoring, endpoint
graph isomorphisms translate those reference indices, preserving elements,
charges, isotopes, hydrogen counts and exact bond orders. Stereo labels are
outside MAPPA's represented graph. Both mappings and the complete index
translations are retained. These cases require a separately corrected audit
before interpreting their previous exact-index misses as algorithmic failures.

All 392 selected witnesses pass joint-choice replay and independent scalar event
counting from the original matrices. Their physical counts agree with the
decoder model counts. The decoder establishes its retained model minimum in
every direction; this is not an exhaustive optimum over all possible searches.
No task failed or exhausted its allowance.

The 16-worker pool ran for 6.088 seconds, including pool startup and shutdown;
the allocation elapsed time was 33 seconds. Summed stage CPU times were 6.014
seconds for matching, including reduction, 12.560 seconds for full decoding,
including minimum/window preparation, and 7.846 seconds for witness retrieval
and audit. The longest direction's decoding took 0.638 CPU seconds. Input
loading, reference-index translation and result publication are outside those
stage timings. Each direction had a 180-second allowance and a 4 GiB memory
guard; actual completion was well below the time allowance.

[Audited completion results](anchored-completion-results.json.gz) retain all
392 mappings, event lists, translations, input checksums, settings, source
identity and timings. This completes supplied references without increasing independent benchmark
recovery. These concrete counts also enter the corrected offset figure.

### Corrected figure: saved GT-consistent mapping counts

[Manuscript Figure 7](../../manuscript/figs/fig7_reference_offsets.png) includes
all 1,851 Golden and 47,974 EnzymeMap reactions, with one color per dataset and
no separate completion legend. Counts at zero are 1,630 (88.06%) and 44,086
(91.90%), respectively. The histograms are:

- Golden: `{0:1630, 1:90, 2:88, 3:24, 4:11, 5:6, 7:1, 12:1}`.
- EnzymeMap: `{-2:9, 0:44086, 1:1384, 2:2190, 3:122, 4:171, 5:9, 6:3}`.

The nine negative EnzymeMap offsets are attained anchored completions with two
fewer events than the original detected minimum. The old two lower bounds at -2
and one at -1 no longer enter the figure. Golden has no negative ticks; the empty
-1 tick is also omitted from EnzymeMap.

[Plot evidence](../../manuscript/evidence/reference_offsets.json.gz) records each
reaction's selected count, unchanged baseline, direction and certificate identity.
[Anchored evidence](../../manuscript/evidence/anchored_offsets.json) retains both
completion directions and index-translation details. The selected hydrogen
assignment is a saved feasible completion, not a new optimized assignment.
The build validates full scope, disjoint original/anchored membership, saved
baseline consistency and Golden certificate counts before rendering the figure.
