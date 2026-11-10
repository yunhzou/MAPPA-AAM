# Ground-truth-containing branch event offsets

The requested metric uses the saved decoder model on both sides:

`n = minimum of GT-containing retained branches - overall detected minimum`.

Ground truth identifies branch membership; its bond events are not counted from
labels. A membership witness's event count is an upper bound on its branch's
minimum, not the plotted branch minimum. Both directions contribute, with the
smallest qualifying count retained once per reaction.

Scope: one seed, discovery tolerance 2, branch cap 100, and the saved decoded
minimum + 5 responses for 1,851 Golden and 47,974 EnzymeMap reactions. Membership
uses the exact annotated heavy-atom indices and unmatched-atom policy of the
independent decoded audit. The event counts themselves use the same original
all-atom decoder model, including its documented boundary approximation.

Existing replayed certificates at the overall minimum immediately establish
zero offset. For other recovered cases, existing context bounds discard branches
that cannot improve the result. Membership identifies a GT-containing context,
then an unconditioned event query locates its minimum. Queries consume only saved
responses; matching and decoding outputs remain unchanged. Unknown queries retain
intervals and stay outside exact bins. Missing GT remains a separate category.

The previous mixed-scope figure was withdrawn. It subtracted heavy-only GT counts
from all-atom minima and did not measure missed minima. It must not be cited.

Analysis revision: `946fc6a`, `bench/analytics.py --positions`. Fifteen affected
tests passed, including a GT witness with two events in a zero-minimum branch,
and a lower branch that does not contain GT.

## Published distribution

[Figure](event-offsets.png) · [PDF](event-offsets.pdf) ·
[Numeric data](event-offsets.json.gz) · [Histogram bounds](histogram.csv)

The later saved analysis processed all 49,825 reactions. The manuscript now uses
that result, combined with the final Golden membership closure, rather than the
earlier in-progress snapshot.

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
extend below zero, but no exact measured offset is negative. The plotted caps
retain this uncertainty without inventing bin assignments.

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
