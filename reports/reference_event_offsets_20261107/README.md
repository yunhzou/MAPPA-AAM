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
