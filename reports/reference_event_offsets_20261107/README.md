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

## Preliminary figure

[Figure](event-offsets.png) · [PDF](event-offsets.pdf) ·
[Numeric data](event-offsets.partial.json.gz) · [Histogram](histogram.partial.csv)

Snapshot updated after the 10 November verification closure: 45,738 exact offsets, with 3,891 further
branch-minimum checks queued. All 49,825 reactions are accounted for once.

| Dataset | Exact offsets | GT not represented | Membership unknown | Offset pending | Total |
|---|---:|---:|---:|---:|---:|
| Golden | 1,621 | 18 | 0 | 212 | 1,851 |
| EnzymeMap | 44,117 | 178 | 0 | 3,679 | 47,974 |

The established Golden offsets are 1,621 at zero. Established EnzymeMap offsets
are 44,085 at zero, 22 at one, 6 at two, 1 at three and 3 at four. Existing
certificates at the detected minimum establish zero immediately, so these
completed results deliberately overrepresent zero. The pending cases must be
resolved before interpreting the shape of the full distribution or reporting
its percentages. Blank bins are not evidence that no pending case belongs there.

The continuation retains 114 completed branch lookups and uses the unchanged
saved responses, with a 30-second query allowance per reaction. Resource-limited
queries will remain explicit intervals outside exact bins. Molecular data,
membership certificates and operational launch evidence remain private.

The [final Golden membership closure](../golden_native_seeds_20261106/decoded-verification-closure.json)
settles all three previously unknown memberships. Case 852 has a replayed
27-event reference witness; both saved directional minima are 27, proving zero
branch offset without another minimum search. Cases 1568 and 1606 are absent.
This updates membership and one zero-offset certificate only; the 3,891 pending
branch-minimum lookups remain incomplete, outside the plotted bins.
