# Remaining 2 GiB failures

The five investigated cuts fit below 2 GiB without changing their search results.
Two storage problems remain in the accepted native implementation: an oversized
canonical-certificate cache and repeated deferred-bond lists in transition records.
This is a diagnosis and controlled experiment, not a deployed production change.

Production revision: `3b43327462d11b587a51e61305436fd71374207d`.
Reaction 5735 forward, 387 atoms per endpoint, original cut/seed plan, branch
limit 100, one seed, random seed 42, matching tolerance 1.0. Measurements below
are from the user's Apple M4 Pro Mac, Python 3.12. Each process uses one compute
thread. The cluster failures identified the cuts; these are local replay peaks.

## Cause 1: the canonical-certificate cache

`DenseGraph::memo` in `src/mappa/mappa_core/native/engine.cpp:76` permits
500,000 entries. It has no byte budget. This is separate from the 64 MiB
`FragmentRepair` cache and the bounded native final-symmetry cache.

For this reaction the exact weighted target graph has 776 vertices after edge
subdivision. Each dense canonical certificate occupies 80,704 bytes. Each cache
entry stores its own certificate string, including when another key has the same
certificate value. These are symmetry-comparison certificates, not enumerated
fragment bijections.

Fresh cut 147 retains 17,171 entries: 1,385,905,752 bytes of certificate storage
and 72,813,248 bytes of keys, **1.359 GiB before map/allocator overhead**.
The cache records **zero hits** in that run. Only 7,359 certificate values are
distinct. Retaining the entries is optional reuse; the certificate computation
and exact comparisons remain necessary. The 500,000-entry bound could permit
37.6 GiB of certificate payload alone for this target.

The native matching cache is correctly near 64 MiB, and the native final-symmetry
cache is only about 3 MiB for this cut. They do not explain the large peak.

## Cause 2: repeated deferred-bond lists

`fragment_scheduler_binding.h:99` builds a new list and new pair lists for every
transition's `match['deferred_edges']`. Sharing mapping/island blocks did not
remove this other expansion.

Cut 147 has 61,755 transitions and only **37 distinct deferred-bond list values**.
The exporter nevertheless creates **3,745,341 pair-list occurrences**. Their
lists alone occupy 366,543,968 bytes, excluding integer objects. The associated
transition-export phase increases RSS by roughly another 0.55 GiB. Saving the
result adds about 0.50 GiB more. The writer already streams through an 8 MiB
spool; that spool does not bound the object graph or pickle's object-tracking
memory. The final compressed cut record is only 1.66 MiB.

## Controlled experiments

The diagnostic engine adds cache counters and a retention-limit setter to the
current source. It does not replace the installed engine. A constructor wrapper
can also share deferred-bond lists by exact value while retaining their existing
record shape. No branch, seed, matching constraint, transition or stop is removed.

For cut 147 alone, reducing certificate retention to 512 entries (about 41 MiB
including keys) reduces peak RSS from **3.136 to 1.785 GiB**. The complete finalized
search graph compares equal. The 512-entry setting is a diagnostic control, not
a proposed new scientific cap or final cache design.

Combining that control with shared deferred-bond lists gives these measurements:

| Cut | Current version, fresh process | Diagnostic changes, reused process | Complete finalized graph |
|---|---:|---:|---|
| 104 | 2.279 GiB | 0.939 GiB | Equal |
| 147 | 3.136 GiB | 1.039 GiB | Equal |
| 228 | 2.539 GiB | 1.056 GiB | Equal |
| 321 | 2.315 GiB | 1.019 GiB | Equal |
| 335 | 2.839 GiB | 1.059 GiB | Equal |

The diagnostic worker ran cuts **111, 104, 147, 228, 321, 335 consecutively**, using
one workspace and releasing each output after saving it. It did not explicitly
collect garbage or clear its caches between cuts. Peak RSS across the whole
sequence was **1.059 GiB**. Full graph equality includes contexts, states,
transitions, stops and compressed fragment records; all five retain 100 terminals.
This validation does not rerun decoding or reference verification.

The unchanged implementation also failed in fresh-process replays, so previous
work is not required to trigger the problem. Reuse aggravates it: unchanged cuts
111, 104 and 147 in sequence reached about 4.30 GiB. Explicit garbage collection
after the sequence did not remove the retained native cache. Allocator retention
also means releasing objects need not immediately return every byte of RSS.

Parent-process RSS sampling runs every 20 ms. The diagnostic guard is 8 GiB so
the baseline can finish; the production limit remains 2 GiB. Constructor probes
and graph-field inspection add overhead. These runs establish memory causality
and record equality, not a throughput guarantee or a universal bound for every
EnzymeMap reaction.

## Production repair indicated by the evidence

Give the canonical-certificate cache a real byte budget and an explicit reuse/
eviction policy. Share equal deferred-bond records during native export, as the
existing Python graph builder already does for bond records. Both changes belong
to storage/reuse code and must preserve search decisions and exact comparisons.
Validate repeated-worker throughput as well as complete-result equality before
deployment; limiting cache residency can cause additional recomputation.

No production source, installed engine, benchmark configuration, running job or
2 GiB worker limit was changed during this diagnosis.

Evidence: [memory-followup.json](memory-followup.json), including raw phase
measurements, cache contents and all six full-graph comparison assertions.
Reproduction uses [memory_probe.py](memory_probe.py),
[instrument_memory.py](instrument_memory.py), and the existing
`mac-point-input.pkl.gz` fixture. Build the instrumented engine into a temporary
directory, then supply its path explicitly with `--instrumented-engine`.
Use `--memo-entries 512 --share-deferred` only for the controlled experiment.
