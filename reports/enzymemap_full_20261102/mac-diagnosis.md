# Reaction 37270: local reproduction of the 2 GiB failures

Both sampled failures reproduce in fresh processes on the Mac using current
revision `beed5b5b90e8989cc2bbff4fecf178cbff51e6b2`. Cut 432 exceeds the limit
while exporting matching history to Python, before symmetry finalization.
Cut 410 finishes computation below the limit and exceeds it while saving.
No decoding is involved.

The exact endpoint graphs, policies and seed orders came from the archived
tasks. Both endpoints contain 661 atoms. Regenerated plans matched all 28 failed
tasks on cut, number, seeds, indices and orders; only these two forward cuts were
executed locally. Each run used a fresh process with one compute thread.

| Forward cut | History states | RSS after finalization | Peak through saving | Matching + finalization | Saving |
|---|---:|---:|---:|---:|---:|
| 432 | 120,139 | 3.24 GiB | 4.23 GiB | 17.67 s | 2.57 s |
| 410 | 67,124 | 1.81 GiB | 2.30 GiB | 12.71 s | 1.15 s |

The table uses completed runs with an **8 GiB local diagnostic ceiling**. Separate
2 GiB runs reproduced both guard failures. The cluster remains at 2 GiB; no jobs
or production code changed. These are point measurements, not a speed comparison.

## What occupies the memory

Cut 432 has 100 completed terminals and 111,373 `combined_live_leaves` cap stops.
The scheduler records each proposed child before admitting it to the next live
frontier. Rejected children retain their history and stop evidence. Thus the
100-branch cap does not bound the total number of historical records.

An isolated observational build measured **1.09 GiB of native history storage**,
including 952,939,904 bytes of full-length state mapping/label and candidate-image
arrays. Native history remains present while its Python representation is built.

Python state fields already share data. Transition deferred-edge lists do not:
120,138 separate lists contain 9,910,765 pair occurrences, but only **47 distinct
sequences**. In the separately reloaded checkpoint, their list containers occupy
**1,276,613,768 bytes (1.19 GiB)** before counting atom integers. Reloaded lists
can have different allocation sizes from the original export, so this is not an
exact allocation measurement of that export. These are repeated deferred-edge
records, not millions of fragment bijections. Independently, the measured
transition-export phase adds approximately 1.46 GiB RSS.

Saving raises RSS again by about 0.99 GiB for cut 432. The compressed output is
only 3,864,653 bytes: a small compressed file does not bound the memory needed to
traverse and serialize its object graph. The writer's 8 MiB spool bounds its byte
buffer, not pickle's object bookkeeping. The measurement isolates the saving
interval; it does not separately measure every serializer allocation.

The completed-fragment cache stayed within its **64 MiB** budget. It cannot
account for the multi-GiB increase identified above. Shared snapshot loading and
saving were omitted, so neither is required to trigger these failures. Native
conditioned-symmetry cache counters are also bounded, but do not measure every
Python finalization object. RSS retained after releasing objects is not, by
itself, proof of a live-object leak.

The concrete repair target is history representation and transport: share
repeated immutable data while retaining scientific choices and cap evidence.
This diagnosis does not change matching, establish full-reaction accuracy, or
show that all 28 failed cuts have the same allocation profile.

## Evidence and reproduction

[Measurements and provenance](mac-diagnosis.json) contain all five phase traces,
public resource metrics, fixture hashes, native allocation accounting and graph
storage counts. Graph inspection ran separately from the timed measurements.
The isolated build preserved state/transition counts and cache statistics; it did
not replace the installed engine. Its timing is not used for the main table.

The existing [phase probe](https://github.com/yunhzou/MAPPA-AAM/blob/f10999c75d7b232e8996cedf996f4d9b3a5e94ef/reports/native_transport_20261031/memory_probe.py) was reused,
removing allocating graph inspections from its timed path. The existing
[instrumentation builder](https://github.com/yunhzou/MAPPA-AAM/blob/f10999c75d7b232e8996cedf996f4d9b3a5e94ef/reports/native_transport_20261031/instrument_memory.py) now
counts current native history rather than retired certificate-cache fields.
No new production diagnostic mode or algorithm was introduced.

From the checkout, with the project environment and the archived local fixtures:

```sh
PYTHONPATH=src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
# Historical command; use the linked script and measured source revision.
python reports/native_transport_20261031/memory_probe.py \
  out/enzymemap_37270_mac_20261102/37270-forward-compute.pkl.gz \
  out/enzymemap_37270_mac_20261102/reproduce432 \
  --cuts 432 --guard-gib 8 --seconds 300
```

Use `37270-forward-publication410.pkl.gz` and `--cuts 410` for the other cut;
use `--guard-gib 2` to reproduce guard termination. Fixtures are trusted internal
checkpoints, identified by hashes and original bundle references in the evidence.
The parent samples child RSS every 20 ms; sampled guard values can overshoot.
Diagnostic scripts received syntax and diff checks. No core source changed, so
the scientific test suite and Golden benchmark were not rerun.

Raw operational records and full artifacts are retained privately.
