# Investigation of the 2 GiB worker failures

Two distributed diagnostic runs replayed exact archived failed tasks against
production algorithm revision 56c4b6c. Only the diagnostic script was installed;
production algorithm files and benchmark parameters were not changed.
Each selected task ran in a fresh spawned process, then repeated in that process.
The parent sampled RSS every 50 ms, with a 12 GiB diagnostic guard. The probe
records growth, native export, Python graph construction, finalization and
publication boundaries. It explicitly collects Python garbage after releasing
each output; first-execution measurements do not depend on this cleanup.

| Reaction/cut | Atoms per endpoint | States | Mapping pairs | Island pairs | First-execution peak |
|---|---:|---:|---:|---:|---:|
| 5735 forward / 111 | 387 | 24,758 | 9,020,201 | 9,020,201 | 2,454,454,272 bytes (2.286 GiB) |
| 6711 forward / 190 | 378 | 25,686 | 9,006,324 | 9,006,324 | 2,441,605,120 bytes (2.274 GiB) |

Both started near 50 MiB and crossed 2 GiB on their first cut. Both completed
successfully with the diagnostic guard. Accumulation from prior cuts is therefore
not necessary to reproduce these failures. This does not establish the cause of
every other failed cut or exclude retained-memory effects in long-lived workers.

For 5735, RSS was 2.254 GiB after native graph export and 2.286 GiB after Python
graph construction. It dropped to approximately 0.939 GiB when the temporary
export record was released. Finalization and disk publication did not create
the offending peak. The published compressed cut artifact was 1,125,601 bytes.

`scheduled_graph_record` in `fragment_scheduler.h` builds fresh Python pairs
for every full mapping/island snapshot in every state. `AAMSearchGraph.from_record`
then interns repeated pairs while the expanded input record is still present.
This is a demonstrated memory-expensive graph representation boundary. A live
branch limit of 100 does not cap the number of historical states in the DAG.
The exporter lines originate in commit e923b248 (September 9), preceding the
persistent-pool refactor. The refactor's 2 GiB guard exposes this existing peak.

Next optimization to validate: share repeated immutable pairs during native
export and reduce overlapping temporary graph representations, preserving every
state, transition and scientific alternative. Do not lower branch/search caps or
silently discard provenance to reduce memory. No such optimization was applied
as part of this investigation.

The optimized transport was subsequently validated and accepted; see
[the Mac comparison](../native_transport_20261031/mac-point.md). The old-export
diagnostic `probe.py` is retained in Git at `246aa3d` and removed from the current
tree because its expanded-record interface no longer exists.

Raw operational records and full artifacts are retained privately.
