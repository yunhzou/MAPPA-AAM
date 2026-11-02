# Golden 1561 forward: decoder memory

This is the pre-refactor diagnosis. The implemented correction and comparison
are recorded in [decoder-refactor.md](decoder-refactor.md).

The Mac replay reproduced the publication failure under a 2 GiB process guard.
The dominant growth occurs while materializing all family event tables and their
plain output records. Serialization then adds enough memory to cross the guard.
This diagnostic did not change the scientific operations or production results.

The saved catalogue is 819,715 compressed bytes. It contains 14,132 families,
11,899 branch groups and 14,132 input paths across 192 cuts, with 186 reactant
atoms and 73 product atoms. The hierarchy branch cap does not cap the combined
catalogue across all cuts. These counts are distinct from symmetry-related
bijections inside a family.

| Measured point on Mac | Process RSS |
|---|---:|
| Catalogue loaded | 208.1 MiB |
| Action annotations prepared | 209.9 MiB |
| 14,000 family tables and their records built | 1.317 GiB |
| Complete record, with live decoding tables | 1.515 GiB |
| Decoding operation returned | 1.508 GiB |
| Last completed serializer sample | 1.753 GiB |
| Guard terminated serializer | 2.0004 GiB |

The output contains **14,132 family/partition tables**, **1,006,967 bond terms**,
and **1,512,124 atom-domain entries**. Only **13 bond-response objects** are
stored in the shared response section. Each exported bond term is a Python
dictionary; those dictionaries alone occupy **273,895,024 bytes / 261.2 MiB**,
excluding their values and all other containers.

`EventInterpretation.tables()` retains every constructed table in `_tables`.
Its key includes the complete family's mapping, actions, required edges, fragment
partition and cut annotations. It shares identical complete table contexts and
bond responses; it does not store every repeated local bond term once globally.
For this case, all 14,132 family table contexts remain distinct.

`decoded_events_record()` walks every family and builds lists/dictionaries for
the whole saved record while the interpretation still retains the live tables.
It also creates a representative event description for every family. At this
point both the working tables and exported bond-term dictionaries exist. Returning
from the operation reduces RSS only slightly; RSS does not measure which freed
objects the allocator has returned to the operating system.

The serializer writes into the same bounded 8 MiB spool used by production, with
gzip level 1 and pickle protocol 5. The output byte buffer is therefore bounded,
but the serializer's object bookkeeping grows with the entire Python object
graph. Its measured `Pickler.__sizeof__()` grew from 0.5 MiB through 2, 4, 8, 16,
32, 64 and **128 MiB** before termination. Process RSS rose from 1.508 GiB to
1.753 GiB by that last completed sample and then crossed 2 GiB. The process was
killed before the next serializer sample; this measurement does not assign all
additional RSS to the currently live pickle bookkeeping table.

The diagnostic counted 344,944 projected endpoint-pair states across the tables.
These are local bond-endpoint projections, not enumerated complete mappings.
The observed footprint is consistent with retained per-family table structures,
their object-heavy export and serialization bookkeeping. It is not evidence that
the decoder expands all fragment bijections.

An inventory-only scan of the same saved catalogue found **252,442 fragment
occurrences across family partitions**, but **2,441 distinct structural fragment
pairs**, keyed by matching policy and the reactant/product atom sets. The most
reused pair appears 9,314 times. There are no coupled-branch entries in this
catalogue. These are acquisition-stage counts, not a proof that every occurrence
has identical conditional symmetry choices. The required collect-first decoder
ordering is now explicit in [design.md](../../design.md).

The next storage design should preserve every family, action correlation and
event response while avoiding simultaneous materialization of the complete live
tables and dictionary export, and bounding serialization bookkeeping. This report
does not implement a new storage format or lower a scientific cap. Streaming
compressed bytes alone cannot solve the object-materialization cost.

Reproduce with [profile_decoder.py](profile_decoder.py) and the source range/hash
in [decoder-memory-source.json](decoder-memory-source.json). The profile calls the
current decoder and export function, wrapping table creation and serializer writes
only to record RSS and counts. A parent samples RSS every 20 ms and kills the
child above 2 GiB. It does not retain the serialized result. Measurements are in
[decoder-memory-mac.jsonl](decoder-memory-mac.jsonl). Python/OS and instrumentation
differ from the distributed run; this is a memory diagnostic, not a matched performance comparison.

Raw operational records and full artifacts are retained privately.
