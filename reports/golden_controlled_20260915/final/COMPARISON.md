# Controlled Golden comparison

CPU model: AMD EPYC 9J14 96-Core Processor. One pinned thread per search. All recovery percentages use 1,851 reactions. The paired completed-search timing cohort contains 1403 identical reactions for all configurations.

Search workflow CPU includes graph preparation, mapping and native output persistence, and excludes reference verification and initial process imports. MAPPA retains raw compressed graphs; SLAP retains native label information and explicit candidates. Separate selected-attempt CPU means include resource-limited calls and describe spent computation, not completed-search latency. See CSV for denominators and completed-case counts.

| Configuration | Recovered | Unknown | Paired mean CPU s/reaction | Selected-attempt mean CPU s/reaction |
|---|---:|---:|---:|---:|
| mappa_c100_s1_uncut | 1489/1851 (80.44%) | 0 | 0.040 | 0.085 (n=1851) |
| mappa_c100_s1_sweep | 1834/1851 (99.08%) | 0 | 0.963 | 3.873 (n=1851) |
| mappa_c100_s2_sweep | 1834/1851 (99.08%) | 0 | 1.619 | 7.260 (n=1851) |
| mappa_c100_s3_sweep | 1837/1851 (99.24%) | 0 | 2.265 | 10.556 (n=1851) |
| mappa_c100_s10_sweep | 1840/1851 (99.41%) | 0 | 6.993 | 33.211 (n=1851) |
| mappa_c2000_s1_uncut | 1508/1851 (81.47%) | 5 | 0.086 | 0.517 (n=1851) |
| mappa_c2000_s1_sweep | 1834/1851 (99.08%) | 2 | 4.317 | 30.055 (n=1851) |
| mappa_c2000_s2_sweep | 1834/1851 (99.08%) | 4 | 8.283 | 48.549 (n=1851) |
| mappa_c2000_s3_sweep | 1837/1851 (99.24%) | 6 | 12.201 | 61.890 (n=1851) |
| mappa_c2000_s10_sweep | 1835/1851 (99.14%) | 12 | 43.843 | 113.777 (n=1851) |
| slap_uncut | 1661/1851 (89.74%) | 1 | 0.105 | 0.659 (n=1851) |
| slap_sweep | 1794/1851 (96.92%) | 5 | 2.793 | 9.895 (n=1851) |

Selected attempts exclude 16 superseded checkpoint-conflict attempts and the initial pilot. They measure retained-attempt search expenditure, not total campaign CPU. MAPPA timings precede fragment competition. Unknown verification is not evidence that a reference is absent.
