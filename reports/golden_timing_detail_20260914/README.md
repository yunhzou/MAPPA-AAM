# Golden timing detail

No mapping searches or reference checks rerun. This audit reads scalar timers from the saved1821common-Mac reactions, checks each MAPPA direction against its saved search execution, and exactly reproduces published bidirectional CPU totals. Direction host fields can refer to later reference verification; timing uses search-host provenance, not scorer location. Endpoint sizes/ties follow the published size-orientation evidence.

SLAP saved records are checked for complete, uniquely ordered cut schedules and successful mapping statuses; graph/mapping/export CPU sums reproduce the published per-case totals. Each source records file is hashed. Its uncut timing uses ordinal0 in both modes/directions. MAPPA uncut search was not independently instrumented, and rescoring time is not substituted. No mixed-host3/10seed timings are compared.

Table4 and Figure2 report mean, median and linearly interpolated95th-percentile CPU cost per reaction, with full1851-reference recovery. Original-cluster default-comparator timings are separately derived from successful-call aggregates, including invalid outputs but excluding exceptions. They exclude startup, queueing, scoring and saving; Chython uses its default8ONNXthreads, so recordedCPUcan exceedwall. Do not compare them to Mac sweep timings as a controlled hardware ranking.

Run derive.py after adapting its repository and saved-run roots to reconstruct the audit.

Coordinate per-reaction decoding times include all140first-pass records and both continuation records; theirCPU sum exactly matches1542.019335s in the published campaign total. Run coordinate_timing.py after derive.py to add this audit. Each result file is hashed. The distribution is over140reaction totals, not142attempts.
