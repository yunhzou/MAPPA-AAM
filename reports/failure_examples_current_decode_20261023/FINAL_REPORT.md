# MAPPA supplied-example decoding — final report

Both supplied examples completed using the current MAPPA search and signed-event decoder. All retained families were decoded with no event-count cutoff, and every emitted witness passed independent raw-WBO rescoring and a producing-family membership audit. The legacy orientation-selection workflow was not used.

| Supplied case | Atoms | Completed families | Decoded event classes | Minimum bond-change events |
|---|---:|---:|---:|---:|
| PR22 | 40 | 7 / 7 | 1 | 0 |
| PMC11894599_10 | 80 | 27,295 / 27,295 | 23,956 | 24 |

PR22 has a complete zero-event mapping. PMC11894599_10 has 12 tied minimum event classes, each with **13 WBO-decrease and 11 WBO-increase events**. There is no zero-event class in the retained search. These events use changes in the supplied WBO values, with thresholds of 0.5 for ordinary pairs and 0.3 for metal-containing pairs, including explicit hydrogens.

## Configuration and scope

The reference and prediction NPZ inputs were unchanged and verified by SHA256. Search used three seeds, branch limit 100, isomorphism tolerance 1.0, and graph/cut floors of 0.15. PR22 used the original uncut configuration. The 80-atom case used all 91 cut-sweep checkpoints, consumed individually to avoid constructing a combined search-history graph.

The 80-atom search was branch capped. Its minimum of 24 is therefore proven over the retained families under these settings; it is not a global impossibility claim over all atom mappings or other search configurations. This report establishes bond-event matching results, without a 3D-orientation or stereochemical verdict.

## Verification and execution

- All 23,957 output witnesses were independently rescored against the original inputs.
- Every witness was checked in its recorded producing family. For the 80-atom case, 19,550 witnesses matched the validated family representative, and 4,406 recorded symmetry-action certificates passed exact replay and allowed-action checks. All 12 minimum witnesses were verified.
- Serial and partitioned decoding agreed on all 23,956 event-class IDs, event counts, and supporting-family IDs. The final partitioned-run witnesses passed the independent audits above; representative choices can differ within the same event class.
- Search and decoding workers retained the original **5 GiB virtual-address-space limit**. Parallel decoding used up to 16 workers; this is a per-process limit, not a total-memory claim.
- PR22 completed on the cluster in 2.30 seconds with about 108 MiB peak resident memory.
- The 80-atom search took 139.38 seconds. Its unrestricted serial decode took 771.21 seconds. The 32 partitioned tasks each completed within 74.64 seconds, including catalogue loading and report writing; at most 16 ran together. Peak resident memory per decoder worker was 3843 MiB.

## Evidence

The frozen source is recorded in [source_sha256.json](source_sha256.json). Machine-readable final results are in [final-results.json](final-results.json). The [decoding evidence archive](decoding-evidence.tar.gz) contains the cut checkpoints, catalogues, decoder outputs, per-family reports, and execution provenance. The final witness audits are retained separately with the extracted results.

Final 80-atom audits: [independent WBO rescore](cluster-evidence/results/search_memory_partitioned/independent-audit.json), [producing-family membership](cluster-evidence/results/search_memory_partitioned/constructive-membership-audit.json), and [serial/parallel comparison](cluster-evidence/serial-parallel-equality.json).

## Viewer and additional configuration

The [single-mapping viewer](viewer_80atom/view.html) distinguishes bond disappearance/appearance from weakening/strengthening using the drawing's WBO > 0.2 bond-presence cutoff. For the displayed witness, the 24 signed events comprise 11 breaking, 10 forming, 2 weakening, and 1 strengthening event. Its display filter leaves the full molecules visible. This does not change the decoder's signed-event score.

The separately requested [cap-10 / threshold-1.0 measurement](cap10_threshold1/REPORT.md) took 1.49 seconds but retained no complete mappings; it is not a successful faster decode.

The exact benchmark engine is preserved in the source archive; the publication branch retains GitHub main's newer TS scorer. The [original input arrays](inputs/) and [publication checks](publication-validation.json) are included.
