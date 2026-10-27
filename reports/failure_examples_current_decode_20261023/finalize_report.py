"""Publish final results only after complete decoding and witness audits pass."""
from pathlib import Path
import gzip,hashlib,json

root=Path(__file__).resolve().parent
forty=root/'cluster/rp_timeout'
eighty=root/'cluster-evidence/results/search_memory_partitioned'
a=json.loads((forty/'independent-audit.json').read_text())
b=json.loads((eighty/'independent-audit.json').read_text())
qa=json.loads((forty/'producing-family-audit-v3.json').read_text())
qb=json.loads((eighty/'constructive-membership-audit.json').read_text())
eq=json.loads((root/'cluster-evidence/serial-parallel-equality.json').read_text())
assert all(x['complete'] and x['all_witnesses_independently_rescored'] and x['unrestricted_event_window'] for x in (a,b))
assert qa['status']==qb['status']==eq['status']=='passed'
assert qa['witnesses']==a['candidate_count'] and qb['witnesses']==b['candidate_count']
assert all(q['status']=='allowed' for q in qa['queries'])
assert qb['identity_action_witnesses']+qb['constructive_action_certificates']==qb['witnesses']
assert qb['minimum_witnesses_verified']==12
data=json.load(gzip.open(eighty/'decoded.json.gz','rt'))
mins=[c for c in data['candidates'] if c['total']==b['minimum_events']]
breakdown=sorted(set((len(c['events']['broken']),len(c['events']['formed'])) for c in mins))
assert breakdown==[(13,11)]
perf=json.loads((root/'cluster-evidence/results/search_memory/failure.json').read_text())['timings']
metrics=json.loads((eighty/'merge-validation.json').read_text())
summary=dict(status='complete',cases=[a,b],family_membership_audits=[{k:v for k,v in x.items() if k!='queries'} for x in (qa,qb)],serial_parallel_comparison=eq,minimum_80_atom_classes=len(mins),minimum_80_atom_broken_formed=breakdown[0],geometry_selection=False,scope='All event counts in all retained complete families under the supplied search settings. The 80-atom search is branch capped, so its minimum is not a universal lower bound over all atom bijections.',search_seconds_80=perf['search'],serial_decode_seconds_80=perf['decode'],parallel_worker_max_seconds=max(x['seconds'] for x in metrics['shards']),parallel_worker_peak_rss_mib=max(x['peak_rss_mib'] for x in metrics['shards']))
(root/'final-results.json').write_text(json.dumps(summary,indent=2)+'\n')
report=f'''# MAPPA supplied-example decoding — final report

Both supplied examples completed using the current MAPPA search and signed-event decoder. All retained families were decoded with no event-count cutoff, and every emitted witness passed independent raw-WBO rescoring and a producing-family membership audit. The legacy orientation-selection workflow was not used.

| Supplied case | Atoms | Completed families | Decoded event classes | Minimum bond-change events |
|---|---:|---:|---:|---:|
| PR22 | 40 | {a['complete_family_count']:,} / {a['family_count']:,} | {a['candidate_count']:,} | {a['minimum_events']} |
| PMC11894599_10 | 80 | {b['complete_family_count']:,} / {b['family_count']:,} | {b['candidate_count']:,} | {b['minimum_events']} |

PR22 has a complete zero-event mapping. PMC11894599_10 has {len(mins)} tied minimum event classes, each with **13 WBO-decrease and 11 WBO-increase events**. There is no zero-event class in the retained search. These events use changes in the supplied WBO values, with thresholds of 0.5 for ordinary pairs and 0.3 for metal-containing pairs, including explicit hydrogens.

## Configuration and scope

The reference and prediction NPZ inputs were unchanged and verified by SHA256. Search used three seeds, branch limit 100, isomorphism tolerance 1.0, and graph/cut floors of 0.15. PR22 used the original uncut configuration. The 80-atom case used all 91 cut-sweep checkpoints, consumed individually to avoid constructing a combined search-history graph.

The 80-atom search was branch capped. Its minimum of 24 is therefore proven over the retained families under these settings; it is not a global impossibility claim over all atom mappings or other search configurations. This report establishes bond-event matching results, without a 3D-orientation or stereochemical verdict.

## Verification and execution

- All {a['candidate_count']+b['candidate_count']:,} output witnesses were independently rescored against the original inputs.
- Every witness was checked in its recorded producing family. For the 80-atom case, {qb['identity_action_witnesses']:,} witnesses matched the validated family representative, and {qb['constructive_action_certificates']:,} recorded symmetry-action certificates passed exact replay and allowed-action checks. All 12 minimum witnesses were verified.
- Serial and partitioned decoding agreed on all 23,956 event-class IDs, event counts, and supporting-family IDs. The final partitioned-run witnesses passed the independent audits above; representative choices can differ within the same event class.
- Search and decoding workers retained the original **5 GiB virtual-address-space limit**. Parallel decoding used up to 16 workers; this is a per-process limit, not a total-memory claim.
- PR22 completed on the cluster in {json.loads((forty/'result.json').read_text())['wall_seconds']:.2f} seconds with about {json.loads((forty/'result.json').read_text())['peak_rss_mib']:.0f} MiB peak resident memory.
- The 80-atom search took {perf['search']:.2f} seconds. Its unrestricted serial decode took {perf['decode']:.2f} seconds. The 32 partitioned tasks each completed within {summary['parallel_worker_max_seconds']:.2f} seconds, including catalogue loading and report writing; at most 16 ran together. Peak resident memory per decoder worker was {summary['parallel_worker_peak_rss_mib']:.0f} MiB.

## Evidence

The frozen source is recorded in [source_sha256.json](source_sha256.json). Machine-readable final results are in [final-results.json](final-results.json). The [decoding evidence archive](decoding-evidence.tar.gz) contains the cut checkpoints, catalogues, decoder outputs, per-family reports, and execution provenance. The final witness audits are retained separately with the extracted results.

Final 80-atom audits: [independent WBO rescore](cluster-evidence/results/search_memory_partitioned/independent-audit.json), [producing-family membership](cluster-evidence/results/search_memory_partitioned/constructive-membership-audit.json), and [serial/parallel comparison](cluster-evidence/serial-parallel-equality.json).
'''
(root/'FINAL_REPORT.md').write_text(report)
manifest=json.loads((root/'manifest.json').read_text());manifest.update(status='complete',final_report='FINAL_REPORT.md',final_results='final-results.json')
(root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Final report and machine-readable results written after all audit gates passed')
