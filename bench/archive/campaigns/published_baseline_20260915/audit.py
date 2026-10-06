"""Audit the published baseline projection; no matching or decoding runs."""
from pathlib import Path
import json,gzip,hashlib,collections
S=Path(__file__).resolve().parent;R=S.parents[1]
read=lambda p:json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
manifest=read(S/'manifest.json')
for name,sha in manifest['sources'].items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==sha,name
windows=read(S/'baseline-windows.json.gz');coord=read(S/'coordinate_minima.json');certs=read(S/'final_dedup.json');old=read(R/'reports/current_validation_20260912/final_dedup.json')
for w,row,cert,original in zip(windows,coord['per_case'],certs['per_case'],old['per_case']):
 assert w['case']==row['case']==cert['case']==original['case']
 assert cert['complete'] and original['complete'] and w['complete_saved_window']
 assert set(w['patterns'])==set(cert['class_ids'])==set(original['baseline_class_ids'])
 minimum=min(p['total'] for p in w['patterns'].values());assert minimum==row['mappa_minimum']
 assert set(row['mappa_minimum_ids'])=={k for k,p in w['patterns'].items() if p['total']==minimum}
 for comparator in row['comparators'].values():
  assert set(comparator['covered_in_mappa_window'])==set(w['patterns'])&set(comparator['ids'])
  assert set(comparator['shared_minimum_ids'])==set(row['mappa_minimum_ids'])&set(comparator['ids'])
for name,total in coord['comparisons'].items():
 rows=[r['comparators'][name] for r in coord['per_case']]
 assert total['count_relations']==dict(collections.Counter(r['count_relation'] for r in rows))
 assert total['covered_in_mappa_window']==sum(len(r['covered_in_mappa_window']) for r in rows)
 assert total['shared_minimum_patterns']==sum(len(r['shared_minimum_ids']) for r in rows)
assert len(windows)==140 and sum(len(w['patterns']) for w in windows)==300
assert sum(len(r['mappa_minimum_ids']) for r in coord['per_case'])==166
assert 'mappa_competition' not in coord['timing'] and 'mappa_decoding' not in coord['timing']
print('Passed: 140 baseline windows, 300 patterns, 166 minimum patterns; no competition.')
