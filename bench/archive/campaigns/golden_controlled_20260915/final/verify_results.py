"""Verify published totals and paired CPU means from individual records."""
from pathlib import Path
from collections import Counter
from statistics import mean
import gzip
import hashlib
import json
import math

base = Path(__file__).resolve().parent
data = json.loads(gzip.decompress((base / 'controlled-results.json.gz').read_bytes()))
common = set(data['common_completed_cases'])
methods = data['methods']
assert len(common) == 1403 and len(methods) == 12
computed_common = set(range(1851))
for key, method in methods.items():
    rows = method['per_case']
    assert len(rows) == 1851 and {r['case'] for r in rows} == set(range(1851))
    assert dict(Counter(r['reference_recovery'] for r in rows)) == method['outcomes']
    completed = {r['case'] for r in rows if r['search_complete']}
    computed_common &= completed
    assert len(completed) == method['completed_search_cases']
    paired = mean(r['api_cpu_seconds'] for r in rows if r['case'] in common)
    spent = mean(r['api_cpu_seconds'] for r in rows)
    assert math.isclose(paired, method['paired_completed_api_cpu']['mean_cpu_seconds'], rel_tol=1e-10)
    assert math.isclose(spent, method['all_attempted_api_cpu']['mean_cpu_seconds'], rel_tol=1e-10)
    for r in rows:
        assert math.isclose(r['api_cpu_seconds'], sum(d['api_cpu_seconds'] for d in r['directions'].values()), abs_tol=1e-8)
assert computed_common == common
for audit in json.loads((base / 'paired-cap-audit.json').read_text()):
    seed = audit['seed']
    a = {r['case']: r['reference_recovery'] for r in methods[f'mappa_c100_s{seed}_sweep']['per_case']}
    b = {r['case']: r['reference_recovery'] for r in methods[f'mappa_c2000_s{seed}_sweep']['per_case']}
    resolved = [c for c in a if a[c] != 'unknown' and b[c] != 'unknown']
    assert len(resolved) == audit['resolved_pairs']
    assert all(a[c] == b[c] for c in resolved)
    assert [c for c in a if a[c] == 'recovered' and b[c] == 'unknown'] == audit['cap100_recovered_now_unknown']
for name, digest in json.loads((base / 'provenance.json').read_text())['published_sha256'].items():
    assert hashlib.sha256((base / name).read_bytes()).hexdigest() == digest, name
print('Passed: 12 configurations, 1,851 cases each, 1,403 paired cases; counts, CPU means, cap pairs and hashes agree.')
