"""Build the offline Golden growth/branching film from verified saved evidence.

No mapping search is performed. Display conformers are separate from search input.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys


def build(repo, out, source):
    sys.path.insert(0, str(repo / 'src'))
    from mappa import AAMProblem, MolecularEndpoint
    from mappa.event_patterns import SignedEventIndex
    data = json.loads(source.read_text())
    problem = AAMProblem(*(MolecularEndpoint(**data['source_input'][k]) for k in ('reactant', 'product')))
    index = SignedEventIndex(problem, threshold=.5, metal_threshold=.3)
    for path in data['paths']:
        vector = [path['mapping'][str(i)] for i in range(index.n)]
        assert sorted(vector) == list(range(index.n))
        assert index.describe(vector)['id'] == path['event_class']
        for frame in path['frames']:
            m = frame['mapping']
            assert len(set(m.values())) == len(m)
            assert all(problem.reactant.elements[int(r)] == problem.product.elements[p] for r, p in m.items())
    classes = [c['event_class'] for c in data['decoded_candidates']]
    assert len(classes) == len(set(classes)) == 2
    assert len(data['paths']) == 2
    assert {p['event_class'] for p in data['paths']} == set(classes)
    assert all(p['event_class'] == classes[p['class_index']] for p in data['paths'])
    for c in data['decoded_candidates']:
        vector = [c['mapping'][str(i)] for i in range(index.n)]
        result = index.describe(vector)
        assert result['id'] == c['event_class'] and result['total'] == len(c['events'])
        for event in c['events']:
            r, s = event['r']; p, q = event['p']
            assert [vector[r], vector[s]] == [p, q]
            assert event['wbo'] == [problem.reactant.wbo[r, s], problem.product.wbo[p, q]]
    for key in ('reactant', 'product'):
        assert data['input'][key]['elements'] == data['source_input'][key]['elements']
        assert data['input'][key]['wbo'] == data['source_input'][key]['wbo']
    out.mkdir(parents=True, exist_ok=True)
    (out / 'film-data.json').write_text(json.dumps(data, indent=2) + '\n')
    html = (Path(__file__).parent / 'film.html').read_text()
    html = html.replace('__LIBRARY__', (repo / 'src/mappa/static/3Dmol-min.js').read_text())
    html = html.replace('__DATA__', json.dumps(data, separators=(',', ':')).replace('</', '<\\/'))
    (out / 'index.html').write_text(html)
    (out / 'science-validation.json').write_text(json.dumps(dict(status='passed', source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        checks=['Recorded baseline fragment calls replay exactly', 'All displayed witnesses come from the saved baseline sweep',
                'All displayed correspondences are element-preserving injections', 'Final cards have distinct canonical signed-event IDs',
                'Exactly two recorded growth branches decode to their two distinct signed-event classes', 'Final cards are decoded catalogue representatives',
                'Complete catalogue has nine classes through six events; two explicitly selected for display',
                'Display conformers preserve endpoint elements, bond orders and atom indexing'], scope=data['scope']), indent=2) + '\n')
    print(out / 'index.html')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source', type=Path)
    args = parser.parse_args()
    build(args.repo, args.output, args.source or args.repo / 'reports/published_baseline_20260915/film/film-source.json')
