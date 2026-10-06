"""Validate two frozen MAPPA witnesses and rebuild the standalone gold film."""
from pathlib import Path
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))


def build(repo, output):
    from mappa import AAMProblem, MolecularEndpoint
    from mappa.event_patterns import SignedEventIndex
    from mappa.final_branches import FinalFamily
    from mappa.family_query import query_path
    source = repo / 'manuscript/animations/gold_rearrangement'
    data = json.loads((source / 'film-data.json').read_text())
    families = json.loads((repo / 'examples/gold_rearrangement/data/selected-families.json').read_text())
    problem = AAMProblem(*(MolecularEndpoint(d['elements'], d['coordinates'], d['wbo'])
                           for d in (data['input']['reactant'], data['input']['product'])))
    import numpy as np
    for endpoint in data['input'].values():
        original = np.asarray(endpoint['coordinates'])
        display = np.asarray(endpoint['display_coordinates'])
        assert np.allclose(np.linalg.norm(original[:, None] - original[None, :], axis=2),
                           np.linalg.norm(display[:, None] - display[None, :], axis=2), atol=1e-10)
    index = SignedEventIndex(problem, threshold=.5, metal_threshold=.3)
    import gzip
    import hashlib
    trace_bytes = gzip.decompress((source / 'trajectory.json.gz').read_bytes())
    trace = json.loads(trace_bytes)
    tree = data['actual_tree']
    assert hashlib.sha256(trace_bytes).hexdigest() == tree['trace_sha256']
    run = trace['runs'][0]
    assert tree['states'] == run['graph']['states']
    assert tree['transitions'] == run['graph']['transitions']
    assert len(tree['states']) == 15 and len(tree['transitions']) == 14
    assert len(tree['terminals']) == 4
    assert set(tree['featured_terminals']) | set(tree['collapsed_subtree']['terminal_ids']) == {p['terminal'] for p in tree['terminals']}
    assert [p['context'] for p in data['paths']] == [21, 21]
    shared = [[s['transition'] for s in p['stages'][:4]] for p in data['paths']]
    assert shared[0] == shared[1] == tree['shared_prefix_transitions']
    checked_frames = 0
    for path in data['paths']:
        raw = next(p for p in run['paths'] if p['terminal'] == path['terminal'])
        unit_growth = 0
        for frame in path['frames']:
            event = raw['frames'][frame['source_event_index']]
            candidate = event['candidates'][event['preferred']] if event['candidates'] else None
            expected = {**event['locked'], **(candidate['witness'] if candidate else {})}
            assert frame['mapping'] == expected
            assert frame['kind'] == event['kind'] and frame['active'] == event['active']
            if frame['kind'] == 'commit':
                previous = raw['frames'][frame['source_event_index'] - 1]
                assert previous['kind'] == 'pop' and len(frame['active']) == len(previous['active']) + 1
                unit_growth += 1
            checked_frames += 1
        assert unit_growth == 59
    assert all(c['archived_fragment_matches'] == 'exact' for c in trace['checks'])
    canonical = []
    for key, path in zip(('a', 'b'), data['paths']):
        mapping = {int(i): j for i, j in path['mapping'].items()}
        assert sorted(mapping) == sorted(mapping.values()) == list(range(65))
        result = index.describe([mapping[i] for i in range(65)])
        assert result['id'] == path['event_class']
        assert result['total'] == path['event_count']
        family = FinalFamily.from_record(families[key]['family'])
        status, _ = query_path(family.as_path(problem), problem, mapping,
                               source_atoms=range(65), timeout_ms=10000)
        assert status == 'recovered', status
        for stage in path['stages']:
            m = {int(i): j for i, j in stage['mapping'].items()}
            assert len(m) == len(set(m.values())) == stage['assigned']
            assert all(problem.reactant.elements[i] == problem.product.elements[j] for i, j in m.items())
        canonical.append(result['id'])
    assert len(set(canonical)) == 2
    output.mkdir(parents=True, exist_ok=True)
    html = (Path(__file__).parent / 'film.html').read_text()
    html = html.replace('__LIBRARY__', (repo / 'src/mappa/static/3Dmol-min.js').read_text())
    html = html.replace('__DATA__', json.dumps(data, separators=(',', ':')).replace('</', '<\\/'))
    (output / 'index.html').write_text(html)
    (output / 'film-data.json').write_text(json.dumps(data, indent=2) + '\n')
    (output / 'rebuild-validation.json').write_text(json.dumps(dict(status='passed',
        complete_atom_count=65,distinct_classes=canonical,
        family_membership_queries=['recovered', 'recovered'], recorded_previews_checked=checked_frames, unit_growth_events_per_path=59, actual_context=21, actual_states=15, actual_transitions=14, actual_terminal_branches=4, collapsed_terminal_branches=2, shared_prefix_transitions=shared[0]), indent=2) + '\n')
    print(output / 'index.html')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    build(args.repo, args.output or args.repo / 'manuscript/animations/gold_rearrangement')
