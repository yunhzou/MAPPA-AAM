"""Golden keeps chemical endpoint symmetry and explicit unmatched-atom semantics."""
from dataclasses import asdict
from itertools import permutations
import json
from pathlib import Path
import sqlite3
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'bench'))
from golden import REFERENCE_POLICY, endpoint_generators
from enzymemap import run
from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob
from mappa.execution import BundleStore
from test_enzymemap_campaign import prepared_index


def test_endpoint_group_matches_independent_colored_graph_permutations():
    heavy = [0, 2, 4, 5]
    feature = dict(heavy=heavy, colors=[['C'], ['O'], ['C'], ['O']],
        bonds=[[0, 1, [1]], [1, 2, [1]], [2, 3, [1]], [3, 0, [1]]])
    generators = endpoint_generators(feature, 6)
    identity = tuple(range(6))
    group, pending = {identity}, [identity]
    while pending:
        current = pending.pop()
        for generator in generators:
            moved = tuple(generator[a] for a in current)
            if moved not in group:
                group.add(moved); pending.append(moved)
    edges = {frozenset((a, b)) for a, b, _ in feature['bonds']}
    expected = set()
    for permutation in permutations(range(4)):
        if any(feature['colors'][a] != feature['colors'][permutation[a]] for a in range(4)):
            continue
        if {frozenset((permutation[a], permutation[b])) for a, b, _ in feature['bonds']} != edges:
            continue
        lifted = list(identity)
        for a, b in enumerate(permutation):
            lifted[heavy[a]] = heavy[b]
        expected.add(tuple(lifted))
    assert group == expected
    assert all(g[1] == 1 and g[3] == 3 for g in group)


@pytest.mark.parametrize('omit_atom', [False, True])
def test_direct_inputs_and_golden_unmatched_policy_across_complete_pipeline(tmp_path, omit_atom):
    index = prepared_index(tmp_path)
    data = json.loads(index.read_text())
    store = BundleStore(tmp_path/'prepared')
    mapping = [(i, i) for i in range(1 if omit_atom else 0, 4)]
    reference = store.put(dict(mapping=mapping, features=[dict(heavy=list(range(4)))]*2,
                              generators=[(), ()]))
    for row in data['directions']:
        preparation = read_blob(Blob(**row.pop('blob')), None)
        row['input'] = asdict(preparation[0].payload['input'])
        row['reference'] = asdict(reference)
    data.update(verification='golden:verify', reference_policy=REFERENCE_POLICY)
    index.write_text(json.dumps(data))
    output = tmp_path/'result'
    assert run(index, output, local=True, workers=2, memory_bytes=2*2**30)['finished']
    with sqlite3.connect(output/'ledger.sqlite') as db:
        rows = db.execute("SELECT reference FROM results WHERE stage='verification'").fetchall()
    for raw, in rows:
        result = store.load(Blob(**json.loads(raw)))
        assert result['status'] == ('not_recovered' if omit_atom else 'recovered')
        assert result['policy'] == REFERENCE_POLICY
