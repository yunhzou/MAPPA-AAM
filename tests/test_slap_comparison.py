"""The paired adapter preserves upstream graphs and strict heavy-atom coverage."""
from collections import Counter
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bench'))
from compare_slap import cut_graphs, heavy_mapping, reference_certificate, slap_graphs, slap_source
from mappa import AAMProblem, MolecularEndpoint


@pytest.mark.parametrize('reaction', ['CCO>>CC=O', 'c1ccccc1>>C1CCCCC1', 'CO.O>>CCO'])
def test_prepared_graph_adapter_matches_upstream_including_balance_and_cuts(reaction):
    pytest.importorskip('slapmapper')
    Chem = pytest.importorskip('rdkit.Chem')
    from slapmapper.aam._smiles import smiles2lgp
    from slapmapper.core import SlapMapper
    endpoints = []
    for side in reaction.split('>>'):
        mol = Chem.AddHs(Chem.MolFromSmiles(side))
        weights = np.zeros((mol.GetNumAtoms(), mol.GetNumAtoms()))
        for bond in mol.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            weights[a, b] = weights[b, a] = bond.GetBondTypeAsDouble()
        endpoints.append(MolecularEndpoint(tuple(a.GetSymbol() for a in mol.GetAtoms()),
                                           np.zeros((len(weights), 3)), weights))
    problem = AAMProblem(*endpoints)
    ours, upstream = slap_graphs(problem), smiles2lgp(reaction)
    assert Counter(ours[0].labels) == Counter(ours[1].labels)
    for actual, expected in zip(ours, upstream):
        assert actual.graph == expected.graph
        assert actual.labels == expected.labels
        assert actual._ini_wl_labels == expected._ini_wl_labels
    edges = [(a, b) for a, neighbors in ours[0].graph.items() for b in neighbors if a < b]
    cut = (edges[0],)
    changed = cut_graphs(ours, cut)
    assert cut[0][1] not in changed[0].graph[cut[0][0]]
    assert cut[0][1] in ours[0].graph[cut[0][0]]
    for mode in (False, True):
        answers = []
        for base in (ours, upstream):
            mapper = SlapMapper(binary=mode)
            np.random.seed(20260910)
            mapper.get_maps(cut_graphs(base, ()),
                            break_sym_targets=[i for i, z in enumerate(base[0].labels) if z > 1])
            answers.append({heavy_mapping([g.labels for g in found['lgp']], problem)
                            for found in mapper.results})
        assert answers[0] == answers[1]
    assert slap_source()['core.py']


def test_reference_certificate_preserves_joint_connectivity_and_unmatched_atoms():
    feature = dict(heavy=[0, 1, 2], colors=[['C', 0, 0, 0, '']]*3,
                   bonds=[(0, 1, (1., '')), (1, 2, (1., ''))])
    features = [feature, feature]
    identity = reference_certificate(features, [(0, 0), (1, 1), (2, 2)])
    assert reference_certificate(features, [(0, 2), (1, 1), (2, 0)]) == identity
    assert reference_certificate(features, [(0, 1), (1, 0), (2, 2)]) != identity
    assert reference_certificate(features, [(0, 0), (1, 1)]) != identity


def test_repeated_heavy_labels_are_rejected_not_expanded():
    endpoint = MolecularEndpoint(('C', 'C'), np.zeros((2, 3)), np.zeros((2, 2)))
    with pytest.raises(ValueError, match='repeated heavy'):
        heavy_mapping(([1, 1], [1, 2]), AAMProblem(endpoint, endpoint))


def test_both_methods_publish_then_decode_and_verify_without_rematching(tmp_path):
    pytest.importorskip('slapmapper')
    from dataclasses import asdict
    from compare_slap import match, decode, verify
    from golden import endpoint_generators
    from mappa.execution import BundleStore
    from mappa.mappa_core.core.scheduler.runtime import Workspace
    from mappa.mappa_core.utils.io.records import endpoint_record

    endpoint = MolecularEndpoint(('C', 'C', 'O'), np.zeros((3, 3)),
                                 np.array([[0., 1., 0.], [1., 0., 1.], [0., 1., 0.]]))
    feature = dict(heavy=[0, 1, 2], colors=[['C', 0, 0, 0, ''], ['C', 0, 0, 0, ''],
                                         ['O', 0, 0, 0, '']],
                   bonds=[(0, 1, (1., '')), (1, 2, (1., ''))])
    reference = dict(mapping=[(0, 0), (1, 1), (2, 2)], features=[feature, feature],
                     generators=[endpoint_generators(feature, 3)]*2)
    store = BundleStore(tmp_path/'bundles')
    workspace = Workspace(tmp_path, store=store)
    inputs = asdict(store.put(dict(name='identity', reactant=endpoint_record(endpoint),
                                  product=endpoint_record(endpoint))))
    ref = asdict(store.put(reference))
    for method in ('mappa', 'slap'):
        payload = dict(name=method+'/0/forward', method=method, input=inputs)
        saved = match(payload, workspace)
        assert saved['search']['completed_cuts'] == 3
        payload['matching'] = store.put(saved)
        if method == 'mappa':
            decoded = decode(payload, workspace)
            assert decoded['window'] == (0, 5)
            payload['decoding'] = store.put(decoded)
        payload.update(reference=ref, reversed=False)
        checked = verify(payload, workspace)
        assert checked['verification']['status'] == 'recovered'
        if method == 'mappa':
            assert checked['decoded_verification']['status'] == 'recovered'
        else:
            from slap_audit import audit_saved, native_call
            audited = audit_saved(dict(payload,
                matching=asdict(payload['matching']),
                verification=asdict(store.put(checked))), workspace)
            assert audited['status'] == 'recovered'
            assert audited['calls'] == saved['search']['calls']
            native = native_call(payload, workspace)
            assert native['status'] == audited['binary_uncut_status'] == 'recovered'
            assert native['witness'] in native['candidates']
            assert native['endpoint_correspondence']


@pytest.mark.parametrize('weights,colors', [
    ([1., 1., 1.], ['C', 'C', 'C', 'C']),
    ([1., 2., 1.], ['C', 'C', 'C', 'C']),
    ([1., 1., 1.], ['O', 'C', 'C', 'O']),
])
def test_two_reference_oracles_agree_with_exhaustive_endpoint_symmetries(weights, colors):
    from itertools import permutations, combinations
    from slap_audit import paired_graph, equivalent
    feature = dict(heavy=list(range(4)), colors=[[c, 0, 0, 0, ''] for c in colors],
                   bonds=[(i, i+1, (w, '')) for i, w in enumerate(weights)])
    matrix = np.zeros((4, 4))
    for i, w in enumerate(weights):
        matrix[i, i+1] = matrix[i+1, i] = w
    symmetries = [p for p in permutations(range(4))
                  if [colors[i] for i in p] == colors and np.array_equal(matrix, matrix[np.ix_(p, p)])]
    features = [feature, feature]
    for truth in [((0, 0), (1, 1), (2, 2), (3, 3)), ((0, 0), (1, 1))]:
        allowed = {tuple(sorted((a[x], b[y]) for x, y in truth)) for a in symmetries for b in symmetries}
        certificate = reference_certificate(features, truth)
        graph = paired_graph(features, truth)
        for source in combinations(range(4), len(truth)):
            for target in permutations(range(4), len(truth)):
                mapping = tuple(zip(source, target))
                expected = mapping in allowed
                assert (reference_certificate(features, mapping) == certificate) == expected
                assert (equivalent(graph, paired_graph(features, mapping)) is not None) == expected


def test_native_label_audit_preserves_unmatched_atoms_and_rejects_invalid_groups():
    from slap_audit import label_mapping
    assert label_mapping([[10, 20], [10, 20]], [['C', 'O'], ['C']]) == ((0, 0),)
    assert label_mapping([[20, 10], [10, 20]], [['C', 'C'], ['C']]) == ((1, 0),)
    with pytest.raises(ValueError, match='element'):
        label_mapping([[10, 20], [20, 10]], [['C', 'O'], ['C', 'O']])
    with pytest.raises(ValueError, match='unique'):
        label_mapping([[10, 10], [10, 10]], [['C', 'C'], ['C', 'C']])
