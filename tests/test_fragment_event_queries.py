"""Independent small-graph oracles for retained symmetry/event selections."""
from collections import Counter
from itertools import permutations, product
import json
import pickle

import numpy as np
import pytest

from mappa import decode_events, decoded_events_record, decoded_events_from_record
from mappa.mappa_core.core.symmetry.minimum import FragmentSymmetryEvents
from test_fragment_optimization import fragment, weights
from test_event_interpretation import saved, endpoint


def automorphs(matrix, elements, anchors=()):
    binary = matrix > .2
    return [p for p in permutations(range(len(matrix)))
            if all(elements[a] == elements[p[a]] for a in range(len(p)))
            and all(p[a] == a for a in anchors)
            and np.array_equal(binary, binary[np.ix_(p, p)])]


def oracle(value):
    left, right = value.source_symmetry, value.target_symmetry
    result = {}
    start = dict(value.witness)
    for source, target in product(automorphs(left.weights, left.elements, left.anchors),
                                  automorphs(right.weights, right.elements, right.anchors)):
        images = [None] * len(source)
        for a in range(len(source)):
            images[source[a]] = target[start[a]]
        count = int(np.count_nonzero(np.triu(abs(left.weights - right.weights[np.ix_(images, images)]) > .5, 1)))
        result[tuple(images)] = count
    return result


@pytest.mark.parametrize('seed', range(24))
def test_target_queries_and_endpoint_conditions_against_complete_oracle(seed):
    rng = np.random.default_rng(seed)
    left, right = [weights(4, [(a, b, float(rng.choice([0., .2, .8, 1., 1.6, 2.])))
                              for a in range(4) for b in range(a + 1, 4)]) for _ in range(2)]
    value = fragment(left, right, (), range(4))
    expected = oracle(value)
    data = FragmentSymmetryEvents(value)
    for count in range(7):
        for image in range(4):
            answer = data.select(count).find(assignments={0: image})
            allowed = {p for p, cost in expected.items() if cost == count and p[0] == image}
            assert answer.status == ('found' if allowed else 'absent')
            if allowed:
                assert tuple(dict(answer.witness)[a] for a in range(4)) in allowed


def arms():
    # C(CH2)4: four exchangeable arms, each carrying a pair of H atoms.
    elements = tuple('CCCCC' + 'H' * 8)
    left, right = [], []
    for arm in range(4):
        carbon, h = arm + 1, 5 + arm * 2
        left.extend(((0, carbon, 1.), (carbon, h, 1.), (carbon, h+1, 2.)))
        right.extend(((0, carbon, 1.), (carbon, h, 1.), (carbon, h+1, 2.)))
    return fragment(weights(13, left), weights(13, right), [(a,b) for a,b,_ in left],
                    range(13), elements=elements)


def test_arm_occupancy_internal_exchanges_and_query_storage():
    data = FragmentSymmetryEvents(arms())
    counts = Counter()
    # Independent construction: choose an arm assignment and one flip per arm.
    mappings = []
    for positions in permutations(range(4)):
        for flips in product((0, 1), repeat=4):
            p = [0] * 13
            for arm, target in enumerate(positions):
                p[arm+1] = target+1
                p[5+2*arm] = 5+2*target+flips[arm]
                p[6+2*arm] = 6+2*target-flips[arm]
            count = 2 * sum(flips)
            counts[count] += 1
            mappings.append((p,count))
    assert counts == {0:24, 2:96, 4:144, 6:96, 8:24}
    for count in range(10):
        answer = data.select(count).find()
        assert answer.status == ('found' if count in counts else 'absent')
    selected = data.select(2)
    for images, count in mappings:
        answer = selected.find(assignments=dict(enumerate(images)))
        assert answer.status == ('found' if count == 2 else 'absent')
    assert data.select(0,8).find(assignments={1:1,2:1}).status == 'absent'
    restored = pickle.loads(pickle.dumps((data.select(2),data.select(0,2))))
    assert restored[0].response is restored[1].response
    assert restored[0].find().event_count == 2


def test_explicit_limit_returns_unknown_not_absence():
    data = FragmentSymmetryEvents(arms())
    answer = data.select(2).find(max_nodes=0)
    assert answer.status == 'unknown' and answer.work.cap_reached


def test_discovery_improvement_keeps_catalogue_and_selects_one_decoded_orbit():
    from test_fragment_optimization import triangle
    value = triangle()
    matrix = value.source_symmetry.weights
    r = endpoint('CCCC', [(i,j,matrix[i,j]) for i in range(4) for j in range(i+1,4) if matrix[i,j]])
    catalogue = saved(r,r,dict(value.witness),fragments=(value,))
    before = catalogue.to_record()
    decoded = decode_events(catalogue)
    selected = decoded.fragments[0]
    assert selected.optimization.event_count == 0
    assert selected.optimization.discovery_witness != value.witness
    assert catalogue.to_record() == before
    # The original discovery remains inspectable, but lies outside this selected
    # decoded orbit. The two API operations make different, explicit promises.
    assert decoded.realize(0,[dict(discovery=dict(value.witness))]) == dict(value.witness)
    assert decoded.select(0,0,6).find(assignments=dict(value.witness)).status == 'absent'
    answer = decoded.select(0,0).find()
    assert answer.status == 'found' and decoded.contributions(0,answer.choices)['total'] == 0


def test_saved_response_restores_without_symmetry_or_starting_search(monkeypatch):
    from mappa.mappa_core.post_processing.decoder import fragments
    value = arms()
    r = endpoint(value.source_symmetry.elements, [(i,j,value.source_symmetry.weights[i,j])
                 for i,j in zip(*np.where(np.triu(value.source_symmetry.weights > 0,1)))])
    decoded = decode_events(saved(r,r,dict(value.witness),fragments=(value,)))
    raw = json.loads(json.dumps(decoded_events_record(decoded)))
    monkeypatch.setattr(fragments, 'optimize_fragment', lambda *a, **k: pytest.fail('reoptimized saved fragment'))
    monkeypatch.setattr(FragmentSymmetryEvents, '__init__', lambda *a, **k: pytest.fail('reprepared saved response'))
    restored = decoded_events_from_record(raw)
    assert restored.minimum_events == 0 and restored.minimum_proven
    for count in range(9):
        answer = restored.select(0,count).find()
        assert answer.status == ('found' if count % 2 == 0 else 'absent')
        if answer.status == 'found':
            assert restored.contributions(0,answer.choices)['total'] == count


def test_restored_window_keeps_only_prepared_fragments():
    r = endpoint('CCC', [(0,1,1.), (1,2,1.)])
    catalogue = saved(r,r,{0:0,1:1,2:2})
    extra = saved(r,r,{0:2,1:1,2:0}, pairs=(((0,),(2,)),((1,),(1,)),((2,),(0,))))
    catalogue.families.extend(extra.families)
    catalogue.branches.update({key:{1} for key in extra.branches})
    result = decode_events(catalogue,event_slack=0)
    restored = decoded_events_from_record(json.loads(json.dumps(decoded_events_record(result))))
    assert len(restored.inventory.fragments) == 2
    assert len(restored.fragments) == 1
    assert restored.minimum_events == 0
    assert restored.select(0,0).find().status == 'found'
