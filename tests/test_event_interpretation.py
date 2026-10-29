"""Independent small oracles for event functions over saved compressed actions."""
from dataclasses import replace
from itertools import permutations
import json
import time
import tracemalloc
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from mappa import AAMProblem, interpret_events
from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex, _response_bounds
from mappa.mappa_core.utils.bonds import is_metal_element
from test_fragment_conversion import endpoint, saved


def dense_counts(problem, mapping, threshold=.5, metal_threshold=.3):
    result = []
    for ep, other, images in ((problem.reactant, problem.product, mapping),
                             (problem.product, problem.reactant, {p:r for r,p in mapping.items()})):
        count = 0
        for a in range(ep.atom_count):
            for b in range(a + 1, ep.atom_count):
                if a not in images and b not in images:
                    continue
                tau = (metal_threshold if metal_threshold is not None
                       and any(is_metal_element(ep.elements[i]) for i in (a,b)) else threshold)
                target = other.wbo[images[a], images[b]] if a in images and b in images else 0.
                count += ep.wbo[a,b] - target > tau
        result.append(int(count))
    return result


@pytest.mark.parametrize('elements,edges,generators,words,expected', [
    ('CCCCCC', [(0,1,1),(0,2,1.8),(3,4,1),(3,5,1.8)],
     ((0,2,1,3,4,5),(0,1,2,3,5,4)), [(),(1,),(2,),(1,2)], [0,2,2,4]),
    ('CCCC', [(0,2,1),(1,3,1),(0,3,.4),(1,2,.4)],
     ((1,0,2,3),(0,1,3,2)), [(),(1,),(2,),(1,2)], [0,4,4,0]),
    ('CCCC', [(0,2,1),(1,3,1),(0,3,.4),(1,2,.4)],
     ((1,0,3,2),), [(),(1,)], [0,0]),
    ('CCC', [(0,1,.6),(0,2,1),(1,2,1.4)],
     ((1,0,2),(0,2,1)), [(),(1,),(2,),(1,2),(2,1),(1,2,1)], [0,0,0,1,1,2]),
    ('CCC', [(0,1,1),(0,2,1.4)], ((0,2,1),), [(),(1,)], [0,0]),
    (('Fe','C','C'), [(0,1,1),(0,2,1.4)], ((0,2,1),), [(),(1,)], [0,2]),
])
def test_original_correlations_and_event_effects(elements, edges, generators, words, expected):
    ep = endpoint(elements, edges)
    catalogue = saved(ep, ep, dict(enumerate(range(len(elements)))), actions=(('group',generators),))
    before = catalogue.to_record()
    interpreted = interpret_events(catalogue)
    assert interpreted.catalogue is catalogue and catalogue.to_record() == before
    actual = []
    for word in words:
        mapping = interpreted.realize(0, [word])
        counts = dense_counts(catalogue.problem, mapping)
        result = interpreted.evaluate(0, [word])
        assert list(result['counts'].values()) == counts
        actual.append(result['total'])
    assert actual == expected
    assert [int(row.sum()) for row in interpreted.counts(0, ([w] for w in words), batch_size=2)] == expected
    # Only the joint-swap and the ordinary 0.4-difference cases are certified.
    assert interpreted.family_actions[0][0].response_invariant == (max(expected) == 0)


@pytest.mark.parametrize('tau', [.3,.5,1.])
def test_vector_response_bounds_preserve_float_edges(tau):
    weights = np.unique([x for v in (0.,.3,.5,1.,1.3,1.5,1.8)
                         for x in (np.nextafter(v,-np.inf),v,np.nextafter(v,np.inf))])
    lows, highs = _response_bounds(weights, weights, tau)
    for weight, low, high in zip(weights, lows, highs):
        signs = np.where(weights-weight < -tau,-1,np.where(weights-weight > tau,1,0))
        assert np.array_equal(signs, [-1]*low + [0]*(high-low) + [1]*(len(weights)-high))


def test_batched_generator_checks_equal_dense_response_signatures():
    rng = np.random.default_rng(83)
    elements = ('Fe','C','C','C','H')
    weights = (0.,.3,np.nextafter(.3,1.),.5,np.nextafter(.5,1.),1.,1.4)
    endpoints = [endpoint(elements, [(a,b,float(rng.choice(weights)))
                 for a in range(5) for b in range(a+1,5)]) for _ in range(2)]
    problem = AAMProblem(*endpoints)
    generators = tuple((0,)+p+(4,) for p in permutations((1,2,3)))
    def signature(a,b):
        pair = sorted((elements[a],elements[b]))
        opposite = {0.} | {float(problem.reactant.wbo[i,j])
            for i in range(5) for j in range(i+1,5) if sorted((elements[i],elements[j])) == pair}
        tau = .3 if 'Fe' in pair else .5
        return tuple(-1 if w-problem.product.wbo[a,b] < -tau else
                     1 if w-problem.product.wbo[a,b] > tau else 0 for w in sorted(opposite))
    expected = []
    for g in generators:
        count = sum(signature(a,b) != signature(g[a],g[b]) for a in range(5) for b in range(a+1,5))
        expected.append((count == 0,count))
    for size in (1,7,65536):
        assert SignedEventIndex(problem).target_action_responses(generators, batch_entries=size) == tuple(expected)


def test_pool_order_original_constraints_and_inverse_words():
    ep = endpoint('CCC', [(0,1,1),(0,2,2)])
    catalogue = saved(ep,ep,{0:0,1:1,2:2},actions=(('pool',(1,2)),('group',((1,2,0),))))
    interpreted = interpret_events(catalogue)
    assert interpreted.realize(0, [(2,1),(-1,)]) == {0:2,1:1,2:0}
    with pytest.raises(ValueError, match='pool choice'):
        interpreted.realize(0, [(1,1),None])
    # The original matching constraint restricts the composed choices even
    # though the stored pool supplies both permutations.
    catalogue.families[0] = replace(catalogue.families[0], actions=(('pool',(1,2)),),
                                    required_edges=((0,1),), policy=(.2,.5))
    interpreted = interpret_events(catalogue)
    with pytest.raises(ValueError, match='original fragment constraints'):
        interpreted.realize(0, [(2,1)])
    assert interpreted.realize(0) == {0:0,1:1,2:2}


def test_invariant_action_is_not_removed_before_overlapping_action():
    ep = endpoint('CCC', [(0,1,.8),(0,2,.8),(1,2,1.6)])
    catalogue = saved(ep,ep,{0:0,1:1,2:2},
        actions=(('group',((0,2,1),)),('group',((1,0,2),))))
    result = interpret_events(catalogue)
    assert result.family_actions[0][0].response_invariant
    assert len(result.family_actions[0]) == 2
    assert result.evaluate(0,[(1,),(1,)])['events'] != result.evaluate(0,[None,(1,)])['events']


def test_partial_boundaries_and_empty_mapping():
    r, p = endpoint('CCC',[(0,1,1),(1,2,2)]), endpoint('CCCC',[(0,1,1),(2,3,2)])
    for mapping in ({0:0},{0:0,1:1},{}):
        source = saved(r,p,mapping,actions=(('pool',(0,1,2,3)),))
        result = interpret_events(source)
        for choice in ((0,1,2,3),(2,3,0,1)):
            actual = result.realize(0,[choice])
            assert list(result.evaluate(0,[choice])['counts'].values()) == dense_counts(source.problem,actual)


def test_count_batches_bound_pair_scratch(monkeypatch):
    ep = endpoint('C'*300,[(0,1,1.)])
    result = interpret_events(saved(ep,ep,{0:0,1:1}))
    original, sizes = result.index.counts, []
    def checked(vectors):
        sizes.append(len(vectors))
        return original(vectors)
    monkeypatch.setattr(result.index,'counts',checked)
    assert [int(row.sum()) for row in result.counts(0,[None]*3)] == [0,0,0]
    assert sizes == [1,1,1]


def test_response_certificate_does_not_freeze_unmatched_target_identity():
    source = saved(endpoint('C',[]), endpoint('CC',[(0,1,1)]), {0:0}, actions=(('pool',(0,1)),))
    result = interpret_events(source)
    assert result.family_actions[0][0].response_invariant
    first, second = result.evaluate(0), result.evaluate(0,[(1,0)])
    assert first['total'] == second['total'] == 1
    assert first['events'] != second['events']


def test_saved_search_graph_uses_same_families_at_each_event_tolerance():
    from mappa import AAMSearchGraph
    from mappa.mappa_core.utils.io.records import problem_from_record
    record = json.loads((Path(__file__).parent / 'fixtures/package_rename/legacy-aam.json').read_text())
    value = SimpleNamespace(problem=problem_from_record(record), graph=AAMSearchGraph.from_record(record['graph']))
    original = json.dumps(value.graph.to_record(),sort_keys=True)
    loose = interpret_events(value,threshold=1.,metal_threshold=1.)
    strict = interpret_events(loose.catalogue)
    assert strict.catalogue is loose.catalogue
    assert json.dumps(value.graph.to_record(),sort_keys=True) == original
    for fid, family in enumerate(strict.catalogue.families):
        assert list(strict.evaluate(fid)['counts'].values()) == dense_counts(value.problem,dict(family.mapping))
        assert strict.realize(fid) == loose.realize(fid) == dict(family.mapping)


def test_large_compressed_space_never_enumerated(monkeypatch):
    import mappa.mappa_core.post_processing.decoder.api as api
    import sympy.combinatorics as groups
    import z3
    def forbidden(*args, **kwargs):
        pytest.fail('interpretation invoked rematching, a solver, or group enumeration')
    monkeypatch.setattr(api,'convert_fragments',forbidden)
    monkeypatch.setattr(groups.PermutationGroup,'generate_schreier_sims',forbidden)
    monkeypatch.setattr(z3,'Solver',forbidden)
    n = 60
    ep = endpoint('C'*n, [(3*i,3*i+1,1.) for i in range(20)] +
                        [(3*i,3*i+2,1.8) for i in range(20)])
    actions = tuple(('pool',(3*i+1,3*i+2)) for i in range(20))
    source = saved(ep,ep,dict(enumerate(range(n))),actions=actions)
    source.families.append(replace(source.families[0]))
    started = time.process_time()
    tracemalloc.start()
    result = interpret_events(source)
    _,peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    elapsed = time.process_time()-started
    assert len(result.family_actions[0]) == 20
    assert all(a is b for a,b in zip(*result.family_actions))
    assert all(not a.response_invariant for a in result.family_actions[0])
    # Sample requested outcomes, never materialize the 2**20 combinations.
    assert sum(next(result.counts(0,[None]))) == 0
    assert sum(next(result.counts(0,[[tuple(reversed(a[1])) for a in actions]]))) == 40
    print(json.dumps(dict(represented_choices=2**20, unique_actions=20,
        cpu_seconds=elapsed,python_peak_bytes=peak)),flush=True)
