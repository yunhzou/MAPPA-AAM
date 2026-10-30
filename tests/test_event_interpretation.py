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


def dense_categories(problem, mapping, threshold=.5, metal_threshold=.3):
    """Independent scalar oracle, including original endpoint pair identities."""
    result = []
    r,p = problem.reactant,problem.product
    for a in range(r.atom_count):
        for b in range(a+1,r.atom_count):
            if a not in mapping and b not in mapping:
                continue
            tau = metal_threshold if metal_threshold is not None and any(
                is_metal_element(r.elements[k]) for k in (a,b)) else threshold
            x,y = mapping.get(a),mapping.get(b)
            old = float(r.wbo[a,b]); new = float(p.wbo[x,y]) if x is not None and y is not None else 0.
            if old-new > tau:
                kind = 'weakening' if new > tau else 'breaking'
            elif new-old > tau:
                kind = 'strengthening' if old > tau else 'forming'
            else:
                continue
            result.append((kind,(a,b),tuple(sorted((x,y))) if x is not None and y is not None else None))
    inverse = {v:k for k,v in mapping.items()}
    for a in range(p.atom_count):
        for b in range(a+1,p.atom_count):
            tau = metal_threshold if metal_threshold is not None and any(
                is_metal_element(p.elements[k]) for k in (a,b)) else threshold
            if (a in inverse) != (b in inverse) and p.wbo[a,b] > tau:
                result.append(('forming',None,(a,b)))
    return sorted(result,key=repr)


def check_contributions(interpreted, choices):
    from collections import Counter
    from mappa.mappa_core.post_processing.decoder.events import EVENT_CATEGORIES
    expected_counts = []
    for choice in choices:
        mapping = interpreted.realize(0,choice)
        expected = dense_categories(interpreted.catalogue.problem,mapping,
            interpreted.index.threshold,interpreted.index.metal_threshold)
        result = interpreted.contributions(0,choice,include_events=True)
        actual = sorted([(e['category'],e['reactant_pair'],e['product_pair']) for e in result['events']],key=repr)
        assert actual == expected
        counts = Counter(e[0] for e in expected)
        assert result['counts'] == {k:counts[k] for k in EVENT_CATEGORIES}
        assert result['total'] == len(expected)
        for kind,count in result['counts'].items():
            assert sum(row[kind] for row in (*result['internal'].values(),*result['boundary'].values())) == count
        expected_counts.append([counts[k] for k in EVENT_CATEGORIES])
    assert list(map(list,interpreted.counts(0,choices,categories=True,batch_size=2))) == expected_counts


def test_sparse_boundaries_keep_shared_choices_and_manual_cut_bonds():
    ep = endpoint('CCCC',[(0,2,1),(1,3,1),(0,3,.4),(1,2,.4)])
    actions = (('pool',(0,1)),('pool',(2,3)))
    catalogue = saved(ep,ep,dict(enumerate(range(4))),actions=actions,
                      pairs=(((0,1),(0,1)),((2,3),(2,3))))
    interpreted = interpret_events(catalogue)
    choices = [[None,None],[(1,0),None],[None,(3,2)],[(1,0),(3,2)]]
    check_contributions(interpreted,choices)
    table = interpreted.tables(0)
    assert not table.internal
    assert set(table.boundary) == {(0,1)}
    assert all(term.action_ids == (0,1) for term in table.boundary[0,1])
    assert [interpreted.contributions(0,c)['total'] for c in choices] == [0,4,4,0]
    # No independent minima of boundary terms or fragment choices are combined.
    catalogue.families[0] = replace(catalogue.families[0],actions=(('group',((1,0,3,2),)),))
    check_contributions(interpret_events(catalogue),[[None],[(1,)]])


@pytest.mark.parametrize('partial',[False,True])
def test_boundary_categories_against_exhaustive_small_permutations(partial):
    r = endpoint('CCCC',[(0,1,1.8),(0,2,.2),(1,2,.8),(2,3,1.4)])
    p = endpoint('CCCCC',[(0,1,.7),(1,2,1.8),(2,3,.4),(0,4,1.2),(3,4,1)])
    mapping = {0:0,1:1} if partial else {0:0,1:1,2:2,3:3}
    pairs = tuple(((a,),(b,)) for a,b in mapping.items())
    interpreted = interpret_events(saved(r,p,mapping,actions=(('pool',tuple(range(5))),),pairs=pairs))
    check_contributions(interpreted,[[perm] for perm in permutations(range(5))])


def test_boundary_categorization_float_thresholds_and_metals():
    for elements in ('CCC',('Fe','C','C')):
        tau = .3 if elements[0] == 'Fe' else .5
        for weight in (tau,np.nextafter(tau,1.),.8,1.8):
            r = endpoint(elements,[(0,1,weight)])
            p = endpoint(elements,[(0,1,tau),(0,2,1.8)])
            result = interpret_events(saved(r,p,{0:0,1:1,2:2},actions=(('pool',(1,2)),),
                pairs=(((0,),(0,)),((1,2),(1,2)))))
            check_contributions(result,[[None],[(2,1)]])


def test_tables_follow_ordered_overlapping_actions_and_required_edges():
    ep = endpoint('CCC',[(0,1,.6),(0,2,1),(1,2,1.4)])
    source = saved(ep,ep,{0:0,1:1,2:2},actions=(('pool',(0,1)),('pool',(1,2))),
        pairs=(((0,),(0,)),((1,2),(1,2))))
    result = interpret_events(source)
    assert result.tables(0).atom_domains == {0:(0,1,2),1:(0,1,2),2:(1,2)}
    check_contributions(result,[[None,None],[(1,0),None],[None,(2,1)],[(1,0),(2,1)]])
    source.families[0] = replace(source.families[0],required_edges=((0,1),),policy=(.2,.1))
    with pytest.raises(ValueError,match='original fragment constraints'):
        interpret_events(source).contributions(0,[None,(2,1)])


def test_fixed_internal_and_boundary_responses_share_across_families():
    ep = endpoint('CCCC',[(0,1,1),(2,3,1)])
    target = endpoint('CCCC',[(0,1,1.8),(2,3,1.8),(0,3,1)])
    source = saved(ep,target,dict(enumerate(range(4))),pairs=(((0,1),(0,1)),((2,3),(2,3))))
    source.families.append(replace(source.families[0]))
    next(iter(source.branches.values())).add(1)
    result = interpret_events(source)
    check_contributions(result,[None])
    assert result.tables(0) is result.tables(1)
    assert result.tables(0).statistics['variable_terms'] == 0
    assert result.contributions(0)['counts'] == dict(breaking=0,forming=1,strengthening=2,weakening=0)


def test_fragment_partitions_are_alternative_views_not_double_counts():
    ep = endpoint('CCC',[(0,1,1),(1,2,1.8)])
    source = saved(ep,ep,{0:0,1:1,2:2},actions=(('pool',(0,2)),))
    source.branches[(source.families[0].policy,(((0,),(0,)),((1,2),(1,2))))] = {0}
    interpreted = interpret_events(source)
    assert len(interpreted.fragment_partitions[0]) == 2
    counts = [interpreted.contributions(0,[(2,0)],partition=i)['counts'] for i in range(2)]
    assert counts[0] == counts[1]


def test_empty_mapping_has_no_internal_or_boundary_events():
    result = interpret_events(saved(endpoint('CC',[(0,1,1)]),endpoint('CCC',[(0,1,1)]),{}))
    check_contributions(result,[None])
    assert result.tables(0).statistics['boundary_terms'] == 0


def test_large_boundary_tables_do_not_expand_combinations(monkeypatch):
    import sympy.combinatorics as groups
    def forbidden(*args,**kwargs):
        pytest.fail('boundary tables expanded group elements')
    monkeypatch.setattr(groups.PermutationGroup,'generate_schreier_sims',forbidden)
    n = 40
    ep = endpoint('C'*n,[(i,i+1,1.8) for i in range(n-1)])
    actions = tuple(('pool',(i,i+1)) for i in range(0,n,2))
    pairs = tuple(((i,i+1),(i,i+1)) for i in range(0,n,2))
    result = interpret_events(saved(ep,ep,dict(enumerate(range(n))),actions=actions,pairs=pairs))
    before = time.process_time()
    table = result.tables(0)
    elapsed = time.process_time()-before
    choices = [[None]*20,[tuple(reversed(a[1])) for a in actions]]
    check_contributions(result,choices)
    assert len(table.boundary) == 19
    assert table.statistics['boundary_terms'] == 76
    assert table.statistics['unique_responses'] == 2
    print(json.dumps(dict(boundary_table_cpu_seconds=elapsed,represented_choices=2**20,
                          statistics=table.statistics)),flush=True)


def test_single_atom_mapping_has_no_bond_tables():
    ep = endpoint('C',[])
    for mapping in ({}, {0:0}):
        result = interpret_events(saved(ep,ep,mapping))
        check_contributions(result,[None])
        assert result.contributions(0)['total'] == 0


@pytest.mark.parametrize('actions,expected_case', [
    ((), 'fixed_fixed'),
    ((('pool',(2,3)),), 'fixed_shuffle'),
    ((('pool',(0,1)),('pool',(2,3))), 'shuffle_shuffle'),
])
def test_cut_endpoint_cases_use_original_bonds(actions, expected_case):
    r = endpoint('CCCC',[(0,2,1)])
    p = endpoint('CCCC',[(0,2,1),(1,3,1)])
    source = saved(r,p,dict(enumerate(range(4))),actions=actions,
                   pairs=(((0,1),(0,1)),((2,3),(2,3))))
    source.families[0].provenance.append(dict(cuts=((0,2),)))
    result = interpret_events(source)
    table = result.tables(0)
    term = next(t for t in table.boundary[0,1] if t.atoms == (0,2))
    assert term.manual_cut and term.endpoint_case == expected_case
    assert len(term.target_pairs) == 2**len(actions)
    choices = [[None]*len(actions),[tuple(reversed(a[1])) for a in actions]]
    check_contributions(result,choices)


def test_joint_cut_endpoint_projection_certifies_constant_boundary():
    ep = endpoint('CCCC',[(0,2,1),(1,3,1)])
    source = saved(ep,ep,dict(enumerate(range(4))),
        actions=(('group',((1,0,3,2),)),),pairs=(((0,1),(0,1)),((2,3),(2,3))))
    source.families[0].provenance.append(dict(cuts=((0,2),)))
    result = interpret_events(source)
    table = result.tables(0)
    term = next(t for t in table.boundary[0,1] if t.atoms == (0,2))
    assert term.target_pairs == ((0,2),(1,3))
    assert term.response.constant == 0
    assert term.endpoint_case == 'shuffle_shuffle'
    assert table.statistics['variable_boundary_terms'] == 0
    check_contributions(result,[[None],[(1,)]])


def test_pair_projection_matches_independent_small_action_composition():
    ep = endpoint('CCCC',[])
    index = SignedEventIndex(AAMProblem(ep,ep))
    identity = tuple(range(4))
    cycle = tuple(tuple((i+step)%4 for i in range(4)) for step in range(4))
    swap = (identity,(0,2,1,3))
    for actions, first, second in (
        ((('group',((1,2,3,0),)),('pool',(1,2))),cycle,swap),
        ((('pool',(1,2)),('group',((1,2,3,0),))),swap,cycle),
    ):
        complete = {tuple(b[a[i]] for i in range(4)) for a in first for b in second}
        for x in range(-1,4):
            for y in range(4):
                if x == y:continue
                expected = {(g[x] if x >= 0 else -1,g[y]) for g in complete}
                assert set(index.endpoint_pairs((x,y),actions)) == expected
    # A same-pool pair is injective and represents 12! choices with 132 pairs.
    ep = endpoint('C'*12,[])
    index = SignedEventIndex(AAMProblem(ep,ep))
    assert len(index.endpoint_pairs((0,1),(('pool',tuple(range(12))),))) == 132


def test_local_contribution_components_retain_shared_actions_and_constraints():
    ep = endpoint('CCCCCC',[(0,1,1),(0,2,1.8),(3,4,1),(3,5,1.8)])
    source = saved(ep,ep,dict(enumerate(range(6))),
        actions=(('pool',(1,2)),('pool',(4,5))),
        pairs=(((0,1,2),(0,1,2)),((3,4,5),(3,4,5))))
    result = interpret_events(source)
    assert result.tables(0).components == (((0,),(0,)),((1,),(1,)))
    check_contributions(result,[[None,None],[(2,1),None],[None,(5,4)],[(2,1),(5,4)]])
    # Equal responses with different feasibility constraints need distinct
    # dependency components, including within the same interpretation cache.
    source.families.append(replace(source.families[0],required_edges=((1,4),)))
    next(iter(source.branches.values())).add(1)
    result = interpret_events(source)
    assert result.tables(0).components == (((0,),(0,)),((1,),(1,)))
    assert result.tables(1).components == (((0,1),(0,1)),)
    source.families[0] = replace(source.families[0],required_edges=(),
                                actions=(('group',((0,2,1,3,5,4),)),))
    result = interpret_events(source)
    assert result.tables(0).components == (((0,1),(0,)),)
    check_contributions(result,[[None],[(1,)]])


def test_sweep_cut_provenance_survives_catalogue_roundtrip():
    from mappa import AAMSearchConfig, search_aam
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    ep = endpoint('CCC',[(0,1,1),(1,2,1)])
    aam = search_aam(AAMProblem(ep,ep),AAMSearchConfig(seed_count=1,sweep_cuts=True))
    source = aam.final_catalogue()
    restored = FinalBranchCatalogue.from_record(source.problem,json.loads(json.dumps(source.to_record())))
    cuts = {tuple(map(tuple,p['cuts'])) for f in restored.families for p in f.provenance}
    assert cuts == {c.cuts for c in aam.graph.contexts}
    assert cuts == {(),((0,1),),((1,2),)}
    result = interpret_events(restored)
    for fid,family in enumerate(restored.families):
        assert result.tables(fid).cut_bonds == frozenset(tuple(p) for row in family.provenance for p in row['cuts'])
        assert result.contributions(fid)['total'] == sum(dense_counts(source.problem,dict(family.mapping)))
