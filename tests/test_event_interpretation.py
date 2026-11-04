"""Independent bond-event and preparation checks for saved fragment relations."""
from dataclasses import replace
from itertools import permutations, product
import json
import time

import numpy as np
import pytest

from mappa import AAMProblem, MolecularEndpoint, decode_events
from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex
from mappa.mappa_core.core.symmetry.responses import response_bounds
from mappa.mappa_core.utils.bonds import is_metal_element
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue, FinalFamily
from growth_cases import saved_fragment


def endpoint(elements, edges):
    matrix = np.zeros((len(elements), len(elements)))
    for a, b, weight in edges:
        matrix[a, b] = matrix[b, a] = weight
    return MolecularEndpoint(tuple(elements), np.zeros((len(elements), 3)), matrix)


def saved(reactant, product, mapping, *, fragments=None, pairs=None):
    """Store explicit fixed-image fragments; no action-only catalogue fixtures."""
    problem = AAMProblem(reactant, product)
    if fragments is None:
        fragments = ((saved_fragment(problem, source_atoms=tuple(sorted(mapping)),
            target_atoms=tuple(sorted(mapping.values())), witness=mapping,
            anchors=tuple(sorted(mapping.items()))),) if mapping else ())
    actions = tuple(('relation', fragment) for fragment in fragments)
    required = tuple(sorted({edge for fragment in fragments for edge in fragment.accepted_edges}))
    catalogue = FinalBranchCatalogue(problem)
    catalogue.families = [FinalFamily(tuple(sorted(mapping.items())), required, actions, (.2, 1.))]
    if pairs is None:
        pairs = ((tuple(sorted(mapping)), tuple(sorted(mapping.values()))),) if mapping else ()
    catalogue.branches[((.2, 1.), pairs)] = {0}
    catalogue.path_count = 1
    return catalogue


def swap(fragment):
    return dict(discovery=dict(zip(fragment.source_atoms, reversed(fragment.target_atoms))))


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


def test_neighbor_domains_and_boundary_filter_match_independent_dense_pairs():
    from itertools import combinations
    rng = np.random.default_rng(73514)
    elements = ('C', 'O', 'Fe', 'C', 'H', 'C', 'O')
    for _ in range(32):
        endpoints = [endpoint(elements, [(a, b, float(rng.choice((0., .1, .3, .5, .8, 1.4))))
                     for a, b in combinations(range(7), 2)]) for _ in range(2)]
        value = AAMProblem(*endpoints)
        index = SignedEventIndex(value)
        domains = {a: tuple(p for p in range(7) if rng.random() < .5)
                   for a in range(7) if rng.random() < .7}
        owners = [int(rng.integers(3)) if a in domains else -1 for a in range(7)]
        def threshold(a, b):
            return .3 if any(is_metal_element(elements[i]) for i in (a, b)) else .5
        expected = set()
        for a, b in combinations(range(7), 2):
            if a not in domains and b not in domains:
                continue
            old = value.reactant.wbo[a, b] > threshold(a, b)
            reachable = any(p != q and value.product.wbo[p, q] > threshold(p, q)
                            for p in domains.get(a, ()) for q in domains.get(b, ()))
            if old or reachable:
                expected.add((a, b))
        assert index.neighboring_pairs(domains) == expected
        assert index.neighboring_pairs(domains, owners=owners) == {
            pair for pair in expected if owners[pair[0]] != owners[pair[1]]}

    # Atom masks retain sparse NumPy indices beyond a machine integer's width.
    ep = endpoint('C' * 131, [(0, 1, 1.)])
    index = SignedEventIndex(AAMProblem(ep, ep))
    domains = {np.int64(65): (0, 1), np.int64(130): (0, 1)}
    assert index.neighboring_pairs(domains) == {(65, 130)}



@pytest.mark.parametrize('tau', [.3,.5,1.])
def test_vector_response_bounds_preserve_float_edges(tau):
    weights = np.unique([x for v in (0.,.3,.5,1.,1.3,1.5,1.8)
                         for x in (np.nextafter(v,-np.inf),v,np.nextafter(v,np.inf))])
    lows, highs = response_bounds(weights, weights, tau)
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


def test_exact_weight_placements_do_not_prepare_opposite_response_labels():
    source = endpoint('COCO', [(0, 1, 1.4), (2, 3, 2.1)])
    target = endpoint('COCO', [(0, 1, 1.2), (2, 3, 1.2)])
    problem = AAMProblem(source, target)
    index = SignedEventIndex(problem)
    exchange = (2, 3, 0, 1)
    assert index.target_action_responses((exchange,), batch_entries=1) == ((True, 0),)
    assert 'response_data' not in index.__dict__
    for images in permutations(range(4)):
        mapping = dict(enumerate(images))
        moved = {r: exchange[p] for r, p in mapping.items()}
        assert dense_counts(problem, mapping) == dense_counts(problem, moved)



def test_count_batches_bound_pair_scratch(monkeypatch):
    ep = endpoint('C'*300,[(0,1,1.)])
    result = decode_events(saved(ep,ep,{0:0,1:1}))
    original, sizes = result.index.counts, []
    def checked(vectors):
        sizes.append(len(vectors))
        return original(vectors)
    monkeypatch.setattr(result.index,'counts',checked)
    assert [int(row.sum()) for row in result.counts(0,[None]*3)] == [0,0,0]
    assert sizes == [1,1,1]



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
    from mappa.mappa_core.utils.io.records import decoded_events_record
    record = json.loads(json.dumps(decoded_events_record(interpreted)))
    expected_counts = []
    for choice in choices:
        mapping = interpreted.realize(0,choice)
        expected = dense_categories(interpreted.index.problem,mapping,
            interpreted.index.threshold,interpreted.index.metal_threshold)
        result = interpreted.contributions(0,choice,include_events=True)
        actual = sorted([(e['category'],e['reactant_pair'],e['product_pair']) for e in result['events']],key=repr)
        assert actual == expected
        # The saved local responses suffice without the derived per-atom maps.
        # Check every saved partition against the independent scalar oracle.
        for number in record['families'][0]['tables']:
            table, exported = record['tables'][number], []
            for region in (*table['internal'], *table['boundary']):
                for term in region['terms']:
                    images = tuple(mapping.get(a) for a in term['atoms'])
                    weight = (0. if None in images else
                              float(interpreted.index.problem.product.wbo[images]))
                    response = record['responses'][term['response']]
                    code = response['categories'][response['target_weights'].index(weight)]
                    if code:
                        exported.append((EVENT_CATEGORIES[code - 1], tuple(term['atoms']),
                                         None if None in images else tuple(sorted(images))))
            occupied = set(mapping.values())
            exported.extend(('forming', None, tuple(pair))
                            for pair in table['product_boundary_pairs']
                            if (pair[0] in occupied) != (pair[1] in occupied))
            assert sorted(exported, key=repr) == expected
        counts = Counter(e[0] for e in expected)
        assert result['counts'] == {k:counts[k] for k in EVENT_CATEGORIES}
        assert result['total'] == len(expected)
        for kind,count in result['counts'].items():
            assert sum(row[kind] for row in (*result['internal'].values(),*result['boundary'].values())) == count
        expected_counts.append([counts[k] for k in EVENT_CATEGORIES])
    assert list(map(list,interpreted.counts(0,choices,categories=True,batch_size=2))) == expected_counts



def test_boundary_categorization_float_thresholds_and_metals():
    for elements in ('CCC',('Fe','C','C')):
        tau = .3 if elements[0] == 'Fe' else .5
        for weight in (tau,np.nextafter(tau,1.),.8,1.8):
            r = endpoint(elements,[(0,1,weight)])
            p = endpoint(elements,[(0,1,tau),(0,2,1.8)])
            fragment = saved_fragment(AAMProblem(r, p), source_atoms=(1, 2))
            result = decode_events(saved(r,p,{0:0,1:1,2:2},fragments=(fragment,),
                pairs=(((0,),(0,)),((1,2),(1,2)))))
            check_contributions(result,[[None],[swap(fragment)]])



def test_fixed_internal_and_boundary_responses_share_across_families():
    ep = endpoint('CCCC',[(0,1,1),(2,3,1)])
    target = endpoint('CCCC',[(0,1,1.8),(2,3,1.8),(0,3,1)])
    source = saved(ep,target,dict(enumerate(range(4))),pairs=(((0,1),(0,1)),((2,3),(2,3))))
    source.families.append(replace(source.families[0]))
    next(iter(source.branches.values())).add(1)
    result = decode_events(source)
    check_contributions(result,[None])
    assert result.tables(0) is result.tables(1)
    assert result.tables(0).statistics['variable_terms'] == 0
    assert result.contributions(0)['counts'] == dict(breaking=0,forming=1,strengthening=2,weakening=0)



@pytest.mark.parametrize('mismatch', ['elements', 'coordinates'])
def test_event_problem_cannot_change_the_saved_atom_indexing(mismatch):
    ep = endpoint('CO', [(0, 1, 1.)])
    catalogue = saved(ep, ep, {0: 0, 1: 1})
    changed = replace(ep, elements=('O', 'C')) if mismatch == 'elements' else replace(
        ep, coordinates=np.asarray([[1., 0., 0.], [0., 0., 0.]]))
    with pytest.raises(ValueError):
        decode_events(catalogue, event_problem=AAMProblem(ep, changed))



def test_empty_mapping_has_no_internal_or_boundary_events():
    result = decode_events(saved(endpoint('CC',[(0,1,1)]),endpoint('CCC',[(0,1,1)]),{}))
    check_contributions(result,[None])
    assert result.tables(0).statistics['boundary_terms'] == 0



def test_single_atom_mapping_has_no_bond_tables():
    ep = endpoint('C',[])
    for mapping in ({}, {0:0}):
        result = decode_events(saved(ep,ep,mapping))
        check_contributions(result,[None])
        assert result.contributions(0)['total'] == 0



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



def test_compiled_pair_program_matches_full_permutations_with_overlapping_actions():
    # An independent full-permutation oracle on four atoms checks correlations,
    # action order and unmatched endpoints. Production only projects pairs.
    identity = tuple(range(4))
    def closure(generators):
        found = {identity}
        frontier = [identity]
        while frontier:
            current = frontier.pop()
            for generator in generators:
                image = tuple(generator[i] for i in current)
                if image not in found:
                    found.add(image)
                    frontier.append(image)
        return found
    pool = ('pool', (0, 1, 2))
    group = ('group', ((1, 0, 2, 3), (0, 2, 3, 1)))
    stages = {
        pool: {tuple(p) + (3,) for p in permutations(range(3))},
        group: closure(group[1]),
        ('pool', (1, 3)): {identity, (0, 3, 2, 1)},
        ('group', (identity,)): {identity},
        ('pool', ()): {identity},
    }
    ep = endpoint('CCCC', [])
    index = SignedEventIndex(AAMProblem(ep, ep))
    programs = ((pool, group), (group, pool),
                (pool, ('pool', (1, 3)), group),
                (('pool', ()), ('group', (identity,))))
    for actions in programs:
        complete = {identity}
        for action in actions:
            complete = {tuple(b[a[i]] for i in range(4))
                        for a in complete for b in stages[action]}
        for x in range(-1, 4):
            for y in range(-1, 4):
                expected = {(g[x] if x >= 0 else -1, g[y] if y >= 0 else -1)
                            for g in complete}
                assert set(index.endpoint_pairs((x, y), actions)) == expected
    assert len(index._pair_programs) == len(programs)
    assert len(index._compiled_pair_actions) == len(stages)



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
    result = decode_events(restored)
    for fid,family in enumerate(restored.families):
        assert result.tables(fid).cut_bonds == frozenset(tuple(p) for row in family.provenance for p in row['cuts'])
        assert result.contributions(fid)['total'] == sum(dense_counts(source.problem,dict(family.mapping)))



def test_preparation_visits_only_bond_neighbors_and_recorded_cut(monkeypatch):
    n = 80
    source_edges = [(i, i+1, 1.) for i in range(n-1)]
    source = saved(endpoint('C'*n, source_edges),
        endpoint('C'*n, source_edges + [(17,42,1.)]), dict(enumerate(range(n))),
        pairs=tuple(((a,),(a,)) for a in range(n)))
    source.families[0].provenance.append(dict(cuts=((0,79),)))
    result = decode_events(source)
    original = result.index.domain_response
    visits = []
    def record(weight, tau, xs, ys):
        visits.append((xs, ys))
        return original(weight, tau, xs, ys)
    monkeypatch.setattr(result.index, 'domain_response', record)
    table = result.tables(0)
    assert len(visits) == table.statistics['neighboring_pairs'] == 81
    assert table.statistics['cut_terms'] == 1
    assert set(table.boundary) == {(0,79),(17,42)}
    check_contributions(result, [None])
    assert result.contributions(0)['counts'] == dict(breaking=0, forming=1, strengthening=0, weakening=0)



def test_decoder_rejects_rematching_configuration():
    from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig
    with pytest.raises(TypeError):
        EventDecodeConfig(branch_limit=100)


@pytest.mark.parametrize('threshold', [.3, .5, 1.])
def test_decoder_session_reuses_relation_preparation(threshold, monkeypatch):
    from mappa.mappa_core.post_processing.decoder.api import EventDecodeSession
    from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig
    from mappa.mappa_core.utils.io.records import decoded_events_record
    ep = endpoint('CCC', [(0, 1, 1.), (0, 2, 1.4)])
    fragment = saved_fragment(AAMProblem(ep, ep), source_atoms=(1, 2))
    source = saved(ep, ep, dict(enumerate(range(3))), fragments=(fragment,))
    session = EventDecodeSession(source.problem, EventDecodeConfig(threshold, None))
    first, second = session.decode(source), session.decode(source)
    assert first.family_actions[0][0] is second.family_actions[0][0]
    assert decoded_events_record(first) == decoded_events_record(second)
    assert first.family_actions[0][0].generators is None
    other = saved(endpoint('CC', []), endpoint('CC', []), {0: 0, 1: 1})
    with pytest.raises(ValueError, match='different endpoints'):
        session.decode(other)


def test_ordered_overlapping_relations_preserve_inverse_words_and_joint_choices():
    from test_relaxed_automorphism import broad_words
    ep = endpoint('CCC', [(0, 1, .8), (0, 2, .8), (1, 2, 1.6)])
    value = AAMProblem(ep, ep)
    pair = saved_fragment(value, source_atoms=(1, 2))
    whole = saved_fragment(value)
    source = saved(ep, ep, dict(enumerate(range(3))), fragments=(pair, whole))
    decoded = decode_events(source)
    words = broad_words(whole.relation.target.generators, 3)
    cycle = words[(1, 2, 0)]
    inverse = tuple(-number for number in reversed(cycle))
    selections = [[None, None], [swap(pair), None],
                  [None, dict(target_word=inverse)], [swap(pair), dict(target_word=inverse)]]
    check_contributions(decoded, selections)
    assert len(decoded.family_actions[0]) == 2
    assert decoded.realize(0, selections[-1]) == {0: 2, 1: 1, 2: 0}
    with pytest.raises(ValueError, match='biject'):
        decoded.realize(0, [dict(discovery={1: 1, 2: 1}), None])
    with pytest.raises(ValueError, match='one choice'):
        decoded.realize(0, [])


@pytest.mark.parametrize('partial', [False, True])
def test_boundary_categories_follow_all_fixed_image_choices(partial):
    ep = endpoint('CCCCC', [(0, 2, 1.), (1, 3, 1.8), (2, 4, .8)])
    atoms = (0, 1, 2, 3) if partial else tuple(range(5))
    fragment = saved_fragment(AAMProblem(ep, ep), source_atoms=atoms)
    source = saved(ep, ep, dict(zip(atoms, atoms)), fragments=(fragment,),
                   pairs=tuple(((a,), (a,)) for a in atoms))
    decoded = decode_events(source)
    selections = [[dict(discovery=dict(zip(atoms, order)))] for order in permutations(atoms)]
    check_contributions(decoded, selections)
    assert decoded.tables(0).components == ((tuple(range(len(atoms))), (0,)),)


@pytest.mark.parametrize('moving,expected_case', [((), 'fixed_fixed'),
    (((2, 3),), 'fixed_shuffle'), (((0, 1), (2, 3)), 'shuffle_shuffle')])
def test_cut_endpoint_cases_retain_original_bonds(moving, expected_case):
    r, p = endpoint('CCCC', [(0, 2, 1.)]), endpoint('CCCC', [(0, 2, 1.), (1, 3, 1.)])
    value = AAMProblem(r, p)
    fragments = tuple(saved_fragment(value, source_atoms=atoms) for atoms in moving)
    source = saved(r, p, dict(enumerate(range(4))), fragments=fragments,
                   pairs=(((0, 1), (0, 1)), ((2, 3), (2, 3))))
    source.families[0].provenance.append(dict(cuts=((0, 2),)))
    decoded = decode_events(source)
    term = next(term for term in decoded.tables(0).boundary[0, 1] if term.atoms == (0, 2))
    assert term.manual_cut and term.endpoint_case == expected_case
    assert len(tuple(term.target_pairs)) == 2 ** len(moving)
    check_contributions(decoded, [[None] * len(moving), [swap(fragment) for fragment in fragments]])


def test_binary_matching_uses_original_wbo_and_preserves_complete_relation_on_roundtrip():
    from mappa.mappa_core.post_processing.decoder.api import EventDecodeSession
    from mappa.mappa_core.utils.io.records import decoded_events_record, problem_from_record
    binary = endpoint(('V', 'O', 'O'), [(0, 1, 1.), (0, 2, 1.)])
    raw = AAMProblem(endpoint(('V', 'O', 'O'), [(0, 1, 2.2), (0, 2, .7)]),
                     endpoint(('V', 'O', 'O'), [(0, 1, .8), (0, 2, 2.1)]))
    fragment = saved_fragment(AAMProblem(binary, binary), accepted_edges=((0, 1), (0, 2)))
    source = saved(binary, binary, dict(enumerate(range(3))), fragments=(fragment,))
    selections = [[None], [dict(discovery={0: 0, 1: 2, 2: 1})]]
    decoded = EventDecodeSession(source.problem, event_problem=raw).decode(source)
    assert [decoded.contributions(0, choice)['total'] for choice in selections] == [2, 0]
    check_contributions(decoded, selections)
    record = decoded_events_record(decoded)
    stored = json.loads(json.dumps(record))
    restored_catalogue = FinalBranchCatalogue.from_record(problem_from_record(stored['problem']), stored['catalogue'])
    restored = decode_events(restored_catalogue, event_problem=problem_from_record(stored['event_problem']))
    assert decoded_events_record(restored) == record
    assert record['matching_problem_sha256'] != record['event_problem_sha256']
    assert source.to_record() == restored_catalogue.to_record()
    check_contributions(restored, selections)


def test_binary_matching_keeps_triple_to_single_event():
    binary = endpoint('CO', [(0, 1, 1.)])
    raw = AAMProblem(endpoint('CO', [(0, 1, 3.)]), endpoint('CO', [(0, 1, 1.)]))
    fragment = saved_fragment(AAMProblem(binary, binary), accepted_edges=((0, 1),))
    decoded = decode_events(saved(binary, binary, {0: 0, 1: 1}, fragments=(fragment,)), event_problem=raw)
    assert decoded.contributions(0)['counts'] == dict(breaking=0, forming=0, strengthening=0, weakening=1)
    check_contributions(decoded, [None])


@pytest.mark.parametrize('tolerance', [1., 10.])
def test_discovery_connectivity_is_enforced_before_relaxed_self_actions(tolerance):
    ep = endpoint('COO', [(0, 1, 1.)])
    fragment = saved_fragment(AAMProblem(ep, ep), accepted_edges=((0, 1),),
                              iso_tolerance=tolerance)
    decoded = decode_events(saved(ep, ep, dict(enumerate(range(3))), fragments=(fragment,)))
    with pytest.raises(ValueError, match='accepted edge'):
        decoded.realize(0, [dict(discovery={0: 0, 1: 2, 2: 1})])


def test_partition_views_do_not_duplicate_counts_or_split_shared_choices():
    ep = endpoint('CCCCCC', [(0, 1, 1.), (0, 2, 1.8), (3, 4, 1.), (3, 5, 1.8)])
    value = AAMProblem(ep, ep)
    fragments = tuple(saved_fragment(value, source_atoms=atoms) for atoms in ((1, 2), (4, 5)))
    source = saved(ep, ep, dict(enumerate(range(6))), fragments=fragments,
        pairs=(((0, 1, 2), (0, 1, 2)), ((3, 4, 5), (3, 4, 5))))
    source.coupled_branches[0] = {tuple(((a,), (a,)) for a in range(6))}
    decoded = decode_events(source)
    separate = next(decoded.tables(0, i) for i, parts in enumerate(decoded.fragment_partitions[0])
                    if len(parts) == 2)
    assert separate.components == (((0,), (0,)), ((1,), (1,)))
    for selections in product(*[(None, swap(fragment)) for fragment in fragments]):
        totals = [decoded.contributions(0, selections, partition=i)['counts'] for i in range(2)]
        assert totals[0] == totals[1]
    joint = saved_fragment(value, accepted_edges=((0, 1), (0, 2), (3, 4), (3, 5)))
    source.families[0].actions = (('relation', joint),)
    decoded = decode_events(source)
    split = next(decoded.tables(0, i) for i, parts in enumerate(decoded.fragment_partitions[0]) if len(parts) == 2)
    assert split.components == (((0, 1), (0,)),)


def test_large_boundary_tables_do_not_expand_complete_choices(monkeypatch):
    import sympy.combinatorics as groups
    import z3
    count = 20
    ep = endpoint('C' * (2 * count), [(i, i + 1, 1.8) for i in range(2 * count - 1)])
    value = AAMProblem(ep, ep)
    atoms = tuple((i, i + 1) for i in range(0, 2 * count, 2))
    fragments = tuple(saved_fragment(value, source_atoms=part) for part in atoms)
    source = saved(ep, ep, dict(enumerate(range(2 * count))), fragments=fragments,
                   pairs=tuple((part, part) for part in atoms))
    def forbidden(*args, **kwargs):
        pytest.fail('local decoding expanded complete choices or called a solver')
    monkeypatch.setattr(groups.PermutationGroup, 'generate_schreier_sims', forbidden)
    monkeypatch.setattr(z3, 'Solver', forbidden)
    decoded = decode_events(source)
    monkeypatch.setattr(decoded.index, 'endpoint_pairs', forbidden)
    table = decoded.tables(0)
    assert len(table.boundary) == count - 1
    assert table.statistics['boundary_terms'] == 4 * (count - 1)
    assert table.statistics['projected_pair_states'] == 0
    check_contributions(decoded, [[None] * count, [swap(fragment) for fragment in fragments]])


def test_constant_boundaries_and_possible_product_neighbors_use_bounded_domains(monkeypatch):
    complete = endpoint('C' * 12, [(a, b, 1.) for a in range(12) for b in range(a + 1, 12)])
    fragment = saved_fragment(AAMProblem(complete, complete))
    source = saved(complete, complete, dict(enumerate(range(12))), fragments=(fragment,),
                   pairs=tuple(((a,), (a,)) for a in range(12)))
    source.families[0].provenance.append(dict(cuts=((0, 1),)))
    decoded = decode_events(source)
    monkeypatch.setattr(decoded.index, 'endpoint_pairs', lambda *args: pytest.fail('expanded constant pair bounds'))
    assert decoded.tables(0).statistics['cut_terms'] == 1
    assert decoded.contributions(0)['total'] == 0
    ep = endpoint('CCCCCC', [(0, 2, 1.), (3, 4, .5)])
    fragment = saved_fragment(AAMProblem(ep, ep), source_atoms=(0, 1))
    decoded = decode_events(saved(ep, ep, dict(enumerate(range(6))), fragments=(fragment,),
        pairs=(((0, 1), (0, 1)), ((2,), (2,)), ((3, 4, 5), (3, 4, 5)))))
    assert decoded.tables(0).statistics['neighboring_pairs'] == 2
    assert {term.atoms for rows in decoded.tables(0).boundary.values() for term in rows} == {(0, 2), (1, 2)}
    check_contributions(decoded, [[None], [swap(fragment)]])


def test_normal_decoding_export_keeps_relations_without_whole_reaction_reports(monkeypatch):
    from mappa.mappa_core.post_processing.decoder import events
    from mappa.mappa_core.post_processing.decoder.records import EventInterpretation
    from mappa.mappa_core.utils.io.records import decoded_events_record
    ep = endpoint('CCC', [(0, 1, 1.), (0, 2, 1.8)])
    fragment = saved_fragment(AAMProblem(ep, ep), source_atoms=(1, 2))
    source = saved(ep, ep, dict(enumerate(range(3))), fragments=(fragment,))
    def forbidden(*args, **kwargs):
        pytest.fail('ordinary decoding invoked whole-reaction reporting')
    monkeypatch.setattr(EventInterpretation, 'evaluate', forbidden)
    monkeypatch.setattr(events, 'classify_bonds', forbidden)
    for name in ('describe', '_describe_events', 'counts'):
        monkeypatch.setattr(SignedEventIndex, name, forbidden)
    decoded = decode_events(source)
    record = decoded_events_record(decoded)
    assert record['catalogue'] == source.to_record() and record['mapping_space_preserved']
    assert not record['exhaustive_event_patterns'] and not record['global_minimum_proven']
    assert 'representatives' not in record and record['tables'] and record['responses']
    assert decoded.contributions(0)['total'] == 0
    assert decoded.contributions(0, [swap(fragment)])['counts'] == dict(
        breaking=0, forming=0, strengthening=1, weakening=1)
