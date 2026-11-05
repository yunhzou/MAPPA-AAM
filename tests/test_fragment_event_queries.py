"""Independent small-graph oracles for retained symmetry/event selections."""
from collections import Counter
from itertools import permutations, product
from dataclasses import replace
import json
import pickle

import numpy as np
import pytest

from mappa import decode_events, decoded_events_record, decoded_events_from_record
from mappa.mappa_core.core.symmetry.minimum import FragmentSymmetryEvents
from mappa.mappa_core.core.symmetry.program import OrderedSymmetryEvents
from test_fragment_optimization import fragment, weights
from test_event_interpretation import saved, endpoint
from growth_cases import saved_fragment
from mappa import AAMProblem


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
    bound = OrderedSymmetryEvents(value.witness, ((('relation', value), data),), right, (),
                                  projections=((data, (0,)),))
    for response in (data, bound):
        for count in range(7):
            for image in range(4):
                answer = response.select(count).find(assignments={0: image})
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


def test_composed_arm_dependence_uses_the_prepared_graph(monkeypatch):
    from sympy.combinatorics import PermutationGroup
    from mappa.mappa_core.core.symmetry.groups import independent_factors

    value = arms()
    prepared = FragmentSymmetryEvents(value)
    response = OrderedSymmetryEvents(value.witness,
        ((('relation', value), prepared),), value.target_symmetry.weights, (),
        projections=((prepared, (0,)),))
    for generator in response.group.generators:
        image = generator.array_form
        assert np.array_equal(response.factor_graph, response.factor_graph[np.ix_(image, image)])

    def unnecessary_stabilizer(*args, **kwargs):
        pytest.fail('the arm graph already proves the carbon/hydrogen orbits are correlated')
    with monkeypatch.context() as patch:
        patch.setattr(PermutationGroup, 'pointwise_stabilizer', unnecessary_stabilizer)
        factors = independent_factors(response.group, response.factor_graph)
    assert len(factors) == 1
    assert response.select(2).find().event_count == 2


def test_composed_placement_reuses_only_a_preserved_graph(monkeypatch):
    from sympy.combinatorics import PermutationGroup
    from mappa.mappa_core.core.symmetry.groups import independent_factors

    value = arms()
    graph = (value.target_symmetry.weights > .2).astype(float)
    action = ('group', value.target_symmetry.generators)
    response = OrderedSymmetryEvents(value.witness, ((action, None),), graph, ())
    def unnecessary_stabilizer(*args, **kwargs):
        pytest.fail('the preserved placement graph already proves the arm dependence')
    with monkeypatch.context() as patch:
        patch.setattr(PermutationGroup, 'pointwise_stabilizer', unnecessary_stabilizer)
        assert len(independent_factors(response.group, response.factor_graph)) == 1
    assert response.binary is None  # Preserve the established branch search ordering.
    for generator in response.group.generators:
        image = generator.array_form
        assert np.array_equal(response.factor_graph, response.factor_graph[np.ix_(image, image)])

    # Variable bond weights are not a placement certificate, but their complete
    # original permutation group remains available to the same event solver.
    variable = OrderedSymmetryEvents(value.witness, ((action, None),),
                                     value.target_symmetry.weights, ())
    assert not np.any(variable.factor_graph)
    assert variable.group.equals(response.group)


@pytest.mark.parametrize('seed', range(24))
def test_pair_minimum_reuses_ordered_binary_choices_without_group_search(seed, monkeypatch):
    from mappa.mappa_core.core.symmetry.actions import apply_action
    from mappa.mappa_core.core.symmetry.minimum import _Search

    rng = np.random.default_rng(seed)
    left, right = [weights(4, [(a,b,float(rng.choice((0., 1., 2.))))
                              for a in range(4) for b in range(a+1,4)]) for _ in range(2)]
    saved_value = fragment(left, right, (), range(4))
    # Preparation may choose another discovery. It is a fixed change of the
    # starting map, not an extra group whose closure is permitted.
    prepared_value = replace(saved_value, relation=replace(saved_value.relation,
        witness=tuple(enumerate(map(int, rng.permutation(4))))))
    prepared = FragmentSymmetryEvents(prepared_value)
    moves = ((1,0,2,3), (0,2,1,3))
    actions = ((('relation', saved_value), prepared),) + tuple(
        (('group', (move,)), None) for move in moves)
    excluded = frozenset(((0,1), (2,3)))
    response = OrderedSymmetryEvents(saved_value.witness, actions, right, (),
        fixed_events=3, sweep=((0,2),1.,.5,1,excluded))
    expected = {}
    for images, first, second in product(oracle(prepared_value), (False,True), (False,True)):
        for selected, move in zip((first,second), moves):
            if selected:
                images = tuple(move[p] for p in images)
        pair = tuple(sorted((images[0], images[2])))
        expected[images] = 3 + int(abs(1. - right[pair]) > .5) + int(pair not in excluded)

    monkeypatch.setattr(_Search, 'run', lambda *args: pytest.fail('sweep rebuilt a full symmetry search'))
    answer = response.solve()
    images = tuple(dict(answer.witness)[a] for a in range(4))
    assert answer.complete and answer.event_count == answer.lower_bound == min(expected.values())
    assert expected[images] == answer.event_count
    assert 0 < answer.work.projected_pair_states <= 5 * 4 * 3 and answer.work.nodes == 0
    realized = np.arange(4)
    for (action, _), choice in zip(actions, answer.choices):
        realized = apply_action(realized, action, choice, 4)
    assert tuple(realized) == images

    # A component reached after other choices were fixed must retain those
    # choices. Compare its exact coset with the independently enumerated group.
    operations = list(response.group.generate_schreier_sims())
    prefix = np.asarray(operations[-1].array_form)
    subgroup = response.group.stabilizer(seed % response.size)
    cost, selected, _ = response.project_pair(prefix, response.pairs, subgroup)
    allowed = [prefix[np.asarray(g.array_form)] for g in subgroup.generate_schreier_sims()]
    counts = [response.score(p, response.pairs) for p in allowed]
    assert cost == min(counts) == response.score(selected, response.pairs)
    assert any(np.array_equal(selected, p) for p in allowed)


def test_pair_projection_keeps_a_cross_action_constraint_joint():
    from mappa.mappa_core.core.symmetry.groups import permutation_group
    target = weights(3, [(0,1,1.)])
    actions = ((('group', ((1,0,2),)), None), (('group', ((0,2,1),)), None))
    response = OrderedSymmetryEvents(tuple(enumerate(range(3))), actions, target,
                                    (((0,1), 1., .5, None),))
    # Choosing both swaps together is a subgroup, but not a direct product of
    # the action frames. Pair projection must not introduce either swap alone.
    joined = permutation_group([(1,0,3,2)],4)
    assert response.project_pair(np.arange(4), response.pairs, joined) is None


def test_native_pair_minimum_preserves_pool_occupancy_and_fixed_action_order():
    from mappa._group_ops import PairAction, EndpointPairProgram
    pool = (0,1,2)
    fixed = (1,2,3,0)
    group = (0,1,3,2)
    program = EndpointPairProgram(4, [PairAction.pool(4,pool),
        PairAction.permutation(fixed), PairAction.group(4,[group])])
    costs = np.array([[3,2,1,0], [0,3,2,1], [1,0,3,2], [2,1,0,3]], dtype=np.int32)
    for a,b in product(range(4), repeat=2):
        expected = set()
        for permutation, selected in product(permutations(pool), (False,True)):
            p = tuple(permutation) + (3,)
            images = [fixed[p[a]], fixed[p[b]]]
            if selected:
                images = [group[p] for p in images]
            expected.add(tuple(images))
        assert set(program.project((a,b))) == expected
        cost, choices, visited = program.minimum((a,b),costs)
        images = (a,b)
        for choice in choices:
            assert sorted(choice) == list(range(4))
            images = tuple(choice[p] for p in images)
        assert images in expected and cost == costs[images] == min(costs[p] for p in expected)
        assert visited <= 3 * 4 * 4


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


@pytest.mark.parametrize('seed', range(12))
def test_ordered_fragment_projection_and_sweep_against_all_small_choices(seed):
    # Empty matching connectivity permits all shuffles. Event weights are
    # deliberately different, so a placement cannot reuse constant costs.
    rng = np.random.default_rng(seed)
    matching = endpoint('CCCC', [])
    value = saved_fragment(AAMProblem(matching, matching), source_atoms=(0, 1, 2))
    catalogue = saved(matching, matching, dict(enumerate(range(4))), fragments=(value,),
                      pairs=(((0, 1), (0, 1)), ((2, 3), (2, 3))))
    shifts = ((2, 3, 0, 1), (0, 2, 1, 3))
    catalogue.families[0] = replace(catalogue.families[0],
        actions=catalogue.families[0].actions + tuple(('group', (g,)) for g in shifts),
        provenance=[dict(cuts=((0, 2),))])
    endpoints = [endpoint('CCCC', [(a, b, float(rng.choice((0., .2, .8, 1.6))))
                 for a in range(4) for b in range(a + 1, 4)]) for _ in range(2)]
    decoded = decode_events(catalogue, event_problem=AAMProblem(*endpoints), event_slack=20)
    table, expected = decoded.tables(0), {}
    for local, move1, move2 in product(permutations(range(3)), (False, True), (False, True)):
        images = list(local) + [3]
        for move, shift in zip((move1, move2), shifts):
            if move:
                images = [shift[p] for p in images]
        # Independent scalar oracle for the agreed ordinary-boundary model.
        total = table.boundary_count
        for a, b in ((0, 1), (2, 3), (0, 2)):
            total += abs(endpoints[0].wbo[a, b] - endpoints[1].wbo[images[a], images[b]]) > .5
        total -= int(endpoints[0].wbo[0, 2] > .5)
        total -= int(tuple(sorted((images[0], images[2]))) in table.product_boundary_set)
        expected[tuple(images)] = int(total)
    program = table.program
    for operation in program.group.generate_schreier_sims():
        prefix = np.asarray(operation.array_form)
        images = tuple(program.mapping(prefix)[a] for a in range(4))
        assert program.fixed_events + program.score(prefix, program.pairs) == expected[images]
    assert decoded.minimum_events == min(expected.values()) and decoded.minimum_proven
    for count in range(max(expected.values()) + 2):
        for image in range(4):
            answer = decoded.select(0, count).find(assignments={0: image})
            allowed = {p for p, cost in expected.items() if cost == count and p[0] == image}
            assert answer.status == ('found' if allowed else 'absent')
            if allowed:
                assert tuple(dict(answer.witness)[a] for a in range(4)) in allowed
