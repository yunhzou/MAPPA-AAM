"""Dense event oracles over tiny complete, two-sided fragment relationships."""
from collections import Counter
from dataclasses import replace
from itertools import product

import numpy as np
import pytest

from mappa import AAMProblem, MolecularEndpoint, decode_events
from mappa.mappa_core.core.matching.representation.records import Fragment
from mappa.mappa_core.core.symmetry import DiscoveryEdge, RelaxedIsomorphism, relaxed_automorphism
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue, FinalFamily
from mappa.mappa_core.utils.io.records import decoded_events_record
from test_event_interpretation import model_events
from test_relaxed_isomorphism import discovery_oracle, relation, self_oracle, weights


def problem(source, target):
    def endpoint(matrix):
        return MolecularEndpoint(('C',) * len(matrix), np.zeros((len(matrix), 3)), matrix)
    return AAMProblem(endpoint(source), endpoint(target))


def catalogue(endpoints, fragments, *, mapping=None, partitions=None, cuts=()):
    if mapping is None:
        mapping = dict(fragments[0].witness)
    mapping = tuple(sorted(mapping.items()))
    actions = tuple(('relation', fragment) for fragment in fragments)
    required = tuple(sorted({edge for fragment in fragments for edge in fragment.accepted_edges}))
    result = FinalBranchCatalogue(endpoints)
    result.families.append(FinalFamily(mapping, required, actions, (.2, 0.),
                                      [dict(cuts=cuts)] if cuts else []))
    if partitions is None:
        partitions = (tuple(r for r, _ in mapping),)
    pairs = tuple((atoms, tuple(sorted(dict(mapping)[r] for r in atoms))) for atoms in partitions)
    result.branches[((.2, 0.), pairs)] = {0}
    return result


def choices(value):
    for discovery, source, target in product(discovery_oracle(value), self_oracle(value.source),
                                             self_oracle(value.target)):
        yield dict(discovery=dict(zip(value.source_atoms, discovery)),
                   source_permutation=source, target_permutation=target)


def dense_events(endpoints, mapping, threshold=.5):
    """Literal pair classification, independent of response-table preparation."""
    source, target = endpoints.reactant.wbo, endpoints.product.wbo
    events = []
    for a in range(len(source)):
        for b in range(a + 1, len(source)):
            if a not in mapping and b not in mapping:
                continue
            pair = tuple(sorted((mapping[a], mapping[b]))) if a in mapping and b in mapping else None
            other = target[pair] if pair is not None else 0.
            if source[a, b] - other > threshold:
                category = 'weakening' if other > threshold else 'breaking'
            elif other - source[a, b] > threshold:
                category = 'strengthening' if source[a, b] > threshold else 'forming'
            else:
                continue
            events.append((category, (a, b), pair))
    occupied = set(mapping.values())
    for p in range(len(target)):
        for q in range(p + 1, len(target)):
            if (p in occupied) != (q in occupied) and target[p, q] > threshold:
                events.append(('forming', None, (p, q)))
    return sorted(events, key=repr)


def assert_counts(decoded, selection, *, partition=0, cut=0):
    mapping = decoded.realize(0, selection)
    physical = dense_events(decoded.event_problem, mapping)
    evaluated = decoded.evaluate(0, selection)
    assert evaluated['total'] == len(physical)
    assert evaluated['counts'] == {
        'broken': sum(event[0] in ('breaking', 'weakening') for event in physical),
        'formed': sum(event[0] in ('forming', 'strengthening') for event in physical),
    }
    expected = model_events(decoded.event_problem, dict(decoded.catalogue.families[0].mapping),
        mapping, decoded.fragment_partitions[0][partition], decoded.cut_contexts[0][partition][cut])
    result = decoded.contributions(0, selection, partition=partition, cut=cut, include_events=True)
    actual = [(event['category'], event['reactant_pair'], event['product_pair'])
              for event in result['events']]
    assert sorted(actual, key=repr) == expected
    assert result['total'] == len(expected)
    counts = Counter(event[0] for event in expected)
    assert result['counts'] == {name: counts[name] for name in result['counts']}
    assert sum(sum(value.values()) for collection in ('internal', 'boundary')
               for value in result[collection].values()) == result['total']
    return tuple(mapping[r] for r in sorted(mapping))


def test_global_placement_composes_local_choices_and_scores_original_boundaries():
    from test_final_branches import placed_fragment_path
    endpoints, path = placed_fragment_path()
    saved = FinalBranchCatalogue(endpoints).add_graph(path.graph)
    # Matching uses equivalent components; scoring can distinguish their WBOs.
    scoring = problem(endpoints.reactant.wbo,
        weights(4, [(0, 1, 1.), (2, 3, 1.8), (1, 2, .7)]))
    decoded = decode_events(saved, event_problem=scoring)
    assert len(decoded.inventory.fragments) == 1
    assert decoded.family_actions[0][0].generators is None
    assert not decoded.family_actions[0][1].response_invariant
    seen = set()
    for image, word in product(((0, 1), (1, 0)), ((), (1,))):
        choice = dict(discovery=dict(zip((0, 1), image)))
        seen.add(assert_counts(decoded, (choice, word)))
    assert seen == {(0, 1), (1, 0), (2, 3), (3, 2)}
    ordinary = decode_events(saved)
    assert ordinary.family_actions[0][1].response_invariant
    assert all(ordinary.contributions(0, (None, word))['total'] == 0 for word in ((), (1,)))


def test_earlier_global_placement_moves_later_anchor_continuation_jointly():
    from growth_cases import saved_fragment
    from test_event_patterns import saved_path
    triangle = [(0, 1, 1.), (1, 2, 1.), (0, 2, 1.)]
    matrix = weights(6, triangle + [(a + 3, b + 3, w) for a, b, w in triangle])
    endpoints = problem(matrix, matrix)
    first = saved_fragment(endpoints, source_atoms=(0, 1), accepted_edges=((0, 1),),
                           iso_tolerance=0., self_tolerance=0.)
    later = saved_fragment(endpoints, source_atoms=(1, 2), anchors=((1, 1),),
                           accepted_edges=((1, 2),), iso_tolerance=0., self_tolerance=0.)
    path = saved_path(endpoints, (first, later))
    graph = replace(path.graph, transitions=(replace(path.graph.transitions[0],
                    placement_generators=((3, 4, 5, 0, 1, 2),)), path.graph.transitions[1]))
    path = next(graph.paths())
    decoded = decode_events(FinalBranchCatalogue(endpoints).add_graph(graph))
    assert tuple(kind for kind, _ in decoded.catalogue.families[0].actions) == (
        'relation', 'relation', 'group')
    observed = set()
    for image, word in product(((0, 1), (1, 0)), ((), (1,))):
        selection = (None, dict(discovery=dict(zip((0, 1), image))), word)
        result = assert_counts(decoded, selection)
        expected = tuple(atom + (3 if word else 0) for atom in (*image, 2))
        assert result == expected
        realized = dict(path.realize({0: selection[1]}, {0: word}).mapping)
        assert tuple(realized[r] for r in range(3)) == expected
        observed.add(result)
    assert len(observed) == 4


def test_all_six_relaxed_triangle_tail_choices_decode_without_discovery_refilter(monkeypatch):
    matrix = weights(4, [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.), (0, 3, 1.)])
    value = relation(matrix, matrix, ((0, 1), (1, 2), (0, 3)))
    fragment = Fragment(value)
    saved = catalogue(problem(matrix, matrix), (fragment,))
    saved.families[0].provenance.extend(
        dict(cuts=(cut,), fragment_partition=((0, 1, 2, 3),))
        for cut in ((0, 2), (2, 3)))
    decoded = decode_events(saved)
    monkeypatch.setattr(decoded.index, 'endpoint_pairs',
        lambda *args: pytest.fail('relation decoding expanded pair projections'))
    assert decoded.cut_contexts[0][0] == (((0, 2),), ((2, 3),))
    seen = {assert_counts(decoded, [choice], cut=0) for choice in choices(value)}
    assert {assert_counts(decoded, [choice], cut=1) for choice in choices(value)} == seen
    assert len(seen) == 6 and len(seen - discovery_oracle(value)) == 2
    table = decoded.tables(0)
    # Both histories cut inside the same fragment: its internal responses already
    # count the original bonds, so neither adds a separate boundary correction.
    assert all(not decoded.tables(0, cut=i).boundary for i in range(2))
    assert all(decoded.tables(0, cut=i).statistics['cut_terms'] == 1 for i in range(2))
    assert table.statistics['projected_pair_states'] == 0
    terms = [term.atoms for rows in (*table.internal.values(), *table.boundary.values()) for term in rows]
    assert len(terms) == len(set(terms))
    assert decoded.inventory.fragments == (fragment,)
    assert decoded.inventory.fragments[0] is fragment
    assert decoded.family_actions[0][0].generators is None
    assert not decoded.family_actions[0][0].response_invariant
    assert decoded_events_record(decoded)['families'][0]['action_responses'] == [None]
    selections = [[choice] for choice in choices(value)]
    actual = list(decoded.counts(0, selections, categories=True, batch_size=3))
    expected = [list(decoded.contributions(0, selection)['counts'].values()) for selection in selections]
    assert np.array_equal(actual, expected)


def test_self_tolerance_remains_joint_and_uses_matching_weights_with_separate_scoring():
    matrix = weights(4, [(0, 1, 1.), (0, 2, 1.6), (0, 3, 2.2)])
    value = relation(matrix, matrix, ((0, 1), (0, 2), (0, 3)), self_tolerance=.7)
    binary = (matrix > .2).astype(float)
    decoded = decode_events(catalogue(problem(matrix, matrix), (Fragment(value),)),
                            event_problem=problem(binary, binary))
    seen = {assert_counts(decoded, [choice]) for choice in choices(value)}
    assert len(seen) == 5 and (0, 3, 2, 1) not in seen
    with pytest.raises(ValueError, match='self tolerance'):
        decoded.realize(0, [dict(target_permutation=(0, 3, 2, 1))])
    with pytest.raises(ValueError, match='biject'):
        decoded.realize(0, [dict(discovery={0: 0, 1: 1, 2: 1, 3: 3})])


def test_shared_history_keeps_final_partition_bond_ownership_and_joint_occupancy():
    matrix = weights(4, [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.), (0, 3, 1.)])
    whole = Fragment(relation(matrix, matrix, ((0, 1), (1, 2), (0, 3))))
    local = relaxed_automorphism('CCCC', matrix, atoms=(0, 1), tolerance=0.)
    previous = Fragment(RelaxedIsomorphism(local, local, ((0, 0), (1, 1)), (),
        ((0, (0, 1)), (1, (0, 1))), (DiscoveryEdge((0, 1), 1., 0.),), ((0, 1),), .2))
    saved = catalogue(problem(matrix, matrix), (whole, previous), cuts=((0, 2),))
    saved.coupled_branches[0] = {(((0, 1), (0, 1)), ((2, 3), (2, 3)))}
    saved.families.append(replace(saved.families[0]))
    next(iter(saved.branches.values())).add(1)
    decoded = decode_events(saved)
    assert decoded.inventory.fragments == (previous, whole)
    assert decoded.family_actions[0][0] is decoded.family_actions[1][0]
    common = decoded.fragment_partitions[0].index(decoded.fragment_partitions[1][0])
    assert decoded.tables(0, common) is decoded.tables(1, 0)
    assert len(decoded.fragment_partitions[0]) == 2
    for first, second in product(choices(whole.relation), choices(previous.relation)):
        one = assert_counts(decoded, [first, second], partition=0)
        assert_counts(decoded, [first, second], partition=1)
        assert len(set(one)) == 4
    split = next(decoded.tables(0, i) for i, part in enumerate(decoded.fragment_partitions[0])
                 if len(part) == 2)
    assert split.components == (((0, 1), (0, 1)),)
    assert split.statistics['cut_terms'] == 1


def test_sparse_fixed_image_preserves_anchors_and_counts_unmatched_boundaries():
    source = weights(5, [(0, 2, 1.), (2, 4, 1.4), (1, 3, 2.), (0, 1, .8)])
    target = weights(6, [(1, 3, 1.), (3, 5, 1.6), (0, 2, 2.), (3, 4, 1.1)])
    left = relaxed_automorphism('C' * 5, source, atoms=(0, 2, 4), anchors=(0,), tolerance=.7)
    right = relaxed_automorphism('C' * 6, target, atoms=(1, 3, 5), anchors=(1,), tolerance=.7)
    value = RelaxedIsomorphism(left, right, ((0, 1), (2, 3), (4, 5)), ((0, 1),),
        ((0, (1,)), (2, (3, 5)), (4, (3, 5))),
        (DiscoveryEdge((0, 2), 1., .5),), ((1, 3), (3, 5)), .2)
    decoded = decode_events(catalogue(problem(source, target), (Fragment(value),),
        partitions=((0,), (2, 4)), cuts=((0, 2),)))
    for choice in choices(value):
        assert_counts(decoded, [choice])
        assert decoded.realize(0, [choice])[0] == 1
    assert decoded.tables(0).statistics['cut_terms'] == 1
    assert decoded.contributions(0)['total'] == 2
    with pytest.raises(ValueError, match='atom domain|anchor'):
        decoded.realize(0, [dict(discovery={0: 3, 2: 1, 4: 5})])


def test_relation_pair_inspection_is_lazy_and_does_not_claim_independent_permissions():
    matrix = weights(4, [(0, 1, 1.), (0, 2, 1.6), (0, 3, 2.2)])
    value = relation(matrix, matrix, ((0, 1), (0, 2), (0, 3)), self_tolerance=.7)
    decoded = decode_events(catalogue(problem(matrix, matrix), (Fragment(value),)))
    term = next(iter(decoded.tables(0).internal[0]))
    pairs = term.target_pairs
    assert iter(pairs) is pairs
    assert len(tuple(pairs)) > 0
    # Conservative lookup domains include an atom motion that the joint self
    # constraint forbids; no independent shuffle permission is inferred.
    assert 3 in decoded.tables(0).atom_domains[1]
    with pytest.raises(ValueError, match='self tolerance'):
        decoded.realize(0, [dict(source_permutation=(0, 3, 2, 1))])


def test_large_symmetric_fragment_prepares_without_solver_or_permutation_expansion(monkeypatch):
    import sympy.combinatorics as groups
    import z3
    matrix = weights(33, [(0, i, 1.) for i in range(1, 33)])
    domains = ((0, (0,)),) + tuple((i, tuple(range(1, 33))) for i in range(1, 33))
    value = relation(matrix, matrix, tuple((0, i) for i in range(1, 33)), domains=domains)
    saved = catalogue(problem(matrix, matrix), (Fragment(value),))
    def forbidden(*args, **kwargs):
        pytest.fail('local event preparation expanded permutations or invoked a solver')
    monkeypatch.setattr(groups.PermutationGroup, 'generate_schreier_sims', forbidden)
    monkeypatch.setattr(groups.PermutationGroup, 'generate_dimino', forbidden)
    monkeypatch.setattr(z3, 'Solver', forbidden)
    decoded = decode_events(saved)
    monkeypatch.setattr(decoded.index, 'endpoint_pairs', forbidden)
    assert decoded.contributions(0)['total'] == 0
    assert decoded.tables(0).statistics['neighboring_pairs'] == 32
    assert decoded.tables(0).statistics['projected_pair_states'] == 0


def test_saved_path_replays_compiled_two_sided_choices_after_json_roundtrip():
    import json
    from mappa.mappa_core.core.search.graph import AAMSearchGraph
    from mappa.mappa_workflow.mapping_analytics.family import query_path
    from test_event_patterns import saved_path

    matrix = weights(4, [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.), (0, 3, 1.)])
    value = relation(matrix, matrix, ((0, 1), (1, 2), (0, 3)))
    endpoints = problem(matrix, matrix)
    path = saved_path(endpoints, (Fragment(value),))
    restored = next(AAMSearchGraph.from_record(json.loads(json.dumps(path.graph.to_record()))).paths())
    reference = {0: 1, 1: 3, 2: 0, 3: 2}
    assert not value.allows_discovery(reference)
    status, result = query_path(restored, endpoints, reference, source_atoms=tuple(range(4)))
    assert status == 'recovered'
    selections = {action['transition']: {name: action[name] for name in (
        'discovery', 'source_permutation', 'target_permutation')} for action in result['actions']}
    assert dict(restored.realize(selections).mapping) == reference
    assert dict(path.realize(selections).mapping) == reference


def test_search_graph_restores_shared_source_preparation_across_distinct_fragments():
    import json
    from mappa.mappa_core.core.search.graph import AAMSearchGraph
    from test_event_patterns import saved_path

    matrix = weights(3, [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.)])
    value = relation(matrix, matrix, ((0, 1), (1, 2)))
    other = replace(value, witness=((0, 1), (1, 2), (2, 0)))
    path = saved_path(problem(matrix, matrix), (Fragment(value), Fragment(other)))
    graph = AAMSearchGraph.from_record(json.loads(json.dumps(path.graph.to_record())))
    first, second = (edge.match for edge in graph.transitions)
    assert first.witness != second.witness
    assert first.source_self_automorph is second.source_self_automorph
    assert first.target_self_automorph is second.target_self_automorph
