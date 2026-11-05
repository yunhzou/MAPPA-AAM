from mappa.mappa_workflow.mapping_analytics.queries import query_assignments
"""Anchored conversion preserves fixed atoms, compressed choices and counts."""
from dataclasses import replace

import numpy as np
import pytest
from growth_cases import saved_fragment

from mappa import AAMProblem, MolecularEndpoint, convert_fragments, decode_events
from mappa.mappa_core.utils.io.records import refined_mapping_record
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue, FinalFamily
from mappa.mappa_core.post_processing.hierarchical_mapping import RefinementConfig
from mappa.mappa_core.post_processing.hierarchical_mapping import refine_mapping


def endpoint(elements, edges):
    matrix = np.zeros((len(elements), len(elements)))
    for a, b, weight in edges:
        matrix[a, b] = matrix[b, a] = weight
    return MolecularEndpoint(tuple(elements), np.zeros((len(elements), 3)), matrix)


def saved(reactant, product, mapping, *, anchors=None, pairs=None):
    catalogue = FinalBranchCatalogue(AAMProblem(reactant, product))
    mapping = tuple(sorted(mapping.items()))
    anchors = mapping if anchors is None else tuple(anchors)
    fragment = saved_fragment(catalogue.problem,
        source_atoms=tuple(r for r, _ in mapping), target_atoms=tuple(sorted(p for _, p in mapping)),
        witness=mapping, anchors=anchors)
    catalogue.families.append(FinalFamily(mapping, (), (('relation', fragment),), (.2, 1.)))
    pairs = pairs if pairs is not None else ((tuple(r for r, _ in mapping),
                                            tuple(sorted(p for _, p in mapping))),)
    catalogue.branches[((.2, 1.), pairs)] = {0}
    catalogue.path_count = 1
    return catalogue


@pytest.fixture(params=['0', '1'])
def backend(request, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', request.param)


@pytest.mark.smoke
def test_free_mapping_is_rediscovered_and_input_untouched(backend):
    ep = endpoint('OCC', [(0, 1, 1.), (0, 2, 2.)])
    catalogue = saved(ep, ep, {0: 0, 1: 2, 2: 1}, anchors=((0, 0),))
    before = catalogue.to_record()
    decoded = refine_mapping(catalogue)
    assert decoded.complete and decoded.observed_minimum == 0
    assert decoded.candidates[0].mapping == {0: 0, 1: 1, 2: 2}
    assert decoded.conversion.rematches[0].anchors == ((0, 0),)
    assert catalogue.to_record() == before


@pytest.mark.smoke
def test_direct_anchors_keep_a_bond_at_the_exact_threshold(backend):
    r = endpoint('CO', [(0, 1, 1.)])
    p = endpoint('CO', [(0, 1, 1.5)])
    decoded = refine_mapping(saved(r, p, {0: 0, 1: 1}))
    assert decoded.complete and decoded.observed_minimum == 0
    assert decoded.candidates[0].endpoint_events['formed'] == ()
    owners = dict(decoded.conversion.branches[0].islands)
    assert owners[0] == owners[1]
    assert decoded.conversion.rematches[0].anchors == ((0, 0), (1, 1))


@pytest.mark.smoke
def test_native_conditional_metal_rule_is_used(backend):
    r = endpoint(('Co', 'C', 'O'), [(0, 1, 1.), (1, 2, 1.)])
    p = endpoint(('Co', 'C', 'O'), [(0, 1, 1.4), (1, 2, 1.4)])
    decoded = refine_mapping(saved(r, p, dict(enumerate(range(3)))))
    assert decoded.complete and decoded.observed_minimum == 1
    owners = dict(decoded.conversion.branches[0].islands)
    assert owners[0] != owners[1] == owners[2]


def test_free_hydrogens_stay_compressed_and_queryable(backend):
    ep = endpoint('CHHH', [(0, h, 1.) for h in (1, 2, 3)])
    decoded = refine_mapping(saved(ep, ep, {0: 0, 1: 2, 2: 3, 3: 1},
                                  anchors=((0, 0),)))
    assert decoded.complete and len(decoded.conversion.branches) == 1
    assert decoded.catalogue.families[0].actions
    candidate = decoded.candidates[0]
    assert query_assignments(decoded, candidate, {1: 3, 2: 1, 3: 2}, same_events=False).status == 'allowed'
    assert query_assignments(decoded, candidate, {1: 3, 2: 3}, same_events=False).status == 'forbidden'


def test_jointly_moving_atoms_are_not_fixed(backend):
    ep = endpoint('CCC', [(0, 1, 1.), (1, 2, 1.), (2, 0, 1.)])
    decoded = refine_mapping(saved(ep, ep, dict(enumerate(range(3))),
        anchors=()))
    assert decoded.complete
    assert decoded.conversion.rematches[0].anchors == ()


def test_grouped_pairs_are_reused_but_different_anchors_are_separate(backend):
    ep = endpoint('OCC', [(0, 1, 1.), (0, 2, 2.)])
    catalogue = saved(ep, ep, {0: 0, 1: 2, 2: 1}, anchors=((0, 0),))
    first = catalogue.families[0]
    catalogue.families += [replace(first, mapping=((0, 0), (1, 1), (2, 2))), replace(first, actions=())]
    next(iter(catalogue.branches.values())).update((1, 2))
    result = convert_fragments(catalogue)
    assert len(result.rematches) == 2 and result.reused_pairs == 0
    assert len(result.branches) == 2
    assert [b.source_family_ids for b in result.branches] == [(0, 1), (2,)]
    assert result.catalogue.families[0].provenance[0]['source_family_ids'] == (0, 1)


def test_conversion_retains_global_placement_coverage_and_program_identity(backend):
    from itertools import permutations
    from test_final_branches import placed_fragment_path
    from mappa.mappa_core.core.search.queries.family import query_path
    value, path = placed_fragment_path()
    catalogue = FinalBranchCatalogue(value).add_graph(path.graph)
    first = catalogue.families[0]
    catalogue.families.append(replace(first, actions=first.actions[:1]))
    next(iter(catalogue.branches.values())).add(1)
    converted = convert_fragments(catalogue)
    assert len(converted.catalogue.families) == 2
    assert [branch.source_family_ids for branch in converted.branches] == [(0,), (1,)]
    images = []
    for family in converted.catalogue.families:
        restored = family.as_path(value)
        images.append({image for image in permutations(range(4), 2)
            if query_path(restored, value, dict(zip((0, 1), image)),
                          source_atoms=(0, 1))[0] == 'recovered'})
        assert all(len(generator) == 4 for kind, group in family.actions if kind == 'group'
                   for generator in group)
    assert images == [{(0, 1), (1, 0), (2, 3), (3, 2)}, {(0, 1), (1, 0)}]
    assert all(decode_events(converted.catalogue).contributions(i)['total'] == 0
               for i in range(2))


@pytest.mark.smoke
def test_partial_boundary_and_wholly_unmatched_material(backend):
    r = endpoint('COHH', [(0, 1, 1), (2, 3, 1)])
    p = endpoint('C', [])
    for left, right in ((r, p), (p, r)):
        decoded = refine_mapping(saved(left, right, {0: 0}))
        assert decoded.complete and decoded.observed_minimum == 1
        assert decoded.candidates[0].mapping == {0: 0}


@pytest.mark.smoke
def test_rematching_stays_inside_original_pairs_and_recombines(backend):
    ep = endpoint('HHHH', [])
    pairs = (((0, 1), (2, 3)), ((2, 3), (0, 1)))
    source = saved(ep, ep, {0: 2, 1: 3, 2: 0, 3: 1}, pairs=pairs,
        anchors=())
    result = convert_fragments(source)
    assert result.complete
    for family in result.catalogue.families:
        mapping = dict(family.mapping)
        assert {mapping[0], mapping[1]} == {2, 3}
        assert {mapping[2], mapping[3]} == {0, 1}
        assert len(set(mapping.values())) == 4


def test_filtered_minimum_is_over_converted_branches(backend):
    ep = endpoint('OCC', [(0, 1, 1.), (0, 2, 2.)])
    catalogue = saved(ep, ep, {0: 0, 1: 2, 2: 1}, anchors=((0, 0),))
    catalogue.families.append(replace(catalogue.families[0], actions=()))
    next(iter(catalogue.branches.values())).add(1)
    decoded = refine_mapping(catalogue, RefinementConfig(event_slack=0))
    assert decoded.complete and decoded.minimum_proven and decoded.event_window == 0
    assert [c.total for c in decoded.candidates] == [0]
    assert sorted(r['total'] for r in decoded.family_reports) == [0, 2]


def test_empty_and_missing_fragment_catalogues_are_explicit():
    ep = endpoint('C', [])
    empty = refine_mapping(FinalBranchCatalogue(AAMProblem(ep, ep)))
    assert not empty.candidates and not empty.minimum_proven
    catalogue = saved(ep, ep, {0: 0})
    catalogue.branches.clear()
    with pytest.raises(ValueError, match='fragment atom sets'):
        refine_mapping(catalogue)


def test_no_decoder_solver_or_old_budget_option(backend, monkeypatch):
    import z3
    def forbidden(*args, **kwargs):
        pytest.fail('conversion invoked a solver')
    monkeypatch.setattr(z3, 'Solver', forbidden)
    ep = endpoint('CH', [(0, 1, 1.)])
    decoded = refine_mapping(saved(ep, ep, {0: 0, 1: 1}))
    assert decoded.complete
    assert all(row['solver_queries'] == 0 for row in decoded.family_reports)
    with pytest.raises(TypeError):
        RefinementConfig(seconds_per_family=10)


def test_artifact_preserves_converted_compression_and_scope(backend):
    ep = endpoint('CHH', [(0, 1, 1.), (0, 2, 1.)])
    decoded = refine_mapping(saved(ep, ep, dict(enumerate(range(3))), anchors=((0, 0),)))
    record = refined_mapping_record(decoded)
    assert record['schema'] == 'mappa.decoded_events/v2'
    assert record['scope'] == 'Direct counts across retained converted branches'
    restored = FinalBranchCatalogue.from_record(decoded.catalogue.problem, record['converted_catalogue'])
    assert restored.families[0].actions == decoded.catalogue.families[0].actions
    assert record['rematches'][0]['search']


@pytest.mark.smoke
def test_incomplete_rematch_does_not_prove_the_minimum(backend, monkeypatch):
    import mappa.mappa_core.post_processing.hierarchical_mapping as module
    from mappa.mappa_core.core.search.graph import SearchStop
    original = module.find_islands
    def capped(*args, **kwargs):
        graph = original(*args, **kwargs)
        return replace(graph, stops=graph.stops + (SearchStop(graph.roots[0], 'capped'),))
    monkeypatch.setattr(module, 'find_islands', capped)
    ep = endpoint('C', [])
    decoded = refine_mapping(saved(ep, ep, {0: 0}))
    assert decoded.observed_minimum == 0 and not decoded.minimum_proven
    assert not decoded.complete


def test_grouping_precedes_recombination_and_result_return(backend, monkeypatch):
    import mappa.mappa_core.post_processing.hierarchical_mapping as module
    ep = endpoint('CHH', [(0, 1, 1.), (0, 2, 1.)])
    catalogue = saved(ep, ep, {0: 0, 1: 1, 2: 2}, anchors=((0, 0),))
    baseline = refine_mapping(catalogue)
    catalogue.families.extend(replace(catalogue.families[0]) for _ in range(99))
    next(iter(catalogue.branches.values())).update(range(100))
    original = module.product
    calls = []
    def counted(*args):
        calls.append(1)
        return original(*args)
    monkeypatch.setattr(module, 'product', counted)
    result = refine_mapping(catalogue)
    assert len(calls) == 1
    assert len(result.conversion.branches) == len(result.catalogue.families) == 1
    assert result.conversion.branches[0].source_family_ids == tuple(range(100))
    assert result.complete == baseline.complete and result.observed_minimum == baseline.observed_minimum
    assert [(c.id, c.mapping, c.total) for c in result.candidates] == [(c.id, c.mapping, c.total) for c in baseline.candidates]
    record = refined_mapping_record(result)
    assert record['converted_branches'][0]['source_family_ids'] == tuple(range(100))
    assert len(record['family_reports']) == 1


def test_fragment_order_is_irrelevant_but_atom_pools_and_policy_are_not(backend):
    ep = endpoint('HHHH', [])
    pairs = (((0, 1), (0, 1)), ((2, 3), (2, 3)))
    catalogue = saved(ep, ep, dict(enumerate(range(4))), pairs=pairs,
                      anchors=())
    catalogue.families.append(replace(catalogue.families[0]))
    catalogue.branches[((.2, 1.), tuple(reversed(pairs)))] = {1}
    catalogue.families.append(replace(catalogue.families[0]))
    other_pairs = (((0, 2), (0, 2)), ((1, 3), (1, 3)))
    catalogue.branches[((.2, 1.), other_pairs)] = {2}
    catalogue.families.append(replace(catalogue.families[0], policy=(.4, 1.)))
    catalogue.branches[((.4, 1.), pairs)] = {3}
    result = convert_fragments(catalogue)
    assert {b.source_family_ids for b in result.branches} == {(0, 1), (2,), (3,)}


def test_same_pairs_in_coupled_and_regular_inputs_share_conversion(backend):
    ep = endpoint('CHH', [(0, 1, 1.), (0, 2, 1.)])
    catalogue = saved(ep, ep, {0: 0, 1: 1, 2: 2}, anchors=((0, 0),))
    pairs = next(iter(catalogue.branches))[1]
    catalogue.families.append(replace(catalogue.families[0]))
    catalogue.coupled_branches[1] = {pairs}
    result = convert_fragments(catalogue)
    assert len(result.branches) == 1
    assert result.branches[0].source_family_ids == (0, 1)


def test_equal_discovery_tolerance_does_not_prove_relaxed_choices_need_no_rematch():
    # D and the two self relationships have separate meanings. Matching a
    # numerical tolerance cannot justify reusing this larger retained relation.
    ep = endpoint('CCC', [(0, 1, 1.)])
    source = saved(ep, ep, dict(enumerate(range(3))), anchors=())
    before = source.to_record()
    result = convert_fragments(source, threshold=1., metal_threshold=1.)
    assert result.rematches and result.catalogue is not source
    assert source.to_record() == before


@pytest.mark.parametrize('threshold,metal_threshold', [(1., .3), (.5, 1.)])
def test_both_tolerances_must_match_to_reuse(backend, threshold, metal_threshold):
    ep = endpoint(('Co', 'O'), [(0, 1, 1.)])
    result = convert_fragments(saved(ep, ep, {0: 0, 1: 1}),
                               threshold=threshold, metal_threshold=metal_threshold)
    assert result.rematches


@pytest.mark.parametrize('elements,limit', [('CO', .5), (('Co', 'O'), .3)])
@pytest.mark.parametrize('offset', [-1, 0, 1])
def test_inclusive_conversion_boundary(backend, elements, limit, offset):
    # Exactly representable difference at equality; adjacent floats stay distinct.
    lower = limit
    upper = 2 * limit
    if offset:
        upper = np.nextafter(upper, -np.inf if offset < 0 else np.inf)
    r = endpoint(elements, [(0, 1, lower)])
    p = endpoint(elements, [(0, 1, upper)])
    result = refine_mapping(saved(r, p, {0: 0, 1: 1}))
    expected = int(abs(upper - lower) > limit)
    assert result.complete and result.observed_minimum == expected
    owners = dict(result.conversion.branches[0].islands)
    assert (owners[0] == owners[1]) == (not expected)
    assert result.index.policy['event_comparison'] == 'absolute difference > threshold'


def test_alternative_partitions_remain_distinct_under_explicit_rematching():
    ep = endpoint('HH', [])
    source = saved(ep, ep, {0: 0, 1: 1}, anchors=())
    source.coupled_branches[0] = {(((0,), (0,)), ((1,), (1,)))}
    result = refine_mapping(source, RefinementConfig(threshold=1., metal_threshold=1.))
    # The two-atom rematch keeps both placements in one global orbit. The
    # separate singleton partition contributes its own fixed placement.
    assert len(result.conversion.branches) == len(result.family_reports) == 2
    assert result.candidates[0].family_ids == [0, 1]
    assert result.catalogue is not source
    assert len({branch.islands for branch in result.conversion.branches}) == 2
    from mappa.mappa_core.core.search.queries.family import query_path
    represented = [{image for image in ((0, 1), (1, 0)) if query_path(
        family.as_path(source.problem), source.problem, dict(enumerate(image)),
        source_atoms=(0, 1))[0] == 'recovered'} for family in result.catalogue.families]
    assert represented == [{(0, 1), (1, 0)}, {(0, 1)}]
