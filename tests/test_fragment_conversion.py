"""Anchored conversion preserves fixed atoms, compressed choices and counts."""
from dataclasses import replace

import numpy as np
import pytest

from mappa import AAMProblem, MolecularEndpoint, convert_fragments
from mappa.artifacts import decoded_events_record
from mappa.final_branches import FinalBranchCatalogue, FinalFamily
from mappa.postprocessing import EventDecodeConfig, decode_events


def endpoint(elements, edges):
    matrix = np.zeros((len(elements), len(elements)))
    for a, b, weight in edges:
        matrix[a, b] = matrix[b, a] = weight
    return MolecularEndpoint(tuple(elements), np.zeros((len(elements), 3)), matrix)


def saved(reactant, product, mapping, *, actions=(), pairs=None):
    catalogue = FinalBranchCatalogue(AAMProblem(reactant, product))
    mapping = tuple(sorted(mapping.items()))
    catalogue.families.append(FinalFamily(mapping, (), actions, (.2, 1.)))
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
    catalogue = saved(ep, ep, {0: 0, 1: 2, 2: 1}, actions=(('pool', (1, 2)),))
    before = catalogue.to_record()
    decoded = decode_events(catalogue)
    assert decoded.complete and decoded.observed_minimum == 0
    assert decoded.candidates[0].mapping == {0: 0, 1: 1, 2: 2}
    assert decoded.conversion.rematches[0].anchors == ((0, 0),)
    assert catalogue.to_record() == before


@pytest.mark.smoke
def test_direct_anchors_split_a_bond_at_the_exact_threshold(backend):
    r = endpoint('CO', [(0, 1, 1.)])
    p = endpoint('CO', [(0, 1, 1.5)])
    decoded = decode_events(saved(r, p, {0: 0, 1: 1}))
    assert decoded.complete and decoded.observed_minimum == 1
    assert decoded.candidates[0].endpoint_events['formed'] == ((0, 1),)
    owners = dict(decoded.conversion.branches[0].islands)
    assert owners[0] != owners[1]
    assert decoded.conversion.rematches[0].anchors == ((0, 0), (1, 1))


@pytest.mark.smoke
def test_native_conditional_metal_rule_is_used(backend):
    r = endpoint(('Co', 'C', 'O'), [(0, 1, 1.), (1, 2, 1.)])
    p = endpoint(('Co', 'C', 'O'), [(0, 1, 1.4), (1, 2, 1.4)])
    decoded = decode_events(saved(r, p, dict(enumerate(range(3)))))
    assert decoded.complete and decoded.observed_minimum == 1
    owners = dict(decoded.conversion.branches[0].islands)
    assert owners[0] != owners[1] == owners[2]


def test_free_hydrogens_stay_compressed_and_queryable(backend):
    ep = endpoint('CHHH', [(0, h, 1.) for h in (1, 2, 3)])
    decoded = decode_events(saved(ep, ep, {0: 0, 1: 2, 2: 3, 3: 1},
                                  actions=(('pool', (1, 2, 3)),)))
    assert decoded.complete and len(decoded.conversion.branches) == 1
    assert decoded.catalogue.families[0].actions
    candidate = decoded.candidates[0]
    assert decoded.query(candidate, {1: 3, 2: 1, 3: 2}, same_events=False).status == 'allowed'
    assert decoded.query(candidate, {1: 3, 2: 3}, same_events=False).status == 'forbidden'


def test_jointly_moving_atoms_are_not_fixed(backend):
    ep = endpoint('CCC', [(0, 1, 1.), (1, 2, 1.), (2, 0, 1.)])
    decoded = decode_events(saved(ep, ep, dict(enumerate(range(3))),
        actions=(('group', ((1, 2, 0),)),)))
    assert decoded.complete
    assert decoded.conversion.rematches[0].anchors == ()


def test_grouped_pairs_are_reused_but_different_anchors_are_separate(backend):
    ep = endpoint('OCC', [(0, 1, 1.), (0, 2, 2.)])
    catalogue = saved(ep, ep, {0: 0, 1: 2, 2: 1}, actions=(('pool', (1, 2)),))
    first = catalogue.families[0]
    catalogue.families += [replace(first, mapping=((0, 0), (1, 1), (2, 2))), replace(first, actions=())]
    next(iter(catalogue.branches.values())).update((1, 2))
    result = convert_fragments(catalogue)
    assert len(result.rematches) == 2 and result.reused_pairs == 0
    assert len(result.branches) == 2
    assert [b.source_family_ids for b in result.branches] == [(0, 1), (2,)]
    assert result.catalogue.families[0].provenance[0]['source_family_ids'] == (0, 1)


@pytest.mark.smoke
def test_partial_boundary_and_wholly_unmatched_material(backend):
    r = endpoint('COHH', [(0, 1, 1), (2, 3, 1)])
    p = endpoint('C', [])
    for left, right in ((r, p), (p, r)):
        decoded = decode_events(saved(left, right, {0: 0}))
        assert decoded.complete and decoded.observed_minimum == 1
        assert decoded.candidates[0].mapping == {0: 0}


@pytest.mark.smoke
def test_rematching_stays_inside_original_pairs_and_recombines(backend):
    ep = endpoint('HHHH', [])
    pairs = (((0, 1), (2, 3)), ((2, 3), (0, 1)))
    source = saved(ep, ep, {0: 2, 1: 3, 2: 0, 3: 1}, pairs=pairs,
        actions=(('pool', (0, 1)), ('pool', (2, 3))))
    result = convert_fragments(source)
    assert result.complete
    for family in result.catalogue.families:
        mapping = dict(family.mapping)
        assert {mapping[0], mapping[1]} == {2, 3}
        assert {mapping[2], mapping[3]} == {0, 1}
        assert len(set(mapping.values())) == 4


def test_filtered_minimum_is_over_converted_branches(backend):
    ep = endpoint('OCC', [(0, 1, 1.), (0, 2, 2.)])
    catalogue = saved(ep, ep, {0: 0, 1: 2, 2: 1}, actions=(('pool', (1, 2)),))
    catalogue.families.append(replace(catalogue.families[0], actions=()))
    next(iter(catalogue.branches.values())).add(1)
    decoded = decode_events(catalogue, EventDecodeConfig(event_slack=0))
    assert decoded.complete and decoded.minimum_proven and decoded.event_window == 0
    assert [c.total for c in decoded.candidates] == [0]
    assert sorted(r['total'] for r in decoded.family_reports) == [0, 2]


def test_empty_and_missing_fragment_catalogues_are_explicit():
    ep = endpoint('C', [])
    empty = decode_events(FinalBranchCatalogue(AAMProblem(ep, ep)))
    assert not empty.candidates and not empty.minimum_proven
    catalogue = saved(ep, ep, {0: 0})
    catalogue.branches.clear()
    with pytest.raises(ValueError, match='fragment atom sets'):
        decode_events(catalogue)


def test_no_decoder_solver_or_old_budget_option(backend, monkeypatch):
    import z3
    def forbidden(*args, **kwargs):
        pytest.fail('conversion invoked a solver')
    monkeypatch.setattr(z3, 'Solver', forbidden)
    ep = endpoint('CH', [(0, 1, 1.)])
    decoded = decode_events(saved(ep, ep, {0: 0, 1: 1}))
    assert decoded.complete
    assert all(row['solver_queries'] == 0 for row in decoded.family_reports)
    with pytest.raises(TypeError):
        EventDecodeConfig(seconds_per_family=10)


def test_artifact_preserves_converted_compression_and_scope(backend):
    ep = endpoint('CHH', [(0, 1, 1.), (0, 2, 1.)])
    decoded = decode_events(saved(ep, ep, dict(enumerate(range(3))), actions=(('pool', (1, 2)),)))
    record = decoded_events_record(decoded)
    assert record['schema'] == 'mappa.decoded_events/v2'
    assert record['scope'] == 'Direct counts across retained converted branches'
    restored = FinalBranchCatalogue.from_record(decoded.catalogue.problem, record['converted_catalogue'])
    assert restored.families[0].actions == decoded.catalogue.families[0].actions
    assert record['rematches'][0]['search']


@pytest.mark.smoke
def test_incomplete_rematch_does_not_prove_the_minimum(backend, monkeypatch):
    import mappa.fragment_conversion as module
    from mappa.search_graph import SearchStop
    original = module.find_islands
    def capped(*args, **kwargs):
        graph = original(*args, **kwargs)
        return replace(graph, stops=graph.stops + (SearchStop(graph.roots[0], 'capped'),))
    monkeypatch.setattr(module, 'find_islands', capped)
    ep = endpoint('C', [])
    decoded = decode_events(saved(ep, ep, {0: 0}))
    assert decoded.observed_minimum == 0 and not decoded.minimum_proven
    assert not decoded.complete


def test_grouping_precedes_recombination_and_result_return(backend, monkeypatch):
    import mappa.fragment_conversion as module
    ep = endpoint('CHH', [(0, 1, 1.), (0, 2, 1.)])
    catalogue = saved(ep, ep, {0: 0, 1: 1, 2: 2}, actions=(('pool', (1, 2)),))
    baseline = decode_events(catalogue)
    catalogue.families.extend(replace(catalogue.families[0]) for _ in range(99))
    next(iter(catalogue.branches.values())).update(range(100))
    original = module.product
    calls = []
    def counted(*args):
        calls.append(1)
        return original(*args)
    monkeypatch.setattr(module, 'product', counted)
    result = decode_events(catalogue)
    assert len(calls) == 1
    assert len(result.conversion.branches) == len(result.catalogue.families) == 1
    assert result.conversion.branches[0].source_family_ids == tuple(range(100))
    assert result.complete == baseline.complete and result.observed_minimum == baseline.observed_minimum
    assert [(c.id, c.mapping, c.total) for c in result.candidates] == [(c.id, c.mapping, c.total) for c in baseline.candidates]
    record = decoded_events_record(result)
    assert record['converted_branches'][0]['source_family_ids'] == tuple(range(100))
    assert len(record['family_reports']) == 1


def test_fragment_order_is_irrelevant_but_atom_pools_and_policy_are_not(backend):
    ep = endpoint('HHHH', [])
    pairs = (((0, 1), (0, 1)), ((2, 3), (2, 3)))
    catalogue = saved(ep, ep, dict(enumerate(range(4))), pairs=pairs,
                      actions=(('pool', (0, 1, 2, 3)),))
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
    catalogue = saved(ep, ep, {0: 0, 1: 1, 2: 2}, actions=(('pool', (1, 2)),))
    pairs = next(iter(catalogue.branches))[1]
    catalogue.families.append(replace(catalogue.families[0]))
    catalogue.coupled_branches[1] = {pairs}
    result = convert_fragments(catalogue)
    assert len(result.branches) == 1
    assert result.branches[0].source_family_ids == (0, 1)
