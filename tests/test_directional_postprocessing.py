"""Batch stage calls preserve the existing directed scientific results."""
from dataclasses import replace
import json

import numpy as np
import pytest

from mappa import (
    AAMProblem, AAMSearchConfig, DirectionalResults, decode_events, deduplicate,
    decoded_events_from_record, decoded_events_record, extract_reaction_cores,
    search_aam_directions,
)
from mappa.mappa_core.post_processing.decoder.api import decode_events as single_decode
from mappa.mappa_core.post_processing.final_branches import deduplicate_final_branches
from mappa.mappa_core.utils.io.records import decoded_events_record as single_record
from test_event_interpretation import endpoint


@pytest.fixture
def matched():
    base = endpoint('COHH', [(0, 1, 1.), (0, 2, 1.), (0, 3, 1.)])
    target = endpoint('COHH', [(0, 1, 1.), (0, 2, 1.), (1, 3, 1.)])
    return search_aam_directions(AAMProblem(base, target),
        AAMSearchConfig(iso_tolerance=2., branch_limit=100), direction='both')


@pytest.mark.smoke
def test_four_stage_calls_preserve_single_direction_results(matched):
    catalogues = deduplicate(matched)
    decoded = decode_events(catalogues, event_slack=5)
    cores = extract_reaction_cores(decoded)
    minimum = extract_reaction_cores(decoded, event_count='minimum')
    assert tuple(catalogues) == tuple(decoded) == tuple(cores) == ('forward', 'reverse')
    assert decoded.forward is decoded['forward'] and decoded.reverse is decoded['reverse']
    assert cores.forward and cores.reverse and minimum.forward and minimum.reverse
    for run, direction in zip(matched, ('forward', 'reverse')):
        expected_catalogue = deduplicate_final_branches(run.aam)
        expected = single_decode(expected_catalogue, event_slack=5)
        assert catalogues[direction].to_record() == expected_catalogue.to_record()
        assert single_record(decoded[direction]) == single_record(expected)
        assert cores[direction] == extract_reaction_cores(expected)
        assert minimum[direction] == extract_reaction_cores(expected, event_count=expected.minimum_events)
        core = cores[direction][0]
        assert cores.to_input_mapping(core.mapping, direction=direction) == run.to_input_mapping(core.mapping)
    # The documented direct route can also consume raw directional searches.
    direct = decode_events(matched, event_slack=5)
    assert decoded_events_record(direct) == decoded_events_record(decoded)


def test_directional_save_restore_preserves_windows_and_never_searches(matched, monkeypatch):
    decoded = decode_events(deduplicate(matched), event_slack=3)
    expected = extract_reaction_cores(decoded)
    record = json.loads(json.dumps(decoded_events_record(decoded)))
    from mappa.mappa_core.core.symmetry.minimum import _Search
    from mappa.mappa_core.core.symmetry.program import OrderedSymmetryEvents
    from mappa.mappa_core.post_processing.decoder import fragments

    def forbidden(*args, **kwargs):
        pytest.fail('restoring or inspecting a directional record repeated search')
    monkeypatch.setattr(_Search, 'run', forbidden)
    monkeypatch.setattr(OrderedSymmetryEvents, 'solve', forbidden)
    monkeypatch.setattr(fragments, 'optimize_fragment', forbidden)
    restored = decoded_events_from_record(record)
    assert extract_reaction_cores(restored) == expected
    assert json.loads(json.dumps(decoded_events_record(restored))) == record
    assert restored.forward.event_window == decoded.forward.event_window
    assert restored.reverse.event_window == decoded.reverse.event_window


@pytest.mark.parametrize('direction', ['forward', 'reverse'])
def test_one_requested_direction_stays_single_and_preserves_partial_anchors(direction):
    base, target = endpoint('COH', []), endpoint('OC', [])
    matched = search_aam_directions(AAMProblem(base, target),
        AAMSearchConfig(anchors=((0, 1),), sweep_cuts=False), direction=direction)
    decoded = decode_events(deduplicate(matched))
    cores = extract_reaction_cores(decoded, event_count='minimum')
    assert len(cores) == 1 and tuple(cores) == (direction,)
    assert cores[direction]
    assert all(cores.to_input_mapping(core.mapping, direction=direction)[0] == 1
               for core in cores[direction])
    missing = 'reverse' if direction == 'forward' else 'forward'
    assert getattr(cores, missing) is None
    with pytest.raises(KeyError):
        cores.to_input_mapping({}, direction=missing)
    assert tuple(decoded_events_from_record(decoded_events_record(decoded))) == (direction,)


def test_original_event_weights_are_oriented_for_reverse_decoding():
    base = endpoint('CO', [(0, 1, 1.)])
    target = endpoint('OC', [(0, 1, 1.)])
    matched = search_aam_directions(AAMProblem(base, target), AAMSearchConfig(sweep_cuts=False))
    raw = AAMProblem(endpoint('CO', [(0, 1, 3.)]), endpoint('OC', [(0, 1, 1.)]))
    decoded = decode_events(deduplicate(matched), event_problem=raw, event_slack=5)
    cores = extract_reaction_cores(decoded, event_count='minimum')
    assert decoded.forward.event_problem is raw
    np.testing.assert_array_equal(decoded.reverse.event_problem.reactant.wbo, raw.product.wbo)
    assert decoded.reverse.event_problem.reactant.elements == ('O', 'C')
    assert cores.forward[0].decreased_bonds == ((0, 1),)
    assert cores.forward[0].increased_bonds == ()
    assert cores.reverse[0].increased_bonds == ((0, 1),)
    assert cores.reverse[0].decreased_bonds == ()
    assert cores.forward[0].physical_event_count == cores.reverse[0].physical_event_count == 1


def test_directional_boundary_rejects_duplicate_or_mixed_searches(matched):
    with pytest.raises(ValueError, match='duplicate forward'):
        deduplicate((matched[0], matched[0]))
    other = AAMProblem(endpoint('N', []), endpoint('N', []))
    mixed = replace(matched[1], plan=replace(matched[1].plan, input_problem=other))
    with pytest.raises(ValueError, match='same endpoint problem'):
        decode_events((matched[0], mixed))
    with pytest.raises(TypeError, match='search_aam_directions'):
        deduplicate(())
    with pytest.raises(ValueError, match='at least one'):
        DirectionalResults()


def test_minimum_selection_is_per_direction_and_does_not_merge_event_ids():
    from test_reaction_cores import decoded_mappings
    ep = endpoint('CCCCON', [(0, 1, 1.), (1, 2, 1.), (2, 3, 1.), (0, 4, 1.), (3, 5, 1.)])
    decoded = DirectionalResults(
        forward=decoded_mappings(ep, ep, [dict(enumerate(range(6)))]),
        reverse=decoded_mappings(ep, ep, [dict(enumerate((1, 0, 2, 3, 4, 5)))]))
    cores = extract_reaction_cores(decoded, event_count='minimum')
    assert cores.forward[0].model_event_counts == (0,)
    assert cores.reverse[0].model_event_counts == (4,)
    exact = extract_reaction_cores(decoded, event_count=0)
    assert exact.forward and exact.reverse == ()
