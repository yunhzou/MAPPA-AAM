"""Small feasibility checks for reusing the matcher in event refinement.

These tests exercise existing APIs, not a completed conversion algorithm.
Only the tiny independent test oracle expands the two allowed assignments.
"""
import json
import time

import pytest

from mappa import FragmentMatchConfig, FragmentMatchContext, match_fragment
from mappa.mappa_core.core.search.branch import find_islands
from mappa.frag import bond_event_threshold, classify_bonds
from mappa.mappa_core.core.matching.representation.state import candidate_from_record
from test_growth_compression_exactness import graph, represented


def anchored_match(source, target, tolerance):
    return match_fragment(
        source, target, seed=0,
        context=FragmentMatchContext(locked_mapping={0: 0}, islands={0: 0}),
        config=FragmentMatchConfig(iso_tolerance=tolerance,
                                   allow_mapped_seed=True, branch_limit=100))


def event_count(source, target, mapping, threshold=.5, metal_threshold=.3):
    broken, formed, _, _ = classify_bonds(
        mapping, source.graph['wbo_matrix'], target.graph['wbo_matrix'],
        dwbo_threshold=threshold,
        elements_R=[source.nodes[a]['element'] for a in source],
        elements_P=[target.nodes[a]['element'] for a in target],
        metal_dwbo_threshold=metal_threshold)
    return len(broken) + len(formed)


@pytest.mark.parametrize('native', ['0', '1'])
def test_anchored_three_atom_rematch(native, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', native)
    source = graph(['C'] * 3, [(0, 1, 1.), (0, 2, 2.)])
    target = graph(['C'] * 3, [(0, 1, 1.), (0, 2, 2.)])
    identity = ((0, 0), (1, 1), (2, 2))
    swapped = ((0, 0), (1, 2), (2, 1))
    outcomes = []
    for tolerance, expected in [(1., {identity, swapped}), (.5, {identity})]:
        cpu = time.process_time()
        result = anchored_match(source, target, tolerance)
        elapsed_cpu = time.process_time() - cpu
        assert not result.capped
        assert all(match.fragment == frozenset(range(3)) for match in result.matches)
        actual = represented([candidate_from_record(m.symmetry)
                              for m in result.matches], target, tolerance)
        assert actual == expected
        assert all(dict(mapping)[0] == 0 for mapping in actual)
        outcomes.append(dict(tolerance=tolerance, branches=len(result.matches),
                             mappings=sorted(actual),
                             events=sorted(event_count(source, target, dict(m)) for m in actual),
                             wall_seconds=result.elapsed_seconds, cpu_seconds=elapsed_cpu))
    # The full island scheduler also takes the intact tight match. It does not
    # automatically preserve the looser swapped alternative as split islands.
    tight = find_islands(source, target, [0, 1, 2], iso_tol=.5,
                         anchor_map={0: 0}, max_branches=100)
    tight_paths = list(tight.paths())
    assert tight_paths
    assert {tuple(sorted(path.mapping.items())) for path in tight_paths} == {identity}
    assert event_count(source, target, dict(swapped)) == 2
    print(json.dumps(dict(native=native, matches=outcomes,
                          tight_scheduler_paths=len(tight_paths),
                          conversion_still_needs_split_alternative=swapped)), flush=True)


@pytest.mark.parametrize('native', ['0', '1'])
def test_matching_threshold_equality_is_not_event_invariance(native, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', native)
    source = graph(['C', 'C'], [(0, 1, 1.)])
    target = graph(['C', 'C'], [(0, 1, 1.5)])
    result = anchored_match(source, target, .5)
    assert not result.capped and len(result.matches) == 1
    assert result.matches[0].fragment == frozenset((0, 1))
    assert event_count(source, target, result.matches[0].mapping) == 1
    print(json.dumps(dict(native=native, delta=.5,
                          accepted_internal_bond=True, events=1)), flush=True)


@pytest.mark.parametrize('elements', [('C', 'O'), ('Co', 'O')])
@pytest.mark.parametrize('delta', [.25, .4, .6])
def test_ordinary_and_metal_anchored_pair_thresholds(elements, delta, monkeypatch):
    # Select the existing pair-specific event policy, then reuse the public
    # matcher. This two-atom check does not implement mixed-fragment policy.
    monkeypatch.setenv('MAPPA_NATIVE', '0')
    source = graph(elements, [(0, 1, 1.)])
    target = graph(elements, [(0, 1, 1. + delta)])
    threshold = bond_event_threshold(elements, 0, 1,
                                     default_threshold=.5, metal_threshold=.3)
    result = anchored_match(source, target, threshold)
    assert not result.capped
    intact = any(match.fragment == frozenset((0, 1)) for match in result.matches)
    events = event_count(source, target, {0: 0, 1: 1})
    assert intact == (delta < threshold)
    assert events == int(delta >= threshold)
    print(json.dumps(dict(elements=elements, delta=delta, threshold=threshold,
                          intact_match=intact, events=events)), flush=True)


def test_mixed_fragment_needs_pair_specific_matching_policy(monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', '0')
    source = graph(['Co', 'C', 'O'], [(0, 1, 1.), (1, 2, 1.)])
    target = graph(['Co', 'C', 'O'], [(0, 1, 1.4), (1, 2, 1.4)])
    result = anchored_match(source, target, .5)
    assert not result.capped and len(result.matches) == 1
    assert result.matches[0].fragment == frozenset((0, 1, 2))
    broken, formed, _, _ = classify_bonds(
        {0: 0, 1: 1, 2: 2}, source.graph['wbo_matrix'], target.graph['wbo_matrix'],
        dwbo_threshold=.5, metal_dwbo_threshold=.3,
        elements_R=['Co', 'C', 'O'], elements_P=['Co', 'C', 'O'])
    assert not broken
    assert [(a, b) for a, b, _, _ in formed] == [(0, 1)]
    print(json.dumps(dict(mixed_fragment=['Co', 'C', 'O'], delta=.4,
                          uniform_iso_tol=.5, intact_match=True,
                          event_pairs=[[0, 1]],
                          needs_pair_specific_matching_policy=True)), flush=True)
