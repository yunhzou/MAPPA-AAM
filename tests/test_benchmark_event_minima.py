"""Exact experiment minima retain matching constraints and count only heavy bonds."""
from dataclasses import replace
from itertools import permutations, product
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'bench'))
import analytics
from analytics import heavy_event_index, minimum_events, saved_response_bounds
from mappa import AAMProblem, decode_events
from mappa.mappa_core.utils.io.records import decoded_events_record
from test_fragment_conversion import endpoint, saved


def analyze(catalogue, reference, **kwargs):
    return minimum_events(decoded_events_record(decode_events(catalogue)), reference, **kwargs)


def rearrangement(*, hydrogen_constraint=False):
    reactant = endpoint('CCCH', [(0, 1, 1.), (1, 3, 1.)])
    product = endpoint('CCCH', [(0, 2, 1.), (1, 3, 1.)])
    catalogue = saved(reactant, product, dict(enumerate(range(4))),
                      actions=(('pool', (1, 2)),))
    if hydrogen_constraint:
        catalogue.families[0] = replace(catalogue.families[0], required_edges=((1, 3),),
                                         policy=(.2, .1))
    return catalogue


def test_minimum_is_over_correlated_choices_not_only_the_saved_witness():
    catalogue = rearrangement()
    before = catalogue.to_record()
    result = analyze(catalogue, {0: 0, 1: 1, 2: 2}, seconds=5)
    assert result['status'] == 'exact'
    assert result['reference_events'] == 2
    assert result['minimum_events'] == 0 and result['difference'] == 2
    assert result['solver_queries'] > 0
    assert dict(result['witness']) == {0: 0, 1: 2, 2: 1, 3: 3}
    assert catalogue.to_record() == before


def test_heavy_event_projection_preserves_hydrogen_matching_constraints():
    catalogue = rearrangement(hydrogen_constraint=True)
    before = catalogue.to_record()
    result = analyze(catalogue, {0: 0, 1: 1, 2: 2}, seconds=5)
    # Swapping C1/C2 would remove both heavy events, but violates the saved C-H bond.
    assert result['status'] == 'exact'
    assert result['minimum_events'] == result['reference_events'] == 2
    assert result['difference'] == 0 and result['same_heavy_coverage']
    assert catalogue.to_record() == before


@pytest.mark.parametrize('reference,expected', [({0: 0, 1: 1, 2: 2}, 0),
                                               ({0: 0, 1: 2, 2: 1}, -2)])
def test_negative_and_zero_differences_are_retained(reference, expected):
    catalogue = rearrangement(hydrogen_constraint=True)
    result = analyze(catalogue, reference, seconds=5)
    assert result['status'] == 'exact' and result['difference'] == expected


def test_event_scope_counts_partial_heavy_boundaries_and_strict_metal_thresholds():
    reactant = endpoint(('Fe', 'C', 'H'), [(0, 1, 1.), (1, 2, 1.)])
    product = endpoint(('Fe', 'C', 'H', 'C'),
                       [(0, 1, 1.3), (1, 2, 2.), (1, 3, 1.)])
    catalogue = saved(reactant, product, {0: 0, 1: 1, 2: 2})
    index = heavy_event_index(catalogue.problem)
    # The H change is excluded; the new bond to unmatched C3 is included. The
    # Fe-C comparison follows the original subtraction, including its float edge.
    expected = int(abs(1.3-1.) > .3) + 1
    result = analyze(catalogue, {0: 0, 1: 1}, seconds=5)
    assert result['status'] == 'exact'
    assert result['reference_events'] == result['minimum_events'] == expected
    assert not np.any(index.r[2]) and not np.any(index.p[2])
    assert catalogue.problem.reactant.wbo[1, 2] == 1.
    assert catalogue.problem.product.wbo[1, 2] == 2.


def test_unfinished_proof_does_not_become_an_exact_minimum_or_zero_bin():
    result = analyze(rearrangement(), {0: 0, 1: 1, 2: 2}, seconds=0)
    assert result['status'] == 'unknown'
    assert result['reference_events'] == 2
    assert result['minimum_events'] is None and result['difference'] is None
    assert result['unresolved_family_count'] == 1


def test_symmetry_invariance_certifies_positive_minimum_without_solver(monkeypatch):
    reactant = endpoint('C', [])
    target = endpoint('CCC', [(0, 1, 1.), (0, 2, 1.)])
    catalogue = saved(reactant, target, {0: 0}, actions=(('pool', (1, 2)),))
    def unexpected(*args, **kwargs):
        pytest.fail('invariant event response needs no constraint model')
    monkeypatch.setattr(analytics, 'compile_path', unexpected)
    result = analyze(catalogue, {0: 0}, seconds=5)
    assert result['status'] == 'exact' and result['minimum_events'] == 2
    assert result['proof_counts'] == {'event_invariance': 1}
    assert result['solver_queries'] == 0


def test_good_witness_does_not_hide_an_unresolved_competing_family(monkeypatch):
    catalogue = rearrangement()
    variable = catalogue.families[0]
    catalogue.families[:] = [replace(variable, actions=()), variable]
    next(iter(catalogue.branches.values())).add(1)
    original = analytics._event_model
    def exhausted(*args, **kwargs):
        raise TimeoutError('controlled encoding budget')
    monkeypatch.setattr(analytics, '_event_model', exhausted)
    result = analyze(catalogue, {0: 0, 1: 1, 2: 2}, seconds=5)
    assert result['status'] == 'unknown'
    assert (result['minimum_lower_bound'], result['minimum_upper_bound']) == (0, 2)
    assert result['minimum_events'] is None and result['difference'] is None
    assert result['unresolved_family_count'] == 1
    monkeypatch.setattr(analytics, '_event_model', original)
    solved = analyze(catalogue, {0: 0, 1: 1, 2: 2}, seconds=5)
    assert solved['status'] == 'exact' and solved['minimum_events'] == 0
    assert solved['witness_family'] == 1


def test_saved_bounds_do_not_add_alternative_partitions_or_hydrogen_events():
    reactant = endpoint('CCCH', [(0, 1, 2.), (0, 2, 2.), (1, 2, 2.), (0, 3, 1.)])
    target = endpoint('CCCH', [(0, 1, 1.), (0, 2, 1.), (1, 2, 1.), (0, 3, 2.)])
    catalogue = saved(reactant, target, dict(enumerate(range(4))))
    catalogue.branches[((.2, 1.), (((0, 1), (0, 1)), ((2, 3), (2, 3))))] = {0}
    record = decoded_events_record(decode_events(catalogue))
    assert len(record['families'][0]['tables']) == 2
    bounds, checked = saved_response_bounds(record, catalogue.problem, float('inf'))
    assert checked == 2 and bounds == [3]
    result = minimum_events(record, {0: 0, 1: 1, 2: 2}, seconds=5)
    assert result['status'] == 'exact' and result['minimum_events'] == 3
    assert result['saved_response_certified_families'] == 1


@pytest.mark.parametrize('partial', [False, True])
def test_saved_bounds_are_admissible_for_ordered_correlated_choices(partial):
    reactant = endpoint('CCCCH', [(0, 1, 1.), (1, 2, 2.), (0, 4, 1.)])
    target = endpoint('CCCCH', [(0, 2, 1.), (1, 3, 2.), (0, 4, 2.)])
    mapping = {0: 0, 1: 1, 2: 2} if partial else dict(enumerate(range(5)))
    catalogue = saved(reactant, target, mapping,
                      actions=(('pool', (1, 2)), ('pool', (2, 3))))
    interpretation = decode_events(catalogue)
    record = decoded_events_record(interpretation)
    index = heavy_event_index(catalogue.problem)
    counts = []
    for choices in product(permutations((1, 2)), permutations((2, 3))):
        mapping = interpretation.realize(0, choices)
        # Independent dense oracle, including exactly-one-mapped target boundaries.
        r, p = catalogue.problem.reactant, catalogue.problem.product
        total = 0
        for a in range(4):
            for b in range(a + 1, 4):
                if a not in mapping and b not in mapping:
                    continue
                target_weight = p.wbo[mapping[a], mapping[b]] if a in mapping and b in mapping else 0.
                total += abs(r.wbo[a, b] - target_weight) > .5
        occupied = set(mapping.values())
        total += sum(p.wbo[a, b] > .5 for a in range(4) for b in range(a + 1, 4)
                     if (a in occupied) != (b in occupied))
        assert int(index.counts([index.vector(mapping)])[0].sum()) == total
        counts.append(int(total))
    bounds, _ = saved_response_bounds(record, catalogue.problem, float('inf'))
    assert bounds[0] <= min(counts)
    result = minimum_events(record, {a: a for a in range(4)}, seconds=5)
    assert result['status'] == 'exact' and result['minimum_events'] == min(counts)


def test_minimum_rejects_rematched_outputs():
    catalogue = rearrangement()
    record = decoded_events_record(decode_events(catalogue))
    record['schema'] = 'mappa.decoded_events/v2'
    with pytest.raises(ValueError, match='original-space'):
        minimum_events(record, {0: 0, 1: 1, 2: 2})


def test_saved_binary_constraints_keep_original_weighted_minimum():
    binary = endpoint('CO', [(0, 1, 1.)])
    catalogue = saved(binary, binary, {0: 0, 1: 1})
    catalogue.families[0] = replace(catalogue.families[0], required_edges=((0, 1),))
    original = AAMProblem(endpoint('CO', [(0, 1, 3.)]), binary)
    record = decoded_events_record(decode_events(catalogue, event_problem=original))
    result = minimum_events(record, {0: 0, 1: 1})
    assert result['status'] == 'exact'
    assert result['minimum_events'] == result['reference_events'] == 1
    assert result['event_problem_sha256'] != result['problem_sha256']
    record['event_problem']['reactant']['wbo'][0][1] = 1.
    record['event_problem']['reactant']['wbo'][1][0] = 1.
    with pytest.raises(ValueError, match='identities'):
        minimum_events(record, {0: 0, 1: 1})
