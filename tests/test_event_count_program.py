"""Independent finite oracles for correlated fragment event contraction."""
from itertools import product

import pytest

from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints
from mappa.mappa_core.core.symmetry.event_counts import EventCountProgram


def test_correlated_cycle_does_not_sum_independent_zero_bounds():
    system = FiniteConstraints()
    choice = system.symbol(range(3))
    images = system.apply_factor([system.constant(i) for i in range(3)],
        ((1, 0, 2), (2, 1, 0), (0, 2, 1)), choice[0])
    costs = [(image[0], i, {(p, i): int(p != i) for p in image[1]})
             for i, image in enumerate(images)]
    response = EventCountProgram(system, costs)
    assert response.project() == {(): 2}
    assert response.project((images[0][0],)) == {(0,): 2, (1,): 2, (2,): 2}


def test_long_chain_contracts_shared_choices_without_full_assignments():
    system = FiniteConstraints()
    values = [system.symbol((0, 1))[0] for _ in range(60)]
    costs = []
    for left, right in zip(values, values[1:]):
        costs.append((left, right, {(a, b): int(a == b) for a, b in product((0, 1), repeat=2)}))
    response = EventCountProgram(system, costs)
    assert response.project() == {(): 0}
    assert response.project((values[0], values[-1])) == {
        (0, 0): 1, (0, 1): 0, (1, 0): 0, (1, 1): 1}
    assert response.project((values[0], values[-1])) is response.project((values[0], values[-1]))


def test_conditional_lookup_keeps_only_the_selected_entry_correlated():
    system = FiniteConstraints()
    selector = system.symbol((0, 1))
    left, right = system.symbol((0, 1)), system.symbol((0, 1))
    system.distinct((left[0], right[0]))
    output = system.lookup(selector, {0: left, 1: right})
    costs = [(selector[0], output[0], {(s, x): -1 if s == x else 2
                                    for s, x in product((0, 1), repeat=2)})]
    response = EventCountProgram(system, costs)
    expected = {}
    for s, a, b in product((0, 1), repeat=3):
        if a == b:
            continue
        x = (a, b)[s]
        key = s, a, b, x
        expected[key] = -1 if s == x else 2
    assert response.project((selector[0], left[0], right[0], output[0])) == expected
    assert response.project() == {(): -1}


def test_repeated_and_fixed_outputs_keep_their_public_order():
    system = FiniteConstraints()
    variable = system.symbol((2, 4))[0]
    response = EventCountProgram(system, [(variable, 7, {(2, 7): 3, (4, 7): 1})])
    assert response.project((7, variable, variable, 7)) == {(7, 2, 2, 7): 3, (7, 4, 4, 7): 1}
    assert response.project() == {(): 1}


def test_infeasible_joint_tables_have_no_attainable_count():
    system = FiniteConstraints()
    a, b = [system.symbol((0, 1))[0] for _ in range(2)]
    system.relation((a, b), ((0, 0), (1, 1)))
    system.distinct((a, b))
    response = EventCountProgram(system)
    assert response.project() == {}
    assert response.project((a,)) == {}


def test_duplicate_rows_keep_minimum_and_incomplete_costs_fail():
    system = FiniteConstraints()
    a, b = [system.symbol((0, 1))[0] for _ in range(2)]
    system.relation((a, b), ((0, 0), (0, 0), (1, 1)))
    with pytest.raises(ValueError, match='cover both endpoint domains'):
        EventCountProgram(system, [(a, b, {(0, 0): 0})])
    response = EventCountProgram(system, [(a, b, {(x, y): x + y for x, y in product((0, 1), repeat=2)})])
    assert response.project((a, b)) == {(0, 0): 0, (1, 1): 2}
