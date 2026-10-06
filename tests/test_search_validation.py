from dataclasses import replace
import numpy as np
import pytest
from mappa import AAMProblem, MolecularEndpoint
from mappa.search_graph import (
    AAMSearchGraph,
    SearchContext,
    SearchState,
    FragmentTransition,
    SearchStop,
)
from mappa.search_validation import RepresentativeValidationWorkspace
from mappa.family_scoring import validate_representative


def example(target_weight=1.0, cuts=(), deferred=()):
    def ep(weight):
        return MolecularEndpoint(
            ("C", "C"), np.zeros((2, 3)), np.array([[0.0, weight], [weight, 0.0]])
        )

    problem = AAMProblem(ep(1.0), ep(target_weight))
    context = SearchContext(
        (0, 1), (0, 1), (0, 1), cuts=cuts, graph_floor=0.2, iso_tolerance=0.5
    )
    mapping = ((0, 0), (1, 1))
    states = (
        SearchState(0, 0, (), (), ()),
        *(SearchState(i, 0, mapping, ((0, 1), (1, 1)), ()) for i in range(1, 4)),
    )

    def match(atoms):
        return dict(
            fragment=atoms,
            deferred_edges=deferred,
            symmetry=dict(witness=dict(mapping), blocks=[], exact_fixed=[0, 1]),
        )

    edges = (
        FragmentTransition(0, 0, 1, 0, (1, 0), match([0, 1]), ((0, 1),)),
        FragmentTransition(1, 0, 2, 0, (1, 0), match([0]), ()),
        FragmentTransition(2, 1, 3, 0, (1, 1), None, ()),
        FragmentTransition(3, 2, 3, 0, (1, 1), None, ()),
    )
    return problem, AAMSearchGraph(
        (context,), (0,), states, edges, (SearchStop(3, "objective_met"),)
    )


def oracle(graph, problem):
    valid = 0
    invalid = []
    for path in graph.paths():
        try:
            validate_representative(path, problem)
            valid += 1
        except AssertionError:
            invalid.append(path.terminal)
    return valid, invalid


@pytest.mark.parametrize("weight", [0.0, 0.2, 0.5, 1.0, 1.5, 1.500000002])
@pytest.mark.parametrize("excluded", ["none", "cut", "deferred"])
def test_shared_dag_checks_match_every_path(weight, excluded):
    p, g = example(
        weight,
        cuts=((0, 1),) if excluded == "cut" else (),
        deferred=((0, 1),) if excluded == "deferred" else (),
    )
    ws = RepresentativeValidationWorkspace(p, cache_entries=2)
    assert ws.validate_graph(g) == oracle(g, p)
    for path in g.paths():
        try:
            validate_representative(path, p)
        except AssertionError:
            with pytest.raises(AssertionError):
                ws.validate_path(path)
        else:
            ws.validate_path(path)


def test_join_certifies_both_histories_with_one_terminal_check():
    p, g = example()
    ws = RepresentativeValidationWorkspace(p)
    assert ws.validate_graph(g) == (2, [])
    assert ws.stats["terminal_checks"] == 1 and ws.stats["certified_paths"] == 2
    assert ws.stats["fallback_path_checks"] == 0


def test_one_invalid_history_does_not_mark_the_other_invalid():
    p, g = example(0.0)
    ws = RepresentativeValidationWorkspace(p)
    assert ws.validate_graph(g) == (1, [3])
    assert ws.stats["fallback_path_checks"] == 2


def test_wrong_elements_and_duplicate_images_fail():
    p, g = example()
    for mapping in [((0, 0), (1, 0)), ((0, 1), (1, 0))]:
        pp = replace(
            p, product=MolecularEndpoint(("C", "O"), np.zeros((2, 3)), p.product.wbo)
        )
        states = (*g.states[:-1], replace(g.states[-1], mapping=mapping))
        gg = replace(g, states=states)
        assert RepresentativeValidationWorkspace(pp).validate_graph(gg) == oracle(
            gg, pp
        )
