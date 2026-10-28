"""Conversion retains saved relations and never certifies a varying family."""
import pytest

from mappa.final_branches import FinalBranchCatalogue
from mappa.fragment_conversion import convert_aam_fragments
from test_event_patterns import problem, saved_path


@pytest.mark.parametrize('native', ['0', '1'])
def test_mixed_threshold_split_and_relation_retention(native, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', native)
    value = problem(['Co', 'C', 'O'], [(0, 1, 1.), (1, 2, 1.)],
                    [(0, 1, 1.4), (1, 2, 1.4)])
    saved = FinalBranchCatalogue(value).add_graph(saved_path(value).graph)
    before = saved.to_record()
    result = convert_aam_fragments(saved)
    assert result.complete and result.minimum == 1
    assert result.branches[0].fragment_pairs == (((0,), (0,)), ((1, 2), (1, 2)))
    assert result.to_record()['original_catalogue'] == before == saved.to_record()
    assert result.statistics['solver_queries'] == result.statistics['expanded_permutations'] == 0


def test_event_dependent_symmetry_is_retained_unresolved():
    value = problem('CCC', [(0, 1, 1.), (0, 2, 2.)], [(0, 1, 1.), (0, 2, 2.)])
    saved = FinalBranchCatalogue(value).add_graph(saved_path(value, [(0, 2, 1)]).graph)
    result = convert_aam_fragments(saved)
    assert not result.complete and result.minimum is None
    assert result.catalogue.families[0].actions == saved.families[0].actions
    assert 'event-dependent' in result.branches[0].reason


def test_compressed_equivalent_alternatives_survive():
    value = problem('COO', [(0, 1, 1.), (0, 2, 1.)], [(0, 1, 1.), (0, 2, 1.)])
    saved = FinalBranchCatalogue(value).add_graph(saved_path(value, [(0, 2, 1)]).graph)
    result = convert_aam_fragments(saved)
    assert result.complete and result.minimum == 0
    assert result.catalogue.families[0].actions


def test_target_only_internal_event_is_not_certified():
    value = problem('CCC', [(0, 1, 1.), (1, 2, 1.)],
                    [(0, 1, 1.), (1, 2, 1.), (0, 2, 1.)])
    saved = FinalBranchCatalogue(value).add_graph(saved_path(value).graph)
    result = convert_aam_fragments(saved)
    assert not result.complete and result.minimum is None
    assert 'target-only' in result.branches[0].reason


def test_missing_paths_preclude_complete_claim():
    value = problem('CO', [(0, 1, 1.)], [(0, 1, 1.)])
    saved = FinalBranchCatalogue(value).add_graph(saved_path(value).graph)
    saved.incomplete_paths = 1
    assert not convert_aam_fragments(saved).complete
