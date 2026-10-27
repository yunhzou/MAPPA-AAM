import itertools

import numpy as np
import pytest

from mappa import AAMProblem, MolecularEndpoint
from mappa.event_patterns import SignedEventIndex
from mappa.final_branches import FinalBranchCatalogue, FinalFamily
from mappa.postprocessing import EventDecodeConfig, decode_events


def catalogue():
    w = np.zeros((4, 4))
    for i, weight in enumerate((.1, 1., 1.8), 1):
        w[0, i] = w[i, 0] = weight
    endpoint = MolecularEndpoint(('C', 'O', 'O', 'O'), np.zeros((4, 3)), w)
    problem = AAMProblem(endpoint, endpoint)
    cat = FinalBranchCatalogue(problem)
    mapping = tuple(enumerate((0, 2, 3, 1)))
    cat.families = [
        FinalFamily(mapping, (), (), (.2, 1.)),
        FinalFamily(mapping, (), (('pool', (1, 2, 3)),), (.2, 1.)),
    ]
    return cat


@pytest.mark.parametrize('slack', [0, 1, 2, 5])
def test_relative_window_matches_exhaustive_permutation_oracle(slack):
    cat = catalogue()
    index = SignedEventIndex(cat.problem)
    patterns = [index.describe((0,) + p) for p in itertools.permutations((1, 2, 3))]
    minimum = min(p['total'] for p in patterns)
    expected = {p['id'] for p in patterns if p['total'] <= minimum + slack}
    decoded = decode_events(cat, EventDecodeConfig(event_slack=slack, seconds_per_family=20))
    assert decoded.complete
    assert decoded.minimum_proven
    assert decoded.observed_minimum == minimum
    assert decoded.event_window == minimum + slack
    assert {c.id for c in decoded.candidates} == expected
    # The nonminimal representative in the first family must disappear when a
    # hidden lower minimum is discovered in the second family.
    assert all(c.total <= decoded.event_window for c in decoded.candidates)


def test_window_is_sent_to_solver_and_tightens_between_families(monkeypatch):
    import mappa.postprocessing as pp
    cat = catalogue()
    cat.families.append(cat.families[1])
    original = pp.extract_path_events
    windows = []
    def capture(*args, **kwargs):
        windows.append(kwargs['max_events'])
        return original(*args, **kwargs)
    monkeypatch.setattr(pp, 'extract_path_events', capture)
    result = decode_events(cat, EventDecodeConfig(event_slack=0, seconds_per_family=20))
    assert result.complete
    assert windows[0] > windows[-1] == 0


def test_interrupted_minimum_is_not_claimed_proven():
    decoded = decode_events(catalogue(), EventDecodeConfig(event_slack=5, seconds_per_family=0))
    assert not decoded.complete
    assert decoded.observed_minimum > 0
    assert not decoded.minimum_proven


def test_invalid_or_ambiguous_budgets_rejected():
    for kwargs in ({'event_slack': -1}, {'event_slack': 1.5}, {'event_slack': 5, 'max_events': 9}):
        with pytest.raises(ValueError):
            EventDecodeConfig(**kwargs)


def test_default_unrestricted_behavior_and_empty_catalogue():
    cat = catalogue()
    full = decode_events(cat, EventDecodeConfig(seconds_per_family=20))
    wide = decode_events(cat, EventDecodeConfig(event_slack=5, seconds_per_family=20))
    assert {c.id for c in full.candidates} == {c.id for c in wide.candidates}
    empty = decode_events(FinalBranchCatalogue(cat.problem), EventDecodeConfig(event_slack=5))
    assert empty.observed_minimum is None and not empty.minimum_proven
