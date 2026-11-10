"""Saved-witness branch counting must not become a new symmetry search."""
from dataclasses import replace
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bench'))
from minimum_branches import composition, saved_branches
from mappa import AAMProblem, decode_events, decoded_events_record
from growth_cases import saved_fragment
from test_event_interpretation import endpoint, saved
from test_fragment_event_queries import arms


def test_internal_shuffles_and_repeated_branches_count_once_without_search(monkeypatch):
    from mappa.mappa_core.core.symmetry.minimum import _Search
    from mappa.mappa_core.core.symmetry.queries import EventSelection
    from mappa.mappa_core.core.symmetry.program import OrderedSymmetryEvents
    value = arms()
    w = value.source_symmetry.weights
    ep = endpoint(value.source_symmetry.elements, [(i, j, w[i, j])
                  for i, j in zip(*np.where(np.triu(w > 0, 1)))])
    catalogue = saved(ep, ep, dict(value.witness), fragments=(value,))
    catalogue.families *= 2
    next(iter(catalogue.branches.values())).add(1)
    record = decoded_events_record(decode_events(catalogue))
    record['window'].pop('branch_minima')
    record['window'].pop('minimum_choices')

    def forbidden(*args, **kwargs):
        pytest.fail('saved-result analysis started a search')
    monkeypatch.setattr(_Search, 'run', forbidden)
    monkeypatch.setattr(EventSelection, 'find', forbidden)
    monkeypatch.setattr(OrderedSymmetryEvents, 'solve', forbidden)
    counts, candidates, evaluated = saved_branches(record, target_count=0)
    assert list(counts.values()) == [{0}]  # 24 minimum mappings, one composition.
    assert len(candidates) == 1 and evaluated == 2
    assert saved_branches(record, reverse=True, target_count=0)[0] == counts


def test_recount_preserves_sweep_and_does_not_claim_an_unstored_tie():
    r = endpoint('CCCO', [(0, 1, 1.), (1, 2, 1.), (0, 3, 1.)])
    p = endpoint('CCCO', [(0, 1, 1.), (1, 2, 1.), (2, 3, 1.)])
    value = saved_fragment(AAMProblem(r, p), source_atoms=(0, 1, 2))
    catalogue = saved(r, p, dict(enumerate(range(4))), fragments=(value,),
                      pairs=(((0, 1, 2), (0, 1, 2)), ((3,), (3,))))
    catalogue.families[0] = replace(catalogue.families[0], provenance=[dict(cuts=((0, 3),))])
    decoded = decode_events(catalogue)
    record = decoded_events_record(decoded)
    counts, _, _ = saved_branches(record, target_count=decoded.minimum_events)
    assert decoded.minimum_events == 0
    assert list(counts.values()) == [{2}]  # The saved local witness does not turn the endpoint.
    assert decoded.contributions(0)['total'] == 2


def test_fragment_atom_sets_define_identity_and_reverse_normalization():
    pieces = ((0, 1), (2, 3))
    original = composition(pieces, {0: 4, 1: 5, 2: 6, 3: 7})
    assert composition(pieces[::-1], {0: 5, 1: 4, 2: 7, 3: 6}) == original
    assert composition(((4, 5), (6, 7)), {4: 0, 5: 1, 6: 2, 7: 3}, True) == original
    assert composition(pieces, {0: 4, 1: 6, 2: 5, 3: 7}) != original
