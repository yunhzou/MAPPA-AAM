"""Reaction-core extraction uses saved choices without starting new searches."""
import json

import pytest

from mappa import (AAMProblem, decode_events, decoded_events_from_record,
                   decoded_events_record, extract_reaction_cores)
from test_event_interpretation import endpoint, saved


def decoded_mappings(base, target, mappings, *, event_slack=5):
    catalogue = saved(base, target, mappings[0], fragments=())
    for mapping in mappings[1:]:
        addition = saved(base, target, mapping, fragments=())
        number = len(catalogue.families)
        catalogue.families.extend(addition.families)
        for key in addition.branches:
            catalogue.branches.setdefault(key, set()).add(number)
    decoded = decode_events(catalogue, event_slack=event_slack)
    _ = decoded.window
    return decoded


def test_different_bond_changes_with_equal_counts_stay_separate_and_replay(monkeypatch):
    ep = endpoint('CCCCON', [(0, 1, 1.), (1, 2, 1.), (2, 3, 1.), (0, 4, 1.), (3, 5, 1.)])
    mappings = [dict(enumerate(order)) for order in
                ((0, 1, 2, 3, 4, 5), (1, 0, 2, 3, 4, 5), (0, 1, 3, 2, 4, 5))]
    decoded = decoded_mappings(ep, ep, mappings)
    record = json.loads(json.dumps(decoded_events_record(decoded)))
    restored = decoded_events_from_record(record)

    def forbidden(*args, **kwargs):
        pytest.fail('reaction-core extraction started an optimization or query')
    from mappa.mappa_core.core.symmetry.minimum import _Search
    from mappa.mappa_core.core.symmetry.queries import EventSelection
    from mappa.mappa_core.core.symmetry.program import OrderedSymmetryEvents
    monkeypatch.setattr(_Search, 'run', forbidden)
    monkeypatch.setattr(EventSelection, 'find', forbidden)
    monkeypatch.setattr(OrderedSymmetryEvents, 'solve', forbidden)
    monkeypatch.setattr(decoded, 'minimum', forbidden)
    monkeypatch.setattr(decoded, 'tables', forbidden)
    cores = extract_reaction_cores(decoded)
    assert [core.physical_event_count for core in cores] == [0, 4, 4]
    assert len({core.id for core in cores}) == 3
    assert extract_reaction_cores(restored) == cores
    assert len(extract_reaction_cores(decoded, event_count=4)) == 2
    assert len(extract_reaction_cores(decoded, event_count=0)) == 1
    assert extract_reaction_cores(decoded, event_count=5) == ()
    for core in cores:
        for family, partition, cut in core.contexts:
            witness = decoded.saved_minimum(family, partition=partition, cut=cut)
            assert decoded.index.describe(decoded.index.vector(dict(witness.witness)))['id'] == core.id
    assert json.loads(json.dumps(decoded_events_record(restored))) == record


def test_window_excludes_saved_witnesses_above_requested_window():
    ep = endpoint('CCCCON', [(0, 1, 1.), (1, 2, 1.), (2, 3, 1.), (0, 4, 1.), (3, 5, 1.)])
    decoded = decoded_mappings(ep, ep, [dict(enumerate(range(6))),
                                       dict(enumerate((1, 0, 2, 3, 4, 5)))], event_slack=3)
    assert len(extract_reaction_cores(decoded)) == 1
    assert extract_reaction_cores(decoded, event_count=4) == ()


def test_equivalent_hydrogen_changes_share_one_core_and_retain_both_contexts():
    base = endpoint('COHH', [(0, 1, 1.), (0, 2, 1.), (0, 3, 1.)])
    target = endpoint('COHH', [(0, 1, 1.), (0, 2, 1.), (1, 3, 1.)])
    decoded = decoded_mappings(base, target, [dict(enumerate((0, 1, 2, 3))),
                                             dict(enumerate((0, 1, 3, 2)))])
    core, = extract_reaction_cores(decoded)
    assert core.physical_event_count == 2 and core.model_event_counts == (2,)
    assert core.contexts == ((0, 0, 0), (1, 0, 0))
    assert core.decreased_bonds == ((0, 3),) and core.increased_bonds == ((1, 3),)
    assert core.base_atoms == core.target_atoms == (0, 1, 3)


def test_fragment_partitions_do_not_define_core_identity_or_replace_physical_count():
    ep = endpoint('CO', [(0, 1, 1.)])
    catalogue = saved(ep, ep, {0: 0, 1: 1}, fragments=())
    catalogue.coupled_branches[0] = {(((0,), (0,)), ((1,), (1,)))}
    decoded = decode_events(catalogue)
    core, = extract_reaction_cores(decoded)
    assert core.physical_event_count == 0
    assert core.model_event_counts == (0, 2) and len(core.contexts) == 2
    assert core.base_atoms == core.target_atoms == ()
    selected, = extract_reaction_cores(decoded, event_count=2)
    assert selected.id == core.id and selected.model_event_counts == (2,)


def test_partial_maps_keep_unmatched_boundary_atoms_in_the_correct_endpoint():
    base, target = endpoint('CO', [(0, 1, 1.)]), endpoint('CN', [(0, 1, 1.)])
    core, = extract_reaction_cores(decoded_mappings(base, target, [{0: 0}]))
    assert core.mapping == ((0, 0),)
    assert core.decreased_bonds == core.increased_bonds == ((0, 1),)
    assert core.base_atoms == core.target_atoms == (0, 1)
    assert core.physical_event_count == 2


def test_original_weights_and_metal_thresholds_define_the_reaction_core():
    binary = endpoint(('Fe', 'C', 'O'), [(0, 1, 1.), (1, 2, 1.)])
    base = endpoint(('Fe', 'C', 'O'), [(0, 1, 1.), (1, 2, 3.)])
    target = endpoint(('Fe', 'C', 'O'), [(0, 1, 1.4), (1, 2, 1.)])
    decoded = decode_events(saved(binary, binary, dict(enumerate(range(3))), fragments=()),
                            event_problem=AAMProblem(base, target))
    core, = extract_reaction_cores(decoded)
    assert core.decreased_bonds == ((1, 2),)  # triple to single: weakening
    assert core.increased_bonds == ((0, 1),)  # 0.4 exceeds the metal threshold
    assert core.physical_event_count == 2


def test_missing_saved_counts_are_rejected_without_calculating_replacements():
    ep = endpoint('CO', [(0, 1, 1.)])
    decoded = decoded_mappings(ep, ep, [{0: 0, 1: 1}])
    record = decoded_events_record(decoded)
    del record['window']['branch_minima']
    restored = decoded_events_from_record(record)
    with pytest.raises(ValueError, match='lacks saved branch counts'):
        extract_reaction_cores(restored)
    for count in (-1, 1.5, True):
        with pytest.raises(ValueError, match='nonnegative integer'):
            extract_reaction_cores(decoded, event_count=count)
