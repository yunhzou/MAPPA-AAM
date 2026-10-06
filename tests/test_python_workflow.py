"""Public notebook workflow: event classes, certified shuffles, orientations."""

import numpy as np
import pytest
from mappa import (
    AAMProblem,
    MolecularEndpoint,
    AAMSearchConfig,
    search_aam,
    search_aam_directions,
)
from mappa.postprocessing import EventDecodeConfig, decode_events


@pytest.fixture
def reaction():
    def endpoint(bonds):
        matrix = np.zeros((6, 6))
        for a, b, w in bonds:
            matrix[a, b] = matrix[b, a] = w
        return MolecularEndpoint(
            ("C", "O", "H", "H", "H", "H"), np.zeros((6, 3)), matrix
        )

    return AAMProblem(
        endpoint([(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 5, 1)]),
        endpoint([(0, 1, 2), (0, 4, 1), (0, 5, 1), (2, 3, 1)]),
    )


@pytest.fixture
def searched(reaction):
    return search_aam(reaction, AAMSearchConfig(branch_limit=2000), workers=1)


def test_event_classes_and_conditional_shuffles(searched):
    decoded = decode_events(searched)
    assert decoded.complete
    assert [c.total for c in decoded.candidates] == [4, 6]
    assert len({c.id for c in decoded.candidates}) == 2
    candidate, other = decoded.candidates
    mapping = dict(candidate.mapping)
    mapping[2], mapping[3] = mapping[3], mapping[2]
    assert decoded.query(candidate, mapping).mapping == mapping
    assert decoded.query(candidate, {2: mapping[2]}).status == "allowed"
    assert decoded.query(candidate, other.mapping).status == "forbidden"
    assert (
        decoded.query(candidate, other.mapping, same_events=False).status == "allowed"
    )
    assert decoded.query(candidate, {2: 2, 3: 2}).status == "forbidden"
    assert decoded.query(candidate, {}, seconds=0).status == "unknown"
    assert decoded.symmetry(candidate)[0]["provenance"]
    assert [c.total for c in decoded.minimum_candidates] == [4]


def test_thresholds_are_separate_and_partial_results_are_marked(searched):
    before = searched.graph.to_record()
    assert [
        c.total
        for c in decode_events(searched, EventDecodeConfig(threshold=1.1)).candidates
    ] == [0]
    incomplete = decode_events(searched, EventDecodeConfig(seconds_per_family=0))
    assert not incomplete.complete
    assert incomplete.candidates
    assert before == searched.graph.to_record()
    with pytest.raises(ValueError):
        incomplete.query(decode_events(searched).candidates[0], {})


def test_bidirectional_anchors_and_unequal_sizes():
    r = MolecularEndpoint(("C", "O", "H"), np.zeros((3, 3)), np.zeros((3, 3)))
    p = MolecularEndpoint(("O", "C"), np.zeros((2, 3)), np.zeros((2, 2)))
    config = AAMSearchConfig(anchors=((0, 1),), sweep_cuts=False)
    runs = search_aam_directions(AAMProblem(r, p), config)
    assert [v.plan.direction for v in runs] == ["R_to_P", "P_to_R"]
    assert runs[1].aam.config.anchors == ((1, 0),)
    for run in runs:
        assert all(
            run.to_input_mapping(path.mapping)[0] == 1 for path in run.aam.graph.paths()
        )
    assert search_aam_directions(AAMProblem(r, p), config, direction="smaller_first")[
        0
    ].plan.reversed
    with pytest.raises(ValueError):
        decode_events(runs[0].aam)


def test_seed_and_cut_controls_preserve_checkpoint_identity(reaction):
    from mappa.aam import checkpoint_manifest
    from dataclasses import replace

    cfg = AAMSearchConfig(sweep_cuts=False, random_seed=9)
    aam = search_aam(reaction, cfg)
    assert aam.metrics.cut_count == 1
    assert len(aam.graph.contexts) == 1
    assert checkpoint_manifest(reaction, cfg) != checkpoint_manifest(
        reaction, replace(cfg, random_seed=10)
    )
    original = checkpoint_manifest(reaction, AAMSearchConfig())["config"]
    assert "random_seed" not in original and "sweep_cuts" not in original


def test_competition_retains_baseline_and_hard_anchors(reaction):
    from mappa.competition import compete_fragments, CompetitionConfig
    from mappa.growth.native import available

    if not available():
        pytest.skip("optional native engine not installed")
    cfg = AAMSearchConfig(branch_limit=2000, anchors=((0, 0),))
    aam = search_aam(reaction, cfg)
    result = compete_fragments(aam, CompetitionConfig(operation_budget=4, seconds=5))
    assert result.repairs  # exercise an accepted repair, not just an empty union
    assert result.baseline is aam
    assert result.counts.get("growth_calls", 0) <= 4
    for repair in result.repairs:
        assert all(
            dict(repair.graph.states[t].mapping)[0] == 0 for t in repair.graph.terminals
        )
    ids = {c.id for c in decode_events(aam).candidates}
    assert ids <= {c.id for c in decode_events(result.final_catalogue()).candidates}


def test_nondefault_seed_and_uncut_agree_across_backends(reaction):
    from mappa.growth.native import available

    if not available():
        pytest.skip("optional native engine not installed")
    cfg = AAMSearchConfig(random_seed=17, sweep_cuts=False, seed_count=2)
    expected = search_aam(reaction, cfg).graph.to_record()
    assert (
        search_aam(reaction, cfg, execution="reused_native").graph.to_record()
        == expected
    )
    # Shared policies intentionally store a different DAG with shared contexts.
    reference_ids = {c.id for c in decode_events(search_aam(reaction, cfg)).candidates}
    shared = search_aam(reaction, cfg, execution="shared_policies")
    assert {c.id for c in decode_events(shared).candidates} == reference_ids
