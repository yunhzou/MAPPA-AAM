import json

import numpy as np

from mappa.cli import main
from mappa import decode_events
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
from mappa.mappa_core.utils.io.records import problem_from_record


def _write_endpoint(path, elements, coordinates, edges, *, modes=False):
    wbo = np.zeros((len(elements), len(elements)))
    for left, right, value in edges:
        wbo[left, right] = wbo[right, left] = value
    payload = {
        "elements": np.asarray(elements),
        "coordinates": np.asarray(coordinates),
        "wbo": wbo,
    }
    if modes:
        payload["frequencies"] = np.asarray([-500.0])
        payload["modes"] = np.ones((1, len(elements), 3))
    np.savez(path, **payload)


def test_typed_cli_decodes_then_uses_explicit_reaction_for_ts(tmp_path):
    elements = ("C", "N", "O", "H")
    xyz = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0], [1.0, 1.0, 0.2],
    ])
    reactant, product, target = (
        tmp_path / "R.npz", tmp_path / "P.npz", tmp_path / "TS.npz")
    _write_endpoint(
        reactant, elements, xyz, ((0, 1, 1.0), (2, 3, 1.0)))
    _write_endpoint(
        product, elements, xyz + 0.1, ((0, 2, 1.0), (1, 3, 1.0)))
    _write_endpoint(
        target, elements, xyz + 0.05,
        ((0, 1, 0.5), (2, 3, 0.5), (0, 2, 0.5), (1, 3, 0.5)),
        modes=True)
    output = tmp_path / "result"

    main([
        "--reactant-npz", str(reactant),
        "--product-npz", str(product),
        "--output", str(output),
        "--seed-count", "1",
        "--workers", "1",
    ])

    decoded = json.loads((output / "events.json").read_text())
    assert decoded['schema'] == 'mappa.decoded_events/v10'
    assert decoded['catalogue_preserved'] and not decoded["mapping_space_preserved"]
    assert decoded['approximate'] and not decoded['minimum_proven']
    assert decoded['event_window']['minimum'] == 4
    assert decoded['event_window']['limit'] == 10
    assert not decoded["exhaustive_event_patterns"]
    assert not decoded["global_minimum_proven"]
    assert 'representatives' not in decoded
    catalogue = FinalBranchCatalogue.from_record(
        problem_from_record(decoded['problem']), decoded['catalogue'])
    interpretation = decode_events(catalogue)
    for family_id in interpretation.retained_family_ids:
        candidate = interpretation.evaluate(family_id)
        assert candidate["mapping"] == {0: 0, 1: 1, 2: 2, 3: 3}
        assert candidate['counts'] == dict(breaking=2, forming=2, strengthening=0, weakening=0)
        assert sorted(event['reactant_pair'] for event in candidate['events']
                      if event['category'] == 'breaking') == [(0, 1), (2, 3)]
        assert sorted(event['product_pair'] for event in candidate['events']
                      if event['category'] == 'forming') == [(0, 2), (1, 3)]
    candidate = interpretation.evaluate(0)
    view = (output / "search.html").read_text()
    assert "$3Dmol.createViewer" in view
    assert "<script src=" not in view
    assert not (output / "rp.json").exists()
    assert not (output / "reaction.json").exists()
    # Choose a recorded witness explicitly; decoding does not select a mechanism.
    record = {
        "schema": "rxn_core.resolved_reaction/v1", "atom_count": 4,
        "mechanisms": [{
            "mapping": {str(atom): image for atom, image in candidate["mapping"].items()},
            "broken_bonds": [event['reactant_pair'] for event in candidate['events']
                             if event['category'] == 'breaking'],
            "formed_bonds": [event['product_pair'] for event in candidate['events']
                             if event['category'] == 'forming'],
            "core_atoms": [0, 1, 2, 3],
        }],
    }
    (output / "reaction.json").write_text(json.dumps(record))

    ts_only = tmp_path / "ts_only"
    main([
        "--stage", "ts",
        "--reactant-npz", str(reactant),
        "--product-npz", str(product),
        "--reaction-json", str(output / "reaction.json"),
        "--target-npz", str(target),
        "--output", str(ts_only),
        "--seed-count", "1",
    ])
    assert (ts_only / "ts_001.json").exists()
    assert not (ts_only / "rp.json").exists()
    ts = json.loads((ts_only / "ts_001.json").read_text())
    assert ts["schema"] == "rxn_core.ts/v2"
    assert ts["mechanisms"][0]["selected"]["sources"] == [
        "product", "reactant"]
