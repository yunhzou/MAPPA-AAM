"""Namespace relocation preserves values; obsolete scientific schemas are rejected."""
import json
import io
import pickle
from importlib.resources import files
from pathlib import Path

import mappa
import pytest
from mappa.mappa_core.utils.io.records import aam_from_record
from mappa.mappa_core.utils.io.checkpoints import read_aam_checkpoint, read_graph_checkpoint, read_raw_cut
from mappa.mappa_core.core.matching.growth import native

FIXTURES = Path(__file__).parent / "fixtures" / "package_rename"


def test_previous_package_checkpoint_global_resolves_to_mappa():
    from mappa.mappa_core.utils.io.checkpoints import _CheckpointUnpickler
    # A protocol-0 GLOBAL from the previous public package namespace.
    stream = io.BytesIO(b"cgraft.domain\nAAMResult\n.")
    assert _CheckpointUnpickler(stream).load() is mappa.AAMResult


def test_archived_bond_response_keeps_values_and_sharing_after_core_move():
    import numpy as np
    from mappa.mappa_core.core.symmetry.responses import BondResponse
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint
    response = BondResponse(1., .5, np.array([.4, 1., 1.8]), np.array([1, 0, 3]))
    payload = pickle.dumps((response, response), protocol=0)
    current = b'cmappa.mappa_core.core.symmetry.responses\nBondResponse\n'
    archived = b'cmappa.mappa_core.post_processing.decoder.records\nBondEventResponse\n'
    assert payload.count(current) == 1
    first, second = _loads_checkpoint(payload.replace(current, archived))
    assert type(first) is BondResponse and first is second
    assert (first.source_weight, first.threshold) == (1., .5)
    assert first.lookup(np.array([1.8, 1.])).tolist() == [3, 0]
    assert first.constant is None


def test_current_decoder_configuration_keeps_its_type_through_bundle_storage(tmp_path):
    from mappa.execution import BundleStore
    from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig
    from mappa.mappa_core.core.scheduler.results import DecodingResult
    original = DecodingResult('case', EventDecodeConfig(.6, .4, event_slack=3), (), (), {},
                              mappa.AAMSearchConfig(), 'reused_native')
    store = BundleStore(tmp_path)
    loaded = store.load(store.put(original))
    assert type(loaded.config) is EventDecodeConfig
    assert loaded.config.event_slack == 3
    assert loaded == original


def test_current_decoder_configuration_rejects_obsolete_rematching_fields():
    from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint
    obsolete = EventDecodeConfig()
    object.__setattr__(obsolete, 'iso_tolerance', 1.)
    with pytest.raises(TypeError, match='iso_tolerance'):
        _loads_checkpoint(pickle.dumps(obsolete))


def test_archived_execution_options_are_not_silently_accepted():
    archived = json.loads((FIXTURES / "legacy-aam.json").read_text())
    with pytest.raises(ValueError, match='checkpoint schema'):
        read_aam_checkpoint(FIXTURES / "legacy-aam.pkl.gz")
    with pytest.raises(ValueError, match='result schema'):
        aam_from_record(archived)


def test_old_graph_formats_cannot_reinterpret_retired_fragment_relationships():
    archived = json.loads((FIXTURES / "legacy-aam.json").read_text())
    with pytest.raises(ValueError, match='graph schema'):
        mappa.AAMSearchGraph.from_record(archived['graph'])
    with pytest.raises(ValueError, match='cut checkpoint schema'):
        read_graph_checkpoint(FIXTURES / "legacy-finalized.pkl.gz")
    with pytest.raises(ValueError, match='raw cut checkpoint schema'):
        read_raw_cut(FIXTURES / "legacy.raw.pkl.gz")


def test_installed_viewer_assets_are_available_under_mappa():
    assert files("mappa").joinpath("presentation/static/3Dmol-min.js").is_file()
    assert files("mappa").joinpath("presentation/static/aam_search.html").is_file()


def test_native_setting_uses_only_current_name(monkeypatch):
    monkeypatch.setattr(native, "_engine", object())
    monkeypatch.delenv("MAPPA_NATIVE", raising=False)
    monkeypatch.delenv("RXN_CORE_NATIVE", raising=False)
    assert native.available()
    monkeypatch.setenv("RXN_CORE_NATIVE", "0")
    assert native.available()
    monkeypatch.setenv("MAPPA_NATIVE", "1")
    assert native.available()
    monkeypatch.setenv("MAPPA_NATIVE", "0")
    monkeypatch.setenv("RXN_CORE_NATIVE", "1")
    assert not native.available()
