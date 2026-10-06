"""The public namespace can change without invalidating archived chemistry."""
import json
import io
from importlib.resources import files
from pathlib import Path

import mappa
from mappa.artifacts import (
    aam_from_record, aam_record, read_aam_checkpoint, read_graph_checkpoint,
    read_raw_cut, write_aam_checkpoint,
)
from mappa.growth import native
from mappa.postprocessing import decode_events

FIXTURES = Path(__file__).parent / "fixtures" / "package_rename"


def test_previous_package_checkpoint_global_resolves_to_mappa():
    from mappa.artifacts import _CheckpointUnpickler
    # A protocol-0 GLOBAL from the previous public package namespace.
    stream = io.BytesIO(b"cgraft.domain\nAAMResult\n.")
    assert _CheckpointUnpickler(stream).load() is mappa.AAMResult


def test_old_checkpoint_globals_resolve_to_mappa_and_preserve_events(tmp_path):
    archived = json.loads((FIXTURES / "legacy-aam.json").read_text())
    result = read_aam_checkpoint(FIXTURES / "legacy-aam.pkl.gz")
    assert isinstance(result, mappa.AAMResult)
    assert isinstance(result.graph, mappa.AAMSearchGraph)
    assert type(result.graph).__module__.startswith("mappa.")
    # JSON-normalize tuple/list differences, without changing the payload.
    assert json.loads(json.dumps(aam_record(result))) == archived
    assert json.loads(json.dumps(aam_record(aam_from_record(archived)))) == archived
    expected = json.loads((FIXTURES / "legacy-events.json").read_text())
    assert sorted(candidate.id for candidate in decode_events(result).candidates) == expected
    current = tmp_path / "mappa.pkl.gz"
    write_aam_checkpoint(result, current)
    assert read_aam_checkpoint(current).graph == result.graph


def test_both_old_binary_cut_formats_keep_their_graph():
    expected = read_aam_checkpoint(FIXTURES / "legacy-aam.pkl.gz").graph
    assert read_graph_checkpoint(FIXTURES / "legacy-finalized.pkl.gz") == expected
    assert read_raw_cut(FIXTURES / "legacy.raw.pkl.gz") == expected


def test_installed_viewer_assets_are_available_under_mappa():
    assert files("mappa").joinpath("static/3Dmol-min.js").is_file()
    assert files("mappa").joinpath("static/aam_search.html").is_file()


def test_native_setting_accepts_legacy_fallback_and_prefers_new_name(monkeypatch):
    monkeypatch.setattr(native, "_engine", object())
    monkeypatch.delenv("MAPPA_NATIVE", raising=False)
    monkeypatch.delenv("RXN_CORE_NATIVE", raising=False)
    assert native.available()
    monkeypatch.setenv("RXN_CORE_NATIVE", "0")
    assert not native.available()
    monkeypatch.setenv("MAPPA_NATIVE", "1")
    assert native.available()
    monkeypatch.setenv("MAPPA_NATIVE", "0")
    monkeypatch.setenv("RXN_CORE_NATIVE", "1")
    assert not native.available()
