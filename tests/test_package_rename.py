"""Archived graph namespaces remain readable; retired execution options are rejected."""
import json
import io
from importlib.resources import files
from pathlib import Path

import mappa
import pytest
from mappa.artifacts import (
    aam_from_record, read_aam_checkpoint, read_graph_checkpoint, read_raw_cut,
)
from mappa.growth import native

FIXTURES = Path(__file__).parent / "fixtures" / "package_rename"


def test_previous_package_checkpoint_global_resolves_to_mappa():
    from mappa.artifacts import _CheckpointUnpickler
    # A protocol-0 GLOBAL from the previous public package namespace.
    stream = io.BytesIO(b"cgraft.domain\nAAMResult\n.")
    assert _CheckpointUnpickler(stream).load() is mappa.AAMResult


def test_archived_execution_options_are_not_silently_accepted():
    archived = json.loads((FIXTURES / "legacy-aam.json").read_text())
    with pytest.raises(TypeError, match='symmetry_repair'):
        read_aam_checkpoint(FIXTURES / "legacy-aam.pkl.gz")
    with pytest.raises(TypeError, match='symmetry_repair'):
        aam_from_record(archived)


def test_both_old_binary_cut_formats_keep_their_graph():
    archived = json.loads((FIXTURES / "legacy-aam.json").read_text())
    expected = mappa.AAMSearchGraph.from_record(archived['graph'])
    finalized = read_graph_checkpoint(FIXTURES / "legacy-finalized.pkl.gz")
    raw = read_raw_cut(FIXTURES / "legacy.raw.pkl.gz")
    assert finalized == raw == expected
    assert type(finalized).__module__.startswith('mappa.')
    assert all(context.growth_policy == 'legacy_wbo_priority'
               for context in finalized.contexts)


def test_installed_viewer_assets_are_available_under_mappa():
    assert files("mappa").joinpath("static/3Dmol-min.js").is_file()
    assert files("mappa").joinpath("static/aam_search.html").is_file()


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
