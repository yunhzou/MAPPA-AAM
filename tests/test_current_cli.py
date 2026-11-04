"""The default CLI decodes events and cannot dispatch to legacy R/P selection."""
import importlib
import json

import numpy as np
import pytest

import mappa
from mappa.cli import main


def endpoints(tmp_path):
    paths = [tmp_path / 'reactant.npz', tmp_path / 'product.npz']
    for path in paths:
        np.savez(path, elements=np.array(['H', 'H']),
                 coordinates=np.array([[0., 0., 0.], [0., 0., 1.]]),
                 wbo=np.array([[0., 1.], [1., 0.]]))
    return ['--reactant-npz', str(paths[0]), '--product-npz', str(paths[1])]


@pytest.mark.smoke
def test_default_cli_decodes_without_any_geometry_selection(tmp_path, monkeypatch):
    assert not hasattr(mappa, 'select_rp_mappings')
    output = tmp_path / 'events'
    main([*endpoints(tmp_path), '--output', str(output)])
    record = json.loads((output / 'events.json').read_text())
    assert record['schema'] == 'mappa.decoded_events/v10'
    assert record['model'] == 'independent_fragments_with_sweep'
    assert record['approximate'] and record['catalogue_preserved']
    assert not record['mapping_space_preserved']
    assert not record['minimum_proven'] and record['model_minimum_proven']
    assert record['event_window']['minimum'] == 0
    assert record['event_window']['limit'] == record['event_window']['slack'] == 6
    assert not record['global_minimum_proven']
    assert 'representatives' not in record
    assert record['families'] and record['tables']
    assert record['responses'] == []
    assert not record['bin_witnesses_guaranteed']
    for family in record['families']:
        for partition, table_id in zip(family['partitions'], family['tables']):
            table = record['tables'][table_id]
            # The uncut H-H fragment has only zero events. Independent swept
            # endpoints can also choose the same target H, producing model bin2;
            # that bin is explicitly not a valid complete bijection certificate.
            expected = [0] if table['sweep_cut'] is None else [0, 2]
            assert record['event_window']['event_counts'][family['family_id']][partition] == expected
    assert (output / 'aam.json').is_file()
    assert (output / 'search.html').is_file()
    assert not (output / 'rp.json').exists()
    assert not (output / 'reaction.json').exists()


def test_decoder_policy_rejects_ambiguous_historical_configs_and_keyword_overrides():
    import pickle
    from mappa import EventDecodeConfig

    assert EventDecodeConfig().event_slack == 6
    config = EventDecodeConfig(event_slack=None)
    assert pickle.loads(pickle.dumps(config)) == config
    historical = object.__new__(EventDecodeConfig)
    object.__setattr__(historical, 'threshold', .5)
    object.__setattr__(historical, 'metal_threshold', .3)
    with pytest.raises(ValueError, match='producing revision'):
        pickle.loads(pickle.dumps(historical))
    for slack in (None, 6):
        with pytest.raises(ValueError, match='not both'):
            mappa.decode_events(None, config, event_slack=slack)


@pytest.mark.smoke
def test_raw_cli_does_not_decode(tmp_path, monkeypatch):
    import mappa.cli as cli
    def forbidden(*args, **kwargs):
        pytest.fail('raw AAM stage invoked event decoding')
    monkeypatch.setattr(cli, 'decode_events', forbidden)
    output = tmp_path / 'raw'
    main([*endpoints(tmp_path), '--stage', 'aam', '--output', str(output)])
    assert (output / 'aam.json').is_file()
    assert not (output / 'events.json').exists()


@pytest.mark.parametrize('args', [
    ['--stage', 'rp'], ['--stage', 'full'], ['--post-workers', '2'],
    ['--stage', 'ts'], ['--target-npz', 'target.npz'],
])
def test_removed_or_implicit_workflows_are_rejected(args):
    with pytest.raises(SystemExit) as exc:
        main(args)
    assert exc.value.code == 2


def test_legacy_selector_has_no_import_or_compatibility_alias():
    for name in ('align_reaction', 'select_rp_mappings', 'RPResult',
                 'RPMechanism', 'write_rp_bundle', 'reaction_context_from_rp'):
        assert not hasattr(mappa, name)
    for name in ('mappa.rp', 'mappa.alignment'):
        with pytest.raises(ModuleNotFoundError) as exc:
            importlib.import_module(name)
        assert exc.value.name == name


def test_ts_cli_requires_and_consumes_explicit_reaction(tmp_path, monkeypatch):
    from mappa import cli, ReactionContext
    from mappa.mappa_workflow.ts_analytics.records import TSResult
    arguments = endpoints(tmp_path)
    reaction = tmp_path / 'reaction.json'
    reaction.write_text(json.dumps({
        'schema': 'rxn_core.resolved_reaction/v1', 'atom_count': 2,
        'mechanisms': [{'mapping': {'0': 0, '1': 1},
                        'broken_bonds': [], 'formed_bonds': [], 'core_atoms': []}],
    }))
    target = tmp_path / 'target.npz'
    np.savez(target, elements=['H', 'H'], coordinates=np.zeros((2, 3)),
             wbo=np.zeros((2, 2)), frequencies=[-100.], modes=np.zeros((1, 2, 3)))
    def forbidden(*args, **kwargs):
        pytest.fail('TS CLI invoked an endpoint search or decoder')
    monkeypatch.setattr(cli, 'search_aam', forbidden)
    monkeypatch.setattr(cli, 'decode_events', forbidden)
    seen = []
    def analyze(context, target, **kwargs):
        assert isinstance(context, ReactionContext)
        seen.append(context.mechanisms[0].mapping.as_dict())
        return TSResult(context, (), 0.)
    monkeypatch.setattr(cli, 'analyze_transition_state', analyze)
    main([*arguments, '--stage', 'ts', '--reaction-json', str(reaction),
          '--target-npz', str(target), '--output', str(tmp_path / 'ts')])
    assert seen == [{0: 0, 1: 1}]
    assert (tmp_path / 'ts/ts_001.json').is_file()


def test_cached_ts_utility_requires_explicit_records():
    import runpy
    from pathlib import Path
    utility = runpy.run_path(str(
        Path(__file__).resolve().parents[1] / 'tools/run_typed_ts_case.py'))
    with pytest.raises(SystemExit) as exc:
        utility['main']([
            '--manifest', 'unused', '--work-root', 'unused',
            '--benchmark-root', 'unused', '--output-root', 'unused',
        ])
    assert exc.value.code == 2


def test_retired_alignment_and_geometry_selection_are_not_exported():
    import mappa.mappa_core.domain as structures
    for name in ('PostAAMSelectionProblem', 'PostAAMMechanism',
                 'AffineChiralityConstraint', 'FixedMappingRMSD'):
        assert not hasattr(structures, name)
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module('mappa.alignment')


def test_removed_mechanism_pipeline_has_no_api_or_options():
    from mappa import AAMSearchConfig
    for name in ('group_mechanisms', 'compile_mechanism_families', 'AAMMechanism',
                 'MechanismResult', 'AnalyticalMechanism', 'AnalyticalAAMResult'):
        assert not hasattr(mappa, name)
    with pytest.raises(ModuleNotFoundError) as exc:
        importlib.import_module('mappa.mechanisms')
    assert exc.value.name == 'mappa.mechanisms'
    for option in ('symmetry_repair', 'symmetry_repair_min_changes',
                   'symmetry_repair_max_evaluations'):
        with pytest.raises(TypeError, match=option):
            AAMSearchConfig(**{option: 1})
