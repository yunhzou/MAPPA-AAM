"""A producing-revision oracle for explicit saved Golden input conversion."""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'bench'))
from catalogue_migration import convert_blob
from mappa import AAMSearchConfig, decode_events, decoded_events_record, decoded_events_from_record
from mappa.mappa_core.core.scheduler.bundles import Blob, BundleWriter, read_blob


def historical_input():
    path = Path(__file__).parent/'fixtures/fragment_matching/golden-625-fragment-v1.bundle'
    raw = path.read_bytes()
    reference = Blob(str(path), 0, len(raw), hashlib.sha256(raw).hexdigest())
    config = json.loads(json.dumps(asdict(AAMSearchConfig(
        seed_count=1, branch_limit=200, iso_tolerance=2., random_seed=42))))
    config['self_tolerance'] = 1.
    return reference, config


def test_conversion_preserves_every_recorded_field_except_approved_policy(tmp_path):
    reference, config = historical_input()
    saved = convert_blob(reference, config)
    record = saved.catalogue.to_record()
    # Independently exported using producing code d30a935, changing only the
    # fragment schema tag and removing the two self-tolerance fields.
    assert hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest() == (
        'f3c56e3d15daba316148b4a63a3df27711c3a6bb2112aa640d32c6c926fc5711')
    assert (saved.case, len(record['families']), len(record['fragments']), record['input_paths']) == (
        '00625/forward', 27, 36, 27)
    converted = BundleWriter(tmp_path, None).put(saved)
    assert read_blob(converted, None).catalogue.to_record() == record
    decoded = decoded_events_record(decode_events(saved.catalogue, event_slack=5))
    assert decoded['minimum_events'] == 3
    assert decoded_events_record(decoded_events_from_record(decoded)) == decoded
    assert saved.catalogue.to_record() == record


@pytest.mark.parametrize('corruption', ['checksum', 'configuration'])
def test_conversion_rejects_wrong_input(corruption):
    reference, config = historical_input()
    if corruption == 'checksum':
        reference = Blob(reference.path, reference.offset, reference.size, '0'*64)
    else:
        config['iso_tolerance'] = 3.
    with pytest.raises(ValueError, match='differs'):
        convert_blob(reference, config)
