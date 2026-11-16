"""Publication and repository boundaries for saved execution evidence."""
import gzip
import base64
import importlib.util
import io
import json
from pathlib import Path
import pickle
import subprocess
import tarfile
import zipfile

import pytest
from pypdf import PdfWriter
from PIL import Image, PngImagePlugin


CHECKER = Path(__file__).resolve().parents[1]/'manuscript/scripts/publication_privacy.py'
spec = importlib.util.spec_from_file_location('publication_privacy', CHECKER)
privacy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(privacy)


def endpoint():
    return f'login-{12:02d}.hpc.example.invalid'


def test_generic_scheduler_example_is_repository_documentation_not_paper_hardware():
    example = b'#SBATCH --account=YOUR_ACCOUNT\n#SBATCH --cpus-per-task=8\n#SBATCH --mem=20G\n'
    assert privacy.repository_findings('example.sbatch', example) == []
    assert privacy.repository_findings('example.md', b'Use Slurm with 8 logical CPUs.') == []
    assert privacy.repository_findings('measurements.json', json.dumps({
        'cpu_seconds': 12.5, 'memory_bytes': 1024, 'branch_limit': 100,
        'cpu_model': 'Example CPU', 'affiliation': 'NVIDIA',
    }).encode()) == []
    with pytest.raises(ValueError, match='Execution-site information'):
        privacy.check_publication('paper.md', b'Use Slurm with 8 logical CPUs.')
    with pytest.raises(ValueError, match='Execution metadata field'):
        privacy.check_publication('paper.json', b'{"concurrency":8}')
    privacy.check_publication('paper.md', b'CPU: Example model. Affiliation: NVIDIA.')


def test_nested_archive_is_scanned_and_pickle_is_never_loaded(monkeypatch):
    def must_not_unpickle(*args, **kwargs):
        raise AssertionError('Evidence scanning must not execute pickle loading')
    payload = pickle.dumps({'saved_host': endpoint()}, protocol=5)
    monkeypatch.setattr(pickle, 'loads', must_not_unpickle)
    archive = io.BytesIO()
    with tarfile.open(fileobj=archive, mode='w:gz') as saved:
        data = gzip.compress(payload)
        member = tarfile.TarInfo('case.pkl.gz')
        member.size = len(data)
        saved.addfile(member, io.BytesIO(data))
    bundle = io.BytesIO()
    with zipfile.ZipFile(bundle, 'w') as saved:
        saved.writestr('evidence.tar.gz', archive.getvalue())
    findings = privacy.repository_findings('release.zip', bundle.getvalue())
    assert findings == ['hostname']
    assert endpoint() not in json.dumps(findings)
    with pytest.raises(ValueError, match='Execution-site information'):
        privacy.check_publication('release.zip', bundle.getvalue())


def test_notebook_output_and_pdf_metadata_are_checked():
    notebook = {'cells': [{'outputs': [{'text': [endpoint()]}]}]}
    assert privacy.repository_findings('demo.ipynb', json.dumps(notebook).encode()) == ['hostname']
    writer = PdfWriter()
    writer.add_blank_page(width=72, height=72)
    writer.add_metadata({'/Subject': endpoint()})
    buffer = io.BytesIO()
    writer.write(buffer)
    assert privacy.repository_findings('paper.pdf', buffer.getvalue()) == ['hostname']
    with pytest.raises(ValueError):
        privacy.check_publication('paper.pdf', buffer.getvalue())


@pytest.mark.parametrize('record', [
    {'host': {'hostname': 'example.invalid'}},
    {'job_id': 123456},
    {'command': ['sbatch', '--cpus-per-task=8']},
])
def test_structured_operational_records_are_not_scientific_results(record):
    assert privacy.repository_findings('run.json', json.dumps(record).encode())


def test_concrete_account_and_raw_inventory_are_rejected():
    account = 'fictional_research'
    script = f'#SBATCH --account={account}\n'.encode()
    assert privacy.repository_findings('run.sbatch', script) == ['account_directive']
    node = f'cpu9-{42:05d}'
    assert privacy.repository_findings('allocation.txt', node.encode()) == ['node_inventory']
    inventory = 'Node' + 'Name=' + 'example.invalid'
    assert privacy.repository_findings('allocation.txt', inventory.encode()) == ['allocation_record']
    rows = f'{12345678}_2 RUNNING 100 5\n{12345679}_[3-8] PENDING 100 5'
    assert privacy.repository_findings('queue.json', json.dumps({'jobs': rows}).encode()) == ['allocation_record']


def test_private_identifiers_can_stay_outside_the_public_checker():
    value = 'fictional_private_project'
    assert privacy.repository_findings('notes.txt', value.encode(), private_values=[value]) == ['private_identifier']


def test_private_absolute_paths_and_explicit_placeholders():
    for directory in ['home', 'Users', 'project', 'h']:
        value = f'/{directory}/fictional_owner/results.json'
        assert privacy.repository_findings('notes.md', value.encode()) == ['storage_path']
        placeholder = f'/{directory}/YOUR_USER/results.json'
        assert privacy.repository_findings('example.md', placeholder.encode()) == []


def test_chemical_labels_and_encoded_pixels_are_not_paths():
    assert privacy.repository_findings('chemistry.md', b'Compare C/H/N labels.') == []
    buffer = io.BytesIO()
    Image.new('RGB', (8, 8), color='white').save(buffer, format='PNG')
    notebook = {'cells': [{'outputs': [{'data': {
        'image/png': base64.b64encode(buffer.getvalue()).decode(),
        'text/plain': ['C/H/N'],
    }}]}]}
    assert privacy.repository_findings('image.ipynb', json.dumps(notebook).encode()) == []
    metadata = PngImagePlugin.PngInfo()
    metadata.add_text('Description', endpoint())
    buffer = io.BytesIO()
    Image.new('RGB', (8, 8)).save(buffer, format='PNG', pnginfo=metadata)
    assert privacy.repository_findings('image.png', buffer.getvalue()) == ['hostname']


def test_repository_scan_includes_staged_files_not_untracked_private_archives(tmp_path):
    subprocess.run(['git', 'init', '-q', str(tmp_path)], check=True)
    (tmp_path/'results.json').write_text('{"cpu_seconds":12.5}')
    (tmp_path/'tracked.json.gz').write_bytes(gzip.compress(json.dumps({'hostname': endpoint()}).encode()))
    (tmp_path/'untracked.json').write_text(json.dumps({'hostname': endpoint()}))
    subprocess.run(['git', '-C', str(tmp_path), 'add', 'results.json', 'tracked.json.gz'], check=True)
    assert privacy.scan_repository(tmp_path) == [
        {'path': 'tracked.json.gz', 'reasons': ['execution_metadata', 'hostname']},
    ]
    (tmp_path/'tracked.json.gz').unlink()
    assert privacy.scan_repository(tmp_path) == []
