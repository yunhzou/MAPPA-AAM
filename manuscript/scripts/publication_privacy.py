"""Check that publication artifacts contain no execution-site information.

CPU model names, algorithm parameters, timings, and author affiliations remain
publishable. Repository checks also inspect compressed evidence, without loading
executable pickle objects or modifying the evidence.
"""
from pathlib import Path
import argparse
import gzip
import io
import json
import re
import subprocess
import tarfile
import zipfile

from pypdf import PdfReader

# Match concrete infrastructure identifiers, not chemistry partitions or CPUs.
PRIVATE_TEXT = re.compile(
    r"/(?:lustre|project|home|h|Users|scratch|mnt)/[^\s\"']+"
    r"|\b(?:cpu|gpu)-\d{4}\b"
    r"|\boci-[A-Za-z0-9_.-]+"
    r"|\b[A-Za-z0-9_.-]*(?:login-\d+)[A-Za-z0-9_.-]*\b"
    r"|\b(?:squeue|sacct|sbatch|slurm)\b"
    r"|\b(?:concurrency limit|logical CPUs|GiB RAM|CPU nodes)\b",
    re.IGNORECASE,
)
PRIVATE_FIELDS = {
    'hostname', 'node_name', 'nodelist', 'slurm_job_id', 'slurm_partition',
    'cpus_per_task', 'mem_per_cpu', 'execution_host', 'login_host',
    'scheduler', 'allocation_id', 'cluster_name', 'concurrency',
}

# Repository documentation may teach Slurm. Concrete endpoints, account values,
# inventories and saved allocation records are a different kind of information.
SITE_TEXT = re.compile(
    r"(?P<hostname>\b[A-Za-z0-9_.-]+\.hpc\.[A-Za-z0-9.-]+\b"
    r"|\b[A-Za-z0-9_.-]*login-\d+[A-Za-z0-9_.-]*\b"
    r"|\boci-[A-Za-z0-9_.-]+)"
    r"|(?P<storage_path>(?-i:/(?:lustre|scratch)/fsw/[^\s\"'<>]+"
    r"|/(?:Users|home|project|h)/(?!YOUR_[A-Z_]+(?:/|$)|USER(?:/|$)|example(?:/|$)|\$|<)[^\s\"'<>]+))"
    r"|(?P<node_inventory>\b(?:cpu|gpu|pool)\d?-(?:\d{4,}|\[[\d,-]+\]))"
    r"|(?P<allocation_record>\b(?:JobId|ArrayJobId|NodeName|NodeList|CPUAlloc|"
    r"CPUEfctv|AllocMem|MaxTRESPU)=|\bQOSMax(?:Node|CPU|Memory)PerUserLimit\b"
    r"|(?m:^\s*JobID(?:Raw)?[|\t,])"
    r"|(?m:^\s*\d{7,}(?:_[\d,\[\]%-]+)?(?:\.[a-z]+)?[^\n]{0,160}\b"
    r"(?:RUNNING|PENDING|COMPLETED|FAILED|CANCELLED|OUT_OF_MEMORY|TIMEOUT)\b))",
    re.IGNORECASE,
)
SITE_FIELDS = {
    'hostname', 'node_name', 'nodelist', 'node_list', 'slurm_job_id',
    'slurm_partition', 'execution_host', 'login_host', 'allocation_id',
    'cluster_name', 'slurm_account', 'slurm_accounting', 'account',
    'cpus_per_task', 'mem_per_cpu', 'cpus_per_node', 'memory_gib_per_node',
    'workers_per_node',
}
ACCOUNT_DIRECTIVE = re.compile(r'^\s*#SBATCH\s+(?:--account[=\s]+|-A\s+)(\S+)', re.MULTILINE)


def _placeholder(value):
    return value.startswith(('$', '<')) or value.upper() in {
        'YOUR_ACCOUNT', 'EXAMPLE_ACCOUNT', 'PROJECT_ACCOUNT', 'ACCOUNT',
        'YOUR_HOST', 'EXAMPLE_HOST', 'HOSTNAME',
    }


def _contents(name, payload):
    """Read archive members as data; never extract paths or unpickle objects."""
    lower = name.lower()
    if lower.endswith(('.tar.gz', '.tgz', '.tar')):
        with tarfile.open(fileobj=io.BytesIO(payload), mode='r|*') as archive:
            for member in archive:
                if member.isfile():
                    with archive.extractfile(member) as stream:
                        yield from _contents(member.name, stream.read())
    elif lower.endswith(('.zip', '.npz')):
        with zipfile.ZipFile(io.BytesIO(payload)) as archive:
            for member in archive.infolist():
                if not member.is_dir():
                    yield from _contents(member.filename, archive.read(member))
    elif lower.endswith('.gz'):
        yield from _contents(name[:-3], gzip.decompress(payload))
    else:
        yield name, payload


def _pdf_text(payload):
    reader = PdfReader(io.BytesIO(payload))
    return '\n'.join([str(reader.metadata), *(page.extract_text() or '' for page in reader.pages)])


def _notebook_text(value):
    """Keep notebook prose, code and text outputs, excluding encoded media."""
    if isinstance(value, dict):
        return {key: _notebook_text(child) for key, child in value.items()
                if not key.startswith(('image/', 'audio/', 'video/'))
                and key not in {'application/pdf', 'application/octet-stream'}}
    if isinstance(value, list):
        return [_notebook_text(child) for child in value]
    return value


def _artifact_text(name, payload):
    suffix = Path(name).suffix.lower()
    if suffix == '.pdf':
        return _pdf_text(payload)
    if suffix == '.ipynb':
        return json.dumps(_notebook_text(json.loads(payload)))
    if suffix in {'.png', '.jpg', '.jpeg', '.gif', '.webp', '.tif', '.tiff'}:
        from PIL import Image
        with Image.open(io.BytesIO(payload)) as image:
            metadata = [value for value in image.info.values() if isinstance(value, str)]
            metadata.extend(value for value in image.getexif().values() if isinstance(value, str))
        return '\n'.join(metadata)
    # Encoded audiovisual samples are not text. Their captions and sidecars are
    # scanned separately; this checker does not perform OCR or media recognition.
    if suffix in {'.mp4', '.mov', '.avi', '.webm', '.mp3', '.wav', '.ogg', '.flac'}:
        return ''
    return payload.decode('utf-8', errors='replace')


def _metadata_findings(value):
    if isinstance(value, dict):
        for key, child in value.items():
            if key.lower() in SITE_FIELDS and child not in (None, '', [], {}):
                if not isinstance(child, str) or not _placeholder(child):
                    yield 'execution_metadata'
            if key.lower() in {'job', 'job_id', 'jobid'} and str(child).isdigit():
                yield 'job_identifier'
            yield from _metadata_findings(child)
    elif isinstance(value, list):
        if value and value[0] in ('sbatch', 'sacct', 'squeue', 'scontrol'):
            yield 'scheduler_command'
        for child in value:
            yield from _metadata_findings(child)
    elif isinstance(value, str):
        yield from (match.lastgroup for match in SITE_TEXT.finditer(value))


def repository_findings(name, payload, *, private_values=()):
    """Return categories only; diagnostics never print the matched values.

    Optional private_values come from an untracked local file. This catches
    standalone account/project names without publishing the denylist itself.
    """
    findings = set()
    for member, data in _contents(name, payload):
        suffix = Path(member).suffix.lower()
        text = _artifact_text(member, data)
        findings.update(match.lastgroup for match in SITE_TEXT.finditer(member + '\n' + text))
        if any(value in text for value in private_values):
            findings.add('private_identifier')
        if any(not _placeholder(match[1]) for match in ACCOUNT_DIRECTIVE.finditer(text)):
            findings.add('account_directive')
        if suffix in {'.json', '.ipynb'}:
            findings.update(_metadata_findings(json.loads(text)))
        elif suffix == '.jsonl':
            for line in text.splitlines():
                if line.strip():
                    findings.update(_metadata_findings(json.loads(line)))
    return sorted(findings)


def scan_repository(root, *, private_values=()):
    """Inspect the tracked working tree, including staged new files."""
    root = Path(root)
    names = subprocess.check_output(['git', 'ls-files', '-z'], cwd=root).decode().split('\0')
    findings = []
    for name in filter(None, names):
        path = root/name
        if not path.exists():
            continue
        reasons = repository_findings(name, path.read_bytes(), private_values=private_values)
        if reasons:
            findings.append(dict(path=name, reasons=reasons))
    return findings


def check_publication(name, payload):
    """Fail before packaging a PDF, document or evidence file with private fields."""
    if name.lower().endswith(('.gz', '.tgz', '.tar', '.zip', '.npz')):
        for member, data in _contents(name, payload):
            check_publication(member, data)
        return
    if repository_findings(name, payload):
        raise ValueError(f'Execution-site information found in publication artifact: {name}')
    suffix = Path(name).suffix.lower()
    if suffix == '.pdf':
        text = _pdf_text(payload)
    elif suffix in {'.tex', '.json', '.md', '.svg', '.bib', '.html'}:
        text = payload.decode('utf-8')
    elif suffix in {'.py', '.sh'}:
        text = payload.decode('utf-8')
        # Code may name scheduler concepts, but must not embed site identifiers.
        concrete = re.compile(r'/(?:lustre|project|home|h|Users|scratch|mnt)/[^\s\"\']+|\b(?:cpu|gpu)-\d{4}\b|\boci-[A-Za-z0-9_.-]+', re.IGNORECASE)
        if concrete.search(text):
            raise ValueError(f'Execution-site identifier found in publication source: {name}')
        return
    elif suffix == '.sbatch':
        raise ValueError(f'Execution script is not a publication artifact: {name}')
    else:
        return
    if PRIVATE_TEXT.search(text):
        # Do not echo the sensitive matched value to logs.
        raise ValueError(f'Execution-site information found in publication artifact: {name}')
    if suffix == '.json':
        def walk(value):
            if isinstance(value, dict):
                if PRIVATE_FIELDS.intersection(k.lower() for k in value):
                    raise ValueError(f'Execution metadata field found in publication artifact: {name}')
                for child in value.values():
                    walk(child)
            elif isinstance(value, list):
                for child in value:
                    walk(child)
        walk(json.loads(text))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, required=True)
    parser.add_argument('--private-values', type=Path,
                        help='Untracked local JSON list of additional exact private identifiers')
    args = parser.parse_args()
    values = json.loads(args.private_values.read_text()) if args.private_values else []
    if not isinstance(values, list) or any(not isinstance(value, str) or not value for value in values):
        parser.error('Private values must be a JSON list of nonempty strings')
    findings = scan_repository(args.repository, private_values=values)
    print(json.dumps(findings, indent=2))
    raise SystemExit(bool(findings))


if __name__ == '__main__':
    main()
