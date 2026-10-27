"""Atomic records and verified archives for node-local benchmark work."""
import hashlib
import json
import os
from pathlib import Path
import tarfile


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for data in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(data)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(value, separators=(',', ':'), allow_nan=False) + '\n')
    tmp.replace(path)


def publish(folder, archive):
    """Call after child processes retire: every included file is immutable."""
    folder, archive = Path(folder), Path(archive)
    archive.parent.mkdir(parents=True, exist_ok=True)
    files = {str(p.relative_to(folder)): dict(sha256=sha(p), bytes=p.stat().st_size)
             for p in sorted(folder.rglob('*')) if p.is_file() and not p.name.endswith('.tmp')
             and p.name != 'archive-manifest.json'}
    save(folder / 'archive-manifest.json', dict(files=files))
    tmp = archive.with_name(archive.name + '.tmp')
    with tarfile.open(tmp, 'w') as tar:
        for name in [*files, 'archive-manifest.json']:
            tar.add(folder / name, arcname=name, recursive=False)
    with tarfile.open(tmp) as tar:
        for name, entry in files.items():
            h = hashlib.sha256()
            with tar.extractfile(name) as stream:
                for data in iter(lambda: stream.read(1024 * 1024), b''):
                    h.update(data)
            assert h.hexdigest() == entry['sha256'], name
    tmp.replace(archive)
    return dict(path=str(archive), sha256=sha(archive), files=len(files), bytes=archive.stat().st_size)


def restore(archive, folder):
    folder = Path(folder)
    assert not folder.exists() or not any(folder.iterdir())
    folder.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as tar:
        manifest = json.load(tar.extractfile('archive-manifest.json'))
        for name, record in manifest['files'].items():
            p = Path(name)
            assert not p.is_absolute() and '..' not in p.parts
            dest = folder / p
            dest.parent.mkdir(parents=True, exist_ok=True)
            with tar.extractfile(name) as source, dest.open('wb') as target:
                import shutil
                shutil.copyfileobj(source, target)
            assert sha(dest) == record['sha256'], name
    return manifest
