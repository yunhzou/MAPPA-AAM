"""Relocate a saved cut index after an explicit filesystem mount change.

Only path addresses change. Verified member hashes, offsets, task settings and
the original index remain intact; this never searches or decodes molecules.
"""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
import os
from pathlib import Path
import tempfile

from mappa.mappa_core.core.scheduler import Blob, BundleWriter, read_blob


def relocate(value, old, new):
    if isinstance(value, str):
        return new + value[len(old):] if value == old or value.startswith(old + '/') else value
    if isinstance(value, dict):
        return {key: relocate(item, old, new) for key, item in value.items()}
    if isinstance(value, list):
        return [relocate(item, old, new) for item in value]
    if isinstance(value, tuple):
        return tuple(relocate(item, old, new) for item in value)
    return value


def prepare(source, destination, old, new, expected_sha256):
    data = source.read_bytes()
    if hashlib.sha256(data).hexdigest() != expected_sha256:
        raise ValueError('Original cut index hash differs')
    original = json.loads(data)
    destination.mkdir(parents=True, exist_ok=False)
    result = dict(original, tasks=[])
    checked_archives = set()
    with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
        writer = BundleWriter(destination/'bundles', scratch)
        for row in original['tasks']:
            location = relocate(row['blob'], old, new)
            tasks = read_blob(Blob(**location), scratch)
            changed = [replace(task, payload=relocate(task.payload, old, new)) for task in tasks]
            assert len(changed) == row['count']
            for task in changed:
                # Check addresses without rescanning or extracting archive trees.
                for member in [*task.payload['source'].values(), task.payload['raw'], task.payload['graph']]:
                    archive = member['archive']
                    if archive not in checked_archives:
                        if not Path(archive).is_file():
                            raise FileNotFoundError(archive)
                        checked_archives.add(archive)
            result['tasks'].append(dict(row, blob=asdict(writer.put(changed))))
    result['relocation'] = dict(source=str(source), sha256=expected_sha256,
                                old_root=old, new_root=new, scientific_work_repeated=False)
    (destination/'index.json').write_text(json.dumps(result, separators=(',', ':'))+'\n')
    print(json.dumps(dict(directions=len(result['tasks']), cut_tasks=sum(r['count'] for r in result['tasks']),
        index_sha256=hashlib.sha256((destination/'index.json').read_bytes()).hexdigest())), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--old-root', required=True)
    parser.add_argument('--new-root', required=True)
    parser.add_argument('--sha256', required=True)
    args = parser.parse_args()
    prepare(args.source, args.destination, args.old_root, args.new_root, args.sha256)
