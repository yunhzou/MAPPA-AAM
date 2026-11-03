"""One reaction's checkpoint index; cut payloads share worker bundles."""
from dataclasses import asdict
import json
from pathlib import Path
import sqlite3

from .bundles import Blob


class CutArchive:
    def __init__(self, directory, manifest, *, resume=False):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        path = self.directory/'manifest.json'
        index = self.directory/'archive.sqlite'
        if resume:
            if json.loads(path.read_text()) != manifest:
                raise ValueError('Checkpoint input or configuration differs from this search')
            if not index.exists():
                raise ValueError('This run requires an aggregated checkpoint index')
        else:
            if path.exists() or index.exists():
                raise ValueError('Existing checkpoints require explicit resume=True')
            path.write_text(json.dumps(manifest)+'\n')
        self.db = sqlite3.connect(index)
        self.db.execute('CREATE TABLE IF NOT EXISTS cuts (number INTEGER PRIMARY KEY, reference TEXT)')
        self.db.commit()

    def references(self):
        return {number: Blob(**json.loads(raw)) for number, raw in
                self.db.execute('SELECT number,reference FROM cuts ORDER BY number')}

    def record(self, number, reference):
        with self.db:
            self.db.execute('INSERT INTO cuts VALUES (?,?)', (number, json.dumps(asdict(reference))))

    def close(self):
        self.db.close()
