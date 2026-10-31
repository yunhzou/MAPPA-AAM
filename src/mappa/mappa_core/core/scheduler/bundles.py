"""Verified byte ranges in a bounded number of append-only worker bundles."""
from dataclasses import dataclass
import gzip
import hashlib
from itertools import groupby
from functools import lru_cache
import json
import os
from pathlib import Path
import pickle
import tempfile
import uuid


@dataclass(frozen=True)
class Blob:
    path: str
    offset: int
    size: int
    sha256: str


@lru_cache(maxsize=1)
def _bundle_roots():
    return json.loads(os.environ.get('MAPPA_BUNDLE_ROOTS', '{}'))


def _bundle_path(path):
    # Relocation changes the storage root, never the saved byte range or checksum.
    for source, destination in _bundle_roots().items():
        relative = Path(path).relative_to(source) if Path(path).is_relative_to(source) else None
        if relative is not None:
            return Path(destination)/relative
    return path


def read_blob(blob, scratch):
    """Only load trusted internal results; verify bytes before unpickling."""
    return next(read_blobs([blob], scratch))[1]


def read_blobs(blobs, scratch):
    """Read adjacent ranges from each bundle through one shared-file handle.

    Callers can order metadata by path and offset without scanning other bytes.
    The bounded memory spool spills larger payloads to local scratch.
    """
    with tempfile.SpooledTemporaryFile(max_size=8*2**20, dir=scratch) as local:
        for path, group in groupby(blobs, key=lambda blob: blob.path):
            with open(_bundle_path(path), 'rb') as stream:
                for blob in group:
                    yield blob, _read_verified(stream, local, blob)


def _read_verified(stream, local, blob):
    if blob.offset < 0 or blob.size <= 0:
        raise ValueError('Invalid result bundle range')
    digest = hashlib.sha256()
    local.seek(0); local.truncate()
    stream.seek(blob.offset)
    remaining = blob.size
    while remaining:
        block = stream.read(min(2**20, remaining))
        if not block:
            raise ValueError('Truncated result bundle')
        digest.update(block)
        local.write(block)
        remaining -= len(block)
    if digest.hexdigest() != blob.sha256:
        raise ValueError('Result bundle checksum differs')
    local.seek(0)
    with gzip.GzipFile(fileobj=local, mode='rb') as compressed:
        from mappa.mappa_core.utils.io.checkpoints import _load_checkpoint
        return _load_checkpoint(compressed)


class BundleWriter:
    """Each worker owns its append stream; metadata names only committed ranges.

    Serialize/compress in a bounded memory spool, spilling to local scratch.
    Publication appends bounded chunks
    without creating a shared file per task. An unacknowledged partial range is
    unreachable; completed ranges are never overwritten, even after retries.
    """
    def __init__(self, directory, scratch, *, limit_bytes=256*2**20):
        self.directory = Path(directory)
        self.scratch = None if scratch is None else Path(scratch)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.limit_bytes = limit_bytes
        self.prefix = uuid.uuid4().hex
        self.sequence, self.size = 0, 0

    def put(self, value):
        with tempfile.SpooledTemporaryFile(max_size=8*2**20, dir=self.scratch) as local:
            with gzip.GzipFile(fileobj=local, mode='wb', mtime=0, compresslevel=1) as compressed:
                pickle.dump(value, compressed, protocol=5)
            return self._append(local)

    def put_file(self, path):
        """Publish an already bundled evidence file without loading it into RAM."""
        with open(path, 'rb') as local:
            local.seek(0, os.SEEK_END)
            return self._append(local)

    def _append(self, local):
        size = local.tell()
        if self.size and self.size+size > self.limit_bytes:
            self.sequence += 1
            self.size = 0
        path = self.directory/f'{self.prefix}-{self.sequence:05d}.bundle'
        digest = hashlib.sha256()
        local.seek(0)
        with path.open('ab') as destination:
            offset = destination.tell()
            while block := local.read(2**20):
                digest.update(block)
                destination.write(block)
            destination.flush()
            os.fsync(destination.fileno())
        blob = Blob(str(path.resolve()), offset, size, digest.hexdigest())
        # Verify publication without materializing or unpickling the result.
        with path.open('rb') as stream:
            stream.seek(offset)
            observed, remaining = hashlib.sha256(), size
            while remaining:
                block = stream.read(min(2**20, remaining))
                if not block:
                    raise ValueError('Truncated publication')
                observed.update(block)
                remaining -= len(block)
        if observed.hexdigest() != blob.sha256:
            raise ValueError('Publication checksum differs')
        self.size += size
        return blob
