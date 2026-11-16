"""Intermediate transport; scientific operations only receive loaded objects."""
from dataclasses import dataclass
from pathlib import Path
import threading
import uuid

from .bundles import Blob, BundleWriter, read_blob


@dataclass(frozen=True)
class MemoryRef:
    store: str
    key: str


class MemoryStore:
    """Own intermediates in RAM. Remote clients transfer objects through the pool."""
    def __init__(self):
        self.identity = uuid.uuid4().hex
        self._values = {}
        self._lock = threading.RLock()

    def put(self, value):
        reference = MemoryRef(self.identity, uuid.uuid4().hex)
        with self._lock:
            self._values[reference.key] = value
        return reference

    def load(self, reference):
        if reference.store != self.identity:
            raise ValueError('Reference belongs to another memory store')
        with self._lock:
            return self._values[reference.key]


class BundleStore:
    """Append results to worker-owned bundles; never create a file per task."""
    def __init__(self, directory, *, scratch=None, limit_bytes=256*2**20):
        self.directory = Path(directory).resolve()
        self.directory.mkdir(parents=True, exist_ok=True)
        self.scratch = scratch
        self.limit_bytes = limit_bytes
        self._writer = None

    def __getstate__(self):
        # Each compute process owns its own append stream.
        return {**vars(self), '_writer': None}

    def put(self, value):
        if self._writer is None:
            self._writer = BundleWriter(self.directory, self.scratch, limit_bytes=self.limit_bytes)
        return self._writer.put(value)

    def load(self, reference):
        return read_blob(reference, self.scratch)


ArtifactRef = Blob | MemoryRef
