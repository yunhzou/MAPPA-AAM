"""Fingerprint the actual Python modules and compiled engines shipped to workers."""
from functools import lru_cache
import hashlib
from pathlib import Path


@lru_cache(maxsize=1)
def code_identity():
    root = Path(__file__).resolve().parents[1]
    paths = sorted(set(root.rglob('*.py')) | set(root.glob('_*.so')) | set(root.glob('_*.pyd')))
    digest = hashlib.sha256()
    for path in paths:
        digest.update(str(path.relative_to(root)).encode())
        with path.open('rb') as source:
            while block := source.read(2**20):
                digest.update(block)
    return digest.hexdigest()
