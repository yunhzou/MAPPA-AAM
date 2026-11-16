"""One persistent worker service per node in an existing Slurm allocation."""
import argparse
import os
import tempfile

from mappa.mappa_core.core.scheduler.runtime import connect, connect_store, run
from mappa.mappa_core.core.scheduler.storage import BundleStore
from .identity import code_identity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--host', required=True)
    parser.add_argument('--port', required=True, type=int)
    parser.add_argument('--identity', required=True)
    parser.add_argument('--code-identity', required=True)
    parser.add_argument('--workers', required=True, type=int)
    parser.add_argument('--memory-bytes', required=True, type=int)
    parser.add_argument('--bundles')
    args = parser.parse_args()
    if args.code_identity != code_identity():
        raise SystemExit('Worker source or compiled engine differs from the coordinator')
    address = args.host, args.port
    key = bytes.fromhex(os.environ['MAPPA_POOL_AUTH'])
    queue = connect(address, key)
    with tempfile.TemporaryDirectory(prefix='mappa-node-', dir=os.environ.get('SLURM_TMPDIR')) as scratch:
        store = BundleStore(args.bundles, scratch=scratch) if args.bundles else connect_store(address, key)
        states = run(queue, workers=args.workers, identity=args.identity,
            scratch=scratch, store=store, memory_bytes=args.memory_bytes)
    if any(state != 'complete' for state in states):
        raise SystemExit(f'Incomplete pool: {states}')


if __name__ == '__main__':
    main()
