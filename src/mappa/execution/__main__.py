"""Launch the common worker locally or inside an existing Slurm allocation."""
import argparse
import json
import os
from pathlib import Path
import socket
import time

from . import TaskQueue, connect, coordinator, run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ('serve', 'worker'):
        command = sub.add_parser(name)
        command.add_argument('--identity', required=True, help='Frozen code/input manifest identity')
        command.add_argument('--auth-file', type=Path, required=True, help='Private shared key file, outside Git')
        command.add_argument('--host', required=True)
        command.add_argument('--port', type=int, required=True)
        if name == 'serve':
            command.add_argument('--database', type=Path, required=True, help='Coordinator-local queue database')
            command.add_argument('--backup', type=Path, required=True, help='Compact shared recovery checkpoint')
        else:
            command.add_argument('--workers', type=int, default=1)
            command.add_argument('--memory-mib', type=int, required=True, help='Per-worker process-tree guard')
            command.add_argument('--scratch', type=Path, required=True, help='Node-local directory')
            command.add_argument('--bundles', type=Path, required=True, help='Shared append-only result directory')
    args = parser.parse_args()
    key = args.auth_file.read_bytes()
    if len(key) < 32:
        parser.error('Use a randomly generated authentication key of at least 32 bytes')
    if args.command == 'worker':
        for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
            os.environ[name] = '1'
        queue = connect((args.host, args.port), key)
        print(json.dumps(run(queue, workers=args.workers, identity=args.identity,
            scratch=args.scratch, bundles=args.bundles, memory_bytes=args.memory_mib*2**20)))
    else:
        queue = TaskQueue(args.database, args.identity)
        temporary = args.backup.with_suffix('.pending')
        with coordinator(queue, (args.host, args.port), key):
            previous = None
            while True:
                state = queue.snapshot()
                if state != previous:
                    print(json.dumps(state), flush=True)
                    previous = state
                queue.backup(temporary)
                temporary.replace(args.backup)
                if not (state.get('pending', 0) or state.get('running', 0)):
                    break
                time.sleep(10)
        queue.close()


if __name__ == '__main__':
    main()
