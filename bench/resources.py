"""File and process accounting shared by retained benchmark drivers."""
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import threading
import time

def read(path):
    with gzip.open(path, 'rt') if path.suffix == '.gz' else path.open() as stream:
        return json.load(stream)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)

class MemoryGuard:
    """Sample only this runner's child PIDs; terminate oversized work explicitly.

    RSS/available-memory checks are polling safeguards, not kernel hard limits.
    A killed case remains a resource-limited result, never a scientific miss.
    """
    def __init__(self, worker_mib, total_mib, reserve_mib):
        self.worker_kib = worker_mib * 1024
        self.total_kib = total_mib * 1024
        self.reserve_kib = reserve_mib * 1024
        self.lock = threading.Lock()
        self.active, self.stops, self.peaks = {}, {}, {}
        self.last_check = 0.
        self.peak_total_kib = 0

    def add(self, process):
        with self.lock:
            self.active[process.pid] = process

    def remove(self, pid):
        with self.lock:
            self.active.pop(pid, None)
            return self.stops.pop(pid, None), self.peaks.pop(pid, 0)

    def check(self):
        with self.lock:
            now = time.monotonic()
            if now - self.last_check < .5:
                return
            self.last_check = now
            live = {pid: p for pid, p in self.active.items() if p.poll() is None}
            if not live:
                return
            result = subprocess.run(['/bin/ps', '-o', 'pid=,rss=', '-p', ','.join(map(str, live))],
                                    text=True, capture_output=True, timeout=3)
            if result.returncode and result.stderr.strip():
                raise RuntimeError('Cannot inspect benchmark child memory: ' + result.stderr.strip())
            rss = {int(parts[0]): int(parts[1]) for line in result.stdout.splitlines()
                   if len(parts := line.split()) == 2 and int(parts[0]) in live}
            for pid, value in rss.items():
                self.peaks[pid] = max(self.peaks.get(pid, 0), value)
            total = sum(rss.values())
            self.peak_total_kib = max(self.peak_total_kib, total)
            available = None
            if sys.platform == 'darwin':
                vm = subprocess.run(['/usr/bin/vm_stat'], text=True, capture_output=True, timeout=3)
                if vm.returncode:
                    raise RuntimeError('Cannot inspect available-memory estimate')
                import re
                page_size = int(re.search(r'page size of (\d+) bytes', vm.stdout).group(1))
                fields = {key.strip(): int(value.strip().rstrip('.')) for line in vm.stdout.splitlines()[1:]
                          if ':' in line for key, value in [line.split(':', 1)]}
                # Matches the free + inactive estimate; compressed memory is not counted as free.
                available = sum(fields.get(key, 0) for key in ('Pages free', 'Pages inactive')) * page_size / 1024
            elif Path('/proc/meminfo').exists():
                available = next(float(line.split()[1]) for line in Path('/proc/meminfo').read_text().splitlines()
                                 if line.startswith('MemAvailable:'))
            victims = {pid: 'worker_memory_limit' for pid, value in rss.items() if value > self.worker_kib}
            remaining = total - sum(rss[pid] for pid in victims)
            for pid in sorted(rss, key=rss.get, reverse=True):
                if remaining <= self.total_kib:
                    break
                if pid not in victims:
                    victims[pid] = 'campaign_memory_limit'
                    remaining -= rss[pid]
            if available is not None and available < self.reserve_kib:
                victims.update({pid: 'system_memory_pressure' for pid in rss})
            for pid, reason in victims.items():
                if live[pid].poll() is None:
                    try:
                        live[pid].kill()
                        self.stops[pid] = reason
                    except ProcessLookupError:
                        pass

def bounded_process(command, env, log, timeout, guard):
    process = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT)
    guard.add(process)
    deadline = time.monotonic() + timeout
    timed_out = False
    try:
        while process.poll() is None:
            guard.check()
            if time.monotonic() >= deadline:
                process.kill()
                timed_out = True
                break
            time.sleep(.1)
        process.wait()
    except BaseException:
        if process.poll() is None:
            process.kill()
        process.wait()
        guard.remove(process.pid)
        raise
    stopped, peak = guard.remove(process.pid)
    state = stopped or ('timeout' if timed_out else 'passed' if process.returncode == 0 else 'failed')
    return state, peak / 1024
