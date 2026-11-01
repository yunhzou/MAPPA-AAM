"""Reaction cache snapshots. One coordinator publishes; workers only stage updates.

The first snapshot for each direction is retained in baseline.cache. latest.cache
is replaced atomically. Neither file contains search outputs or mutable locks.
"""
from dataclasses import asdict
import hashlib
import fcntl
import json
import os
from pathlib import Path
import pickle
import struct
import time
import zlib


def cache_identity(problem, config, code):
    def endpoint(value):
        return hashlib.sha256(pickle.dumps((tuple(value.elements), value.wbo.shape,
            value.wbo.astype('<f8').tobytes()), protocol=5)).hexdigest()
    source, target = endpoint(problem.reactant), endpoint(problem.product)
    pair = tuple(sorted((source, target)))
    key = hashlib.sha256(pickle.dumps((pair, asdict(config), code), protocol=5)).hexdigest()
    return key, str(int(source > target))


def _header(stream):
    prefix = stream.read(4)
    if len(prefix) != 4:
        raise ValueError('Truncated cache header')
    size, = struct.unpack('<I', prefix)
    if size > 65536:
        raise ValueError('Invalid cache header size')
    header = json.loads(stream.read(size))
    if header['schema'] != 'mappa.reaction-cache/1' or header['codec'] != 'zlib1':
        raise ValueError('Unsupported reaction cache schema')
    return header


def _section(stream, header, direction):
    for name, entry in header['sections'].items():
        if name != direction:
            stream.seek(entry['size'], 1)
            continue
        data = stream.read(entry['size'])
        if len(data) != entry['size'] or hashlib.sha256(data).hexdigest() != entry['sha256']:
            raise ValueError('Reaction cache checksum differs')
        return data
    return None


def read_section(path, key, direction):
    # Opening once keeps header and payload on the same inode during replacement.
    with open(path, 'rb') as stream:
        header = _header(stream)
        if header['key'] != key:
            raise ValueError('Reaction cache identity differs')
        return _section(stream, header, direction), header['version']


def _read_all(path, key, *, exclude=None):
    with open(path, 'rb') as stream:
        header = _header(stream)
        if header['key'] != key:
            raise ValueError('Reaction cache identity differs')
        result = {}
        for name, entry in header['sections'].items():
            if name == exclude:
                stream.seek(entry['size'], 1)
                continue
            data = stream.read(entry['size'])
            if len(data) != entry['size'] or hashlib.sha256(data).hexdigest() != entry['sha256']:
                raise ValueError('Reaction cache checksum differs')
            result[name] = data
        return result


def _write(path, key, version, sections):
    header = dict(schema='mappa.reaction-cache/1', codec='zlib1', key=key, version=version,
                  sections={name: dict(size=len(data), sha256=hashlib.sha256(data).hexdigest())
                            for name, data in sorted(sections.items())})
    encoded = json.dumps(header, separators=(',', ':')).encode()
    with open(path, 'xb') as stream:
        stream.write(struct.pack('<I', len(encoded)))
        stream.write(encoded)
        for name in sorted(sections):
            stream.write(sections[name])
        stream.flush()
        os.fsync(stream.fileno())
    return header


class ReactionCaches:
    """Coordinator-owned metadata; only one publisher per reaction at a time."""
    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.owner = (self.directory/'.owner.lock').open('a+b')
        try:
            fcntl.flock(self.owner, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            self.owner.close()
            raise RuntimeError('Reaction cache directory already has an active coordinator') from None
        for path in self.directory.glob('*.pending'):
            path.unlink()
        self.records, self.pending, self.writers = {}, {}, set()

    def close(self):
        for token in list(self.pending):
            self.finish(token)
        self.owner.close()

    def describe(self, identity):
        key, direction = identity
        if key not in self.records:
            record = dict(key=key, version=0, sections=[], baseline_sections=[], errors=[], baseline_error=False,
                baseline=str(self.directory/f'{key}.baseline.cache'),
                latest=str(self.directory/f'{key}.latest.cache'))
            for field in ('baseline', 'latest'):
                path = Path(record[field])
                if path.exists():
                    try:
                        with path.open('rb') as stream:
                            header = _header(stream)
                        if header['key'] != key:
                            raise ValueError('Reaction cache identity differs')
                    except (OSError, ValueError, KeyError) as error:
                        record['errors'].append(f'{field} header: {error}')
                        record['baseline_error'] |= field == 'baseline'
                        continue
                    record['version'] = header['version']
                    record['sections'] = list(header['sections'])
                    if field == 'baseline':
                        record['baseline_sections'] = list(header['sections'])
            if record['errors'] and not record['version']:
                record['version'] = -1  # Let the worker record an explicit failed load.
            self.records[key] = record
        return dict(self.records[key], direction=direction)

    def reserve(self, identity, version, token):
        record = self.describe(identity)
        if record['version'] != version or record['key'] in self.writers or record['baseline_error']:
            return None
        plan = dict(record, version=version+1, token=token,
                    staging=str(self.directory/f'{record["key"]}.{token}.pending'),
                    first_staging=str(self.directory/f'{record["key"]}.{token}.first.pending'))
        self.pending[token] = plan
        self.writers.add(record['key'])
        return plan

    def finish(self, token, publication=None):
        plan = self.pending.pop(token, None)
        if plan is None:
            return False
        self.writers.discard(plan['key'])
        try:
            if publication is None:
                return False
            if publication['version'] != plan['version']:
                raise ValueError('Stale cache publication')
            # Baseline entries are only added for previously unseen directions.
            if publication['baseline_added']:
                os.replace(plan['first_staging'], plan['baseline'])
                self.records[plan['key']]['baseline_sections'] = publication['sections']
            os.replace(plan['staging'], plan['latest'])
            self.records[plan['key']] = {key: value for key, value in plan.items()
                if key in ('key', 'version', 'baseline', 'latest', 'baseline_sections')}
            self.records[plan['key']]['sections'] = publication['sections']
            if publication['baseline_added']:
                self.records[plan['key']]['baseline_sections'] = publication['sections']
            self.records[plan['key']]['publisher'] = token
            self.records[plan['key']].update(errors=[], baseline_error=False)
            return True
        finally:
            for key in ('staging', 'first_staging'):
                Path(plan[key]).unlink(missing_ok=True)


class WorkerCache:
    """Task-boundary load/save around a worker's existing prepared session."""
    def __init__(self):
        self.clear()

    def clear(self):
        self.identity, self.version, self.repair = None, 0, None
        self.saved_revision = 0
        self.pending = None

    def begin(self, descriptor):
        self.descriptor = descriptor
        self.used = False
        self.metrics = dict(enabled=descriptor is not None, load_seconds=0., save_seconds=0., loaded_bytes=0,
                            saved_bytes=0, baseline_loads=0, errors=[])

    def use(self, repair):
        if repair is None:
            return
        self.used = True
        if self.descriptor is None:
            self.repair, self.before = repair, repair.stats()
            return
        current = self.descriptor
        self.metrics['errors'].extend(current.get('errors', ()))
        identity = current['key'], current['direction']
        if self.repair is not repair or self.identity != identity:
            self.identity, self.version, self.repair = identity, 0, repair
            self.saved_revision = repair.stats()['revision']
            self.pending = None
        if self.pending is not None and current.get('publisher') == self.pending[0]:
            self.version, self.saved_revision = self.pending[1:]
            self.pending = None
        if current['version'] != self.version:
            started = time.perf_counter()
            for field in ('latest', 'baseline'):
                try:
                    data, version = read_section(current[field], *identity)
                    if data is None and field == 'latest' and current['direction'] in current['baseline_sections']:
                        continue
                    if data is not None:
                        repair.loads(zlib.decompress(data))
                        self.metrics['loaded_bytes'] += len(data)
                        self.metrics['baseline_loads'] += int(field == 'baseline')
                    self.version = current['version'] if field == 'baseline' else version
                    self.saved_revision = repair.stats()['revision']
                    break
                except (OSError, ValueError, KeyError, zlib.error) as error:
                    self.metrics['errors'].append(f'{field}: {error}')
            else:
                raise RuntimeError('Neither reaction cache version could be loaded')
            self.metrics['load_seconds'] += time.perf_counter()-started
        self.before = repair.stats()

    def proposal(self):
        if not self.used:
            return None
        stats = self.repair.stats()
        self.metrics.update(fragment_calls=stats['calls']-self.before['calls'],
            fragment_hits=stats['full_hits']-self.before['full_hits'],
            retained_bytes=stats['resident_bytes'], entries=stats['entries'])
        if self.descriptor is not None:
            missing = self.descriptor['direction'] not in self.descriptor['sections'] and stats['entries'] > 0
            if missing or stats['revision'] != self.saved_revision or self.metrics['baseline_loads']:
                return dict(version=self.version)
        return None

    def save(self, plan):
        if plan is None:
            return None
        started = time.perf_counter()
        data = zlib.compress(self.repair.dumps(), 1)
        direction = plan['direction']
        added = direction not in plan['baseline_sections']
        baseline = _read_all(plan['baseline'], plan['key']) if added and plan['baseline_sections'] else {}
        try:
            sections = (_read_all(plan['latest'], plan['key'], exclude=direction)
                        if any(name != direction for name in plan['sections']) else {})
        except (OSError, ValueError, KeyError) as error:
            self.metrics['errors'].append(f'save latest: {error}')
            sections = _read_all(plan['baseline'], plan['key'], exclude=direction)
        if added:
            baseline[direction] = data
            _write(plan['first_staging'], plan['key'], plan['version'], baseline)
        sections[direction] = data
        _write(plan['staging'], plan['key'], plan['version'], sections)
        self.metrics['saved_bytes'] = sum(map(len, sections.values())) + (sum(map(len, baseline.values())) if added else 0)
        self.metrics['save_seconds'] = time.perf_counter()-started
        # The next task's publisher token acknowledges our snapshot without
        # reading back the same bytes this process just wrote.
        self.pending = plan['token'], plan['version'], self.repair.stats()['revision']
        return dict(version=plan['version'], baseline_added=added, sections=list(sections))
