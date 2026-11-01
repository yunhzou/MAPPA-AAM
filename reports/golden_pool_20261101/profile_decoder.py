"""Measure saved-catalogue decoding and pickle memory without changing either."""
import argparse
import gzip
import json
import os
from pathlib import Path
import pickle
import subprocess
import sys
import tempfile
import time

import psutil


def profile(bundle, source):
    from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob
    from mappa.mappa_core.core.scheduler.results import DecodingResult
    from mappa.mappa_core.post_processing.decoder.api import EventDecodeSession
    from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig, EventInterpretation
    from mappa.mappa_core.utils.io.records import decoded_events_record

    process = psutil.Process()
    started = time.monotonic()

    def sample(stage, **values):
        print(json.dumps(dict(stage=stage, seconds=time.monotonic()-started,
                              rss_bytes=process.memory_info().rss, **values)), flush=True)

    original = json.loads(source.read_text())['original']
    saved = read_blob(Blob(str(bundle.resolve()), 0, bundle.stat().st_size,
                           original['sha256']), None)
    sample('catalogue_loaded', families=len(saved.catalogue.families),
           branches=saved.catalogue.branch_count, cuts=len(saved.scope))
    config = EventDecodeConfig()
    session = EventDecodeSession(saved.catalogue.problem, config)
    sample('session_ready')
    original_tables = EventInterpretation.tables

    def observed_tables(interpretation, family_id, partition=0):
        table = original_tables(interpretation, family_id, partition)
        if (family_id+1) % 1000 == 0 and partition == 0:
            sample('building_record', family=family_id+1, tables=len(interpretation._tables))
        return table

    EventInterpretation.tables = observed_tables

    def decode():
        interpreted = session.decode(saved.catalogue)
        sample('actions_ready')
        record = decoded_events_record(interpreted)
        sample('record_built_with_live_tables', tables=len(record['tables']),
               responses=len(record['responses']))
        return DecodingResult(saved.case, config, saved.scope, saved.stops,
                              record, saved.search_config, saved.execution)

    output = decode()
    del saved
    sample('operation_returned')
    terms = sum(t['statistics']['internal_terms'] + t['statistics']['boundary_terms']
                for t in output.record['tables'])
    term_dict_bytes = sum(sys.getsizeof(term) for t in output.record['tables']
                          for kind in ('internal', 'boundary')
                          for region in t[kind] for term in region['terms'])
    domain_entries = sum(sum(map(len, t['atom_domains'].values())) for t in output.record['tables'])
    sample('record_shape', terms=terms, term_dictionary_bytes_excluding_values=term_dict_bytes,
           domain_entries=domain_entries,
           projected_pair_states=sum(t['statistics']['projected_pair_states'] for t in output.record['tables']))

    class ObservedStream:
        def __init__(self, stream):
            self.stream, self.last_size = stream, 0

        def write(self, block):
            size = serializer.__sizeof__()
            if size != self.last_size:
                sample('pickle_allocation', pickler_bytes=size)
                self.last_size = size
            return self.stream.write(block)

    # Same protocol, gzip level and bounded spool as BundleWriter.put().
    with tempfile.SpooledTemporaryFile(max_size=8*2**20) as spool:
        with gzip.GzipFile(fileobj=spool, mode='wb', mtime=0, compresslevel=1) as compressed:
            serializer = pickle.Pickler(ObservedStream(compressed), protocol=5)
            serializer.dump(output)
            sample('pickle_complete', pickler_bytes=serializer.__sizeof__())
            del serializer
        sample('compressed_output', bytes=spool.tell())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bundle', type=Path)
    parser.add_argument('source', type=Path)
    parser.add_argument('--child', action='store_true')
    parser.add_argument('--memory-gib', type=float, default=2.)
    args = parser.parse_args()
    if args.child:
        profile(args.bundle, args.source)
        return
    command = [sys.executable, __file__, str(args.bundle), str(args.source), '--child']
    child = subprocess.Popen(command, env={**os.environ, 'OMP_NUM_THREADS': '1',
                              'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'})
    process = psutil.Process(child.pid)
    peak = 0
    while child.poll() is None:
        try:
            peak = max(peak, process.memory_info().rss)
        except psutil.NoSuchProcess:
            break
        if peak > args.memory_gib*2**30:
            child.kill()
            break
        time.sleep(.02)
    code = child.wait()
    print(json.dumps(dict(stage='guard', memory_bytes=int(args.memory_gib*2**30),
                          peak_rss_bytes=peak, exit_code=code)), flush=True)
    raise SystemExit(code)


if __name__ == '__main__':
    main()
