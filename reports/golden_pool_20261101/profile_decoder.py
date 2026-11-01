"""Measure saved-catalogue decoding and pickle memory without changing either."""
import argparse
import gc
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import subprocess
import sys
import tempfile
import time

import psutil


def profile(bundle, source, digest_only=False, timings=False, compare_record=None):
    from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob
    from mappa.mappa_core.core.scheduler.results import DecodingResult
    from mappa.mappa_core.post_processing.decoder.api import EventDecodeSession
    from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig, EventInterpretation
    from mappa.mappa_core.post_processing.decoder import events
    from mappa.mappa_core.post_processing.decoder.fragments import FragmentInventory, FragmentDecoder
    from mappa.mappa_core.utils.io.records import decoded_events_record

    process = psutil.Process()
    started = time.monotonic()

    def sample(stage, **values):
        print(json.dumps(dict(stage=stage, seconds=time.monotonic()-started,
                              rss_bytes=process.memory_info().rss, **values)), flush=True)

    measured = {}

    def observe(owner, name, label):
        original = getattr(owner, name)
        measured[label] = dict(calls=0, seconds=0.)

        def timed(*args, **kwargs):
            start = time.perf_counter()
            try:
                return original(*args, **kwargs)
            finally:
                measured[label]['seconds'] += time.perf_counter()-start
                measured[label]['calls'] += 1

        setattr(owner, name, timed)

    if timings:
        for owner, name, label in (
            (FragmentInventory, '__init__', 'inventory'),
            (EventDecodeSession, 'decode', 'inventory_and_action_annotations'),
            (EventInterpretation, 'tables', 'family_table_assembly'),
            (FragmentDecoder, 'fragment', 'fragment_lookup_and_decode'),
            (EventInterpretation, 'evaluate', 'representative_descriptions'),
            (EventInterpretation, 'realize', 'representative_validation'),
            (events.SignedEventIndex, 'describe', 'whole_reaction_description'),
            (events, 'classify_bonds', 'whole_reaction_bond_classification'),
            (events.SignedEventIndex, '_describe_events', 'event_identity_and_description'),
            (events.SignedEventIndex, 'counts', 'whole_reaction_count_check'),
        ):
            observe(owner, name, label)

    original = json.loads(source.read_text())['original']
    saved = read_blob(Blob(str(bundle.resolve()), 0, bundle.stat().st_size,
                           original['sha256']), None)
    sample('catalogue_loaded', families=len(saved.catalogue.families),
           branches=saved.catalogue.branch_count, cuts=len(saved.scope))
    if timings:
        families = saved.catalogue.families
        keys = {(f.policy, f.mapping, f.required_edges, f.actions) for f in families}
        sample('catalogue_deduplication', saved_type=type(saved).__name__,
               input_paths=saved.catalogue.path_count, unique_family_keys=len(keys),
               duplicate_family_keys=len(families)-len(keys),
               unique_representative_mappings=len({f.mapping for f in families}))
        del families, keys
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
        if timings:
            sample('timing_breakdown', inclusive=True, functions=measured,
                   unique_fragments=len(interpreted.inventory.fragments),
                   fragment_occurrences=sum(len(partition) for partitions in interpreted.inventory.uses
                                            for partition in partitions),
                   decoded_fragment_contexts=len(interpreted.local.fragments),
                   local_bond_contexts=len(interpreted.local.bonds))
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
    unique_terms = {id(term): sys.getsizeof(term) for t in output.record['tables']
                    for kind in ('internal', 'boundary')
                    for region in t[kind] for term in region['terms']}
    sample('shared_record_shape', unique_term_records=len(unique_terms),
           unique_term_dictionary_bytes_excluding_values=sum(unique_terms.values()))
    del unique_terms

    def record_digest(record):
        digest = hashlib.sha256()
        for block in json.JSONEncoder(sort_keys=True, separators=(',', ':'), default=float).iterencode(record):
            digest.update(block.encode())
        return digest.hexdigest()

    if digest_only or compare_record:
        actual = record_digest(output.record)
        sample('record_digest', sha256=actual, schema=output.record['schema'])
        if compare_record:
            assert output.record['schema'] == 'mappa.decoded_events/v4'
            assert 'representatives' not in output.record
            del output
            gc.collect()
            # Trusted, single-blob published baseline. Load sequentially to bound RAM.
            with gzip.open(compare_record, 'rb') as stream:
                previous = pickle.load(stream).record
            assert previous['schema'] == 'mappa.decoded_events/v3'
            removed = len(previous.pop('representatives'))
            previous['schema'] = 'mappa.decoded_events/v4'
            expected = record_digest(previous)
            sample('baseline_comparison', baseline=str(compare_record),
                   removed_reports=removed, normalized_sha256=expected,
                   all_other_fields_identical=actual == expected)
            assert actual == expected
        return

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
    parser.add_argument('--digest-only', action='store_true')
    parser.add_argument('--timings', action='store_true', help='Measure nested decoder stages; includes instrumentation overhead')
    parser.add_argument('--compare-record', type=Path,
                        help='Compare all retained fields with a trusted single-blob v3 baseline; implies digest-only')
    args = parser.parse_args()
    if args.child:
        profile(args.bundle, args.source, args.digest_only, args.timings, args.compare_record)
        return
    command = [sys.executable, __file__, str(args.bundle), str(args.source), '--child']
    if args.digest_only:
        command.append('--digest-only')
    if args.timings:
        command.append('--timings')
    if args.compare_record:
        command.extend(['--compare-record', str(args.compare_record)])
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
