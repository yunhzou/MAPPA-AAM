"""One saved EnzymeMap reaction through a real, persistent two-node Slurm pool."""
from dataclasses import asdict, is_dataclass
import gzip
import json
import os
from pathlib import Path
import socket
from statistics import median
import subprocess
import sys
import time

import numpy as np
from mappa import AAMSearchConfig, search_aam, deduplicate, decode_events, query_path
from mappa.execution import WorkPool, BundleStore, StageResults
from mappa.execution.identity import code_identity
from mappa.mappa_core.core.search.graph import AAMSearchGraph
from mappa.mappa_core.utils.io.records import problem_from_record, decoded_events_record
from test_work_pool import normalize_provenance


def read(path):
    with gzip.open(path, 'rt') as stream:
        return json.load(stream)


def plain(value):
    return json.loads(json.dumps(value))


def measured(call):
    start = time.perf_counter()
    value = call()
    return value, time.perf_counter()-start


def dense_counts(problem, mapping):
    assert set(mapping) == set(range(problem.source_atom_count))
    assert set(mapping.values()) == set(range(problem.target_atom_count))
    assert set(problem.reactant.elements) <= {'C', 'H', 'N', 'O', 'P'}
    order = [mapping[i] for i in range(problem.source_atom_count)]
    before = problem.reactant.wbo
    after = problem.product.wbo[np.ix_(order, order)]
    upper = np.triu(np.ones(before.shape, dtype=bool), 1)
    decrease, increase = before-after > .5, after-before > .5
    return {name: int(np.count_nonzero(upper & mask)) for name, mask in {
        'breaking': decrease & (after <= .5), 'weakening': decrease & (after > .5),
        'forming': increase & (before <= .5), 'strengthening': increase & (before > .5)}.items()}


def main():
    folder = Path(sys.argv[1]).resolve()
    problem = problem_from_record(read(folder/'input.json.gz'))
    reference = read(folder/'reference.json.gz')
    expected = read(folder/'expected_decoding.json.gz')
    config = AAMSearchConfig(seed_count=1, branch_limit=100, iso_tolerance=1., random_seed=42)
    report = dict(job=os.environ['SLURM_JOB_ID'], nodes=os.environ['SLURM_JOB_NODELIST'],
        partition=os.environ['SLURM_JOB_PARTITION'], coordinator=socket.gethostname(),
        source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        source_identity=code_identity(), config=asdict(config), direct=[], distributed=[])
    for trial in range(3):
        baseline, search = measured(lambda: search_aam(problem, config, execution='reused_native'))
        catalogue, reduction = measured(lambda: deduplicate(baseline))
        record, decoding = measured(lambda: decoded_events_record(decode_events(catalogue)))
        assert plain(record) == expected
        report['direct'].append(dict(search=search, deduplication=reduction, decoding=decoding))
    store = BundleStore(folder/'bundles')
    start = time.perf_counter()
    with WorkPool.slurm(workers=2, store=store, memory_bytes=2*2**30, case_seconds=120) as pool:
        report['pool_creation_seconds'] = time.perf_counter()-start
        report['worker_service_pid'] = pool._launch.pid
        for trial in range(3):
            matched, search = measured(lambda: pool.search({'446': problem}, config).wait_and_validate())
            catalogues, reduction = measured(lambda: pool.deduplicate(matched).wait_and_validate())
            decoded, decoding = measured(lambda: pool.decode(catalogues).wait_and_validate())
            cuts = [store.load(ref) for ref in matched.load('446').cuts]
            assert AAMSearchGraph.combine([cut.graph for cut in cuts]) == baseline.graph
            actual = dict(decoded.load('446').record)
            actual['catalogue'] = normalize_provenance(actual['catalogue'], cuts)
            assert plain(actual) == expected
            report['distributed'].append(dict(search=search, deduplication=reduction, decoding=decoding))
            print('TRIAL', trial, report['distributed'][-1], flush=True)
        catalogues.save_manifest(folder/'catalogues.json')
        decoded.save_manifest(folder/'decoded.json')
        # Reopen the catalogue manifest and decode it through the SAME Slurm services.
        reopened = StageResults.open(folder/'catalogues.json', store=store)
        again = pool.decode(reopened).wait_and_validate()
        assert again.load('446').record == decoded.load('446').record
        # Exercise a genuine seed-batch dependency boundary on the same reaction.
        multi = AAMSearchConfig(seed_count=3, branch_limit=100, iso_tolerance=1., random_seed=42)
        multi_baseline = search_aam(problem, multi, execution='reused_native')
        multi_matched = pool.search({'446-three-seeds': problem}, multi, seed_batch=1).wait_and_validate()
        multi_cuts = [store.load(ref) for ref in multi_matched.load('446-three-seeds').cuts]
        assert AAMSearchGraph.combine([cut.graph for cut in multi_cuts]) == multi_baseline.graph
        multi_catalogues = pool.deduplicate(multi_matched).wait_and_validate()
        multi_decoded = pool.decode(multi_catalogues).wait_and_validate()
        actual = dict(multi_decoded.load('446-three-seeds').record)
        actual['catalogue'] = normalize_provenance(actual['catalogue'], multi_cuts)
        assert plain(actual) == plain(decoded_events_record(decode_events(deduplicate(multi_baseline))))
        rows = pool.queue.results()
        assert all(row['state'] == 'complete' for row in rows.values())
        workers = {}
        for row in rows.values():
            result = row['result']
            workers.setdefault(result['worker'], set()).add(result['compute_pid'])
        assert all(len(pids) == 1 for pids in workers.values()), workers
        assert len({worker.split(':')[0] for worker in workers}) == 2, workers
        assert len(workers) == 4, workers
        report['persistent_workers'] = {worker: sorted(pids) for worker, pids in workers.items()}
        report['tasks'] = len(rows)
        report['states'] = pool.queue.snapshot()
        with gzip.open(folder/'task_metrics.json.gz', 'wt') as stream:
            json.dump(rows, stream, default=lambda value: asdict(value) if is_dataclass(value) else str(value))
        # Verification remains a separate operation after all decoding passes.
        catalogue = catalogues.load('446').catalogue
        events = decode_events(catalogue)
        checks = 0
        for fid in range(len(catalogue.families)):
            counts = dense_counts(problem, events.realize(fid))
            for partition in range(len(events.fragment_partitions[fid])):
                assert events.contributions(fid, partition=partition)['counts'] == counts
                checks += 1
        target = dict(reference['mapping'])
        for fid, family in enumerate(catalogue.families):
            status, evidence = query_path(family.as_path(problem), problem, target,
                source_atoms=tuple(target), complete_reference=True, timeout_ms=10000)
            if status == 'recovered':
                mapping = dict(evidence['mapping'])
                assert all(mapping[a] == b for a, b in target.items())
                assert events.tables(fid).evaluate(events.index.vector(mapping))['counts'] == dense_counts(problem, mapping)
                report['verification'] = dict(status=status, family=fid, heavy_atoms=len(target), partition_checks=checks)
                break
        else:
            raise AssertionError('Exact annotated reference was not recovered')
    report['pool_wall_seconds'] = time.perf_counter()-start
    report['bundle_files'] = len(list((folder/'bundles').glob('*.bundle')))
    assert report['bundle_files'] <= 5
    report['median_seconds'] = {mode: median(sum(row.values()) for row in report[mode])
                                for mode in ('direct', 'distributed')}
    report['full_records_equal'] = True
    report['multi_batch_seed_equivalence'] = True
    report['manifest_reopen_equal'] = True
    report['success'] = True
    (folder/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
