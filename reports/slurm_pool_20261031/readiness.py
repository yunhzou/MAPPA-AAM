"""Same 116-cut case, before and after all four workers report ready."""
from dataclasses import asdict
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from mappa import AAMSearchConfig, search_aam
from mappa.execution import WorkPool, BundleStore
from mappa.mappa_core.core.search.graph import AAMSearchGraph
from mappa.mappa_core.utils.io.records import problem_from_record
from run_case import read


def main():
    folder = Path(sys.argv[1]).resolve()
    folder.mkdir(parents=True, exist_ok=True)
    problem = problem_from_record(read(folder.parent/'slurm_pool_case446_20261031/input.json.gz'))
    config = AAMSearchConfig(seed_count=1, branch_limit=100, iso_tolerance=1., random_seed=42)
    baseline = search_aam(problem, config, execution='reused_native')
    store = BundleStore(folder/'bundles')
    report = dict(job=os.environ['SLURM_JOB_ID'], source=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(), passes=[])
    with WorkPool.slurm(workers=2, store=store, memory_bytes=2*2**30, case_seconds=120) as pool:
        for name in ('immediate', 'all-ready'):
            if name == 'all-ready':
                deadline = time.monotonic()+60
                while pool.queue.db.execute('SELECT count(*) FROM workers').fetchone()[0] < 4:
                    pool.check_workers()
                    if time.monotonic() > deadline:
                        raise TimeoutError('Not all four workers became ready')
                    time.sleep(.02)
            submitted = time.time()
            matched = pool.search({name:problem}, config).wait_and_validate()
            completed = time.time()
            cuts = [store.load(ref) for ref in matched.load(name).cuts]
            assert AAMSearchGraph.combine([cut.graph for cut in cuts]) == baseline.graph
            keys = [key for group in matched.tasks.values() for key in group]
            rows = pool.queue.results(keys)
            attempts = [pool.queue.db.execute('SELECT worker,started,finished FROM attempts WHERE task=?',(key,)).fetchone() for key in keys]
            report['passes'].append(dict(name=name, submitted=submitted, completed=completed,
                tasks=len(keys), attempts=attempts, graph_equal=True,
                cpu_seconds=sum(row['result']['cpu_seconds'] for row in rows.values()),
                workers=sorted({row[0] for row in attempts})))
            print(name, 'complete', completed-submitted, 'workers', report['passes'][-1]['workers'], flush=True)
        report['readiness'] = [dict(zip(('worker','pid','ready','pending','runnable'),row)) for row in
                              pool.queue.db.execute('SELECT * FROM workers ORDER BY ready')]
        report['states'] = pool.queue.snapshot()
    (folder/'readiness.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='passes'},indent=2),flush=True)


if __name__ == '__main__':
    main()
