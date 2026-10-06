"""Time public decoding of trusted saved searches without rerunning AAM.

One serial decoder per process, four processes by default, a five-minute
external watchdog and the same memory guards as the end-to-end driver.
Reference answers are consulted only after decoding, including per-family
class support when a completed reference journal is available.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
import time

from final_end_to_end import REPO, read, save


def child(args):
    from mappa.artifacts import read_aam_checkpoint
    from mappa.final_branches import FinalBranchCatalogue
    import mappa.postprocessing as pp
    folder = args.output / f'case{args.child}'
    reference = next(row for row in read(REPO / 'reports/published_baseline_20260915/baseline-windows.json.gz')
                     if row['case'] == args.child)
    archive = args.archives / f'case{args.child}/cuts/aam.pkl.gz'
    aam = read_aam_checkpoint(archive)
    catalogue = FinalBranchCatalogue(aam.problem).add_aam(aam)
    original = pp.extract_path_events
    rows = []
    with (folder / 'families.jsonl').open('w', buffering=1) as journal:
        def measured(*positional, **kwargs):
            result = original(*positional, **kwargs)
            row = dict(family=len(rows), **{k:v for k,v in result.items() if k != 'patterns'},
                       pattern_ids=sorted(p['id'] for p in result['patterns']))
            rows.append(row)
            journal.write(json.dumps(row, separators=(',', ':')) + '\n')
            return result
        pp.extract_path_events = measured
        start, cpu = time.perf_counter(), time.process_time()
        result = pp.decode_events(catalogue, pp.EventDecodeConfig(
            max_events=reference['window'], seconds_per_family=300, max_patterns_per_family=None))
        timing = dict(wall_seconds=time.perf_counter()-start, cpu_seconds=time.process_time()-cpu)
    # No reference result is used by the decoder or any cache.
    ids = sorted(c.id for c in result.candidates)
    old_journal = args.archives / f'case{args.child}/families.jsonl'
    old_rows = [json.loads(line) for line in old_journal.read_text().splitlines()] if old_journal.exists() else []
    checked = [r for r in old_rows if r['complete']]
    same_support = all(rows[r['family']]['complete'] and
                       rows[r['family']]['pattern_ids'] == sorted(r['pattern_ids']) for r in checked)
    save(folder / 'result.json', dict(case=args.child, complete=result.complete, decode=timing,
         families=len(rows), pattern_ids=ids, same_class_ids=set(ids)==set(reference['patterns']),
         same_family_count=len(rows)==reference['baseline_families'],
         compared_family_supports=len(checked), same_compared_family_supports=same_support,
         archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest()))


def campaign(args):
    import psutil
    cases = list(range(140)) if args.cases == 'all' else [int(x) for x in args.cases.split(',')]
    workers = min(args.workers, 8, max(1, (os.cpu_count() or 2)-1))
    sources = list((REPO/'src/mappa').glob('*.py')) + [Path(__file__).resolve()]
    manifest = dict(cases=cases, workers=workers, watchdog_seconds=300,
        source_sha256={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        scope='Public decoder only, saved final search checkpoints; catalogue/read excluded. '
              'Published event windows, .5/.3 thresholds, no class cap; one numerical thread/process.')
    if (args.output/'manifest.json').exists() and read(args.output/'manifest.json') != manifest:
        raise ValueError('Output directory belongs to different sources/settings')
    save(args.output/'manifest.json', manifest)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1',
               NUMEXPR_NUM_THREADS='1', PYTHONHASHSEED='0')
    lock, active = threading.Lock(), {}
    def run(case):
        folder = args.output / f'case{case}'; folder.mkdir(parents=True, exist_ok=True)
        started, peak, cpu, stopped = time.perf_counter(), 0., 0., None
        with (folder/'run.log').open('w') as log:
            proc = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), '--archives',str(args.archives),
                '--output',str(args.output),'--child',str(case)],env=env,stdout=log,stderr=subprocess.STDOUT)
            q = psutil.Process(proc.pid)
            try:
                while proc.poll() is None:
                    try:
                        rss=q.memory_info().rss/2**20; usage=q.cpu_times(); cpu=usage.user+usage.system
                    except psutil.NoSuchProcess:
                        break
                    peak=max(peak,rss)
                    with lock:
                        active[proc.pid]=rss; total=sum(active.values()); largest=max(active,key=active.get)
                    if time.perf_counter()-started > 300:stopped='watchdog'
                    elif rss>4096 or (total>8192 and largest==proc.pid):stopped='memory_guard'
                    elif psutil.virtual_memory().available<6144*2**20:stopped='system_memory_pressure'
                    if stopped:proc.kill(); break
                    time.sleep(.2)
                code=proc.wait()
            finally:
                if proc.poll() is None:proc.kill(); proc.wait()
                with lock:active.pop(proc.pid,None)
        row=dict(case=case,status=stopped or ('finished' if code==0 else 'error'),
                 peak_mib=peak,sampled_process_cpu_seconds=cpu)
        save(folder/'execution.json',row);print(json.dumps(row),flush=True)
        return row
    start=time.perf_counter()
    with ThreadPoolExecutor(max_workers=workers) as pool:rows=list(pool.map(run,cases))
    save(args.output/'campaign.json',dict(rows=rows,wall_seconds=time.perf_counter()-start))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archives',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--cases',default='all')
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--child',type=int)
    args=parser.parse_args()
    child(args) if args.child is not None else campaign(args)
