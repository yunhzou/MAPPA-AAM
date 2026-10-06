"""Fresh one-seed coordinate search and public event decoding, with stage timing.

No archived search or decoder answers are reused. Published windows define the
output scope; saved class IDs are checked only after measurement. Four serial
workers by default, one numerical thread each, and an external five-minute
watchdog per reaction. An interrupted case is never a complete timing result.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import threading
import time

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'src'))


def read(path):
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_text())


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clock():
    return time.perf_counter(), time.process_time()


def elapsed(start):
    wall, cpu = clock()
    return dict(wall_seconds=wall-start[0], cpu_seconds=cpu-start[1])


def child(args):
    startup = clock()
    import numpy as np
    import random
    from mappa import AAMProblem, MolecularEndpoint, AAMSearchConfig, search_aam
    from mappa.final_branches import FinalBranchCatalogue
    import mappa.postprocessing as pp

    folder = args.output / f'case{args.child}'
    expected = next(row for row in read(REPO / 'reports/published_baseline_20260915/baseline-windows.json.gz')
                    if row['case'] == args.child)
    config = AAMSearchConfig(seed_count=1, branch_limit=2000, random_seed=42)
    decode_config = pp.EventDecodeConfig(max_events=expected['window'], seconds_per_family=300,
                                        max_patterns_per_family=None)
    state = dict(case=args.child, stage='input', complete=False, config=asdict(config),
                 decode_config=asdict(decode_config), competition=False, stages={},
                 setup=elapsed(startup), completed_families=0, processed_families=0)
    start = clock()
    def checkpoint():
        save(folder / 'progress.json', dict(state, elapsed=elapsed(start)))

    t = clock()
    raw = read(args.inputs / str(args.child) / 'input.json')
    problem = AAMProblem(*(MolecularEndpoint(**raw[s]) for s in ('reactant', 'product')),
                         raw.get('name', ''))
    state['stages']['input'] = elapsed(t)
    state['atoms'] = problem.source_atom_count
    random.seed(42)
    np.random.seed(42)
    state['stage'] = 'search'
    checkpoint()
    t = clock()
    aam = search_aam(problem, config, execution='reused_native', workers=1,
                     intermediate_dir=folder / 'cuts', archive_format='checkpoint', resume=False)
    state['stages']['search'] = elapsed(t)
    state['search_metrics'] = asdict(aam.metrics)
    state['search_capped'] = aam.graph.capped
    state['stage'] = 'catalogue'
    checkpoint()
    t = clock()
    catalogue = FinalBranchCatalogue(problem).add_aam(aam)
    state['stages']['catalogue'] = elapsed(t)
    state.update(paths=catalogue.path_count, branches=catalogue.branch_count,
                 families=len(catalogue.families), incomplete_paths=catalogue.incomplete_paths)
    state['stage'] = 'decode'
    checkpoint()
    original = pp.extract_path_events
    observed = {}
    last_checkpoint = time.perf_counter()
    with (folder / 'families.jsonl').open('w', buffering=1) as journal:
        def measured(*positional, **kwargs):
            nonlocal last_checkpoint
            result = original(*positional, **kwargs)
            for p in result['patterns']:
                observed.setdefault(p['id'], p['total'])
            journal.write(json.dumps(dict(family=state['processed_families'],
                **{k:v for k,v in result.items() if k != 'patterns'},
                pattern_ids=[p['id'] for p in result['patterns']]), separators=(',', ':')) + '\n')
            state['processed_families'] += 1
            state['completed_families'] += int(result['complete'])
            if time.perf_counter() - last_checkpoint > 3:
                state['observed_classes'] = observed
                checkpoint()
                last_checkpoint = time.perf_counter()
            return result
        pp.extract_path_events = measured
        t = clock()
        decoded = pp.decode_events(catalogue, decode_config)
        state['stages']['decode'] = elapsed(t)
        pp.extract_path_events = original
    state['stage'] = 'serialize'
    checkpoint()
    t = clock()
    save(folder / 'candidates.json', dict(complete=decoded.complete,
        candidates=[asdict(c) for c in decoded.candidates], family_reports=decoded.family_reports))
    state['stages']['serialize'] = elapsed(t)
    state.update(stage='done', complete=decoded.complete, end_to_end=elapsed(start),
                 observed_classes={c.id:c.total for c in decoded.candidates})
    # Validation is outside the measured pipeline; no reference affects search or decode.
    state['validation'] = dict(same_class_ids=set(state['observed_classes']) == set(expected['patterns']),
        same_family_count=len(catalogue.families) == expected['baseline_families'],
        expected_classes=len(expected['patterns']), expected_families=expected['baseline_families'])
    save(folder / 'result.json', state)
    checkpoint()
    print(json.dumps(dict(case=args.child, complete=decoded.complete,
                          timing=state['end_to_end'], validation=state['validation'])), flush=True)


def campaign(args):
    import psutil
    workers = min(args.workers, 8, max(1, (os.cpu_count() or 2)-1))
    cases = list(range(140)) if args.cases == 'all' else list(map(int, args.cases.split(',')))
    args.output.mkdir(parents=True, exist_ok=True)
    identity = dict(commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip(),
        driver_sha256=sha(Path(__file__)), inputs={str(c):sha(args.inputs / str(c) / 'input.json') for c in range(140)},
        sources={str(p.relative_to(REPO)):sha(p) for base in ('src/mappa','native/src')
                 for p in sorted((REPO/base).rglob('*')) if p.suffix in ('.py','.cpp','.h','.so')},
        windows_sha256=sha(REPO/'reports/published_baseline_20260915/baseline-windows.json.gz'),
        workers=workers, numerical_threads=1, watchdog_seconds=300,
        memory_limits_mib=dict(worker=4096, total=8192, reserve=6144), python=sys.version,
        cpu=subprocess.check_output(['sysctl','-n','machdep.cpu.brand_string'],text=True).strip()
            if sys.platform=='darwin' else platform.processor(), platform=platform.platform(),
        scope='Fresh coordinate R-to-P, one-seed sweep, cap2000, competition off; same published event windows; public decode_events with no class cap and 300s/family soft budget inside 300s/reaction watchdog. Includes input loading, search checkpoint I/O, catalogue construction, decoding, result serialization and light progress journaling. Excludes interpreter/import setup and post-run reference-set comparison.')
    manifest = args.output / 'manifest.json'
    if manifest.exists():
        if read(manifest) != identity:
            raise ValueError('Run identity changed; choose a new output directory')
    else:
        save(manifest, identity)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1', NUMEXPR_NUM_THREADS='1', PYTHONHASHSEED='0',
               PYTHONDONTWRITEBYTECODE='1', MAPPA_NATIVE='1')
    lock = threading.Lock()
    active = {}
    def run(case):
        folder = args.output / f'case{case}'
        if (folder/'execution.json').exists():
            return read(folder/'execution.json')
        folder.mkdir(parents=True, exist_ok=True)
        started = time.perf_counter()
        peak = cpu = 0.
        stop = None
        with (folder/'run.log').open('w') as log:
            proc = subprocess.Popen([sys.executable, str(Path(__file__)), '--inputs',str(args.inputs),
                '--output',str(args.output),'--child',str(case)],env=env,stdout=log,stderr=subprocess.STDOUT)
            q = psutil.Process(proc.pid)
            with lock:
                active[proc.pid] = 0.
            try:
                while proc.poll() is None:
                    try:
                        rss=q.memory_info().rss/2**20
                        usage=q.cpu_times()
                        cpu=usage.user+usage.system
                    except psutil.NoSuchProcess:
                        break
                    peak=max(peak,rss)
                    with lock:
                        active[proc.pid]=rss
                        total=sum(active.values())
                        largest=max(active,key=active.get)
                    if time.perf_counter()-started > 300:
                        stop='watchdog'
                    elif rss>4096 or (total>8192 and largest==proc.pid):
                        stop='memory_guard'
                    elif psutil.virtual_memory().available < 6144*2**20:
                        stop='system_memory_pressure'
                    if stop:
                        proc.kill()
                        break
                    time.sleep(.2)
                code=proc.wait()
            finally:
                if proc.poll() is None:
                    proc.kill()
                    proc.wait()
                with lock:
                    active.pop(proc.pid,None)
        row=dict(case=case,status=stop or ('finished' if code==0 else 'error'),
                 process_wall_seconds=time.perf_counter()-started, sampled_process_cpu_seconds=cpu,
                 peak_mib=peak, exit_code=code)
        if (folder/'result.json').exists():
            result=read(folder/'result.json')
            row.update(complete=result['complete'], end_to_end=result['end_to_end'],
                       validation=result['validation'])
        else:
            row['complete']=False
        save(folder/'execution.json',row)
        print(json.dumps(row),flush=True)
        return row
    started=time.perf_counter()
    rows=[]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for future in as_completed([pool.submit(run,c) for c in cases]):
            rows.append(future.result())
            save(args.output/'campaign.json',dict(requested_cases=cases,workers=workers,rows=rows,
                 wall_seconds=time.perf_counter()-started))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--cases',default='all')
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--child',type=int)
    args=parser.parse_args()
    child(args) if args.child is not None else campaign(args)
