"""Recompute representation counts and competition gains from fresh archives.

Only compressed paths and existing complete event certificates are traversed.
No permutations, atom mappings, or event classes are newly enumerated.
"""
import os, sys, time, json, statistics
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
ROOT = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT/'optimized/src'), str(ROOT/'optimized/bench')]
from run_event_campaign import read, save, sha, MemoryGuard, bounded_process
OUT = ROOT/'coordinate-audit'

def child(case):
    if case == 0:
        import subprocess
        subprocess.run([sys.executable,str(ROOT/'score_checkpoints.py'),'test'],check=True)
    from mappa.artifacts import read_aam_checkpoint
    from mappa.final_branches import FinalBranchCatalogue
    from mappa.search_graph import frozen_value
    started = time.perf_counter()
    decoded = read(ROOT/'decoded-optimized'/f'case{case}'/'result.json')
    assert decoded['complete_saved_window']
    archives = [ROOT/'runs/coordinate/seed1'/f'case{case}'/'R_to_P/cuts/aam.pkl.gz',
                *sorted((ROOT/'competition'/f'case{case}').glob('takeover*.pkl.gz'))]
    identity = read(ROOT/'decoded-optimized'/f'case{case}'/'identity.json')
    assert identity['archives'] == {str(p): sha(p) for p in archives}
    interner, keys = {}, set()
    baseline_families = None
    catalogue = None
    for ordinal, archive in enumerate(archives):
        aam = read_aam_checkpoint(archive)
        if catalogue is None:
            catalogue = FinalBranchCatalogue(aam.problem)
        catalogue.add_graph(aam.graph, str(archive))
        if ordinal == 0:
            baseline_families = len(catalogue.families)
        edgekeys, local = {}, set()
        for path in aam.graph.paths():
            fragments = []
            for eid in path.transitions:
                edge = aam.graph.transitions[eid]
                if edge.match is None:
                    continue
                if eid not in edgekeys:
                    symmetry = {k: v for k, v in edge.match['symmetry'].items()
                                if k not in {'multiplicity', 'automorph_group_source'}}
                    fragment = (tuple(edge.match['fragment']), edge.preserved_bonds,
                                frozen_value(symmetry),
                                frozen_value(edge.match.get('deferred_edges', ())))
                    edgekeys[eid] = interner.setdefault(fragment, len(interner))
                fragments.append(edgekeys[eid])
            local.add((tuple(sorted(aam.graph.states[path.terminal].mapping)), tuple(fragments)))
        if ordinal == 0 and len(aam.graph.terminals) <= 10:
            assert len(local) == len(aam.graph.literal_branches())
        keys.update(local)
    assert catalogue.branch_count == decoded['branches']
    assert len(catalogue.families) == decoded['families']
    assert catalogue.path_count == decoded['paths']
    baseline_patterns, union, completed = set(), set(), set()
    with (ROOT/'decoded-optimized'/f'case{case}'/'families.jsonl').open() as journal:
        for line in journal:
            row = json.loads(line)
            patterns = {p['id'] for p in row['patterns']}
            union.update(patterns)
            if row['family'] < baseline_families:
                baseline_patterns.update(patterns)
            if row['complete']:
                completed.add(row['family'])
    assert completed == set(range(len(catalogue.families)))
    assert union == set(decoded['patterns'])
    representatives = read(ROOT/'competition'/f'case{case}'/'result.json')['patterns']
    new_patterns = union - baseline_patterns
    row = dict(case=case, old_literal_branches=len(keys), final_branches=catalogue.branch_count,
               flat_saved_families=len(catalogue.families), input_paths=catalogue.path_count,
               baseline_families=baseline_families, baseline_classes=len(baseline_patterns),
               baseline_class_ids=sorted(baseline_patterns),
               window=decoded['window'], classes=len(union), complete=True,
               new_classes=sorted(new_patterns),
               new_classes_absent_from_repair_representatives=sorted(new_patterns-set(representatives)),
               wall_seconds=time.perf_counter()-started)
    save(OUT/f'case{case}.json', row)

def main():
    OUT.mkdir(exist_ok=True)
    guard = MemoryGuard(3072, 8192, 6144)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1', NUMEXPR_NUM_THREADS='1', PYTHONHASHSEED='0')
    def run(case):
        with (OUT/f'case{case}.log').open('w') as log:
            status, peak = bounded_process([sys.executable, __file__, str(case)], env, log, 300, guard)
        return dict(case=case, status=status, peak_mib=peak)
    with ThreadPoolExecutor(max_workers=3) as pool:
        execution = list(pool.map(run, range(140)))
    save(OUT/'execution.json', execution)
    assert all(r['status'] == 'passed' for r in execution)
    rows = [read(OUT/f'case{c}.json') for c in range(140)]
    save(OUT/'summary.json', dict(cases=140, per_case=rows,
         old_branches=sum(r['old_literal_branches'] for r in rows),
         new_branches=sum(r['final_branches'] for r in rows),
         old_median=statistics.median(r['old_literal_branches'] for r in rows),
         new_median=statistics.median(r['final_branches'] for r in rows),
         new_window_class_count=sum(len(r['new_classes']) for r in rows),
         new_window_cases=[r['case'] for r in rows if r['new_classes']],
         scope='Fresh saved archives and full decoder journals; literal ordered grouping and unordered final grouping; baseline families are the initial prefix in the combined catalogue.'))

if __name__ == '__main__':
    if len(sys.argv)>1 and sys.argv[1] == 'wait':
        while not (ROOT/'comparator-campaigns-finished.json').exists():
            time.sleep(2)
        main()
    elif len(sys.argv)>1:
        child(int(sys.argv[1]))
    else:
        main()
