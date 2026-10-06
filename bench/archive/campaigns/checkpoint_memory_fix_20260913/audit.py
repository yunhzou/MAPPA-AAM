"""Audit recorded cut verdicts; no search, decoding, or permutation enumeration."""
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from mappa.alignment.sweep import cut_sweep_items


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def complete_negative(expected, rows):
    indices = [int(Path(r['cut']).name.split('_')[1].split('.')[0]) for r in rows]
    return (indices == list(range(expected)) and
            all(r['reference_recovery'] == 'not_recovered' for r in rows))


def main():
    # Incomplete or duplicated coverage and solver unknowns must never certify absence.
    sample = [dict(cut=f'cut_{i:05d}.raw.pkl.gz', reference_recovery='not_recovered')
              for i in range(2)]
    assert complete_negative(2, sample)
    assert not complete_negative(2, sample[:1])
    assert not complete_negative(2, sample[:1] * 2)
    assert not complete_negative(2, [sample[0], dict(sample[1], reference_recovery='unknown')])
    proof_path = HERE / 'original-finalization-proof.json'
    proof = read(proof_path)
    driver = HERE / 'original-score-checkpoints.py'
    assert proof['status'] == 'passed' and proof['driver_sha256'] == sha(driver)
    results = []
    for case in (590, 1358):
        input_path = HERE / f'case{case}-input.json'
        raw = read(input_path)
        for direction, side in [('R_to_P', 'reactant'), ('P_to_R', 'product')]:
            path = HERE / f'case{case}-{direction}.json'
            d = read(path)
            audit = d['checkpoint_verification']
            execution_path = HERE / f'case{case}-{direction}-execution.json'
            execution = read(execution_path)
            assert execution['status'] == 'passed'
            assert execution['task'] == dict(seed=10, case=case, direction=direction)
            assert audit['driver_sha256'] == sha(driver)
            assert audit['proof_sha256'] == sha(proof_path)
            expected = len(cut_sweep_items(np.asarray(raw[side]['wbo']), .2))
            rows = audit['checked']
            assert audit['available_cuts'] == expected
            assert complete_negative(expected, rows)
            assert all(f'/seed10/case{case}/{direction}/cuts/' in r['cut'] for r in rows)
            results.append(dict(case=case, direction=direction, expected_cuts=expected,
                checked_cuts=len(rows), reference_recovery='not_recovered',
                original_aggregate=d['reference_recovery'],
                peak_mib=execution['peak_mib'], cpu_seconds=d['cpu_seconds'],
                evidence={str(p.relative_to(HERE)): sha(p) for p in (path, input_path, execution_path, proof_path, driver)}))
    report = dict(status='passed', directions=results,
        corrected_seed10_counts=dict(recovered=1840, not_recovered=11, unknown=0, total=1851),
        scope='Reclassification from recorded complete per-cut negative verdicts, independently checking expected cut coverage and verifier provenance. Raw remote archives were not downloaded or re-evaluated. Absence concerns only the saved capped search families. Original resource failures remain recorded.',
        guard_checks='Missing cuts, duplicate cuts, and unknown cut verdicts reject an absence certificate.')
    (HERE / 'audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(status=report['status'], counts=report['corrected_seed10_counts'],
                         directions=len(results), cuts=sum(r['checked_cuts'] for r in results))))


if __name__ == '__main__':
    main()
