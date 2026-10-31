"""Partition a stopped campaign into disjoint batches, reusing its saved results."""
import argparse
import hashlib
import json
from pathlib import Path
import sqlite3


def prepare(index_path, previous, output, count=40, additional_ledgers=()):
    index = json.loads(Path(index_path).read_text())
    previous = Path(previous)
    plan = json.loads((previous/'plan.json').read_text())
    assert plan['input_index_sha256'] == hashlib.sha256(Path(index_path).read_bytes()).hexdigest()
    with sqlite3.connect(f'file:{previous}/ledger.sqlite?mode=ro', uri=True) as source:
        records = source.execute('SELECT * FROM results').fetchall()
        schema = source.execute("SELECT sql FROM sqlite_master WHERE name='results'").fetchone()[0]
    merged = {(row[0], row[1]): row for row in records}
    for ledger in additional_ledgers:
        with sqlite3.connect(f'file:{ledger}?mode=ro', uri=True) as source:
            for row in source.execute('SELECT * FROM results'):
                key = row[0], row[1]
                if key in merged and merged[key] != row:
                    raise ValueError(f'Conflicting checkpoint results for {key}')
                merged[key] = row
    records = list(merged.values())
    for batch in range(count):
        root = Path(output)/str(batch)
        root.mkdir(parents=True, exist_ok=False)
        rows = [row for row in index['directions'] if row['case'] % count == batch]
        # A batch manifest contains only its own scientific input rows.
        subset = dict(config=index['config'], directions=rows)
        path = root/'index.json'
        path.write_text(json.dumps(subset))
        result = root/'result'; result.mkdir()
        names = {f"{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}" for row in rows}
        updated = dict(plan, input_index=str(path.resolve()),
                       input_index_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                       cases=len({row['case'] for row in rows}), directions=len(rows))
        (result/'plan.json').write_text(json.dumps(updated,indent=2)+'\n')
        with sqlite3.connect(result/'ledger.sqlite') as destination:
            destination.execute(schema)
            destination.executemany('INSERT INTO results VALUES (?,?,?,?,?,?)',
                                    [record for record in records if record[0] in names])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('index'); parser.add_argument('previous'); parser.add_argument('output')
    parser.add_argument('--count',type=int,default=40)
    parser.add_argument('--additional-ledger',action='append',default=[])
    args = parser.parse_args()
    prepare(args.index, args.previous, args.output,count=args.count,
            additional_ledgers=args.additional_ledger)
