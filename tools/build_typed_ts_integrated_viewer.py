#!/usr/bin/env python3
"""Pack typed per-case TS viewers into one self-contained batch navigator."""
from __future__ import annotations

from mappa.presentation.viewers import embedded_collection_html

import argparse
import json
from pathlib import Path


def build(results_root, output):
    root = Path(results_root).resolve()
    output = Path(output).resolve()
    paths = sorted((root / "cases").glob("*/ts_scores.json"))

    def documents():
        for score_path in paths:
            document = json.loads(score_path.read_text())
            rankings = document.get('rankings') or ()
            ranked = sum(len(item['ranked_initial_guesses']) for item in rankings)
            best = max((row['score'] for item in rankings
                        for row in item['ranked_initial_guesses']), default=None)
            best_label = 'n/a' if best is None else f'{best:.5f}'
            atoms, case = int(document['atom_count']), document['case']
            yield dict(id=case, atoms=atoms, mechanisms=len(rankings), ranked=ranked, best=best,
                seconds=float(document.get('elapsed_seconds', 0.0)),
                meta=f'{atoms} atoms · {len(rankings)} mech · {ranked} ranked · best {best_label}',
                detail=f'{case} — {atoms} atoms, {len(rankings)} mechanisms, {ranked} ranked rows',
                document=(score_path.parent / 'view.html').read_bytes())

    html = embedded_collection_html(documents(), title='Typed TS analysis',
        summary=f'{len(paths)} cases · R/P alignment · ranked initial guesses · reactive-mode scores')
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html)
    return {"output": str(output), "cases": len(paths),
            "bytes": output.stat().st_size}


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    print(json.dumps(build(args.results_root, args.output), indent=2))


if __name__ == "__main__":
    main()
