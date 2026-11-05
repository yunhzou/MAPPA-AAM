#!/usr/bin/env python3
"""Build one offline HTML navigator containing every batch case viewer."""
from __future__ import annotations

from mappa.presentation.viewers import embedded_collection_html

import argparse
import hashlib
import json
from pathlib import Path


def _read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def _rounded(value, digits=6):
    if isinstance(value, float):
        return round(value, digits)
    return value


def _compact_view_document(document: bytes) -> bytes:
    """Remove non-rendered diagnostics and compact display-only coordinates.

    The standalone batch navigator only needs the chosen final symmetry blocks,
    not the branch witnesses used to derive them.  Interpolation coordinates are
    display data; six decimal places retain sub-micro-Angstrom precision while
    avoiding large JSON representations of binary floating-point values.
    """
    text = document.decode("utf-8")
    marker = "const DATA = "
    start = text.find(marker)
    if start < 0:
        raise ValueError("case viewer has no DATA declaration")
    data_start = start + len(marker)
    data, consumed = json.JSONDecoder().raw_decode(text[data_start:])
    data_end = data_start + consumed

    for mechanism in data.get("mechanisms", []):
        symmetry = mechanism.get("branch_symmetry")
        if isinstance(symmetry, dict):
            # These are the only symmetry fields consulted by rDegMap() and
            # mappedProductDegMap() in the case viewer.
            mechanism["branch_symmetry"] = {
                key: symmetry[key]
                for key in ("blocks", "color_groups")
                if key in symmetry
            }

        interpolation = mechanism.get("endpoint_interpolation")
        if not isinstance(interpolation, dict):
            continue
        for key in ("coordination_constraints", "primitive_counts",
                    "schema_version"):
            interpolation.pop(key, None)
        for frame in interpolation.get("frames", []):
            frame["coords"] = [
                [_rounded(component) for component in point]
                for point in frame.get("coords", [])
            ]
            residuals = frame.get("constraint_residuals")
            if isinstance(residuals, dict):
                frame["constraint_residuals"] = {
                    key: _rounded(value)
                    for key, value in residuals.items()
                }
            clashes = frame.get("clashes")
            if isinstance(clashes, dict):
                clashes["minimum_radius_ratio"] = _rounded(
                    clashes.get("minimum_radius_ratio"))

    compact_data = json.dumps(data, separators=(",", ":"))
    return (text[:data_start] + compact_data + text[data_end:]).encode("utf-8")


def _case_records(batch_root: Path):
    records = []
    for tier in ("small", "medium", "large"):
        manifest = _read_json(batch_root / "manifests" / f"{tier}.json")
        for case in manifest["cases"]:
            records.append((int(case["source_index"]), tier, case))
    return sorted(records)


def build(batch_root: Path, output: Path, compact: bool = False):
    batch_root = batch_root.resolve()
    output = output.resolve()
    inventory = _read_json(batch_root / "inventory.json")
    batch_summary = _read_json(batch_root / "batch_summary.json")
    cases = []

    def documents():
        for _source_index, tier, case in _case_records(batch_root):
            step = str(case["step_id"])
            summary = _read_json(batch_root / "cases" / step / "summary.json")
            rp = _read_json(batch_root / "cases" / step / "rp_stage.json")
            document = (batch_root / "views" / step / "view.html").read_bytes()
            mechanisms = rp.get("mechanisms") or []
            atoms = int(case['atom_count'])
            seconds = round(float(summary.get('elapsed_seconds', 0.0)), 3)
            violations = sum(int((mechanism.get('index_chirality') or {}).get(
                'selected_index_chirality_violation_count', 0)) for mechanism in mechanisms)
            record = dict(id=step, tier=tier, atoms=atoms, mechanisms=len(mechanisms),
                seconds=seconds, violations=violations, sha256=hashlib.sha256(document).hexdigest(),
                badge=tier, meta=f'{atoms} atoms · {len(mechanisms)} mechanisms · {seconds:.1f} s',
                detail=f'{step} — {atoms} atoms, {len(mechanisms)} mechanisms, violations={violations}')
            cases.append(record)
            yield record | dict(document=_compact_view_document(document) if compact else document)

    html = embedded_collection_html(documents(),
        title=f"AAM integrated viewer — {inventory['case_count']} cases",
        summary=(f"{inventory['case_count']} cases · {batch_summary['error_count']} errors · "
                 f"{batch_summary['total_elapsed_case_seconds']:.1f} case-seconds"))
    total_mechanisms = sum(case['mechanisms'] for case in cases)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html, encoding="utf-8")
    return {
        "output": str(output),
        "cases": len(cases),
        "mechanisms": total_mechanisms,
        "bytes": output.stat().st_size,
        "sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--compact", action="store_true",
        help="omit non-rendered diagnostics and compact display coordinates",
    )
    args = parser.parse_args()
    print(json.dumps(build(
        args.batch_root, args.output, compact=args.compact), indent=2))


if __name__ == "__main__":
    main()
