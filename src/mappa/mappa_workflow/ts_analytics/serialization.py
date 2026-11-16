"""TS-specific result formatting, separate from core persistence."""
from __future__ import annotations
from pathlib import Path
from mappa.mappa_workflow.ts_analytics.records import TSResult
from mappa.mappa_core.utils.io.records import _json_dump

def ts_record(result: TSResult):
    """Return a compact JSON boundary representation of one typed TS run."""
    from mappa.mappa_workflow.ts_analytics.scoring import POLICY
    return {
        "scoring_policy": POLICY,
        "schema": "rxn_core.ts/v2",
        "target": result.mechanisms[0].target.molecule.label
        if result.mechanisms else "",
        "elapsed_seconds": result.elapsed_seconds,
        "mechanisms": [{
            "id": index,
            "status": item.status,
            "reason": item.reason,
            "reactant_core_assignments": (
                0 if item.reactant_core_aam is None
                else len(item.reactant_core_aam.assignments)),
            "product_core_assignments": (
                0 if item.product_core_aam is None
                else len(item.product_core_aam.assignments)),
            "candidate_count": len(item.candidates),
            "selected": None if item.selected is None else {
                "assignment": {str(a): int(b)
                               for a, b in item.selected.assignment.pairs},
                "sources": sorted(item.selected.sources),
                "score": item.selected.score,
                "diagnostics": dict(item.selected.diagnostics),
                "wbo_progress": item.selected.wbo_progress,
                "mode_index": item.selected.mode_index,
                "frequency": item.selected.frequency,
                "event_terms": [dict(term)
                                for term in item.selected.event_terms],
            },
        } for index, item in enumerate(result.mechanisms, 1)],
    }


def write_ts_record(result: TSResult, path):
    _json_dump(path, ts_record(result))
    return Path(path)
