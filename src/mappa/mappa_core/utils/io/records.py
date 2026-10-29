"""Core result records and JSON boundaries; no workflow or viewer dependencies."""
from __future__ import annotations
import json, gc
from pathlib import Path
from mappa.mappa_core.domain import AAMProblem, AAMSearchConfig, AAMResult, AAMSearchMetrics, MolecularEndpoint, ReactionContext, ResolvedMechanism, AtomBijection
from mappa.mappa_core.core.search.graph import AAMSearchGraph

def endpoint_record(endpoint, *, array_lists=False):
    return dict(elements=endpoint.elements,
                coordinates=endpoint.coordinates.tolist() if array_lists else endpoint.coordinates,
                wbo=endpoint.wbo.tolist() if array_lists else endpoint.wbo,
                label=endpoint.label, energy=endpoint.energy, metadata=dict(endpoint.metadata))


def problem_from_record(record):
    return AAMProblem(MolecularEndpoint(**record['reactant']), MolecularEndpoint(**record['product']), record['name'])


def _aam_header(result):
    from dataclasses import asdict
    return {'schema': 'rxn_core.aam/v1', 'name': result.problem.name,
            'reactant': endpoint_record(result.problem.reactant, array_lists=True),
            'product': endpoint_record(result.problem.product, array_lists=True),
            'config': asdict(result.config), 'metrics': asdict(result.metrics)}


def aam_record(result: AAMResult, *, copy_graph=True):
    """Persist the search graph before selecting mechanisms or representatives."""
    return {**_aam_header(result),'graph': result.graph.to_record(copy=copy_graph)}


def aam_from_record(record, *, copy_graph=True):
    if record['schema'] != 'rxn_core.aam/v1':
        raise ValueError('unsupported AAM result schema')
    return _aam_result(record,AAMSearchGraph.from_record(record['graph'],copy=copy_graph))


def _aam_result(record,graph):
    def endpoint(raw):
        return MolecularEndpoint(tuple(raw['elements']), raw['coordinates'], raw['wbo'], raw['label'],
                                  energy=raw.get('energy'), metadata=raw.get('metadata', {}))
    return AAMResult(AAMProblem(endpoint(record['reactant']), endpoint(record['product']), record['name']),
                     AAMSearchConfig(**record['config']), graph,
                     AAMSearchMetrics(**record['metrics']))


def aam_json(result):
    """Encode the complete record while borrowing its immutable graph data."""
    enabled=gc.isenabled()
    try:
        gc.disable()
        return json.dumps(aam_record(result,copy_graph=False))+'\n'
    finally:
        if enabled:gc.enable()


def _mapping_record(mapping):
    return {str(source): int(target)
            for source, target in mapping.as_dict().items()}


def decoded_events_record(result):
    """Serialize event witnesses and their explicitly scoped completeness."""
    from mappa.mappa_core.post_processing.decoder.records import DecodedEvents
    if not isinstance(result, DecodedEvents):
        raise TypeError("decoded_events_record requires DecodedEvents")
    return {
        "schema": "mappa.decoded_events/v2",
        "method": "anchored_fragment_conversion",
        "scope": "Direct counts across retained converted branches",
        "converted_catalogue": result.catalogue.to_record(),
        "rematches": [item.to_record() for item in result.conversion.rematches],
        "converted_branches": [vars(branch) for branch in result.conversion.branches],
        "complete": result.complete,
        "minimum_proven": result.minimum_proven,
        "observed_minimum": result.observed_minimum,
        "event_window": result.event_window,
        "config": vars(result.config),
        "incomplete_paths": result.catalogue.incomplete_paths,
        "event_policy": result.index.policy,
        "family_reports": result.family_reports,
        "candidates": [{
            "id": candidate.id,
            "mapping": {str(a): int(b) for a, b in candidate.mapping.items()},
            "events": candidate.events,
            "endpoint_events": candidate.endpoint_events,
            "total": candidate.total,
            "family_ids": candidate.family_ids,
        } for candidate in result.candidates],
    }


def reaction_record(reaction: ReactionContext):
    """Serialize the exact R/P information required by the TS stage."""
    return {
        "schema": "rxn_core.resolved_reaction/v1",
        "name": reaction.problem.name,
        "atom_count": reaction.problem.atom_count,
        "mechanisms": [{
            "mapping": _mapping_record(item.mapping),
            "broken_bonds": [list(bond) for bond in item.broken_bonds],
            "formed_bonds": [list(bond) for bond in item.formed_bonds],
            "core_atoms": list(item.core_atoms),
        } for item in reaction.mechanisms],
    }


def reaction_from_record(record, problem: AAMProblem,
                         config: AAMSearchConfig | None = None):
    """Materialize a TS-stage context without executing R/P search."""
    if record.get("schema") != "rxn_core.resolved_reaction/v1":
        raise ValueError("unsupported resolved reaction schema")
    if int(record.get("atom_count", -1)) != problem.atom_count:
        raise ValueError("resolved reaction atom count differs from endpoints")
    mechanisms = []
    for raw in record.get("mechanisms") or ():
        mechanisms.append(ResolvedMechanism(
            mapping=AtomBijection.from_mapping(
                {int(a): int(b) for a, b in raw["mapping"].items()},
                degree=problem.atom_count),
            broken_bonds=tuple(tuple(bond)
                               for bond in raw.get("broken_bonds") or ()),
            formed_bonds=tuple(tuple(bond)
                               for bond in raw.get("formed_bonds") or ()),
            core_atoms=tuple(raw.get("core_atoms") or ()),
        ))
    if not mechanisms:
        raise ValueError("resolved reaction contains no mechanisms")
    return ReactionContext(
        problem, config or AAMSearchConfig(), tuple(mechanisms))


def _json_dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, default=float))
