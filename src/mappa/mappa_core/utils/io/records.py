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


def refined_mapping_record(result):
    """Serialize event witnesses and their explicitly scoped completeness."""
    from mappa.mappa_core.post_processing.hierarchical_mapping import RefinedMapping
    if not isinstance(result, RefinedMapping):
        raise TypeError("refined_mapping_record requires RefinedMapping")
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
        "counting_statistics": result.counting_statistics,
        "candidates": [{
            "id": candidate.id,
            "mapping": {str(a): int(b) for a, b in candidate.mapping.items()},
            "events": candidate.events,
            "endpoint_events": candidate.endpoint_events,
            "total": candidate.total,
            "family_ids": candidate.family_ids,
        } for candidate in result.candidates],
    }


def decoded_events_record(result):
    """Persist unchanged mappings and local event functions, not an enumeration.

    Tables and responses shared across branches are stored once. Representative
    counts are labelled explicitly; they do not prove a minimum over choices.
    """
    from mappa.mappa_core.post_processing.decoder.records import EventInterpretation
    if not isinstance(result, EventInterpretation):
        raise TypeError('decoded_events_record requires EventInterpretation')
    tables, table_ids, responses, response_ids, families = [], {}, [], {}, []
    def term_record(term):
        response = term.response
        if response not in response_ids:
            response_ids[response] = len(responses)
            responses.append(dict(source_weight=response.source_weight,
                threshold=response.threshold, target_weights=response.target_weights.tolist(),
                categories=response.categories.tolist()))
        return dict(atoms=term.atoms, action_ids=term.action_ids,
            response=response_ids[response], initial_pair=term.initial_pair,
            endpoint_case=term.endpoint_case, manual_cut=term.manual_cut)
    for fid, partitions in enumerate(result.fragment_partitions):
        saved = []
        for partition in range(len(partitions)):
            table = result.tables(fid, partition)
            identity = id(table)
            if identity not in table_ids:
                table_ids[identity] = len(tables)
                tables.append(dict(fragments=table.fragments, cut_bonds=sorted(table.cut_bonds),
                    internal=[dict(fragment=k, terms=[term_record(t) for t in v])
                              for k,v in table.internal.items()],
                    boundary=[dict(fragments=k, terms=[term_record(t) for t in v])
                              for k,v in table.boundary.items()],
                    product_boundary_pairs=table.product_boundary_pairs.tolist(),
                    atom_domains=table.atom_domains, atom_actions=table.atom_actions,
                    components=table.components, statistics=table.statistics))
            saved.append(table_ids[identity])
        families.append(dict(family_id=fid, tables=saved,
            action_responses=[r.generators for r in result.family_actions[fid]]))
    representatives = [dict(family_id=i, **result.evaluate(i))
                       for i in range(len(result.catalogue.families))]
    return dict(schema='mappa.decoded_events/v3', method='compressed_event_interpretation',
        scope='Original compressed mappings with correlated local event responses',
        catalogue=result.catalogue.to_record(), event_policy=result.index.policy,
        problem=dict(name=result.catalogue.problem.name,
            reactant=endpoint_record(result.catalogue.problem.reactant,array_lists=True),
            product=endpoint_record(result.catalogue.problem.product,array_lists=True)),
        families=families, tables=tables, responses=responses,
        representatives=representatives, mapping_space_preserved=True,
        incomplete_paths=result.catalogue.incomplete_paths,
        exhaustive_event_patterns=False, global_minimum_proven=False)


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
