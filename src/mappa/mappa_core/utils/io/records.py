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
    return {'schema': 'mappa.aam/v2', 'name': result.problem.name,
            'reactant': endpoint_record(result.problem.reactant, array_lists=True),
            'product': endpoint_record(result.problem.product, array_lists=True),
            'config': asdict(result.config), 'metrics': asdict(result.metrics)}


def aam_record(result: AAMResult, *, copy_graph=True):
    """Persist the search graph before selecting mechanisms or representatives."""
    return {**_aam_header(result),'graph': result.graph.to_record(copy=copy_graph)}


def aam_from_record(record, *, copy_graph=True):
    if record['schema'] != 'mappa.aam/v2':
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


def decoded_events_record(result, *, array_lists=True):
    """Persist the windowed fragment model and unchanged mapping constraints.

    Tables and responses shared across branches are stored once. Whole-reaction
    reports and canonical event IDs are explicit queries, not part of export.
    Per-atom lookup maps are derived from the saved catalogue, not persisted.
    Binary transport retains shared numeric arrays; JSON export converts them.
    """
    from mappa.mappa_core.post_processing.decoder.records import EventInterpretation
    if not isinstance(result, EventInterpretation):
        raise TypeError('decoded_events_record requires EventInterpretation')
    tables = []
    families = [dict(family_id=fid, tables=[]) for fid in range(len(result.catalogue.families))]
    for fid, number, table in result._export_tables():
        for fragment in result.catalogue.families[fid].fragments:
            result.local.prepare(fragment)
        if table is not None:
            tables.append(dict(fragments=table.fragments, cut_bonds=sorted(table.cut_bonds),
                fixed_boundary=[dict(fragments=region, counts=tuple(map(int, counts)))
                                for region, counts in table.fixed_boundary.items()],
                reactant_boundary_pairs=table.reactant_boundary_pairs.tolist(),
                product_boundary_pairs=table.product_boundary_pairs.tolist(),
                event_window=result.event_window, components=table.components))
        families[fid]['tables'].append(number)
        del table
    from dataclasses import asdict
    fragment_ids = {fragment.key:i for i,fragment in enumerate(result.inventory.fragments)}
    prepared = [dict(fragment=fragment_ids[key], optimization={name:item
                        for name,item in asdict(value.optimization).items()
                        if name not in ('search_seconds', 'elapsed_seconds', 'work')},
                     response=value.response.to_record())
                for key,value in result.local.prepared.items() if key in fragment_ids]
    record = dict(schema='mappa.decoded_events/v14', method='binary_fragment_events',
        scope='Independent ordinary boundaries with selected fragment starting correspondences and binary symmetry',
        catalogue=result.catalogue.to_record(array_lists=array_lists), event_policy=result.index.policy,
        problem=dict(name=result.catalogue.problem.name,
            reactant=endpoint_record(result.catalogue.problem.reactant,array_lists=array_lists),
            product=endpoint_record(result.catalogue.problem.product,array_lists=array_lists)),
        families=families, tables=tables,
        fragment_responses=prepared, fragment_seconds=result.local.fragment_seconds,
        retained_contexts=result.window.contexts,
        window=result.window.to_record(),
        mapping_space_preserved=True,
        decoded_coverage='selected_start_binary_symmetry',
        minimum_proven=result.minimum_proven, minimum_lower_bound=result.minimum_lower_bound,
        event_counts_approximate=True, minimum_events=result.minimum_events,
        event_slack=result.event_slack, event_window=result.event_window,
        minimum_scope='independent_fragment_model',
        incomplete_paths=result.catalogue.incomplete_paths,
        exhaustive_event_patterns=False, global_minimum_proven=False)
    if result.event_problem_sha256 != result.matching_problem_sha256:
        record.update(matching_problem_sha256=result.matching_problem_sha256,
            event_problem_sha256=result.event_problem_sha256,
            event_problem=dict(name=result.event_problem.name,
                reactant=endpoint_record(result.event_problem.reactant, array_lists=array_lists),
                product=endpoint_record(result.event_problem.product, array_lists=array_lists)))
    return record


def decoded_events_from_record(record):
    """Restore saved response decisions without rerunning discovery optimization."""
    from dataclasses import replace
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    from mappa.mappa_core.post_processing.decoder.api import EventDecodeSession
    from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig
    from mappa.mappa_core.post_processing.decoder.fragments import DecodedFragment
    from mappa.mappa_core.core.symmetry.minimum import FragmentMinimum, FragmentSymmetryEvents, SymmetryWork
    if record['schema'] != 'mappa.decoded_events/v14':
        raise ValueError('require a saved binary-fragment event response')
    problem = problem_from_record(record['problem'])
    events = problem_from_record(record['event_problem']) if 'event_problem' in record else problem
    catalogue = FinalBranchCatalogue.from_record(problem,record['catalogue'])
    policy = record['event_policy']
    config = EventDecodeConfig(policy['threshold'],policy['metal_threshold'],record['event_slack'],record['fragment_seconds'])
    result = EventDecodeSession(problem,config,event_problem=events).decode(catalogue)
    for saved in record['fragment_responses']:
        fragment = result.inventory.fragments[saved['fragment']]
        values = dict(saved['optimization'])
        for name in ('witness','discovery_witness'):
            values[name] = tuple(tuple(pair) for pair in values[name])
        for name in ('source_operation','target_operation'):
            values[name] = tuple(values[name])
        values.update(work=SymmetryWork(), search_seconds=0., elapsed_seconds=0.)
        minimum = FragmentMinimum(**values)
        selected = replace(fragment,relation=replace(fragment.relation,witness=minimum.discovery_witness))
        response = FragmentSymmetryEvents.from_record(selected,saved['response'],event_problem=events,
            threshold=config.threshold,metal_threshold=config.metal_threshold)
        result.local.prepared[fragment.key] = DecodedFragment(fragment,minimum,response)
    result.local.frozen = True
    from mappa.mappa_core.post_processing.decoder.window import EventWindow
    result.window = EventWindow.from_record(record['window'])
    return result


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
