"""Trusted saved graphs and explicit historical class relocation."""
from __future__ import annotations
import json, gc, gzip, pickle, io
from pathlib import Path
from mappa.mappa_core.core.search.graph import AAMSearchGraph

from mappa.mappa_core.utils.io.records import _aam_header, _aam_result, aam_from_record

_CHECKPOINT_MODULES = {
    ('mappa.artifacts', '_CheckpointUnpickler'): 'mappa.mappa_core.utils.io.checkpoints',
    ('mappa.frag', 'WeightedNode'): 'mappa.mappa_core.core.matching.graphs',
    ('mappa.frag', 'WeightedGraph'): 'mappa.mappa_core.core.matching.graphs',
    ('mappa.smiles', 'FormalWBOEndpoint'): 'mappa.mappa_core.utils.io.smiles',
    ('mappa.search_validation', 'RepresentativeValidationWorkspace'): 'mappa.mappa_core.core.search.branching.validation',

    ('mappa.domain', 'MolecularEndpoint'): 'mappa.mappa_core.domain',
    ('mappa.domain', 'AAMProblem'): 'mappa.mappa_core.domain',
    ('mappa.domain', 'AAMSearchConfig'): 'mappa.mappa_core.domain',
    ('mappa.domain', 'AAMSearchMetrics'): 'mappa.mappa_core.domain',
    ('mappa.domain', 'AAMResult'): 'mappa.mappa_core.domain',
    ('mappa.domain', 'ResolvedMechanism'): 'mappa.mappa_core.domain',
    ('mappa.domain', 'ReactionContext'): 'mappa.mappa_core.domain',
    ('mappa.domain', 'VibrationalModes'): 'mappa.mappa_core.domain',
    ('mappa.domain', 'TransitionStateTarget'): 'mappa.mappa_workflow.ts_analytics.records',
    ('mappa.domain', 'AtomAssignment'): 'mappa.mappa_core.domain',
    ('mappa.domain', 'CoreAAMBranch'): 'mappa.mappa_core.domain',
    ('mappa.domain', 'CoreAAMResult'): 'mappa.mappa_core.domain',
    ('mappa.domain', 'TSScore'): 'mappa.mappa_workflow.ts_analytics.records',
    ('mappa.domain', 'TSScoringConfig'): 'mappa.mappa_workflow.ts_analytics.records',
    ('mappa.domain', 'TSMechanismResult'): 'mappa.mappa_workflow.ts_analytics.records',
    ('mappa.domain', 'TSResult'): 'mappa.mappa_workflow.ts_analytics.records',

    ('mappa.postprocessing', 'EventDecodeConfig'): 'mappa.mappa_core.post_processing.decoder.records',
    ('mappa.postprocessing', 'EventCandidate'): 'mappa.mappa_core.post_processing.decoder.records',
    ('mappa.postprocessing', 'ShuffleQuery'): 'mappa.mappa_workflow.mapping_analytics.queries',
    ('mappa.postprocessing', 'DecodedEvents'): 'mappa.mappa_core.post_processing.decoder.records',
    ('mappa.fragment_conversion', 'FragmentRematch'): 'mappa.mappa_core.post_processing.decoder.conversion',
    ('mappa.fragment_conversion', 'ConvertedBranch'): 'mappa.mappa_core.post_processing.decoder.conversion',
    ('mappa.fragment_conversion', 'FragmentConversion'): 'mappa.mappa_core.post_processing.decoder.conversion',
    ('mappa.event_patterns', 'SignedEventIndex'): 'mappa.mappa_core.post_processing.decoder.events',
    ('mappa.event_patterns', '_EventTerms'): 'mappa.mappa_workflow.mapping_analytics.event_query',
    ('mappa.final_branches', 'FinalFamily'): 'mappa.mappa_core.post_processing.final_branches',
    ('mappa.final_branches', 'FinalBranchCatalogue'): 'mappa.mappa_core.post_processing.final_branches',
    ('mappa.family_query', 'SymbolicActions'): 'mappa.mappa_workflow.mapping_analytics.family',
    ('mappa.family_query', 'CompiledFamily'): 'mappa.mappa_workflow.mapping_analytics.family',
    ('mappa.reference_constraints', 'Choice'): 'mappa.mappa_workflow.mapping_analytics.constraints',
    ('mappa.reference_constraints', 'ReferenceConstraints'): 'mappa.mappa_workflow.mapping_analytics.constraints',

    ('mappa.execution.queue', 'Task'): 'mappa.mappa_core.core.scheduler.queue',
    ('mappa.execution.queue', 'TaskQueue'): 'mappa.mappa_core.core.scheduler.queue',
    ('mappa.execution.bundles', 'Blob'): 'mappa.mappa_core.core.scheduler.bundles',
    ('mappa.execution.bundles', 'BundleWriter'): 'mappa.mappa_core.core.scheduler.bundles',
    ('mappa.execution.runtime', 'Workspace'): 'mappa.mappa_core.core.scheduler.runtime',
    ('mappa.execution.runtime', '_QueueServer'): 'mappa.mappa_core.core.scheduler.runtime',

    ('mappa.mappa_core.core.search.branch', '_StateKey'): 'mappa.mappa_core.core.search.branching.construction',
    ('mappa.mappa_core.core.search.branch', 'BranchLimitExceeded'): 'mappa.mappa_core.core.search.branching.construction',
    ('mappa.mappa_core.core.search.branch', '_Branch'): 'mappa.mappa_core.core.search.branching.construction',
    ('mappa.mappa_core.core.search.shared', '_SharedGraphBuilder'): 'mappa.mappa_core.core.search.branching.shared',
    ('mappa.mappa_core.core.search.shared', 'SharedSeedSearch'): 'mappa.mappa_core.core.search.branching.shared',
    ('mappa.mappa_core.core.search.shared', 'AdaptiveResult'): 'mappa.mappa_core.core.search.branching.shared',
    ('mappa.mappa_core.core.search.shared', 'AdaptiveSearchCondition'): 'mappa.mappa_core.core.search.branching.shared',
    ('mappa.mappa_core.core.search.reuse', 'CutReplayView'): 'mappa.mappa_core.core.search.branching.reuse',
    ('mappa.mappa_core.core.search.reuse', '_CutSession'): 'mappa.mappa_core.core.search.branching.reuse',
    ('mappa.mappa_core.core.search.reuse', 'CutReplay'): 'mappa.mappa_core.core.search.branching.reuse',
    ('mappa.mappa_core.core.search.reuse', 'FragmentRepair'): 'mappa.mappa_core.core.search.branching.reuse',
    ('mappa.mappa_core.core.search.competition', 'CompetitionConfig'): 'mappa.mappa_core.core.search.branching.competition',
    ('mappa.mappa_core.core.search.competition', 'CompetitionResult'): 'mappa.mappa_core.core.search.branching.competition',
    ('mappa.mappa_core.core.search.subgraph', 'SubgraphMatch'): 'mappa.mappa_core.core.search.queries.subgraph',
    ('mappa.mappa_core.core.search.subgraph', 'SubgraphSearchResult'): 'mappa.mappa_core.core.search.queries.subgraph',
    ('mappa.mappa_core.core.search.planning', 'DirectedAAMResult'): 'mappa.mappa_core.core.search.api',

    ('mappa.mappa_core.core.matching.state', '_SymBlock'): 'mappa.mappa_core.core.matching.representation.state',
    ('mappa.mappa_core.core.matching.state', '_SymCand'): 'mappa.mappa_core.core.matching.representation.state',
    ('mappa.mappa_core.core.matching.records', '_IsoResult'): 'mappa.mappa_core.core.matching.representation.records',
    ('mappa.mappa_core.core.matching.records', 'IslandBranchLimitExceeded'): 'mappa.mappa_core.core.matching.representation.records',
    ('mappa.mappa_core.core.matching.records', 'FragmentPlacement'): 'mappa.mappa_core.core.matching.representation.records',
    ('mappa.mappa_core.core.matching.records', 'FragmentMatchContext'): 'mappa.mappa_core.core.matching.representation.records',
    ('mappa.mappa_core.core.matching.records', 'FragmentMatchConfig'): 'mappa.mappa_core.core.matching.representation.records',
    ('mappa.mappa_core.core.matching.records', 'FragmentMatchResult'): 'mappa.mappa_core.core.matching.representation.records',
    ('mappa.mappa_core.core.matching.canonical', '_PartialMappingCanonicalizer'): 'mappa.mappa_core.core.matching.symmetry.canonical',
    ('mappa.mappa_core.core.matching.canonical', '_CandidateAutomorphismCanonicalizer'): 'mappa.mappa_core.core.matching.symmetry.canonical',
    ('mappa.mappa_core.core.matching.orbits', '_OrbitMap'): 'mappa.mappa_core.core.matching.symmetry.orbits',
    ('mappa.mappa_core.core.matching.symmetry', 'SymmetryWorkspace'): 'mappa.mappa_core.core.matching.symmetry.workspace',
    ('mappa.mappa_core.core.matching.symmetry', 'ConditionedSymmetryWorkspace'): 'mappa.mappa_core.core.matching.symmetry.workspace',
    ('mappa.mappa_core.core.matching.extend', '_ExtensionContext'): 'mappa.mappa_core.core.matching.growth.extend',
    ('mappa.mappa_core.core.matching.frontier', 'GrowthRandom'): 'mappa.mappa_core.core.matching.growth.frontier',
    ('mappa.mappa_core.core.matching.native', '_NativeGraphView'): 'mappa.mappa_core.core.matching.growth.native',

    ('mappa.alignment.post_aam', 'AtomPermutation'): 'mappa.mappa_core.domain',
    ('mappa.alignment.post_aam', 'AtomBijection'): 'mappa.mappa_core.domain',
    ('mappa.alignment.post_aam', 'PermutationGroup'): 'mappa.mappa_core.domain',
    ('mappa.alignment.post_aam', 'SymmetryDomain'): 'mappa.mappa_core.domain',
    ('mappa.alignment.post_aam', 'FragmentMatch'): 'mappa.mappa_core.domain',
    ('mappa.alignment.post_aam', 'AAMHierarchy'): 'mappa.mappa_core.domain',
    ('mappa.alignment.post_aam', 'AAMHierarchyView'): 'mappa.mappa_core.domain',
    ('mappa.alignment.post_aam', 'AAMHierarchyChain'): 'mappa.mappa_core.domain',

    ('mappa.search_symmetry', 'SymmetryWorkspace'): 'mappa.mappa_core.core.matching.symmetry.workspace',
    ('mappa.conditioned_symmetry', 'ConditionedSymmetryWorkspace'): 'mappa.mappa_core.core.matching.symmetry.workspace',
    ('mappa.fragment', 'FragmentMatchConfig'): 'mappa.mappa_core.core.matching.representation.records',
    ('mappa.fragment', 'FragmentMatchContext'): 'mappa.mappa_core.core.matching.representation.records',
    ('mappa.fragment', 'FragmentMatchResult'): 'mappa.mappa_core.core.matching.representation.records',
    ('mappa.fragment', 'FragmentPlacement'): 'mappa.mappa_core.core.matching.representation.records',
    ('mappa.growth.frontier', 'GrowthRandom'): 'mappa.mappa_core.core.matching.growth.frontier',
    ('mappa.growth.native', '_NativeGraphView'): 'mappa.mappa_core.core.matching.growth.native',
    ('mappa.growth.result', 'IslandBranchLimitExceeded'): 'mappa.mappa_core.core.matching.representation.records',
    ('mappa.growth.result', '_IsoResult'): 'mappa.mappa_core.core.matching.representation.records',
    ('mappa.matcher.bond_policy', 'BondTolerance'): 'mappa.mappa_core.core.matching.policy',
    ('mappa.matcher.bond_policy', 'ResolvedBondTolerance'): 'mappa.mappa_core.core.matching.policy',
    ('mappa.matcher.canonical', '_CandidateAutomorphismCanonicalizer'): 'mappa.mappa_core.core.matching.symmetry.canonical',
    ('mappa.matcher.canonical', '_PartialMappingCanonicalizer'): 'mappa.mappa_core.core.matching.symmetry.canonical',
    ('mappa.matcher.extend', '_ExtensionContext'): 'mappa.mappa_core.core.matching.growth.extend',
    ('mappa.matcher.orbits', '_OrbitMap'): 'mappa.mappa_core.core.matching.symmetry.orbits',
    ('mappa.matcher.policy', 'AttributeNodeMatchPolicy'): 'mappa.mappa_core.core.matching.policy',
    ('mappa.matcher.policy', 'CallableNodeMatchPolicy'): 'mappa.mappa_core.core.matching.policy',
    ('mappa.matcher.policy', 'ElementNodeMatchPolicy'): 'mappa.mappa_core.core.matching.policy',
    ('mappa.matcher.state', '_SymBlock'): 'mappa.mappa_core.core.matching.representation.state',
    ('mappa.matcher.state', '_SymCand'): 'mappa.mappa_core.core.matching.representation.state',

    ('mappa.aam', '_GrowthCounts'): 'mappa.mappa_core.core.search.session',
    ('mappa.aam', 'SearchSession'): 'mappa.mappa_core.core.search.session',
    ('mappa.aam', 'AAMCheckpointResult'): 'mappa.mappa_core.core.search.api',
    ('mappa.alignment.branch', '_StateKey'): 'mappa.mappa_core.core.search.branching.construction',
    ('mappa.alignment.branch', 'BranchLimitExceeded'): 'mappa.mappa_core.core.search.branching.construction',
    ('mappa.alignment.branch', '_Branch'): 'mappa.mappa_core.core.search.branching.construction',
    ('mappa.adaptive_search', 'AdaptiveResult'): 'mappa.mappa_core.core.search.branching.shared',
    ('mappa.adaptive_search', 'AdaptiveSearchCondition'): 'mappa.mappa_core.core.search.branching.shared',
    ('mappa.shared_seed_search', '_SharedGraphBuilder'): 'mappa.mappa_core.core.search.branching.shared',
    ('mappa.shared_seed_search', 'SharedSeedSearch'): 'mappa.mappa_core.core.search.branching.shared',
    ('mappa.cut_replay', 'CutReplayView'): 'mappa.mappa_core.core.search.branching.reuse',
    ('mappa.cut_replay', '_CutSession'): 'mappa.mappa_core.core.search.branching.reuse',
    ('mappa.cut_replay', 'CutReplay'): 'mappa.mappa_core.core.search.branching.reuse',
    ('mappa.cut_replay', 'FragmentRepair'): 'mappa.mappa_core.core.search.branching.reuse',
    ('mappa.search_orientation', 'AAMSearchPlan'): 'mappa.mappa_core.core.search.planning',
    ('mappa.search_orientation', 'DirectedAAMResult'): 'mappa.mappa_core.core.search.api',
    ('mappa.search_graph', '_BondRecords'): 'mappa.mappa_core.core.search.graph',
    ('mappa.search_graph', 'SearchContext'): 'mappa.mappa_core.core.search.graph',
    ('mappa.search_graph', 'SearchState'): 'mappa.mappa_core.core.search.graph',
    ('mappa.search_graph', 'FragmentTransition'): 'mappa.mappa_core.core.search.graph',
    ('mappa.search_graph', 'SearchStop'): 'mappa.mappa_core.core.search.graph',
    ('mappa.search_graph', 'SearchPath'): 'mappa.mappa_core.core.search.graph',
    ('mappa.search_graph', 'PathRealization'): 'mappa.mappa_core.core.search.graph',
    ('mappa.search_graph', 'SearchBranch'): 'mappa.mappa_core.core.search.graph',
    ('mappa.search_graph', 'AAMSearchGraph'): 'mappa.mappa_core.core.search.graph',
    ('mappa.search_graph', 'SearchGraphBuilder'): 'mappa.mappa_core.core.search.graph',
    ('mappa.subgraph', 'SubgraphMatch'): 'mappa.mappa_core.core.search.queries.subgraph',
    ('mappa.subgraph', 'SubgraphSearchResult'): 'mappa.mappa_core.core.search.queries.subgraph',
    ('mappa.competition', 'CompetitionConfig'): 'mappa.mappa_core.core.search.branching.competition',
    ('mappa.competition', 'CompetitionResult'): 'mappa.mappa_core.core.search.branching.competition',
}


class _CheckpointUnpickler(pickle.Unpickler):
    """Resolve pre-rename globals in trusted internal checkpoints.

    This is namespace compatibility, not a safe loader for untrusted pickle.
    """

    def find_class(self, module, name):
        for legacy in ("rxn_core", "graft"):
            if module == legacy or module.startswith(legacy + "."):
                module = "mappa" + module[len(legacy):]
                break
        # Explicit class-path relocation for checkpoints saved before ownership consolidation.
        module = _CHECKPOINT_MODULES.get((module, name), module)
        return super().find_class(module, name)


def _load_checkpoint(stream):
    return _CheckpointUnpickler(stream).load()


def _loads_checkpoint(data):
    return _load_checkpoint(io.BytesIO(data))


def write_aam_checkpoint(result,path):
    """Trusted internal snapshot preserving shared Python graph objects.

    Keep the producing engine with this file. JSON remains the interchange
    format; this pickle-based checkpoint must never be loaded from strangers.
    """
    path=Path(path);temporary=path.with_suffix(path.suffix+'.tmp')
    record={**_aam_header(result),'schema':'rxn_core.aam_checkpoint/v1','graph':result.graph}
    enabled=gc.isenabled()
    try:
        gc.disable()
        with gzip.open(temporary,'wb',compresslevel=1) as stream:
            pickle.dump(record,stream,protocol=5)
        temporary.replace(path)
    finally:
        if enabled:gc.enable()


def read_aam_checkpoint(path):
    """Read only trusted internally produced snapshots, never untrusted pickle."""
    enabled=gc.isenabled()
    try:
        gc.disable()
        with gzip.open(path,'rb') as stream:record=_load_checkpoint(stream)
        if record['schema']!='rxn_core.aam_checkpoint/v1':
            raise ValueError('Unsupported AAM checkpoint schema')
        return _aam_result(record,record['graph'])
    finally:
        if enabled:gc.enable()


def write_graph_checkpoint(graph,path):
    """Save a trusted, finalized cut graph without expanding shared values."""
    _write_graph_checkpoint(graph,path,'rxn_core.finalized_cut/v1')


def _write_graph_checkpoint(graph,path,schema):
    path=Path(path);temporary=path.with_suffix(path.suffix+'.tmp')
    with gzip.open(temporary,'wb',compresslevel=1) as stream:
        pickle.dump((schema,graph),stream,protocol=5)
    temporary.replace(path)


def read_graph_checkpoint(path):
    """Load only internally produced finalized-cut checkpoints."""
    with gzip.open(path,'rb') as stream:schema,graph=_load_checkpoint(stream)
    if schema!='rxn_core.finalized_cut/v1':raise ValueError('Unsupported cut checkpoint schema')
    return graph


def read_aam(stream):
    """Load owned, acyclic JSON directly into the typed graph."""
    enabled=gc.isenabled()
    try:
        gc.disable()
        return aam_from_record(json.load(stream),copy_graph=False)
    finally:
        if enabled:gc.enable()


def raw_cut_paths(directory):
    """Completed raw cuts in either supported format, in canonical cut order."""
    directory=Path(directory);paths={}
    for path in (*directory.glob('cut_*.json'), *directory.glob('cut_*.raw.pkl.gz')):
        index=int(path.name.split('_')[1].split('.')[0])
        if index in paths:raise ValueError(f'Duplicate raw cut checkpoint {index}')
        paths[index]=path
    return [paths[index] for index in sorted(paths)]


def write_raw_cut(graph,path):
    """Persist the requested raw-cut format atomically; binary is trusted-only."""
    path=Path(path)
    if path.name.endswith('.raw.pkl.gz'):
        _write_graph_checkpoint(graph,path,'rxn_core.raw_cut/v1')
    elif path.suffix=='.json':
        temporary=path.with_suffix('.json.tmp')
        # Interchange JSON uses the fast C encoder. Large persistent searches
        # select compact checkpoints to avoid expanding shared graph payloads.
        temporary.write_text(json.dumps(graph.to_record(copy=False)))
        temporary.replace(path)
    else:raise ValueError(f'Unsupported raw cut format: {path}')


def read_raw_cut(path, *, tuple_pool=None):
    """Read only trusted internal binary cuts; JSON remains interchange-safe."""
    path=Path(path)
    if path.name.endswith('.raw.pkl.gz'):
        with gzip.open(path,'rb') as stream:schema,graph=_load_checkpoint(stream)
        if schema!='rxn_core.raw_cut/v1':raise ValueError('Unsupported raw cut checkpoint schema')
        return graph
    if path.suffix=='.json':
        return AAMSearchGraph.from_record(json.loads(path.read_bytes()),copy=False,tuple_pool=tuple_pool)
    raise ValueError(f'Unsupported raw cut format: {path}')
