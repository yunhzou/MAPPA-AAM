"""Decode saved catalogues or finalized cuts through the common task pool."""
from dataclasses import asdict, replace
import time

from mappa.mappa_core.core.scheduler.queue import Task
from mappa.mappa_core.utils.io.records import endpoint_record, problem_from_record, decoded_events_record
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig
from mappa.mappa_core.post_processing.decoder.api import decode_events


def plan_event_decoding(case, problem, catalogue, config, writer, *, case_seconds=3600., memory_bytes=2**30):
    """Schedule an already saved catalogue without repeating matching."""
    if catalogue.problem_sha256 != FinalBranchCatalogue(problem).problem_sha256:
        raise ValueError('Catalogue belongs to different endpoints')
    source = writer.put(dict(reactant=endpoint_record(problem.reactant), product=endpoint_record(problem.product), name=problem.name))
    saved = writer.put(catalogue.to_record())
    return [Task(f'{case}/decode', 'mappa.mappa_core.core.scheduler.decoding:decode_catalogue',
        dict(input=source, catalogue=saved, config=asdict(config)), case=case,
        case_seconds=case_seconds, memory_bytes=memory_bytes)]


def plan_cut_decoding(tasks, cuts, config=None):
    """One task per finalized cut; all its seeds share the event index and tables."""
    config = config or EventDecodeConfig()
    cuts = tuple(cuts)
    keys = set(cuts)
    finals = {task.key: task for task in tasks if task.key in keys}
    for key in cuts:
        final = finals[key]
        yield replace(final, key=key.rsplit('/', 1)[0]+'/decode',
            operation='mappa.mappa_core.core.scheduler.decoding:decode_cut', dependencies=(key,),
            payload=dict(input=final.payload['input'], cut=key, number=final.payload['number'],
                         seeds=final.payload['seeds'], config=asdict(config)))


def decode_catalogue(payload, workspace):
    problem = problem_from_record(workspace.load(payload['input']))
    catalogue = FinalBranchCatalogue.from_record(problem, workspace.load(payload['catalogue']))
    return decoded_events_record(decode_events(catalogue, EventDecodeConfig(**payload['config'])))


def decode_graph(problem, graph, config=None):
    """Interpret a cut once, retaining exact seed membership after deduplication."""
    started, cpu = time.perf_counter(), time.process_time()
    catalogue = FinalBranchCatalogue(problem).add_graph(graph)
    families = [set() for _ in graph.contexts]
    for fid, family in enumerate(catalogue.families):
        for origin in family.provenance:
            families[origin['context']].add(fid)
    seeds = []
    for seed, ids in enumerate(families):
        tiers = {fid: (sum(problem.reactant.elements[a] != 'H' for a, _ in catalogue.families[fid].mapping),
                       len(catalogue.families[fid].mapping)) for fid in ids}
        tier = max(tiers.values(), default=None)
        seeds.append(dict(seed=seed, family_ids=sorted(ids), coverage=tier,
                          selected_family_ids=sorted(fid for fid in ids if tiers[fid] == tier)))
    setup = dict(wall_seconds=time.perf_counter()-started, cpu_seconds=time.process_time()-cpu)
    started, cpu = time.perf_counter(), time.process_time()
    result = decode_events(catalogue, config)
    for fid, partitions in enumerate(result.fragment_partitions):
        for partition in range(len(partitions)):
            result.tables(fid, partition)
    decoding = dict(wall_seconds=time.perf_counter()-started, cpu_seconds=time.process_time()-cpu)
    started, cpu = time.perf_counter(), time.process_time()
    record = decoded_events_record(result)
    serialization = dict(wall_seconds=time.perf_counter()-started, cpu_seconds=time.process_time()-cpu)
    return dict(decoded=record, seeds=seeds, timing=dict(setup=setup, decoding=decoding, record=serialization))


def decode_cut(payload, workspace):
    from .search import _session
    started, cpu = time.perf_counter(), time.process_time()
    problem = _session(payload, workspace).context[0]
    finalized = workspace.dependency(payload['cut'])
    if len(finalized['graph'].contexts) != payload['seeds'] or finalized['cut'] != payload['number']:
        raise ValueError('Finalized cut or seed denominator differs')
    loading = dict(wall_seconds=time.perf_counter()-started, cpu_seconds=time.process_time()-cpu)
    value = decode_graph(problem, finalized['graph'], EventDecodeConfig(**payload['config']))
    value['timing']['loading'] = loading
    parent = workspace.dependencies[payload['cut']]['blob']
    return dict(value, cut=payload['number'], input=payload['input'], finalized=parent)
