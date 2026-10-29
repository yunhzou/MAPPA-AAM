"""Run the public fragment conversion through the common portable executor."""
from dataclasses import asdict

from mappa.mappa_core.core.scheduler.queue import Task
from mappa.mappa_core.core.scheduler.search import _endpoint, _restore_problem


def plan_event_decoding(case, problem, catalogue, config, writer, *, case_seconds=3600., memory_bytes=2**30):
    """One saved case/cut catalogue, retaining pair grouping inside conversion."""
    if catalogue.problem_sha256 != catalogue.__class__(problem).problem_sha256:
        raise ValueError('Catalogue belongs to different endpoints')
    source = writer.put(dict(reactant=_endpoint(problem.reactant), product=_endpoint(problem.product), name=problem.name))
    saved = writer.put(catalogue.to_record())
    return [Task(f'{case}/convert', 'mappa.mappa_core.core.scheduler.decoding:decode_catalogue',
        dict(input=source, catalogue=saved, config=asdict(config)), case=case,
        case_seconds=case_seconds, memory_bytes=memory_bytes)]


def decode_catalogue(payload, workspace):
    from mappa.artifacts import decoded_events_record
    from mappa.final_branches import FinalBranchCatalogue
    from mappa.postprocessing import EventDecodeConfig, decode_events
    problem = _restore_problem(workspace.load(payload['input']))
    catalogue = FinalBranchCatalogue.from_record(problem, workspace.load(payload['catalogue']))
    return decoded_events_record(decode_events(catalogue, EventDecodeConfig(**payload['config'])))
