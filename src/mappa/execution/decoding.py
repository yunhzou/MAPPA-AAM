"""Independent production family decoding with explicit absolute event bounds."""
from dataclasses import asdict

from .queue import Task
from .search import _endpoint, _restore_problem


def plan_event_decoding(case, problem, families, config, writer, *, case_seconds=3600., memory_bytes=2**30):
    """Schedule existing signed-event decoder calls, retaining incomplete reports.

    A caller deriving minimum+n windows must first establish a valid absolute
    bound for each requested seed prefix. A bound from only the ten-seed minimum
    is not silently reused for other prefixes. None means unbounded decoding.
    Complete and partial families use the same mapped-boundary event rule.
    """
    if config.event_slack is not None:
        raise ValueError('Supply an absolute window established for the requested prefix')
    from ..final_branches import FinalBranchCatalogue
    FinalBranchCatalogue(problem)
    reference = writer.put(dict(reactant=_endpoint(problem.reactant), product=_endpoint(problem.product),
        name=problem.name))
    tasks = []
    for index, family in enumerate(families):
        saved = writer.put(family.to_record())
        tasks.append(Task(f'{case}/decode/{index:06d}', 'mappa.execution.decoding:decode_family',
            dict(input=reference, family=saved, config=asdict(config), family_index=index),
            case=case, case_seconds=case_seconds, memory_bytes=memory_bytes))
    return tasks


def decode_family(payload, workspace):
    from ..event_patterns import SignedEventIndex, extract_path_events
    from ..final_branches import FinalFamily
    from ..postprocessing import EventDecodeConfig
    config = EventDecodeConfig(**payload['config'])
    policy = (config.threshold, config.metal_threshold)
    def prepare():
        problem = _restore_problem(workspace.load(payload['input']))
        return problem, SignedEventIndex(problem, threshold=config.threshold, metal_threshold=config.metal_threshold)
    problem, index = workspace.prepared(('decode', payload['input'], policy), prepare)
    family = FinalFamily.from_record(workspace.load(payload['family']))
    report = extract_path_events(family.as_path(problem), problem, index,
        max_events=config.max_events, seconds=config.seconds_per_family, max_patterns=config.max_patterns_per_family)
    return dict(family_index=payload['family_index'], report=report, max_events=config.max_events)
