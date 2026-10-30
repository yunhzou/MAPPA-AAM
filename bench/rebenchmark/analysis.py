"""Original Golden scoring, independently retained for every seed and cut."""
from dataclasses import replace, asdict
import gc
import json
from pathlib import Path
import time

from .evidence import read, save, sha
from .search import select_contexts


def analyze_seed(plan, ref, graph, seed, out, number, **provenance):
    """Evaluate one finalized context; both campaign executors call this operation."""
    from mappa import AAMResult
    from mappa.mappa_core.domain import AAMSearchMetrics
    from golden_evaluation import evaluate_planned
    from publication_analysis import rank_archive
    one=select_contexts(graph,[seed])
    result=AAMResult(plan.problem,plan.config,one,AAMSearchMetrics.from_record({},0))
    c,w=time.process_time(),time.perf_counter()
    classes=rank_archive(result,plan)
    for row in classes:
        row.update(cut=number,seed=seed+1)
    save(out/'classes.json',classes)
    ranking=dict(cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w)
    c,w=time.process_time(),time.perf_counter()
    verdict=evaluate_planned(result,plan,ref['features'],ref['mapping'],
        seconds=3600)
    save(out/'evaluation.json',dict(**verdict,cpu_seconds=time.process_time()-c,
         wall_seconds=time.perf_counter()-w,cut=number,seed=seed+1,ranking=ranking,
         **provenance))
    return dict(seed=seed+1, outcome=verdict['reference_recovery'])


def analyze_cut(payload):
    folder, number = payload
    folder = Path(folder)
    from .inputs import load_plan
    from mappa.mappa_core.utils.io.checkpoints import read_raw_cut, write_graph_checkpoint
    from mappa.mappa_core.core.matching.graphs import build_graph
    from mappa.mappa_core.core.matching.symmetry.workspace import ConditionedSymmetryWorkspace
    from mappa.mappa_core.core.search.branching.symmetry import finalize_graph_symmetry
    plan = load_plan(folder)
    ref = read(folder/'reference.json')
    path = folder/'cuts'/f'cut_{number:05d}.raw.pkl.gz'
    start, wall = time.process_time(), time.perf_counter()
    graph = read_raw_cut(path)
    target = build_graph(plan.problem.product.elements,plan.problem.product.wbo,bond_cut=plan.config.graph_floor)
    graph, stats = finalize_graph_symmetry(graph,target,iso_tolerance=plan.config.iso_tolerance,
        workspace=ConditionedSymmetryWorkspace(target,plan.config.iso_tolerance))
    finalization = dict(cpu_seconds=time.process_time()-start,wall_seconds=time.perf_counter()-wall,stats=stats)
    finalized = folder/'finalized'/f'{number:05d}.pkl.gz'
    finalized.parent.mkdir(parents=True,exist_ok=True)
    write_graph_checkpoint(graph,finalized)
    rows=[]
    for seed in range(len(graph.contexts)):
        out=folder/'analysis'/f'{number:05d}'/f'{seed:02d}'
        row = analyze_seed(plan, ref, graph, seed, out, number,
            raw_sha256=sha(path), finalized_sha256=sha(finalized))
        rows.append(row)
    save(folder/'analysis'/f'{number:05d}'/'complete.json',dict(cut=number,seeds=rows,finalization=finalization))
    del graph
    gc.collect()
    return dict(cut=number,seeds=rows)


def analyze(folder,workers):
    import multiprocessing as mp
    from .inputs import load_plan
    from mappa.mappa_core.core.search.planning import cut_sweep_items
    folder=Path(folder);plan=load_plan(folder)
    cuts=cut_sweep_items(plan.problem.reactant.wbo,plan.config.cut_floor)
    payloads=[(str(folder),n) for n in range(len(cuts))]
    start=time.perf_counter()
    with mp.get_context('fork').Pool(min(workers,len(cuts))) as pool:
        for row in pool.imap_unordered(analyze_cut,payloads,chunksize=1):
            save(folder/'analysis-progress'/f'{row["cut"]:05d}.json',row)
    save(folder/'analysis.json',dict(complete=True,cuts=len(cuts),wall_seconds=time.perf_counter()-start))
