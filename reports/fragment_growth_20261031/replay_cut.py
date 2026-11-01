"""Measure installed growth policies on an archived cut, with bounded result reuse.

Use one fresh process per policy; compare complete raw and finalized graphs.
The earlier diagnostic-only runner is retained in Git at 30f0bc9.
"""
import argparse
from collections import Counter
from dataclasses import asdict, replace
import gzip
import hashlib
import json
from pathlib import Path
import pickle
import platform
import resource
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
POLICIES = ("incremental", "existence")


def balance_order(samples):
    """Give both execution orders equal weight, including for an odd sample count."""
    means={str(order):{policy:statistics.mean(row[policy] for i,row in enumerate(samples)
                                            if i%2==order) for policy in POLICIES}
           for order in (0,1)}
    balanced={policy:statistics.mean(row[policy] for row in means.values()) for policy in POLICIES}
    return dict(order_mean_cpu_seconds=means, balanced_cpu_seconds=balanced,
                balanced_ratio=balanced['existence']/balanced['incremental'],
                order_sample_counts=[len(samples[::2]),len(samples[1::2])])


def paired(input_path, output, pairs):
    """Alternate warm policies in one process to control CPU/process variation."""
    from mappa import _engine
    from mappa.mappa_core.core.search.session import SearchSession
    from mappa.mappa_core.core.search.branching.symmetry import finalize_graph_symmetry
    from mappa.mappa_core.core.matching.symmetry.workspace import ConditionedSymmetryWorkspace
    with gzip.open(input_path, "rb") as stream:
        saved=pickle.load(stream)
    source, task=saved['source'], saved['task']
    sessions={policy: SearchSession(source.problem,replace(source.config,growth_policy=policy),
                                   source.execution) for policy in POLICIES}
    symmetry={policy: ConditionedSymmetryWorkspace(session.context[2],source.config.iso_tolerance)
              for policy,session in sessions.items()}
    expected=None

    def run(policy):
        nonlocal expected
        session=sessions[policy]
        started=time.process_time()
        graph,_=session.search_cut(task.payload['cut'],seed_indices=task.payload['indices'])
        graph,_=finalize_graph_symmetry(graph,session.context[2],
            iso_tolerance=source.config.iso_tolerance,workspace=symmetry[policy])
        elapsed=time.process_time()-started
        if expected is None:
            expected=graph
        else:
            assert graph==expected
        return elapsed

    for policy in POLICIES:
        print(json.dumps(dict(phase='prepare',policy=policy)),flush=True)
        run(policy)
    samples=[]
    for pair in range(pairs):
        times={policy:run(policy) for policy in (POLICIES if pair%2==0 else POLICIES[::-1])}
        samples.append(times)
        print(json.dumps(dict(pair=pair,**times)),flush=True)
    record=dict(samples=samples,median_ratio=statistics.median(
        row['existence']/row['incremental'] for row in samples),
        medians={p:statistics.median(row[p] for row in samples) for p in POLICIES},
        caches={p:session.repair.stats() for p,session in sessions.items()},
        input_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),
        engine_sha256=hashlib.sha256(Path(_engine.__file__).read_bytes()).hexdigest(),
        identical_complete_graphs=True,
        scope='Same process and thread, independent prepared sessions; alternating warmed policy order. '
              'CPU time includes full cut search/export and finalization. Equality checks are outside timing.')
    record.update(balance_order(samples))
    output.write_text(json.dumps(record,indent=2)+'\n')


def replay(input_path, output, policy, repeats):
    from mappa import _engine
    from mappa.mappa_core.core.search.session import SearchSession
    from mappa.mappa_core.core.search.branching.symmetry import finalize_graph_symmetry
    from mappa.mappa_core.core.matching.symmetry.workspace import ConditionedSymmetryWorkspace

    with gzip.open(input_path, "rb") as stream:
        saved = pickle.load(stream)
    source, task = saved["source"], saved["task"]
    output.mkdir(parents=True, exist_ok=True)
    config = replace(source.config, growth_policy=policy)
    session = SearchSession(source.problem, config, source.execution)
    target = session.context[2]
    symmetry = ConditionedSymmetryWorkspace(target, source.config.iso_tolerance)
    engine_path = Path(_engine.__file__)
    record = dict(policy=policy, task=task.key, config=asdict(config),
        atoms=[len(source.problem.reactant.elements), len(source.problem.product.elements)],
        platform=platform.platform(), python=sys.version,
        revision=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        source_diff=subprocess.check_output(["git", "diff", "--", "src"], cwd=ROOT, text=True),
        input_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),
        engine_sha256=hashlib.sha256(engine_path.read_bytes()).hexdigest(), runs=[])

    for repetition in range(repeats):
        print(json.dumps(dict(phase="search_started", policy=policy, repetition=repetition)), flush=True)
        wall, cpu = time.perf_counter(), time.process_time()
        graph, counts = session.search_cut(task.payload["cut"], seed_indices=task.payload["indices"])
        search_wall, search_cpu = time.perf_counter()-wall, time.process_time()-cpu
        assert tuple(c.seed_order for c in graph.contexts) == task.payload["orders"]
        wall, cpu = time.perf_counter(), time.process_time()
        raw = graph
        graph, _ = finalize_graph_symmetry(raw, target,
            iso_tolerance=source.config.iso_tolerance, workspace=symmetry)
        finalization_wall, finalization_cpu = time.perf_counter()-wall, time.process_time()-cpu
        key_bytes, value_bytes = session.repair.target_view.graph.memo_storage_bytes()
        hits, misses, entries = session.repair.target_view.graph.memo_stats()
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if sys.platform != "darwin":
            peak *= 1024
        row = dict(repetition=repetition, search_wall_seconds=search_wall,
            search_cpu_seconds=search_cpu, finalization_wall_seconds=finalization_wall,
            finalization_cpu_seconds=finalization_cpu, peak_process_rss_bytes=peak,
            states=len(graph.states), transitions=len(graph.transitions),
            terminals=len(graph.terminals), search_counts=counts,
            repair=session.repair.stats(),
            certificate_cache=dict(key_bytes=key_bytes, value_bytes=value_bytes,
                                   entries=entries, hits=hits, misses=misses))
        # Save compressed records, never enumerate represented atom bijections.
        # Serialization and equality checking are outside search/RSS measurements.
        for name, value in (("raw", raw), ("graph", graph)):
            with gzip.open(output / f"{name}-{repetition}.pkl.gz", "wb") as stream:
                pickle.dump(value, stream, protocol=5)
        record["runs"].append(row)
        (output / "measurements.json").write_text(json.dumps(record, indent=2)+"\n")
        print(json.dumps(row), flush=True)
        del raw, graph, value


def compare(directories, output):
    """Check complete ordered records without expanding mapping families."""
    measurements = [json.loads((p / "measurements.json").read_text()) for p in directories]
    assert len({r["input_sha256"] for r in measurements}) == 1
    assert len({r["engine_sha256"] for r in measurements}) == 1
    expected = {}
    for stage in ("raw", "graph"):
        with gzip.open(directories[0] / f"{stage}-0.pkl.gz", "rb") as stream:
            expected[stage] = pickle.load(stream)
    comparisons = []
    for directory, measurement in zip(directories, measurements):
        for run in measurement["runs"]:
            for stage in ("raw", "graph"):
                with gzip.open(directory / f"{stage}-{run['repetition']}.pkl.gz", "rb") as stream:
                    graph = pickle.load(stream)
                fields = ("contexts", "roots", "states", "transitions", "stops")
                equality = {field: getattr(expected[stage], field) == getattr(graph, field) for field in fields}
                comparisons.append(dict(policy=measurement["policy"], repetition=run["repetition"],
                                        stage=stage, equal=all(equality.values()), fields=equality))
                del graph
    record = dict(measurements=measurements, comparisons=comparisons,
        exact_graph_equal=all(row["equal"] for row in comparisons),
        stops=dict(Counter(stop.reason for stop in expected["graph"].stops)),
        scope="Single archived cut, one compute thread per process. Complete native hierarchy, "
              "matching and symmetry finalization. No decoding or verification. Both policies use the same "
              "bounded complete-fragment cache. Peak RSS is cumulative "
              "process residency measured before saving each graph. Repeated-run peak can include "
              "earlier serialization. Equality checks include complete ordered compressed records, "
              "not expanded atom bijections.")
    record["summary"] = {}
    for measurement in measurements:
        record["summary"][measurement["policy"]] = [
            dict(repetition=run["repetition"],
                 total_cpu_seconds=run["search_cpu_seconds"]+run["finalization_cpu_seconds"],
                 total_wall_seconds=run["search_wall_seconds"]+run["finalization_wall_seconds"],
                 peak_process_gib=run["peak_process_rss_bytes"]/2**30,
                 certificate_cache_gib=(run["certificate_cache"]["key_bytes"]+
                                        run["certificate_cache"]["value_bytes"])/2**30)
            for run in measurement["runs"]]
    output.write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps(dict(exact_graph_equal=record["exact_graph_equal"], comparisons=comparisons)), flush=True)
    if not record["exact_graph_equal"]:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--policy", choices=POLICIES)
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--compare", type=Path, nargs="+", metavar="RUN_DIRECTORY")
    parser.add_argument("--paired", type=int, metavar="PAIRS")
    args = parser.parse_args()
    if args.paired is not None:
        if not args.input or not args.output or args.paired<2:
            parser.error("paired measurement requires --input, --output and at least two pairs")
        paired(args.input,args.output,args.paired)
    elif args.compare:
        if not args.output:
            parser.error("comparison requires --output")
        compare(args.compare, args.output)
    else:
        if not all((args.input, args.output, args.policy)):
            parser.error("replay requires --input, --output and --policy")
        replay(args.input, args.output, args.policy, args.repeats)


if __name__ == "__main__":
    main()
