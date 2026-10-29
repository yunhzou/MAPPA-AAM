"""Command-line composition of the typed mappa API."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from mappa.mappa_core.core.search.api import search_aam
from mappa.mappa_core.utils.io.records import decoded_events_record, reaction_from_record
from mappa.presentation.viewers import write_aam_bundle
from mappa.mappa_workflow.ts_analytics.serialization import write_ts_record
from mappa.mappa_core.utils.io import load_cached_xtb
from mappa.mappa_core.domain import AAMProblem, AAMSearchConfig, MolecularEndpoint, VibrationalModes
from mappa.mappa_workflow.ts_analytics.records import TransitionStateTarget
from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig
from mappa.mappa_core.post_processing.decoder.api import decode_events
from mappa.mappa_workflow.ts_analytics.api import analyze_transition_state


def _endpoint_npz(path, label):
    with np.load(path, allow_pickle=False) as data:
        return MolecularEndpoint(
            tuple(str(item) for item in data["elements"]),
            data["coordinates"], data["wbo"], label=label)


def _endpoint_cache(path, label):
    elements, coordinates, wbo, _xyz = load_cached_xtb(path)
    return MolecularEndpoint(tuple(elements), coordinates, wbo, label=label)


def _endpoint(args, prefix, label):
    npz = getattr(args, f"{prefix}_npz")
    cache = getattr(args, f"{prefix}_cache")
    if bool(npz) == bool(cache):
        raise SystemExit(
            f"provide exactly one --{prefix}-npz or --{prefix}-cache")
    return _endpoint_npz(npz, label) if npz else _endpoint_cache(cache, label)


def _config(args):
    anchors = tuple(tuple(map(int, value.split(":"))) for value in args.anchor)
    return AAMSearchConfig(
        cut_floor=args.cut_floor,
        graph_floor=args.graph_floor,
        iso_tolerance=args.iso_tolerance,
        event_threshold=args.event_threshold,
        metal_event_threshold=args.metal_event_threshold,
        seed_count=args.seed_count,
        branch_limit=args.branch_limit,
        anchors=anchors,
    )


def _parser():
    parser = argparse.ArgumentParser(prog="mappa")
    parser.add_argument("--reactant-npz")
    parser.add_argument("--product-npz")
    parser.add_argument("--reactant-cache")
    parser.add_argument("--product-cache")
    parser.add_argument("--target-npz", action="append", default=[])
    parser.add_argument("--stage", choices=("aam", "events", "ts"),
                        default="events")
    parser.add_argument("--reaction-json",
                        help="explicit resolved-reaction record required by --stage ts")
    parser.add_argument("--name", default="reaction")
    parser.add_argument("--output", default="mappa_output")
    parser.add_argument("--workers", type=int, default=1,
                        help="CPU processes used by the cut sweep")
    parser.add_argument("--max-events", type=int)
    parser.add_argument("--cut-floor", type=float, default=0.2)
    parser.add_argument("--graph-floor", type=float, default=0.2)
    parser.add_argument("--iso-tolerance", type=float, default=1.0)
    parser.add_argument("--event-threshold", type=float, default=0.5)
    parser.add_argument("--metal-event-threshold", type=float, default=0.3)
    parser.add_argument("--seed-count", type=int, default=1)
    parser.add_argument("--branch-limit", type=int, default=100)
    parser.add_argument("--anchor", action="append", default=[], metavar="R:P")
    return parser


def main(argv=None):
    parser = _parser()
    args = parser.parse_args(argv)
    if args.stage == "ts":
        if not args.reaction_json or not args.target_npz:
            parser.error("--stage ts requires --reaction-json and --target-npz")
    elif args.reaction_json or args.target_npz:
        parser.error("--reaction-json and --target-npz require --stage ts")
    reactant = _endpoint(args, "reactant", "R")
    product = _endpoint(args, "product", "P")
    config = _config(args)
    problem = AAMProblem(reactant, product, name=args.name)
    output = Path(args.output)
    if args.stage == "ts":
        reaction = reaction_from_record(
            json.loads(Path(args.reaction_json).read_text()),
            problem, config)
        output.mkdir(parents=True, exist_ok=True)
    else:
        aam = search_aam(
            problem, config,
            workers=max(1, args.workers),
            intermediate_dir=output / 'aam_search')
        output = write_aam_bundle(aam, output)
        if args.stage == "events":
            decoded = decode_events(aam, EventDecodeConfig(
                threshold=args.event_threshold,
                metal_threshold=args.metal_event_threshold,
                max_events=args.max_events,
                branch_limit=args.branch_limit))
            (output / "events.json").write_text(
                json.dumps(decoded_events_record(decoded), indent=2))
        print(output.resolve())
        return
    for index, path in enumerate(args.target_npz, 1):
        with np.load(path, allow_pickle=False) as data:
            target = TransitionStateTarget(
                MolecularEndpoint(
                    tuple(str(item) for item in data["elements"]),
                    data["coordinates"], data["wbo"],
                    label=Path(path).stem),
                VibrationalModes(data["frequencies"], data["modes"]))
        ts = analyze_transition_state(reaction, target, search_config=config)
        write_ts_record(ts, output / f"ts_{index:03d}.json")
    print(output.resolve())


if __name__ == "__main__":
    main()
