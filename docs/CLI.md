# Command line

The `mappa` command runs one reaction in one search direction. For both
directions, reaction cores or many reactions, use the
[Python API](PYTHON_API.md) or the [work pool](shared-execution.md).

## Matching and decoding

The default stage is signed-event decoding. NPZ endpoints contain `elements`,
`coordinates`, and `wbo` arrays:

```bash
mappa --stage events \
  --reactant-npz R.npz --product-npz P.npz \
  --workers 4 --output mapping
```

Use `--stage aam` for raw search only. Both stages write `aam.json` and the
self-contained `search.html` viewer. The events stage also writes `events.json`
with the original catalogue, shared conditional responses and approximate model
minimum/window. This retained scope does not prove that capped search explored
every mapping or that the event minimum is physical.
Existing xTB caches can be supplied with `--reactant-cache` and `--product-cache`.

Search settings are `--iso-tolerance`, `--seed-count`, `--branch-limit`,
`--cut-floor`, `--graph-floor` and repeated `--anchor R:P`. Event settings are
`--event-threshold` and `--metal-event-threshold`.

## Transition-state scoring

TS processing requires an explicitly supplied resolved-reaction record:

```bash
mappa --stage ts \
  --reactant-npz R.npz --product-npz P.npz \
  --reaction-json selected_reaction.json --target-npz guess.npz \
  --output ts_scores
```

Target NPZ files additionally contain `frequencies` and `modes`. Repeat
`--target-npz` to score several targets; each writes `ts_001.json`,
`ts_002.json` and so on. See the [Python API guide](PYTHON_API.md) for supplying
a `ReactionContext` and serializing it with `reaction_record`, and the
[TS scoring guide](TS_SCORING.md) for the score itself.

## Native engine

Build the engine with `python src/mappa/mappa_core/native/build_engine.py`; see the
[native build guide](../src/mappa/mappa_core/native/README.md).

`MAPPA_NATIVE=1` is the default native selection. A missing or broken native
extension raises an error. Set `MAPPA_NATIVE=0` to explicitly select Python.
Custom node policies and growth tracing require explicit Python selection;
native requests reject these unsupported capabilities. The trajectory viewer
calls Python tracing explicitly.
