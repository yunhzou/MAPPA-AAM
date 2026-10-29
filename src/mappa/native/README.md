# Native MAPPA implementations

This directory owns MAPPA's C++ code beside the Python package. The compiled
imports are `mappa._engine` and `mappa._group_ops`.

- `engine.cpp` and the fragment/growth headers implement compressed matching,
  growth, conditioned symmetry and reusable search work. Python calls these
  through `mappa_core/core/matching/native.py` and the native search adapters.
- `group_ops.cpp` implements permutation projection and conjugation used by
  the shared matcher and hierarchy types.
- `autgrp.cpp` implements the automorphism graph used by conditioned symmetry.
- `third_party/nauty/` holds the vendored dependency, separate from MAPPA code.

## Build from the repository root

Install the package and its build dependencies, then build the growth engine:

```bash
python -m pip install -e ".[dev]"
python -m pip install pybind11
python src/mappa/native/build_engine.py
```

Package installation builds `_group_ops` through `setup.py`.
The engine build compiles `engine.cpp`, `autgrp.cpp` and the nauty
core, then places `_engine.<abi>.so` in `src/mappa/`. It requires a C++17 compiler.
Build products live under `build/engine/` and are not committed.

`MAPPA_NATIVE=1` is the default. Missing or broken selected native execution is
an error. `MAPPA_NATIVE=0` explicitly selects Python growth. Custom node policies
and tracing use Python where those capabilities are not supported natively.

## Verify

`tests/test_native_engine.py` compares every growth call and the complete search
graph between Python and native engines. The routine smoke checks also cover
public workflows. A source relocation does not require a dataset benchmark.

## Third-party license

The vendored nauty 2.8.8 sources are by Brendan McKay and Adolfo Piperno.
See [Apache License 2.0](third_party/nauty/LICENSE-2.0.txt) and
[copyright notices](third_party/nauty/COPYRIGHT). These files are preserved
unchanged with the dependency.
