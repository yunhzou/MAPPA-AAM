# Native MAPPA implementations

This directory owns MAPPA's C++ code beside the Python package. The compiled
imports are `mappa._engine` and `mappa._group_ops`.

- `engine.cpp` and the fragment/growth headers implement compressed matching,
  growth, fragment-local relaxed relationships and reusable search work. Python calls these
  through `mappa_core/core/matching/growth/native.py` and the native search adapters.
- `fragment_query.h` owns first-witness matching and witness extension.
  At saturation it discovers one representative per certified target-set orbit,
  using compact orbit exclusions without expanding equivalent sets first.
  `fragment_completion.h` retains the complete discovery constraints and
  prepares anchored relaxed self-relationships for each side. Global placement
  preparation is shared per target reaction. Its certified exchanges preserve
  the full target WBO matrix, active topology and source-node compatibility,
  while fixing prior locks and explicitly fixed target labels. Each result
  retains this placement relationship separately from its fragment-local
  relationships; saved hierarchy transitions transport the entire later
  continuation jointly. This prepared endpoint record is separate from the
  bounded completed-result cache.
  Frontier selection stays in `engine.cpp`; neither module owns hierarchy caps
  or work-pool scheduling.
- `find_fragment_witness` binds `FragmentQuery.find` directly for a supplied
  accepted fragment and its locks. It returns one witness or proven absence,
  together with the visited search-state count. It does not perform full
  compressed completion or select a different solver after failure.
- `group_ops.cpp` implements permutation projection and conjugation used by
  the shared matcher and hierarchy types. It also compiles shared decoder actions
  and projects ordered two-endpoint orbits without Python traversal or complete
  permutation enumeration.
- `set_orbit.h` supplies compact target-set transport and pointwise stabilizers.
  It keeps correlated generator supports together without materializing the
  complete permutation group or target-set orbit.
- `autgrp.h` shares the binary structural graph solver between native fragment
  preparation and the public bindings in `autgrp.cpp`.
- `third_party/nauty/` holds the vendored dependency, separate from MAPPA code.

## Build from the repository root

Install the package and its build dependencies, then build the growth engine:

```bash
python -m pip install -e ".[dev]"
python -m pip install pybind11
python src/mappa/mappa_core/native/build_engine.py
```

Package installation builds `_group_ops` through `setup.py`.
The engine build compiles `engine.cpp`, `autgrp.cpp` and the nauty
core, then places `_engine.<abi>.so` in `src/mappa/`. It requires a C++17 compiler.
Intermediate build products use a temporary directory (`TMPDIR` when set)
and are removed after the extension is published.

`MAPPA_NATIVE=1` is the default. Missing or broken selected native execution is
an error. `MAPPA_NATIVE=0` explicitly selects Python growth. Custom node policies
and tracing require explicit Python selection; native requests reject them.
The trajectory viewer calls Python tracing explicitly. Growth uses
one-witness feasibility followed by exact compressed completion. There is no
growth-mode selector or candidate-prefix cache.

## Verify

`tests/test_native_engine.py` compares every growth call and the complete search
graph between Python and native engines. The routine smoke checks also cover
public workflows. A source relocation does not require a dataset benchmark.

## Third-party license

The vendored nauty 2.8.8 sources are by Brendan McKay and Adolfo Piperno.
See [Apache License 2.0](third_party/nauty/LICENSE-2.0.txt) and
[copyright notices](third_party/nauty/COPYRIGHT). These files are preserved
unchanged with the dependency.
