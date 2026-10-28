# Alignment Architecture

This document describes the module boundaries for the WBO graph alignment
implementation. The main rule is that each layer owns one concept and exports
only the abstraction needed by the layer above it.

Implemented interface: [AAM search graph API](AAM_SEARCH_GRAPH_API.md).
The [original design analysis](AAM_SEARCH_GRAPH_DESIGN.md) records the rationale.
The [retro detection and assembly contract](RETRO_ASSEMBLY.md) describes the
separate matching, occupation, exact assembly, ranking, and viewer layers.

## Public Surface

`mappa.__init__` is the convenience API used by scripts and the pipeline.
It re-exports stable entry points such as:

- `match_fragment`
- `search_aam`
- `group_mechanisms`
- `compile_mapping_families` / `compile_mechanism_families`
- `mappa.postprocessing.decode_events`
- `build_graph`
- `classify_bonds`

Code outside the package should prefer these top-level imports unless it is
debugging internals.

## Typed AAM API and C++ execution boundary

`aam.search_aam(problem: AAMProblem, config: AAMSearchConfig) -> AAMResult`
is the typed Python AAM interface. `domain.py` defines the input, configuration,
mechanism, result, and metrics objects. `search_graph.py` owns states, fragment
transitions, paths, and stops. A live `_Branch` owns only cumulative frontier
state and a graph-node reference, not copied histories. Its ancestry is a DAG.
Optional mechanisms retain Python `AAMBranch` projections with an `AAMHierarchy`.

The C++ backend replaces computation inside `growth.island.grow_island`:
`growth/native.py` prepares native graph views, calls `_engine.grow_island`,
restores original atom indices, and returns Python `_IsoResult` objects.
The kernel owns its internal candidate state, growth loop, extension, and
deduplication; Python owns cut/seed orchestration and completed AAM objects.
The Python growth implementation remains available as the reference engine.

`fragment.match_fragment` is the common Python contract above this boundary;
both the AAM scheduler and fragment detection use it. Its separate
`_native.paired_mapping_invariant` kernel accelerates occupation-family
grouping. Augmentation, assembly, ranking, and visualization remain Python.
See [the native backend guide](../native/README.md) for pseudocode, build
instructions, and the distinction between active and currently unused kernels.

## Molecule Alignment: `mappa.alignment`

Purpose: turn two WBO matrices into scored molecule-level alignments.

Files:

- `alignment/api.py`: public dataclasses and entry points:
  `MatchCandidate`, `MatchResult`, `match_wbo_graphs`,
  `align_from_arrays`, `analyze_alignment`, `cut_edges_above_floor`.
- `alignment/branch.py`: molecule-level branch scheduling, branch dedupe,
  seed-order generation, and explicit symmetry repair.
- `alignment/sweep.py`: R-P sweep-cut mechanism discovery via no-cut plus
  one-edge R cuts, with symmetry-canonical mechanism signatures.  Parallel
  cut-sweep work is grouped by cut, so each worker computes the cut-specific R
  orbit map once and reuses it for all seed orders; product orbits are
  invariant and are computed once in the worker initializer.
- `alignment/ts_core.py`: mechanism-local endpoint->TS/IG core mapping used by
  the ranker after R-P mechanisms are known. The BGCP pipeline runs it from
  both R and P, pulls P-derived mappings back through the R-P witness, and
  scores the union in R-core indexing.

This layer owns whole-molecule decisions: seed orders, branch lifecycle,
mechanism discovery, mechanism-local TS core matching, and final representative
choice. It does not implement single-fragment growth or symmetry-block witness
search.

## Fragment Growth: `mappa.growth`

Purpose: grow one connected R fragment against a target graph.

Files:

- `growth/island.py`: the live island-growth loop. It initializes seed
  candidates, pops frontier edges, calls the matcher extension rule, commits
  successful atoms, records deferred boundary edges, and locks the island only
  after heap saturation.
- `growth/frontier.py`: heap/frontier utilities and final uniqueness checks.
- `growth/result.py`: `_IsoResult`, the locked fragment result.
- `growth/trace.py`: diagnostic event formatting and failure explanations.

This layer owns traversal order and locking. It does not decide whether a
candidate target atom is valid; validity is delegated to `mappa.matcher`.

## Symmetry Matching: `mappa.matcher`

Purpose: represent and extend symmetry-compressed partial mappings.

Files:

- `matcher/state.py`: `_SymCand`, `_SymBlock`, witness/materialization helpers.
- `matcher/extend.py`: one-atom extension for a list of compressed candidates.
- `matcher/support.py`: witness search inside unresolved symmetry blocks.
- `matcher/dedupe.py`: orbit plus one-hop-boundary dedupe signatures.
- `matcher/orbits.py`: exact pynauty orbit grouping on 0.2-tolerance
  WBO-colored graphs, plus the old WL/color-refinement helper for fallback
  and comparison.
- `matcher/primitives.py`: WBO access and tolerance helpers.
- `matcher/chemistry.py`: post-hoc chemistry-relevant symmetry expansion.

This layer owns the central compression model:

```text
_SymCand = witness mapping + interchangeable/correlated _SymBlock pools
```

The witness is one deterministic assignment. The blocks encode the actual
symmetry state. A block with two R atoms and four P atoms represents
`P(4, 2)` injective concrete assignments without enumerating them.

## Chemistry And IO Helpers

- `chemistry_computations/`: file and external-computation utilities.
  - `xyz.py`: XYZ parsing and formatting, including extended XYZ mode output.
  - `xtb.py`: cached xtb single-point/Hessian execution and WBO cache loading.
  - `frames.py`: coordinate reindexing into the R atom frame.
- `frag.py`: WBO graph construction, `classify_bonds`, and local witness
  materialization. It re-exports old XYZ/xtb names only for compatibility.
- `modes.py`: normal-mode parsing and score-feature vectors.
- `cli.py`: explicit `aam`, `events` (default), and `ts` stages. Search and
  event decoding never call a geometry selector. TS processing requires a
  caller-supplied resolved-reaction record.
- `alignment/mapping_families.py`: exact structural and event-conditioned
  family compilation and membership. It performs no 3D witness selection.

The legacy R/P selector, its public wrappers, and its CLI modes are removed.
There are no compatibility aliases or automatic replacements.

These modules should not depend on matcher internals unless they are explicitly
running alignment.

## Dependency Direction

The intended direction is:

```text
mappa.__init__
  -> alignment
      -> growth
          -> matcher
      -> frag
  -> chemistry_computations / modes / align / pipeline helpers
```

Lower layers must not import higher layers. In particular:

- `matcher` must not import `growth` or `alignment`.
- `growth` must not import `alignment`.
- `alignment` may import `growth`, `matcher`, `frag`, and
  `chemistry_computations`.
- scripts may import public top-level names, or internal packages only when
  producing debug traces.

## Naming

The old `pq.py` module and related compatibility facades were removed. The
priority queue is now an implementation detail of `growth/island.py`, not a
package-level abstraction.
