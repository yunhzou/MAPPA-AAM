# High-level MAPPA pseudocode

Added Algorithm1 near the start of Methods, with a short definition of a live branch. It connects sweep conditions and seed orders to conditioned fragment growth, structural placement forks, compressed symmetry/constraints, optional bounded competition, unordered final grouping, and separate event decoding. Figure1a/b/c references explicitly connect growth, competition, and decoding. Isomorphism and group internals remain named operations, with their detailed definitions in the existing subsections.

Reviewed against the fragment loop in `alignment/branch.py`, the `match_fragment` wrapper, and the separate competition module. The outline retains no-growth branches for later seeds, repeats passes until no progress, records capped search status, keeps baseline branches when competition adds repairs, and restricts optional event decoding to complete families. The overview does not treat a single search mask as a one-event restriction. Bidirectional provenance remains orientation-specific.

No algorithm code changes or experiment reruns. Added the existing local `algorithmic.sty` to the source bundle. Tectonic fetched two missing font resources, then built successfully. All21 pages rendered; algorithm and Methods inspected individually and reflowed pages reviewed. Numerical/molecular/reference/overflow checks passed.

PDF SHA256: 512dc62108153525e0f50d3d7d06337ae3c6f3523f684a17937961d59f7d9be5
