# Critical fragment matching cases

Run `python -m pytest -q -s tests/test_fragment_matching_contract.py`.
These acceptance tests check one-witness growth, direct compressed completion,
and exact saved mapping relations in both Python and native execution.

- **C6 → C6:** one compressed candidate representing exactly 12 rotations and
  reflections. Construction must produce one final candidate, without replaying
  smaller candidate fragments or calling candidate deduplication.
- **C6 → two disconnected C6 rings:** two candidates, each representing exactly
  the 12 mappings into its own component. Existence still returns only one witness.
- **Competing destinations:** twelve carbon atoms require eleven reachable carbon
  destinations, while a separate pool leaves the total target count large enough.
  Element and ordinary degree checks pass. Reject before assignment branching;
  removing one competing atom must yield a verified witness through the same solver.
- **14! tolerance pool:** different exact bond weights remove target symmetry,
  while the accepted tolerance still permits every permutation. Completion must
  return one pool without enumerating those permutations.
- **Free hydrogen pools:** five eligible hydrogens in equivalent pairs of H–C
  and H–O components plus H–N form three choices, not five singleton branches.
  Component identity remains part of global symmetry, not a restriction on a
  valid local free-atom pool.

`enzymemap_16175.json` contains exact sparse endpoint WBOs and two queries from
reaction 16175, reverse direction, cut 0, first seed 345, random seed 42, tolerance
1.0. The fixture records the archived input SHA-256 and its source revision.
It contains no worker state, caches, pickle objects or cluster paths.

- **73 atoms, match exists:** the previous Python solver needed 114,795 search
  states and 7.75 s on the Mac to find one witness. This was slow positive matching,
  not evidence of absence. The saved witness is checked directly against every
  accepted constraint. The work budget is at most twice the number of source atoms.
- **317 atoms, no match:** adding edge 134–135 to the 316-atom matched predecessor
  makes this query impossible. The saved predecessor witness is checked directly.
  An accepted 28-atom subtree has no element-compatible subgraph monomorphism even
  ignoring bond weights; NetworkX verifies that independent absence proof. The
  production solver must reject at its initial constraint propagation, before
  branching on assignments. This particular negative query was already quick in
  the previous solver; it is not presented as the cause of the archived timeout.

As a negative control, the archived Python solver was also run on the competing
destinations case: it visited 153,711 search states before a 2 s diagnostic stop.
The current solver rejects it in one state. The archived implementation remains
in Git, not in the test fixtures or an alternate production execution path.

Both Python and native execution use their existing `FragmentQuery.find` kernel
for positive and negative queries. The native binding `find_fragment_witness`
exposes that kernel directly; it does not introduce another solver. Query tests
run in isolated processes: 10 s external deadline, 2 s kernel ceiling and 512 MiB
peak process RSS. Timeouts and crashes are failures, never no-match results.
Search-state limits catch lost pruning without relying only on wall time.

The ring checks validate the complete represented relation against an independent
rotation/reflection formula. Python tests observe actual candidate construction;
native completion constructs candidates only at terminal emission. Both report
one or two final emissions, with no prefix candidates or candidate deduplication.
Further checks cover growth-order independence, internal tolerance pools, custom
policies, subsequent locks, and serialization of the exact conditioned group.
