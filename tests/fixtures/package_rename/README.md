# Package-rename compatibility fixtures

These small trusted files were produced by `rxn_core` at commit `d22f353`, before renaming the package to `mappa`. They contain the four-atom search used by `toy_archive` in `tests/test_search_trajectory.py`: permuted identical C/O/O/H endpoints, one seed, a branch limit of 32, matching tolerance 0.1 and graph floor 0.4.

The full checkpoint, finalized-cut checkpoint and raw-cut checkpoint retain original `rxn_core.*` pickle globals. Their fragment schema predates the two-sided relaxed relationship. Tests require current scientific readers to reject these files explicitly; namespace-only relocation is checked separately. Current checkpoints are round-tripped by the typed API and transport tests. These fixtures are project-generated; pickle readers remain restricted to trusted inputs.
