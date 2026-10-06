# Package-rename compatibility fixtures

These small trusted files were produced by `rxn_core` at commit `d22f353`, before renaming the package to `mappa`. They contain the four-atom search used by `toy_archive` in `tests/test_search_trajectory.py`: permuted identical C/O/O/H endpoints, one seed, a branch limit of 32, matching tolerance 0.1 and graph floor 0.4.

The full checkpoint, finalized-cut checkpoint and raw-cut checkpoint retain original `rxn_core.*` pickle globals. The JSON archive and decoded event IDs establish the expected scientific output. Tests load these files through the public readers and verify that classes resolve to `mappa`, graph contents and event IDs are unchanged, and new checkpoints round-trip. These fixtures are project-generated; pickle readers remain restricted to trusted inputs.
