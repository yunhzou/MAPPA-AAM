# Multi-candidate discussion and Figure 4

The paper now distinguishes alternative event patterns, residual symmetry, and potential downstream geometric uses. It does not equate distinct mappings with distinct experimentally established pathways or claim that MAPPA/SLAP exhaust all mappings.

## Evidence

- `multicandidate_example.json`: current public Python API, one worker, six explicit atoms, one seed with cut sweep, branch cap 2,000, iso tolerance 1, final event thresholds 0.5/0.3. Full decoding of seven saved families returns two event classes (four and six changes). Both fixed figure witnesses are certified members. H3/H4 exchange preserves class A; class B is forbidden under the class-A restriction and allowed without that restriction. Total search/decoding/query check: 0.5264 CPU seconds. This is a functional example, not a timing benchmark.
- `output_multiplicity.json`: all 12,957 archived comparator records, read without rerunning any mapper. RXNMapper/LocalMapper/Indigo/Chython/RDT return at most one explicit record; binary SLAP returns several in 522 reactions and weighted SLAP in 417. Explicit records are not symmetry-expanded bijections or unique event classes.
- The runtime used the already built package in `work/python-api-audit/src`. Every Python file in that package was compared byte for byte with the current repository source before the run; all matched. The optional compiled group operations were available there. The bare source checkout has no built `_group_ops` extension; no source changes were needed.
- The figure builder independently verifies the signed events and all depicted endpoint bonds against the saved input matrices. It does not rerun AAM.

## Literature and scope

Ali, Mizuno, Akiyama, Nagata, and Komatsuzaki, *Enumeration Approach to Atom-to-Atom Mapping Accelerated by Ising Computing*, JCIM 65(4), 1901–1910 (2025), https://doi.org/10.1021/acs.jcim.4c01871, is a counterexample to global exclusivity. It enumerates optimal mappings and clusters them by molecular symmetry. It was not added as an unmeasured performance row.

De, Krummenacher, Schaefer, and Goedecker, *Finding Reaction Pathways with Optimal Atomic Index Mappings*, PRL 123, 206102 (2019), https://doi.org/10.1103/PhysRevLett.123.206102 (open preprint https://arxiv.org/abs/1906.06077), supports the relevance of atomic-index choices to pathway searches. This does not establish transition-state or barrier improvements for MAPPA.

`audit_outputs.py` reproduces the small example using an installed/built package plus the checked-in notebook helper and reads the archived comparator files. Use `--repo` to select a checkout; optionally supply `--package-src` for a built source tree and `--archive-dir` for the rescored outputs. The expensive benchmark campaigns are not rerun.
