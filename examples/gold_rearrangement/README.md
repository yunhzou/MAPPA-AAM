# Gold-catalyzed rearrangement: two oxygen-fate patterns

This 65-atom example tests whether an endpoint-only MAPPA search retains oxygen correspondences compatible with different mechanisms considered by González Pérez et al., *J. Org. Chem.* **2009**, 74, 2982–2991 ([DOI: 10.1021/jo802516k](https://pubs.acs.org/doi/10.1021/jo802516k)).

The endpoints are the supplied AuPPh₃ complexes **2 → 9**, including the catalyst and explicit hydrogens. They are not the catalyst-free structures 1 → 10. The 3D growth view shows all 65 atoms; the final candidate cards enlarge the organic heavy-atom core and follow the original epoxide oxygen, O*.

See the [original published pathway pictures](../../manuscript/animations/gold_rearrangement/pathway-reference/README.md) alongside the two oxygen-pattern interpretations.

## Result

The default search recovered both oxygen-fate patterns, without reference mappings, intermediate structures or anchors as input. Each displayed witness is a complete 65-atom bijection certified against a returned compressed family. Their canonical signed-event IDs differ.

| Reactant oxygen | Pattern A, consistent with route a | Pattern B, consistent with routes b/c |
|---|---|---|
| A: epoxide, atom 5 | Product ester-link O, atom 24 | Product ketone O, atom 5 |
| B: ester-link, atom 16 | Product ester carbonyl O, atom 26 | Product ester carbonyl O, atom 26 |
| C: carbonyl, atom 18 | Product ketone O, atom 5 | Product ester-link O, atom 24 |

Indices in this table are **one-based indices in the supplied XYZ files**. JSON and Python use zero-based indices. The source carbon framework is preserved, allowing permutations of the two equivalent methyl groups.

The paper examines three mechanisms. Route a starts with 1,3-ester migration and is favored energetically. Routes b and c start with 1,2-ester migration and oxirane activation, respectively, and converge at intermediate 14. The endpoint comparison distinguishes **two oxygen-fate patterns**, not three distinct endpoint mechanisms. The correspondences were inferred from the reported schemes and cross-checked by tracing the supplied stationary structures; they are not experimentally verified isotope labels. AAM recovery supports a pathway hypothesis, not its kinetic feasibility or the sequence of intermediate steps.

## Reproduce

Install MAPPA with its native extension as described in the repository README, then run:

```bash
python examples/gold_rearrangement/run_gold.py --output gold-output --workers 4
```

Settings are explicit: one seed, branch cap 100, uncut plus single-edge sweep, `iso_tolerance=1`, graph/cut floor 0.2, competition off. A 300-second watchdog wraps the run. Post-processing uses event thresholds 0.5 for nonmetal pairs and 0.3 for metal pairs. Search and selected-witness checks remain separate.

The recorded native run returned 1,796 deduplicated final fragment branches and 2,795 retained families across 69 contexts, with no incomplete returned paths. The search took about 1.0 seconds wall time with four workers; checkpoint writing and catalogue construction brought the recorded stage to about 1.9 seconds. These are example timings, not a benchmark or exhaustive decoding time. The branch cap was reached in parts of the search; completeness of all possible AAM solutions is not claimed.

The example selects carbon-framework-preserving **representatives already present** in the family catalogue, certifies full family membership, and computes their canonical event signatures. It does not enumerate all bijections or exhaustively decode the union. The selected patterns have 12 and 11 signed threshold events, respectively, including bond-order and metal-coordination changes. Those numbers are not counts of elementary reaction steps. The animation displays two selected event classes, not the full output count.

### Inputs and provenance

`data/reactant.json` and `data/product.json` contain fixed supplied DFT geometries and newly calculated **GFN2-xTB bond weights**, not bond weights reported by the original paper. The single points use tblite 0.7.0, charge +1 and a closed shell. Tiny negative numerical WBO values were clipped to zero after symmetrization; raw minima and counts are recorded. No geometry optimization was performed.

The original endpoint XYZ files and their SHA-256 hashes are included. Optional WBO regeneration:

```bash
OMP_NUM_THREADS=1 python examples/gold_rearrangement/recompute_wbo.py \
  examples/gold_rearrangement/data/2PPh3.xyz /tmp/reactant-wbo.json
```

`selected-families.json` stores the two film families and full witness checks. Both film witnesses are siblings in cut context 21, with a shared 63-atom prefix. The standalone recovery script can select a different representative with the same event class. The film includes verified atom-level replay; the linked detailed trajectory viewer retains all four branches in this cut. `oxygen-provenance-audit.json` records the independent stationary-structure tracing; this was performed after the search and used only for interpretation. Intermediate geometry files are not required to reproduce MAPPA's recovery. Selected pathway schemes are reproduced on the linked reference page with attribution; the complete source PDF and full supplied archive are not redistributed.

[Animations and offline viewer](../../manuscript/animations/gold_rearrangement/README.md)
