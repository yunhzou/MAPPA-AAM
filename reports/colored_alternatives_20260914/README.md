# Colored bond changes and a real multi-candidate example

Figure 4 replaces the constructed methanol example with all three saved minimum-event MAPPA classes for coordinate case 127 (23 explicit atoms). Every class has four signed events, yet the event patterns differ. The final SLAP-sweep comparison includes all three; native SLAP includes one. This is an illustration from the mapping-unverified 140-case development collection, not an independent accuracy result or mechanism assignment.

`coordinate_case127.json` preserves the original WBO matrices, all three full witness mappings, their supporting family constraints and source fragment groups, the certified same-class methyl H10/H11 swap, and source hashes. The three witnesses were checked by exact membership queries against their saved families. No case-127 AAM search or event decoding was repeated. Read-only extraction and membership checks took 1.18 CPU seconds. Comparator class identifiers were cross-checked against the final `coordinate-slap-comparison.json.gz`.

Colors use separate visual channels: pastel regions track actual saved fragments within each candidate; red dashed bond strokes indicate negative WBO events, and blue strokes indicate positive events. Event colors and counts derive independently from original WBO values (four signed events per mapping); the drawing uses conventional single/double bonds only for readability. Hydrogens are included in all calculations and omitted from the top panel for readability. The lower panels use a certified methyl-H symmetry shuffle from candidate A and the retained ethyl fragment from candidate B.

Figure 1's competing-fragment and decoded-outcome drawings now follow the same bond-color convention. Numerical benchmark results and algorithm code are unchanged.

Rebuild figures with `manuscript/scripts/build_figures.py`; `build_alternatives_figure.py` validates the saved-case event lists and depiction bonds. `extract_case.py` records the extraction from the local archived campaign; it does not run search or exhaustive decoding.
