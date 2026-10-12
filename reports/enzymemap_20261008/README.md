# EnzymeMap: completed decoding and reference-miss investigation

**All 47,974 reactions have completed processing under the combined recorded configurations and output policies.** The follow-up resolved all 108 previously unfinished reactions and classified every one of the original 807 strict-reference misses. The final evidence review closed on 12 October 2026 UTC.

Processing completion and reference agreement answer different questions. The completed output union covers **47,335 / 47,974 references (98.6680%)**. The **100% completion** figure does not mean 100% reference agreement, full atom coverage, or exhaustive exploration of all possible mappings.

## Results and denominators

| Measure | Recorded result | Interpretation |
|---|---:|---|
| Completed processing after this investigation | **47,974 / 47,974 (100%)** | Every row has completed decoding under a recorded output policy |
| Previously unfinished reactions resolved | **108 / 108** | Independently audited continuations, including partial mappings and preparation fixes |
| Strict decoded-reference agreement after follow-up | **47,335 / 47,974 (98.6680%)** | Any returned candidate matches; union across recorded configurations |
| Additional unique strict-reference matches | **168** | Added to the original 47,167; repeated configurations are not additional reactions |
| Remaining strict-reference differences | **639** | All classified; none is an unfinished computation |

The preceding benchmark remains unchanged: **47,866 complete decodes**, **47,167 strict decoded-reference hits**, and **47,351 raw-family reference hits**, each over **47,974** rows. Its report and case outcomes are retained byte-for-byte. The earlier baseline/retry report is also a historical record, not the final completion result.

The original 108 unfinished rows and the original 807 strict-reference misses are different selections. Resolving 108 unfinished rows and recovering 168 additional references are therefore compatible results.

## What was measured

MAPPA searches for compressed families of atom correspondences, then decodes those families into signed bond-event classes with explicit representative mappings. This investigation distinguishes:

- **Processing completion:** the recorded search and decoding work is complete for the selected policy. Completeness is relative to the retained families and declared event window.
- **Raw-family reference inclusion:** an exact membership check establishes that a reference-equivalent mapping exists in a retained symbolic family. This need not be the representative returned by decoding.
- **Strict decoded-reference agreement:** at least one actual returned candidate matches the annotated heavy-atom relation, modulo endpoint chemical symmetry. The comparison preserves charge, isotope, CIP and bond-stereo distinctions and explicitly represents unmatched atoms.
- **Event-class agreement:** a returned bond-event class matches even when its chosen representative fails the finer strict mapping comparison.

Candidate-set agreement is not top-1 accuracy. Candidate ordering is not a chemical probability ranking, and reference agreement alone does not establish a physical reaction pathway. Failures were retained in the full denominator throughout the investigation.

## Protocol and changes made

The starting engine is commit [`ea6a7a3`](https://github.com/yunhzou/MAPPA-AAM/commit/ea6a7a3d6bfb4361d582dbc60dcfb87f183d1f8e). The original run used forward search, formal bond orders, explicit hydrogens, iso tolerance 1, one seed with random seed 42, branch cap 100 and no balance repair. The follow-up investigated recorded caps 100, 200 and 1,000 and used minimum-plus-five decoding where specified. Original input molecules and reference mappings were preserved.

The final totals combine historical completed outputs and versioned follow-ups; they are **not a fresh uniform rerun at one cap**. Individual manifests and the case ledger identify each configuration. Full decoding does not make a branch-limited search globally exhaustive.

### Completing the unfinished outputs

The 108 originally unfinished rows comprised 35 incomplete family decodes, 40 decoding timeouts, six search timeouts, 16 explicit-hydrogen composition mismatches and 11 input/reference preparation mismatches. All are now resolved in the completion reconciliation.

Completed families were retained and reused through verified archive handoffs. Long-running cases were diagnosed before continuation. Parallel family work and exact alternative solver encodings addressed difficult residual queries. The last case, index 28666, completed with ten candidates, proven minimum eight, and event window thirteen; its strict reference occurs among the returned candidates.

For the final difficult families, a source-indexed free-pool permutation encoding and endpoint symmetry reduction removed redundant permutation work while preserving physical family membership. Validation included exhaustive small equivalence tests, fresh control outputs, physical-witness membership checks, residual-query reconstruction and independent UNSAT replay. A solver timeout or unknown result was never promoted to a proof of absence.

### Partial mappings: 16 reactions

The full-bijection decoding path was inappropriate for these composition-mismatched inputs. The follow-up used explicit partial mappings with the Golden coverage priority: maximize mapped heavy atoms, then mapped atoms, then enumerate the applicable minimum-plus-five event window. Unmatched atoms remain explicit. **All 16 completed and recovered the strict heavy reference**, without adding atoms or altering molecules. See the partial-output audit.

### Independent input/reference preparation: 11 reactions

An atom-index consistency assumption failed when input and reference graphs were prepared independently. The corrected preparation path keeps each graph's own coordinates and features and compares their resulting mappings appropriately. **All 11 completed and recovered the strict reference**. Neither endpoint molecule nor the reference was rewritten. See the independent-preparation audit.

## Explaining the original 807 strict-reference misses

| Final classification | Reactions |
|---|---:|
| Reference recovered by completed follow-up decoding | 141 |
| Reference recovered after independent input/reference preparation | 11 |
| Reference recovered using explicit partial mappings | 16 |
| Reference absent from the completed recorded search configurations | 472 |
| Reference event class returned, but representative differs under stereo-sensitive scoring | 166 |
| Reference event class returned, but representative differs under charge-sensitive scoring | 1 |
| **Total** | **807** |

The first three categories recover **168** references. The remaining **639** comprise **472 scoped search absences** and **167 representative distinctions**. The 166 stereo-sensitive cases include index 6511; the charge-sensitive case is index 12009, involving Fe²⁺/Fe³⁺ with a zero-event class. See the completed witness diagnoses, earlier representative analysis, and charge check.

The machine-readable ledger retains the historical categories: 166 earlier representative distinctions plus one newly diagnosed output distinction for case 6511. The table above groups those by the final scientific explanation; it does not change their evidence records.

An absence establishes that the reference is outside the completed recorded family union for the tested settings. It is not a claim that the mapping is impossible under every cap, seed or pathway. These results neither require forcing all outputs to equal the annotation nor establish general algorithm correctness.

## Validation

The publication is gated by:

1. **Full decoder evidence audit:** 254 declared task rows, comprising **226 validated** configurations and **28 superseded** attempts. Every required current decoder task has validated evidence. These are task rows, not unique reactions.
2. **Exact-reference proof audit:** **211 / 211 task configurations** passed, with 26 positive and 185 negative task outcomes. These include controls and repeated caps, not 211 additional reactions or the entire 807-case set.
3. **Complete miss ledger:** all **807 reactions** classified; all **464 required search, decoding and additional-output tasks** complete, with no unresolved categories.
4. **Completion reconciliation:** all **108 missing rows** have audited completed outputs; union with the original 47,866 covers exactly indices 0–47,973.
5. **Local evidence preservation:** all **131,552 packaged members** covering the 211 exact-reference proof tasks were checked against their original archive manifests. The final 28666 full archive was also verified locally.

The final accounting check completed successfully after correcting a reporting assumption for case 27963: its reference is proven absent, so a positive reference witness cannot be reconstructed. The corrected check binds all 529 negative cut proofs to the exact decoded cut hashes and matching input/reference records. It did not change a mapping or rerun decoding.

Operational retries and cache reuse make this campaign unsuitable for claiming a controlled cold-start speedup. No new timing comparison is claimed here.

## Publication scope and evidence preservation

This GitHub publication contains the methods and results write-up plus [machine-readable aggregate results](results.json). The repository README and benchmark index link here.

The full case-level dataset, 807-case ledger, 108-row reconciliation, proof/query records, audit manifests and research source snapshot remain preserved in the local evidence package and original cluster archives. They are **not included in this report-only GitHub commit**. The integrity check run before publication recomputed the totals from those records and verified all 131,552 packaged proof-record members against their original audited archive manifests.

The historical baseline reports are preserved locally without rewriting their results. A complete independent molecular replay requires the retained raw graph and decoding archives, original dataset and runtime dependencies; the aggregate JSON alone is not sufficient for that replay.

No recurring monitor or investigation job remains active. Processing completion and the classified reference differences are the conclusions of this completed investigation; no general correctness or universal reference-recovery claim is made.
