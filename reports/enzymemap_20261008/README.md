# MAPPA on EnzymeMap: final benchmark report

## Final results

MAPPA completed processing for **all 47,974 EnzymeMap reactions (100%)** under the recorded search configurations and decoding policies. The returned candidate sets contain a strict reference-equivalent mapping for **47,335 reactions (98.6680%)**. Every remaining reference difference has a completed, evidence-supported classification.

| Measure | Final result |
|---|---:|
| Reactions evaluated | 47,974 |
| Completed processing and decoding | **47,974 / 47,974 (100%)** |
| Strict candidate-set reference agreement | **47,335 / 47,974 (98.6680%)** |
| Classified reference differences | 639 |
| Unresolved processing cases | **0** |

Completion describes enumeration within the retained mapping families and the applicable output policy. Reference agreement describes whether an actual returned candidate matches the annotation. These measures are reported separately.

## Evaluation protocol

The evaluation uses the MAPPA engine based on commit [`ea6a7a3`](https://github.com/yunhzou/MAPPA-AAM/commit/ea6a7a3d6bfb4361d582dbc60dcfb87f183d1f8e), with forward search, formal bond orders, explicit hydrogens, iso tolerance 1, and one seed with random seed 42. Recorded branch caps are 100, 200 and 1,000. Decoding uses the declared output policy for each configuration, including the minimum-plus-five event window where specified. Input molecules and reference mappings are preserved.

The final reference-coverage figure is the **union of returned candidates across the recorded configurations and output policies**. Each reaction contributes once to the denominator and at most once to the match count. The result is therefore a combined-configuration benchmark rather than a uniform single-cap run.

### Reference comparison

A strict match requires a returned candidate to reproduce the annotated heavy-atom relation modulo endpoint chemical symmetry. The comparison includes charge, isotope, CIP and bond-stereo distinctions, with unmatched atoms represented explicitly.

MAPPA decodes compressed mapping families into signed bond-event classes with explicit representative mappings. A reference event class can be present even when its returned representative differs under the finer, stereo- and charge-sensitive mapping comparison. The report distinguishes these outcomes from reference absence in the retained search families.

Coverage is evaluated over the candidate set. It does not measure top-1 chemical ranking or establish the physical reaction pathway.

### Partial mappings and graph preparation

For **16 reactions with unequal explicit-hydrogen composition**, the output policy supports partial mappings. It prioritizes mapped heavy atoms, then mapped atoms, and enumerates the applicable minimum-plus-five event window. Unmatched atoms remain explicit. **All 16 reactions completed and matched the strict heavy reference.** Processing completion therefore accommodates partial atom coverage where the endpoint compositions differ.

For **11 reactions evaluated with independently prepared input and reference graphs**, each graph retains its own coordinates, features and atom indexing. **All 11 reactions completed and matched the strict reference.** Neither endpoint molecule nor its annotation was altered.

## Classification of all reactions

| Final outcome | Reactions |
|---|---:|
| Strict reference-equivalent candidate returned | **47,335** |
| Reference absent from the completed recorded search configurations | 472 |
| Reference event class returned; representative differs under stereo-sensitive comparison | 166 |
| Reference event class returned; representative differs under charge-sensitive comparison | 1 |
| **Total** | **47,974** |

The **639 strict-reference differences** consist of **472 search absences** and **167 representative distinctions**. All have conclusive classifications within the recorded evaluation scope.

For the 472 search absences, the reference is outside the completed retained family unions for the tested configurations. This conclusion is limited to those caps, seeds and search settings.

For the 166 stereo-sensitive cases, the returned representative agrees with the reference when only stereo labels are ignored. The single charge-sensitive case, index 12009, contains a matching zero-event class with Fe²⁺/Fe³⁺ distinctions under the strict mapping comparison. The event-class and strict-reference results remain separate in the reported totals.

## Validation

The final results passed independent checks of decoding completeness, physical mapping membership, reference comparison and evidence integrity.

| Validation component | Result |
|---|---:|
| Current decoder task configurations validated | **226** |
| Exact-reference proof task configurations validated | **211 / 211** |
| Required search, decoding and additional-output tasks completed | **464 / 464** |
| Reactions covered by the detailed reference-miss ledger | **807 / 807 classified** |
| Packaged proof-record members verified against original archive manifests | **131,552** |

Task counts include controls and multiple configurations of the same reaction; they are distinct from reaction counts. The 211 exact-reference checks comprise 26 positive and 185 negative task outcomes.

Validation includes physical-witness membership checks, reconstruction of saved residual queries, independent proof replay, and SHA-256 verification of archived evidence. The exact permutation encodings used for decoding were checked through exhaustive small equivalence tests and full control outputs. The all-row reconciliation covers exactly **47,974 unique reaction indices**, with no unresolved processing cases.

## Interpretation and supporting records

The completed benchmark establishes **100% processing completion** and **98.6680% strict candidate-set reference agreement** across the recorded configurations. The remaining differences are fully classified within that scope. Search completeness is relative to retained families and declared event windows; the results do not establish exhaustive coverage of every possible mapping, seed, cap or chemical pathway.

[Machine-readable aggregate results](results.json) accompany this report. The detailed case ledger, completion reconciliation, proof/query records, audit manifests and research source snapshot are retained in the local evidence package and original cluster archives. This GitHub publication contains the final report and aggregate results; full independent molecular replay requires those retained artifacts and their runtime dependencies.
