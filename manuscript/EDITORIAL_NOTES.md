# Editorial scope and remaining author decisions

The introduction now distinguishes alternative generation from single-answer
prediction and describes mapping as an endpoint correspondence, not a mechanism.
Base and target identify search roles, not fixed reactant/product identities.
Either endpoint can be the base; reverse search exchanges the roles. Prose,
mathematical notation and figure labels use this convention.
Supplied author metadata is retained. The Methods follow the current design in
separate scientific subsections: constraint discovery, inputs, existence
growth, compressed completion and symmetry, hierarchy/cut sweep, atom-set
reduction, local event responses, branch composition/windowing, and reuse/execution.
Reduction is a subphase of matching; the two scientific stages are Matching and
Decoding. Their CPU costs use that same grouping. The verifier appears only in
Evaluation. The case 9 alternatives figure is introduced early to motivate
retaining different correspondences.

## Author decisions remaining before submission

1. Corresponding authors, contributions and competing interests.
2. Selection/provenance, electronic-structure method and redistribution permission
   for the 140 archived coordinate/WBO pairs.
3. Permanent data/archive identifiers if required.
4. Applicability of the expanded acknowledgements and funding statements supplied
   in `f095f30a00d5e01e5ad79da2f47201123b2463b8`. Preserve that supplied text pending
   author confirmation; do not invent missing roles or declarations.

## Results and interpretation

The main evidence now covers Golden one/ten seeds and completed EnzymeMap at
tolerance 2 and cap 100, with decoded min+5 windows. Catalogue and decoded
membership have different representation scopes; Golden also uses different
reference-equivalence policies. Catalogue evaluation is complete with zero
unresolved reactions. Report its coverage as the primary full-output comparison;
label the three and two Golden unknowns as bounded decoded-window audit queries,
not unfinished matching or decoding.
EnzymeMap's conclusive reaction accounting includes separately justified absence
from the containing catalogue. Never turn a timed-out query into direct absence.

The branch-offset histogram is preliminary, with 3,891 checks pending. Existing
zero-offset certificates resolve first, so the subset is biased toward zero.
Include the figure and its definition, but do not infer final tail frequencies
or full-dataset percentages. It compares two counts from the same all-atom model.

Include reduction in matching time; keep decoding and external verification separate.
The saved Golden matching measurements predate indexed containment; the manuscript
reports approximately linear seed scaling and omits the cache-control comparison.
EnzymeMap timings combine successful original work and recoveries. Do not claim
a uniform final-version rerun or uninterrupted elapsed time.

The separately paired current SLAP study compares the same inputs and CPU model.
It supports catalogue coverage of 1,839/1,851 for MAPPA and 1,794/1,851 for expanded
SLAP, with timing on explicitly common completed cohorts. Default SLAP is faster;
MAPPA search is cheaper than the expanded sweep. The expanded protocol changes
both weight mode and search scope, so its gain is not a sweep-only ablation.
SLAP's published 86.9% is inclusion among alternatives on a 1,758-reaction version;
our 1,851-reaction reproduction is similar, not an exact dataset replication.

Minimum-choice counts use the saved one-seed decoded results, not the paired
study. Keep the 22 unfinished counts and two unresolved minimum-set memberships
explicit. Uniform-choice accuracy is a bounded expectation under that selector,
not observed top-1 accuracy or evidence of superiority over another top-1 method.

## Scientific limits to retain

Binary local symmetry is computed on each fragment, not inherited from the
whole endpoint. Certified global placement actions have the stronger invariance
needed to move a later continuation. Reduction keeps one original fragment per
atom-set pair; it does not union discarded discovery relations. Decoding selects
one starting orbit after bounded improvement. Ordinary boundaries use the
accepted fixed model; only the sweep correction retains joint endpoint dependence.
Model-window coverage is not physical event-bin accuracy. Keep these distinctions
in the main paper at a readable level, with mechanics in the supplement.
The new sweep, A-substituted carbon-ring and C--C triple-N diagrams illustrate
constraint discovery and boundary dependence. The latter two are analytically
checked heavy-atom graph examples without hydrogen constraints; they are not
new benchmark observations. Greediness helps reduce coupling but does not prove
fragment independence.

## Molecular figure conventions

Reserve molecule colors for fragment membership. Keep original black single/double bond strokes. Mark bond cleavage/weakening with a vivid red × directly on the affected bond, and formation/strengthening with paired bright green inward arrows. Keep the marks compact so fragment colors and black bond lines remain visible. Match fragment colors across mapped endpoints; different candidate partitions may assign colors differently and must say so. The selected discussion example is Golden case 9, with the verified heavy-atom reference labeled. Ground truth refers to the dataset annotation, not a validated mechanism. Count events from the evaluation matrices, including H, while distinguishing heavy-only depictions and unannotated H identities.

## Publication policy

Report only recorded CPU model names when identifying compute hardware. Exclude
site names, machine identifiers, addresses, accounts, storage paths, resource
inventories and scheduling configuration from the manuscript and bundle.
Author affiliations are separate from hardware disclosure. Build validation
checks the PDF and every bundle member. Detailed source/configuration hashes
belong in repository evidence, not the main scientific narrative.

The main Results compare all evaluated methods by coverage, with a single
prediction treated as a singleton set. Archived default-comparator predictions
retain their re-evaluation scope; their timings and historical MAPPA comparisons
remain in the supplement. Selected molecular examples remain illustrations, not new execution traces.
The coordinate study is not an accuracy benchmark. Cap 100 remains the default
study setting; cap 2,000 is labeled as a larger-budget historical comparison.
