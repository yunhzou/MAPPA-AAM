# Editorial scope and remaining author decisions

The Introduction follows the agreed narrative: AAM's meaning and applications,
ambiguity, taxonomy and prior work, computational challenge, and MAPPA's
contribution and principal findings. The alternatives example precedes the
taxonomy and motivates retaining different correspondences. The taxonomy
organizes related work into alternative generation and prediction, connected
by subsequent selection. Base and target are introduced in Methods as search
roles; either endpoint can be the base. Algorithm prose, mathematical notation
and figure labels use this convention.
Supplied author metadata is retained. The Methods follow the current design in
separate scientific subsections: constraint discovery, inputs, existence
growth, compressed completion and symmetry, hierarchy/cut sweep, atom-set
reduction, local event responses, branch composition/windowing, and reuse/execution.
Reduction is a subphase of matching; the two scientific stages are Matching and
Decoding. Their CPU costs use that same grouping. Growth and saturation matching
have separate subsubsections and a zero/one/many decision table. Multiplicities
in that table are conceptual target-set counts, not work performed during growth;
a rejected extension is distinct from a saturated fragment. Decoding is organized
around local symmetry, internal events, starting correspondence, ordinary
boundaries, sweep correction and window selection. The verifier appears only in
Evaluation. The EnzymeMap alternatives figure is introduced early to motivate
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
use the completed exact-index decoded audit separately,
with every reaction resolved.
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
The sweep, cyclopropanol/propan-2-ol and acetonitrile/ethylamine diagrams illustrate
constraint discovery and boundary dependence. The latter two are analytically
checked heavy-atom graph examples without hydrogen constraints; they are not
new benchmark observations. Greediness helps reduce coupling but does not prove
fragment independence.

Explain the intended mitigation alongside this limitation: random growth,
multiple seeds and cut sweeps explore different fragment decompositions. Another
decomposition can put a problematic attachment inside a fragment and represent
the same correspondence at its correct count. This does not repair the original
branch's fixed boundary model. Do not claim a guaranteed correct bin or a measured
reduction in boundary error without a dedicated comparison.

## Molecular figure conventions

Use chemical structure depictions for molecular examples, with the established vector drawing workflow. The saturation figure shows one carbon ring and two target rings distinguished by F/Cl attachments: two atom-set placements, with internal rotations and reflections retained within each placement. It is a fixed-fragment illustration without hydrogen constraints, not a complete reaction. Direct graph checks accompany the drawing. Do not replace molecular geometry with circles and group names when a chemical depiction can explain the same point.

Use the same pastel palette, black bond strokes and compact red cleavage/weakening and green formation/strengthening marks throughout. In algorithm figures, colors identify matched fragments. In the introductory alternatives figure, fixed colors track atom origins across both mappings; state this in its legend. The accepted example is EnzymeMap case 13175: the reference assigns the formaldehyde carbon to the three-carbon molecule, whereas an alternative assigns it to CO₂. Both retain the complete C4H6O5 composition. Show all nine heavy atoms, number only the four carbons, and label the displayed four versus two events as heavy-atom counts. Optimizing the unannotated H assignments separately for these fixed heavy-atom maps yields six total events for both. The saved witnesses have eight and six events, but the saved reference H completion is not optimal. Use the example to illustrate ambiguity; do not claim a strict all-atom advantage or global optimality. The independent H-completion check does not establish decoded membership of the new H assignments. Ground truth refers to the dataset annotation, not a validated mechanism. Keep molecular components separated and inspect atom labels at full resolution before accepting the figure.

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
