# Agent instructions

Before starting project work, clone or reuse [PersonalCodingMoat](https://github.com/yunhzou/PersonalCodingMoat.git) and read its [README.md](README.md) and this file in full. Apply these principles alongside the target project's instructions. Reload them for each new session or after context loss.

## 1. Explanation

Keep explanations and documentation short, plain, and understandable to researchers. When communicating with the user, skip niche details unless they matter so the user can focus on important decisions. Say what is needed, then stop. Don’t repeat yourself or keep selling an idea. Explain more when the user asks or when a design decision needs discussion.

### 1.1 Paper writing

Keep one manuscript source: the Overleaf-compatible project at
`manuscript/main.tex`, with its shared `includes/`, `assets/` and `figs/`.
Edit that project directly. `mappa_overleaf.zip` is its generated export, not
another editable copy. Do not restore parallel manuscript trees or entry points.

Keep the appendix limited to methods and evidence needed for the current paper.
Do not include superseded benchmarks, old failure analyses, development history,
or operational storage/retry protocols. Keep those in repository records and Git.
Remove their unused tables, figures and build paths from the Overleaf export.
Retain assumptions and measurement scope needed to interpret the current results.

Keep the abstract focused on the problem, contribution, implications and main
results. Mention the central approach only briefly, using accessible concepts
such as shared molecular fragments. Do not introduce detailed methodology,
execution settings or unexplained terms such as seed counts. Retain informative
coverage and CPU-time results without turning the abstract into a Methods summary.

Build the Introduction as a scientific argument: AAM's meaning and applications,
ambiguity, taxonomy and prior work, the computational challenge, then MAPPA's
contribution and principal findings. Organize the related-work discussion around
alternative generation and prediction, with a concise comparison table and an
explanation of how a selector connects them. Use the molecular alternatives
figure beside the ambiguity discussion, before the taxonomy.
Treat author feedback as scientific intent, not sentences to insert. Rewrite
the argument when its structure fails; avoid accumulating qualifications,
terminology definitions and responses to prior discussions.

Use direct scientific prose: state the method, comparison and finding. Avoid
repeated defensive disclaimers. Explain material limitations where they affect
interpretation; keep software-version and run-provenance bookkeeping in the
supplement or reproducibility documentation.

Use the Discussion to interpret findings, practical trade-offs and downstream
possibilities. Keep its tone measured. Avoid repeating implementation details,
audit accounting or lists of untested capabilities; state only limitations
needed to interpret the claims being discussed. Describe the actual input
representation in Evaluation without adding unrelated validation disclaimers.

Describe AAM as correspondence between two molecular endpoints. Introduce the
algorithmic roles **base** and **target** in Methods: either endpoint can be the
base, and reversing the search exchanges these roles. Use base/target consistently
in algorithm prose, notation and figure labels, not fixed reactant/product labels.

Write the main manuscript for researchers who have not read the code. Organize
the algorithm into clear scientific stages: explain each stage's purpose, input,
operation and output, and how it connects to the next stage. Introduce the core
idea in plain language before notation or technical detail. Define essential
terms once and use them consistently. Follow the current approved `design.md`
and verify descriptions against the implementation.

Keep the main text at the level needed to understand the method, why it works
and what the evidence supports. Retain rules, assumptions and limitations that
affect scientific meaning. Put detailed derivations, solver mechanics, data
structures, cache and scheduling protocols, and full reproducibility settings in
the supplement or repository documentation unless essential to a central claim.
Simplify the explanation without changing the algorithm or overstating its
guarantees.

Open Methods with the complete fragment-growth-to-decoding workflow before
introducing equations or solver details. Use explicit stage names in subsection
titles (Matching Stage, Decoding Stage). Reduction is the final subphase of
matching, not a third scientific stage; include its cost in matching time.
Integrate the design
philosophy where it explains each operation. Keep figures on the same or adjacent
page to their first discussion, with the molecular overview on the first page
of the main text (Introduction). Keep the title and abstract on their separate
page; never place a scientific figure beneath the abstract.

Within Matching, give growth matching, saturation matching and conditional
branching their own short subsections. Define the growing fragment, locked
assignments, witness and atom-set placement before using them. Explain the
zero/one/many match transitions with a decision table and a concrete example;
positive multiplicity never triggers branching during growth. A failed edge is
deferred; saturation requires exhaustion of eligible extensions. Apply the same
structure to decoding: local symmetry, internal events, boundary contribution,
sweep correction and the event window. Use notation to expose the decisions,
not long paragraphs listing implementation details.

Use molecular examples, figures and concise pseudocode to clarify the stages.
Depict a matched fragment as a relationship between its base and target atom
sets: show both endpoint views, identify corresponding atoms consistently, and
distinguish the fragment from its surrounding molecular context. A one-sided
molecular drawing cannot stand in for a completed fragment relationship.
Keep one-witness growth, compressed completion, final family reduction and local
decoding distinct. Show these operations through molecular changes and relationships,
not a pipeline of text boxes. Reference verification is external evaluation
logic: explain it in the evaluation, outside the MAPPA pipeline and its overview
figure. Avoid repeating the same
explanation across prose, algorithms and captions, and keep repair history and
retired implementation choices out of the main narrative.

Present results for the implementation and scientific scope actually measured.
Separate historical studies from current-method evidence. During revision, ask
whether each technical detail helps the reader understand or assess a central
claim; move details that interrupt that explanation to supporting material.

Compare methods in the main Results by reference coverage over their returned
outputs, treating a single prediction as a singleton candidate set. Distinguish
coverage from ranking accuracy and state the evaluated output scope. The current
paper uses verified decoded minimum-plus-five coverage as its single MAPPA
coverage metric; separate catalogue audits remain in the benchmark evidence. Use a
molecular alternatives example early in the Introduction to motivate the output
objective, rather than leaving it at the end of Results.

Organize the main Results around coverage across datasets and methods,
computational cost and one-versus-multiple-seed scaling, reference-mapping event
offsets, and ambiguity at the minimum. Keep comparisons with other methods and
the one- and ten-seed results. Report findings rather than experiment history;
keep secondary timing cohorts and detailed audit accounting in supporting material.
Do not publish queued work as a result figure. Use the latest measured analysis
and retain uncertainty in the evidence. Figure 7 subtracts the saved detected
minimum from a verified GT-consistent mapping's event count, using the lowest
count among already saved witnesses across both directions. Do not optimize a
GT-containing branch or substitute a lower bound. Use the same all-atom event
model on both sides. Include saved anchored completions of missing references
in the same distribution, with one color per dataset and no separate completion
legend; these do not change independent coverage. Omit unoccupied negative ticks.

Report ranking accuracy only for a defined selection rule. An unranked candidate
list's first output is not a top-1 prediction; exclude accuracy statistics based
solely on that output order from the manuscript results.

## 2. Coding principles

### 2.1 Coding style

Write short, elegant, generalizable code with clear names. Avoid defensive coding and case-specific patches. Fix the underlying logic so the solution applies across cases.

Avoid abbreviations and unnecessary technical terms in names and comments. Follow the explanation rules in section 1: communicate plainly and efficiently for humans. Comments should explain reasoning that the code cannot show.

Keep code easy to follow. Don’t make it clever or compress it just to reduce line count.

### 2.2 Architectural design

Code structure must coevolve with vibe coding. During a rapid proof of concept, it is okay to prototype quickly, even in a single file with hundreds of lines. Once the user confirms that a feature is useful, refactor it into a coherent, modular architecture.

Discuss the high-level design with the developer. Establish the key concepts, module responsibilities, interfaces, and interactions before focusing on niche functionality. Think about the long-term structure and specify which parts of the code belong in each module.

Keep closely related functionality together, minimize dependencies between unrelated components, and separate reusable core logic from experiment-specific scripts. Modularity means meaningful separation of responsibilities, not simply dividing code into smaller files.

Follow the architecture previously approved by the user. Consider how new functionality fits and whether the design needs to change. Avoid competing structures or abstractions without a concrete need.

Revisit the design with the developer when refactoring reveals conflicting responsibilities, functionality does not fit the core modules, or the necessary architectural changes are still unresolved. Resolve these questions before proceeding with the affected refactoring.

You may provide architectural choices, tradeoffs, and recommendations, but the user makes the final decision. Once agreed, record the design in `design.md` and create diagrams showing module responsibilities, hierarchy, interfaces, and interactions.

## 3. Cleaning

Code must never have fallback options. Duplicate functionality is unacceptable. After refactoring, remove the superseded code and options completely from the repository. Leave no obsolete code behind.

Keep `design.md` as the single current architecture and scientific contract.
Update affected usage guides when interfaces change; remove conflicting plans.
Reusable operations belong in their existing source owner, portable experiment
runners in `bench/`, and input/output utilities in `tools/`. Reports retain
evidence and source revisions, not executable copies of retired diagnostics.
Historical claims must identify their measured revision and scope.

When you observe multiple fragmented pieces of code in the repository, treat this as a clear signal that it needs cleaning and refactoring. Follow principle 2.

For example, if several quickly written scripts exist in the repository and the user has confirmed that their functionality is useful, refactor that functionality into reusable core modules instead of leaving it in standalone scripts. Modularity and reusability must be enforced.

## 4. Version control

When the desired outcome is reached and the user confirms that a feature is good, immediately commit it to Git. Clearly document what experiment was performed, what was accepted, and whether refactoring and validation have been completed.

Version design principles as well as code. When the user explicitly specifies how the design should be structured, record those decisions and their reasoning in `design.md` and commit them.

Version control is not limited to preserving code. It must also support effective information exchange by documenting the reasoning and development process. Proper versioning also supports thorough cleaning.

## 5. Working on a cluster

Development often takes place locally and on a cluster. Whenever development spans multiple machines, always use Git to synchronize the work.

There must be a single source of truth for the project, and communication with the user must consistently refer to it. Each machine must have only one copy of the codebase. Enforce this, and use Git to reconcile changes whenever the copies become unsynchronized.

Keep cluster configuration, accounts, node/host names, paths, quotas and raw
scheduler records in private external configuration or a Git-ignored private
area. Commit only portable launchers and sanitized scientific evidence. Check
generated reports and logs for infrastructure details before publishing.

## 6. Validation scope

Use the test commands in README.md. For routine code refactors, run the smoke
checks and tests for affected modules once. Keep exhaustive chemistry tests for
matching/counting changes and the full suite for releases or changes spanning
core algorithms. Documentation-only changes need content/link checks, not runtime
tests. Do not repeat passing tests unless code changed or a specific unresolved
risk requires it. Do not rerun benchmarks to validate a structural cleanup.
Report the actual test scope. Retain independent scientific oracles; remove tests
that only exercise retired test helpers rather than supported behavior.
