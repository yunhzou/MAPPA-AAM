# Transition-state scoring

MAPPA exposes one TS scorer: `significant_subset_motion_v1` through
`analyze_transition_state` (or `mappa.mappa_workflow.ts_analytics.scoring.subset_score` for an already mapped
geometry). It permits asynchronous event motion, combines substantial event
support with contact selectivity and a bounded synchrony preference, and retains
endpoint-consensus ranking. Every result records `scoring_policy`.

Missing or invalid scoring inputs are errors, not a reason to substitute another
score. Results saved with an earlier scoring formula must be recomputed before
they are presented as current scores. A high score is not TS convergence or IRC
confirmation.

The reaction is supplied explicitly as a `ReactionContext`; no hidden
reactant/product selection is performed. From the command line, use
[`mappa --stage ts`](CLI.md#transition-state-scoring).

## Bond changes against both endpoints

Every TS score includes `diagnostics.endpoint_bond_changes` with separate
`reactant_to_ts` and `product_to_ts` lists. Each entry reports zero-based atom
indices, elements, endpoint/TS WBOs and distances, the directional WBO change,
its existing metal/nonmetal event threshold, and whether it is an intended
reaction pair or touches atoms outside the reaction core. These are descriptive
bond-order changes, **not** IRC validation or proof of an intermediate.

For complete spectator-bond coverage, supply
`TransitionStateTarget(..., reactant_atom_mapping=(...))`, a full R→TS atom
bijection. Identity mapping must be explicitly declared for atom-preserving
workflows; arbitrary file ordering is never assumed. Without that mapping the
report says `mapped_core_pairs_only` and `outside_core_assessed=false`.
Disagreement between the declared identity and a core assignment is reported.
Bond diagnostics do not alter the score or impose a rejection threshold.
