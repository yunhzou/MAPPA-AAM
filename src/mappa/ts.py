"""Typed transition-state processing composed from partial AAM and scoring."""
from __future__ import annotations

import time
from collections import defaultdict

import numpy as np

from mappa.mappa_core.core.search.partial import search_core_assignments
from .domain import (
    AAMSearchConfig,
    AtomAssignment,
    ReactionContext,
    TSMechanismResult,
    TSResult,
    TSScore,
    TSScoringConfig,
    TransitionStateTarget,
)
from .scoring import POLICY, pair_geometry, subset_score
from .bond_diagnostics import endpoint_bond_changes
from dataclasses import replace


def _clip01(value):
    return max(0.0, min(1.0, float(value)))


def _assignment_terms(rp, mechanism, target, assignment):
    reactant = rp.problem.reactant
    product = rp.problem.product
    target_molecule = target.molecule
    mapping_rp = mechanism.mapping.as_dict()
    mapping_rt = assignment.as_dict()
    terms = []
    for kind, bonds in (("broken", mechanism.broken_bonds),
                        ("formed", mechanism.formed_bonds)):
        for left, right in bonds:
            product_left, product_right = mapping_rp[left], mapping_rp[right]
            target_left, target_right = mapping_rt[left], mapping_rt[right]
            w_r = float(reactant.wbo[left, right])
            w_p = float(product.wbo[product_left, product_right])
            w_t = float(target_molecule.wbo[target_left, target_right])
            delta = abs(w_p - w_r)
            if delta < 1e-12:
                progress = 1.0
            elif kind == "formed":
                progress = _clip01((w_t - w_r) / delta)
            else:
                progress = _clip01((w_r - w_t) / delta)
            terms.append({
                "kind": kind,
                "reactant_bond": (left, right),
                "product_bond": (product_left, product_right),
                "target_bond": (target_left, target_right),
                "wbo_reactant": w_r,
                "wbo_product": w_p,
                "wbo_target": w_t,
                "event_weight": delta,
                "target_progress": progress,
            })

    total_weight = sum(term["event_weight"] for term in terms)
    progress = (1.0 if total_weight < 1e-12 else sum(
        term["event_weight"] * term["target_progress"] for term in terms
    ) / total_weight)
    return tuple(terms), float(progress)


def _score_assignment(target, assignment, sources, terms, progress, pairs):
    imaginary = np.flatnonzero(target.vibrations.frequencies < 0)
    if not len(imaginary):
        return None
    scored = [(int(k), subset_score(target.molecule.coordinates,
                target.vibrations.displacements[k], target.molecule.elements,
                terms, pairs)) for k in imaginary]
    k, diagnostic = max(scored, key=lambda item: (item[1]['score'], -item[0]))
    return TSScore(assignment=assignment, sources=frozenset(sources),
                   score=float(diagnostic['score']), wbo_progress=progress,
                   mode_index=k, frequency=float(target.vibrations.frequencies[k]),
                   event_terms=terms, diagnostics=dict(policy=POLICY, **diagnostic))


def analyze_transition_state(
        reaction: ReactionContext, target: TransitionStateTarget, *,
        search_config: AAMSearchConfig | None = None,
        scoring_config: TSScoringConfig | None = None) -> TSResult:
    """Evaluate one TS target under every selected R/P mechanism.

    R->TS and P->TS are independent partial-AAM searches.  P assignments are
    pulled into the R frame through the already selected R/P bijection, then
    exact tuples are unioned with source provenance before mode scoring.
    """
    if (not isinstance(reaction, ReactionContext)
            or not isinstance(target, TransitionStateTarget)):
        raise TypeError(
            "typed ReactionContext and TransitionStateTarget are required")
    search_config = search_config or reaction.config
    scoring_config = scoring_config or TSScoringConfig()
    started = time.perf_counter()
    results = []
    pending = []
    problem = reaction.problem
    for mechanism in reaction.mechanisms:
        core_r = mechanism.core_atoms
        if not core_r:
            results.append(TSMechanismResult(
                mechanism=mechanism,
                target=target,
                reactant_core_aam=None,
                product_core_aam=None,
                candidates=(),
                selected=None,
                status="no_reaction_core",
                reason=(
                    "the selected R/P mechanism has no broken or formed "
                    "bonds, so reactive-mode scoring is undefined"),
            ))
            continue
        core_p = tuple(mechanism.mapping.images[atom] for atom in core_r)
        from_r = search_core_assignments(
            problem.reactant, target.molecule, core_r,
            config=search_config,
            assignment_limit=scoring_config.core_assignment_limit)
        from_p = search_core_assignments(
            problem.product, target.molecule, core_p,
            config=search_config,
            assignment_limit=scoring_config.core_assignment_limit)
        inverse = mechanism.mapping.inverse().as_dict()
        candidates = defaultdict(set)
        objects = {}
        for assignment in from_r.assignments:
            candidates[assignment.pairs].add("reactant")
            objects[assignment.pairs] = assignment
        for assignment in from_p.assignments:
            pulled = AtomAssignment(tuple(
                (inverse[source_p], target_atom)
                for source_p, target_atom in assignment.pairs))
            candidates[pulled.pairs].add("product")
            objects[pulled.pairs] = pulled
        assignments = []
        for key in sorted(candidates):
            terms, progress = _assignment_terms(reaction, mechanism, target, objects[key])
            assignments.append((objects[key], candidates[key], terms, progress))
        pending.append((mechanism, from_r, from_p, assignments))

    # One normalization set across ALL retained mechanism/core assignments,
    # exactly as in significant_subset_motion_v1, never per-candidate selection.
    geometry = pair_geometry(target.molecule.coordinates, target.molecule.elements)
    events = {tuple(sorted(e['target_bond'])) for _, _, _, assignments in pending
              for _, _, terms, _ in assignments for e in terms}
    pairs = [(int(i), int(j)) for k, (i, j) in enumerate(zip(geometry['i'], geometry['j']))
             if geometry['distance'][k] <= 1.2 * geometry['radii_sum'][k]
             or (int(i), int(j)) in events]
    for mechanism, from_r, from_p, assignments in pending:
        scores = tuple(filter(None, (_score_assignment(target, assignment, sources,
                            terms, progress, pairs)
                        for assignment, sources, terms, progress in assignments)))
        scores = tuple(replace(score, diagnostics=dict(score.diagnostics,
            endpoint_bond_changes=endpoint_bond_changes(reaction, mechanism, target, score.assignment)))
            for score in scores)
        eligible = scores
        if scoring_config.prefer_endpoint_consensus:
            consensus = tuple(score for score in scores
                              if len(score.sources) == 2)
            if consensus:
                eligible = consensus
        selected = max(
            eligible,
            key=lambda score: (score.score, -score.mode_index,
                               tuple(-image for _, image in score.assignment.pairs)),
            default=None)
        results.append(TSMechanismResult(
            mechanism=mechanism,
            target=target,
            reactant_core_aam=from_r,
            product_core_aam=from_p,
            candidates=scores,
            selected=selected,
            status="scored",
        ))
    return TSResult(
        reaction=reaction,
        mechanisms=tuple(next(item for item in results if item.mechanism is mechanism)
                         for mechanism in reaction.mechanisms),
        elapsed_seconds=time.perf_counter() - started,
    )
