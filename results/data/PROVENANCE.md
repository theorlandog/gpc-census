# PROVENANCE
Constraint systems ranks 9-10: parsed from M. Altunbulak, PhD thesis (Bilkent, 2008), Appendix A;
counts validated against Altunbulak-Klyachko CMP 282, 287 (2008): 52/60/93/125/161.
(3,9): 52 inequalities, 58 vertices.
Vertices: exact enumeration (lrslib, rational arithmetic). Validation: embedding coherence
(rank r -> r+1), frozen-core lifts, particle-hole self-duality of (4,8) and (5,10) at
inequality/vertex/verdict level, agreement with published rank-8 vertex tables and (4,9)=103.
Census verdicts: two-stage solver certificates (CP-SAT integer stage, ortools 9.15.6755;
CBC real stage). (3,10) vertex 89 resolved with extended real-stage budget (7200 s).
Rank-10 constraint lists are conjecturally complete per AK2008; rank-10 rows conditional.
Verdicts transport to dual systems (4,7),(5,8),(5,9),(6,9),(6,10),(7,10) by particle-hole
complement bijection on determinants.

## states.jsonl (constructed extremal states)

Per-vertex extremal states, produced by the routed state pipeline
(scripts/solve_all.py --all): DESIGN-INT vertices built directly from their
classification witness (exact by construction), DESIGN-REAL and INTERFERENCE
vertices through the weights-first block/clique solver (max_clique=3), certified
by exact characteristic-polynomial identity where the closed form is recognized
(status EXACT), else labeled TIER-C (exact numeric state, closed form not yet
recognized) or FAIL (outside the current ansatz family). Each record carries the
classification verdict, so the dataset cross-checks census_master.csv.

This file is COMPLETE: all 799 vertices of the nine determinate systems
(ranks 6 to 10) carry a certified closed-form extremal state. Regenerated fresh
it is deterministic up to the interference phase gauge, which does not affect
the certified spectrum.

Correction (v0.11): three DESIGN-REAL records, (3,9) index 53, (3,10) index 106,
and (4,10) index 53, shipped valid extremal states on one-hop-connected supports
rather than genuine design witnesses. They were replaced by design witnesses
transported (padding / frozen-core lift) from the (3,8) parent design at
denominator 12, each verified by the exact characteristic-polynomial identity;
their feature_table.csv rows (one-hop-class and max-class-size counts) were
regenerated to match.

## cascade.jsonl, cascade_summary.json (fixed-spectrum sparsification cascade)

NUMERICAL, not certified. One record per vertex from scripts/census_cascade.py,
which runs gpc_census.cascade over every state of states.jsonl: continuation
along the FIXED-SPECTRUM fiber driving the smallest surviving amplitude to zero
against the minimal-polynomial certificate prod_i (rho - lambda_i I) = 0 taken
over distinct eigenvalues, with the orbital-phase gauge pinned and each accepted
round gated on the eigenvalue multiplicities as well as the certificate.

The field support_upper is the smallest support reached. Read it alongside
support_certified: the certified support is an UPPER BOUND on the minimal
support s_Q of the vertex, and support_upper is a better upper bound where the
cascade landed a drop. Neither is a minimality claim, and support_upper is not
certified: support_upper_evidence records "numerical" on every record, and the
endpoint amplitudes are shipped (final_amplitudes_real, final_amplitudes_imag,
final_support_dets) so the endpoint can be rechecked without rerunning the
search. Fiber tangent dimensions are written under fiber-qualified keys
(fixed_spectrum_*, fixed_rho_*); the two are different objects.

## oos_bracket.json (out-of-sample bracket holdout)

Scored pre-registration, scripts/oos_bracket.py against the predictions in
docs/prereg_bracket_3_11_3_12.md, which were committed before the script was
run. Holdout: 19 certified true vertices of (3,11) and their 19 (3,12) face
embeddings. Result 4 PASS, 2 FAIL, both FAILs at (3,11) v43. NUMERICAL.

## reconciliation.json (superseded-claim reconciliation)

scripts/reconcile_claims.py re-derives the four claims this session corrects
(v103 rigidity by fiber, the isolated support-14 reading, the support column as
an upper bound, and the DESIGN-REAL denominator ratio count) so each correction
points at a measurement rather than at prose. The denominator-ratio item is
exact rational arithmetic; the rest are numerical.
