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

Update (v0.12): the (4,9) index 65 record (v_B = (20,14,14,14,14,4,4,4,4)/23)
now ships the REAL extremal state: the theta=0 endpoint of the wall family
k(t) = (1+t, 8-t, 4, 3, 2-t, 2+t, 1, 2)/23 on the same support, at
t0 = 7/17 - 18*sqrt(35)/85 (reality wall 85 t^2 - 70 t - 119), squared
weights in Q(sqrt(35)). Generated and certified by scripts/vb_real_state.py
(exact characteristic-polynomial identity; independently re-verified by
scripts/verify_states_standalone.py). closed_form.weights for this record
holds exact radical strings instead of integers; its feature_table.csv
all_real flag is now 1. The previously shipped phase-carrying fiber point is
the t=0 rational anchor of the same family, recoverable from the family
metadata in the record.

## cascade.jsonl, cascade_summary.json (fixed-spectrum sparsification cascade)

NUMERICAL, not certified. One record per vertex from scripts/census_cascade.py,
which runs gpc_census.cascade over every state of states.jsonl: continuation
along the FIXED-SPECTRUM fiber driving the smallest surviving amplitude to zero
against the minimal-polynomial certificate prod_i (rho - lambda_i I) = 0 taken
over distinct eigenvalues, with the orbital-phase gauge pinned and each accepted
round gated on the eigenvalue multiplicities as well as the certificate.

The field support_upper is the smallest support reached. Read it alongside
support_certified: the certified support is an UPPER BOUND on the minimal
support, and support_upper is a better upper bound where the cascade landed a
drop. The two bound DIFFERENT quantities: support_certified bounds s_Q^NO (its
state has rho = diag(lambda)) while support_upper bounds only s_Q^free, since
no cascade endpoint has a diagonal 1-RDM. See orbit_check.json. Neither is a minimality claim, and support_upper is not
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

## cascade_exact.json, v103_support10_certificate.json (EXACT support bounds)

EXACT, in contrast to cascade.jsonl. scripts/exactify_cascade.py refines each
cascade endpoint in 40-digit arithmetic, recognizes the squared weights by
integer relation detection, gauge-fixes the phases, and verifies the
characteristic-polynomial identity det(rho - x I) = prod_m (n_m/D - x) in exact
rational and radical arithmetic. An entry is marked CERTIFIED only when the
weights are rationals summing to 1, the state is real up to the orbital-phase
gauge, and the exact identity holds; every other status is reported with the
reason it fell short and its bound stays numerical.

69 of the 71 sparsified endpoints are CERTIFIED, so those vertices carry an
exactly certified support bound strictly better than the library state's. The
two exceptions, (5,10) index 144 and index 256, have partly irrational weights.

v103_support10_certificate.json is the same result for (3,10) index 103 by an
independent standalone route, kept separate because it is the instance that
changes a published claim: support 14 in the library, support 10 certified here,
state denominator 46852 = 2^2*13*17*53. The two routes agree on support,
weights and denominator; their sign patterns differ by a gauge choice and both
satisfy the exact identity.

## v103_endpoint.json (high-precision endpoint)

The 60-digit refinement of the v103 support-10 endpoint with its
integer-relation recognition, kept as the numerical input the exact certificate
was derived from.

## gauge_min.jsonl, gauge_min_summary.json (natural-orbital gauge minimization)

Produced by scripts/gauge_census.py (engine src/gpc_census/gauge.py). One record
per vertex: the library support, the exact gauge-invariant lower bound (occupied
block-profile classes, refined by per-class flattening ranks), the support
reached by minimizing over the natural-orbital gauge family
U(m_1) x ... x U(m_k), and the acceptance gate on any state that moved.

The headline is a null result with an exact half: no state of the 799 admits a
support-reducing rotation in its own gauge family, and for 629 of them the
library support attains the lower bound, so gauge-minimality is certified by a
rank computation rather than by a search. The other 170 are open in both
directions. The bounds are exact (integer ranks of exactly built matrices); the
absence of a reduction is a search result, whose power is pinned by the
recovery test in tests/test_gauge.py.

These records bound s_Q^NO, the minimum support over states with rho exactly
diag(lambda). They are NOT comparable with the cascade's free-basis bounds.

## orbit_check.json (2-RDM gate on the cascade endpoints)

Produced by scripts/orbit_check.py. The 2-RDM spectrum is invariant under every
one-body unitary, so it decides whether a cascade endpoint is a new state or the
library state in different orbitals. Over all 71 endpoints: 63 agree with the
library state to ~1e-16 (same state, rotated), 8 differ by 3.4e-3 to 4.7e-2 (no
one-body rotation can do that, so they are genuinely different extremal states),
and ZERO have a diagonal 1-RDM, which is why every cascade bound is a bound on
s_Q^free and not on the natural-orbital quantity. The v103 entry carries the
explicit connecting unitary (residual 7.2e-16, off-block entries 0.196, so a
U(10) rotation and not a gauge one). Every one of the 63 same-state endpoints
ships its rotation matrix, all verified, worst residual 1.7e-15.

## holonomy_fields.json (gauge-invariant loop holonomy fields)

EXACT, with an independent numerical cross-check. scripts/holonomy_field_scan.py
walks every certified INTERFERENCE state whose support carries loops (63 states,
82 independent loops) and records, per loop: the integer kernel vector, the
exact holonomy cosine, its minimal polynomial, its degree, its Galois group
computed with sympy galois_group rather than inferred from the shape of the
polynomial, the extracted square-root radicands, and a PSLQ integer-relation
recognition performed at 120 digits and RE-VERIFIED at 240, with the
verification residual recorded on every row. A relation that fails
re-verification is deleted, not downgraded; the discriminating evidence for that
protocol is reproduced in the precision_lesson block.

Result: all 82 degrees lie in {1, 2, 4} (37 / 32 / 13) and all 82 Galois groups
are 2-elementary (C1 / C2 / V4). The loop-kernel bases are certified saturated,
so no holonomy relation is missed. The channel_hypotheses block verifies the
hypotheses the theorem of docs/holonomy_rationality.md runs on: all 84 one-hop
channels have rational |rho_AB|^2, and for all 63 states the within-channel
difference vectors generate the loop lattice, saturated over Z. The
trisection_audit block shows the pipeline's TRISECTED phase tier is a property
of how the amplitude phase is printed, not of the state: no cubic subextension
occurs anywhere in those 13 states.

## denominator_arithmetic.json (where a state denominator comes from)

EXACT. scripts/denominator_arithmetic.py splits the state denominator into an
incidence stage (solve A_S p = lambda, sum p = 1 over Q, with the Smith
invariant factors of A_S recorded) and a pinning stage (the one-hop channel
conditions restricted to the incidence solution space). Worked examples at v89,
where the silence condition on the incidence line loses its quadratic terms and
the surviving linear equation 130 p9 - 1 = 0 mints the prime 5, and at v103,
where two channel equations have identical quadratic parts and their difference
13 p12 = p9 + p13 carries the minted prime 13.

The census block scores the predictions pre-registered in
docs/prereg_denominator_arithmetic.md: P-DEN-2 FAIL (18 exceptions, mechanism
recorded), P-DEN-3 FAIL as predicted (1 exception, v103, pinned as predicted),
P-DEN-4 PASS (736 states, 0 exceptions). The (4,9) index 65 record is excluded
from the ratio predictions: its weights are not rational, so it has no state
denominator.

## two_block_family.json (the two-block family, and the cross-block rule)

MIXED, labeled per field. Channel splits, block structure and transport
classes are EXACT combinatorics; tangent dimensions are NUMERICAL, inherited
from cascade.jsonl. scripts/two_block_family.py assembles the four two-block
INTERFERENCE vertices in full (channels split into within-block and
cross-block, both fiber tangents, holonomy fields, denominator arithmetic,
gauge-minimality status, cascade behaviour) and scores the predictions
pre-registered in docs/prereg_two_block_inversion.md.

Result: a CROSS-BLOCK one-hop channel is NECESSARY for the fixed-rho and
fixed-spectrum tangent dimensions to differ, on all 799 states across 403
transport classes with zero exceptions. The converse is false at 82 states,
all with incidence kernel 0. The two proposed counting rules (gap = X and
gap = 2X) are falsified. Independent scope is reported in TRANSPORT CLASSES
rather than row counts, per standing rule R2 of docs/prereg_template.md,
because padding and frozen-core siblings carry a state's fiber data unchanged.

The tangent_sector_split block decomposes v103's six fixed-spectrum
directions as 3 weight plus 3 phase, reproducing the silence-rank lemma's
formula term by term; at (3,10) v108 the sectors do not decouple, because the
state is genuinely complex and the lemma's real-state bookkeeping does not
apply there.

## oos_bracket.json: provenance_independence block (added 2026-07)

The scorecard was scored at 38 rows. A later audit traced every measurement to
its root origin and found 34 of them are transports or face embeddings of
rank-10 states drawn from 17 distinct donors, leaving 4 rows resting on 2
states actually constructed at the new rank. The provenance_independence block
records this per measurement, and every verdict now carries independent_scope
and independent_distinct_states beside its row scope, plus a corrected reading
alongside the original verdict. The original verdicts are unchanged; see the
PROVENANCE CORRECTION section of docs/prereg_bracket_3_11_3_12.md.

## attainer_gate_audit.json (diagonality gate over the whole library)

NUMERICAL, cross-checked in exact arithmetic at (3,10) v17.
scripts/attainer_gate_audit.py runs gpc_census.validate.check_vertex_attainer
over all 799 records. 645 pass (all 631 DESIGN-INT, all 12 DESIGN-REAL, and
the two cancellation-regime states v89 and v103); 154 fail, and they are
exactly the remaining INTERFERENCE records.

Every failing record still attains the spectrum (worst error 1.1e-15) at unit
norm, and none is a merely permuted relabelling: the failures are a basis
statement, not a wrong-state statement. The shipped certificate is the
characteristic-polynomial identity, which certifies spec(rho) = lambda and
asserts diagonality only for design states, where one-hop freedom gives it by
inspection. The operative consequence is that the support column of those 154
records bounds s_Q^free and not s_Q^NO.

## rigidity_rationality.json (does rigidity imply rational weights)

MIXED. Rigidity is NUMERICAL, inherited from cascade.jsonl; the weight field
is EXACT, recomputed from the closed form rather than read off the
integer-versus-string encoding. scripts/rigidity_rationality_census.py joins
the two per record and tallies the implication "fixed-rho first-order rigid
implies rational squared weights": 740 of 740 with zero counterexamples, the
single irrational record ((4,9) v65, weights in Q(sqrt 35)) being deformable.

The artifact carries a power_warning block that must travel with the number:
all 645 natural-orbital representatives are rigid and all 59 deformable
records fail the diagonality gate, so the implication has no
rigid-versus-deformable contrast among gate-passing records and its
discriminating power lies entirely in the non-natural-orbital population.
