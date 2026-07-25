# Census-wide sparsification cascade (fiber program Task 1, 2026-07)

Artifacts: `results/data/cascade.jsonl` (one record per vertex),
`results/data/cascade_summary.json`. Engine: `src/gpc_census/cascade.py`.
Driver: `scripts/census_cascade.py`. Tests: `tests/test_cascade.py`.

EVIDENCE LEVEL: everything here is NUMERICAL. A landed round is a continuation
endpoint whose minimal-polynomial certificate and eigenvalue multiplicities hold
to the residuals quoted, not a proof that an exact state exists there. Task 3
(alpha-theory or exactification) is what upgrades these; until then no claim
below may be quoted as certified.

## What was run

The v103 cascade generalized to every certified state of the shipped ledger,
all 799, without manual intervention: parameterize by the squared modulus of
the smallest surviving amplitude, pin the orbital-phase gauge, solve the smooth
minimal-polynomial certificate prod_i (rho - lambda_i I) = 0 with one factor per
DISTINCT eigenvalue, continue the parameter to zero, drop the amplitude, repeat.
Where the greedy smallest-first branch stalls, the run retries from random
fiber points (deform first, then cascade) before recording a floor.

Two additions to the method as handed over, both load bearing:

1. MULTIPLICITY GATE. The minimal polynomial constrains the eigenvalue SET, not
   the multiplicities. Without an explicit check, a round can zero an amplitude,
   satisfy the certificate to 1e-17, and land on a DIFFERENT stratum whose
   1-RDM has the same distinct eigenvalues in different numbers. This is not
   hypothetical: 33 of the 799 states have their cascade stopped by exactly this
   (`blocked_reason: multiplicity_drift`), the first of them v_B. Every accepted
   round is now gated on max |sorted eig(rho) - sorted lambda| < 1e-8.
2. POST-DROP POLISH. Driving |c_j|^2 to zero leaves c_j at solver tolerance, so
   deleting it and renormalizing perturbs the 1-RDM at FIRST order in that
   leftover. Each accepted round is re-solved on the reduced support, which
   moved the worst residual among all drops from 2.7e-10 to 2.2e-16.

## Headline result

71 of the 799 certified states sparsify, removing 89 amplitudes in total. All 71
are INTERFERENCE vertices; no DESIGN-INT or DESIGN-REAL state moves. 69 of the
71 improved bounds are then certified in EXACT arithmetic (see below), so this
is not only a numerical result.

The census support column is therefore an UPPER BOUND on the minimal support
s_Q, and it is now known STRICT at 71 vertices, not just at v103. The new column
`support_upper` in `results/data/cascade.jsonl` records the best support reached,
with `support_upper_provenance` naming the method and `support_upper_evidence`
recording that it is numerical.

Drops by system: (3,8) 4, (3,9) 6, (3,10) 22, (4,9) 8, (4,10) 9, (5,10) 22.
Drop sizes: 55 states drop one amplitude, 15 drop two, one drops four (v103).

## The regularity the run exposes

Cross-tabulating the first-order fixed-spectrum tangent dimension at the
certified point against the number of amplitudes the cascade removes:

| fixed-spectrum tangent dim | states | amplitudes dropped |
|---|---|---|
| 0 | 726 | 0 for all 726 |
| 1 | 11 | 1 for all 11 |
| 2 | 46 | 1 for 44, 0 for 2 |
| 4 | 15 | 2 for all 15 |
| 6 | 1 (v103) | 4 |

Two directions, stated at their true strength:

- FIRST-ORDER RIGID IMPLIES NO SPARSIFICATION, 726 for 726, no exceptions. This
  is the expected direction (the cascade moves along the fiber, and a rigid
  point has no fiber directions to move along), and the run confirms it with no
  counterexample. It is an observation about this cascade, not a theorem: a
  second-order deformation invisible to the first-order test would evade it.
- POSITIVE TANGENT IMPLIES SPARSIFICATION, 71 of 73. The two exceptions,
  (4,10) v62 and (5,10) v178, both have tangent 2, kernel dim 1, one channel,
  and both floor at 2 to 3e-6, which is the regime the handoff explicitly marks
  as BLOCKED ON THIS PATH rather than blocked in general. They are recorded as
  path floors. No minimality claim is made for them, and none should be.

The correspondence is close enough to be worth stating as a candidate law and
far too clean to be assumed: the mined-law mortality rate in this project is
about one in three, and this is a post-hoc pattern in a corpus the cascade
itself explored. It is offered as a hypothesis for out-of-sample testing, not
as an established regularity.

## v103: the certified support is not minimal, and 12 was not the floor either

The previous session's greedy cascade reached support 12 and recorded a floor of
1.4e-6 at the next round. With the continuation schedule and fallback targets in
the current engine, v103 goes further:

    support 14 (certified) -> 13 -> 12 -> 11 -> 10

dropping determinants (2,6,9), (1,3,6), (2,3,6), (1,2,9), every round at
certificate residual below 1e-16 and spectrum error at machine precision,
superseding the <= 12 in `docs/unified_fiber_dimension.md` section 3. The
earlier floor at (2,3,6) was a property of that continuation path, exactly as
the handoff's warning said, and not a property of the fiber.

### s_Q(v103) <= 10 is CERTIFIED, not numerical

The endpoint turned out to be exactly recognizable, so this one bound does not
stay at [N]. Refining it by Gauss-Newton in 60-digit arithmetic
(`scripts/v103_endpoint_precision.py`, residual 1.1e-62) and running integer
relation detection on the squared weights recognizes ALL TEN as rationals, and
its single loop holonomy is pi, so the state is real up to the orbital-phase
gauge. Substituting the exact values back and rebuilding the 1-RDM in exact
rational and radical arithmetic verifies the characteristic-polynomial identity

    det(rho - x I) = (x - 9/17)^4 (x - 5/34)^6

with no floating point anywhere (`scripts/v103_support10_certificate.py`,
`tests/test_v103_support10.py`, artifact
`results/data/v103_support10_certificate.json`). The state:

| determinant | squared weight | sign |
|---|---|---|
| (1,2,3) | 13/34 | + |
| (0,2,7) | 5/34 | + |
| (0,5,6) | 53/442 | + |
| (0,1,8) | 195/1802 | - |
| (0,3,4) | 5/68 | + |
| (0,4,9) | 5/68 | + |
| (3,8,9) | 35/901 | + |
| (1,6,9) | 1/34 | + |
| (3,5,9) | 18/901 | + |
| (0,1,5) | 84/11713 | + |

Its state denominator is 46852 = 2^2 * 13 * 17 * 53, against 195364 = 442^2 =
2^2 * 13^2 * 17^2 for the certified support-14 library state. The prime 53 is
new: it appears in neither the spectrum denominator 34 nor the library state's
denominator, which is one more reason no grid search over multiples of the
natural denominator was ever going to find this state.

So s_Q(v103) <= 10 **\[E\]**. This is an upper bound and NOT a minimality
claim: the endpoint is first-order rigid in both fiber senses, which makes it
isolated, but isolation of one point says nothing about other components of the
support-10 or sparser strata. The rest of the census-wide cascade remains
numerical.

The endpoint is structurally interesting, and is the best Task 3 and Task 4
target this run produced:

- its 1-RDM has eigenvalues (18^4, 5^6)/34 with the correct multiplicities,
- it is first-order RIGID in both senses at the endpoint (fixed-spectrum tangent
  0, fixed-rho tangent 0), so it is an ISOLATED point of the support-10 stratum
  modulo gauge, which is precisely the kind of object alpha-theory certifies and
  exact recognition can pin down,
- it has left the cancellation regime: its 1-RDM off-diagonal reaches 7.4e-2,
  so the walk from support 14 to support 10 carries the state from the
  cancellation regime into the magnitude regime,
- every one of its ten squared weights is an exact rational (table above), and
  the (0,2,7) weight is still exactly 5/34, the fiber-constant invariant the
  handoff flagged, now confirmed to survive four more cascade rounds.

The exactification lead the handoff recorded is therefore closed for this
endpoint: the high-precision re-solve plus integer relation detection worked on
the first attempt, and the recognized values verify exactly. The remaining
exactification targets are the other 70 sparsified states, none of which have
been refined or recognized.

## 69 of the 71 improved bounds are CERTIFIED, not numerical

The v103 recipe generalizes. `scripts/exactify_cascade.py` runs it on every
cascade endpoint: refine to 40 digits by Gauss-Newton, recognize the squared
weights by integer relation detection, gauge-fix the phases, and verify the
characteristic-polynomial identity in exact rational and radical arithmetic.
A state is reported CERTIFIED only if the weights are rationals summing to 1,
the state is real up to the orbital-phase gauge, and the exact charpoly matches.

    endpoints attempted     71
    CERTIFIED               69
    weights not rational     2

So 69 vertices now carry an EXACTLY certified improved support bound, 55 of them
one determinant better than the library and 13 two better, plus v103 at four
better. Artifact: `results/data/cascade_exact.json`, which ships the exact
weights, signs, support and state denominator of every certified endpoint.

Two independent routes agree at v103: the standalone certificate
(`scripts/v103_support10_certificate.py`) and this pipeline produce identical
support, identical squared weights and identical state denominator 46852. Their
sign patterns differ, which is expected and not a discrepancy: signs are
gauge dependent, the two runs fix the gauge differently, and both patterns pass
the exact characteristic-polynomial check.

The two that resisted are (5,10) v144 and (5,10) v256, both dropping 10 to 8 at
refinement residual below 1e-40. Their endpoints are partly rational (four of
eight squared weights at denominator 37) and partly not, which is the signature
of a positive-dimensional stratum whose distinguished exact point the cascade did
not land on. This is the same phenomenon the paper records for v_B, whose real
extremal states sit at quadratic-irrational weights in Q(sqrt 15) and Q(sqrt 35)
rather than on any rational grid. Extending recognition to quadratic surds is
the obvious next step for these two and was not done here; their bounds stay
NUMERICAL.

One methodological note worth keeping. The first version of this pipeline
gauge-fixed the phases by an ordinary least-squares fit, which spreads the
residual over every coordinate instead of concentrating it on the
gauge-invariant directions, and it wrongly reported v103's endpoint as not real
up to gauge. Fitting exactly on a maximal independent set of support rows is the
correct test. The bug was caught only because v103 had already been certified by
hand, which is an argument for keeping at least one instance verified by an
independent route.

## Cancellation-regime census: still exactly two

The regime (1-RDM exactly diagonal with at least one one-hop channel, so every
channel is silent by cancellation) has exactly two members across all 799
certified states: v89 and v103. The run adds none. The silence-rank lemma's
genericity question (Task 6) remains on a sample of size two and is still
unanswerable from census instances alone; a synthetic-support search is the only
route left short of new certified states.

Note the census counts the CERTIFIED representative of each vertex. The regime
is a property of the representative, not of the vertex: v103's own cascade
leaves the regime by round 4. So "exactly two" means two certified
representatives, and a different fiber point of the same vertex may or may not
be in the regime.

## Representative-swap check (v0.12 real v_B record)

Everything above was first computed against the ledger that shipped the
phase-carrying fiber point over v_B = (4,9) index 65, and re-run against the
v0.12 ledger that ships the real theta=0 endpoint of the same wall family.
Nothing moved: the same 71 states sparsify, the same 89 amplitudes come out,
the same 33 states stop on multiplicity drift, the cancellation regime still
has two members, and all 799 per-vertex verdicts and tangent dimensions are
identical. Exact certification is identical too, 69 of 71 with the same
weights.

v_B itself is a useful check on representative dependence. Both
representatives give incidence kernel 1, fixed-spectrum tangent 2 and
fixed-rho tangent 1, and both cascade to support 7 by dropping the same
determinant (0,1,2,4), certifying at the same squared weights (9/23, 4/23,
3/23, 12/115, 8/115, 1/23, 2/23) with state denominator 115. The two endpoints
differ only by a gauge phase on one amplitude: starting from the real
representative the cascade lands on the real form directly, with an
all-positive sign pattern, where the complex representative needed the
gauge-fixing step to get there. Fiber-level quantities are representative
independent here, which is what the theory predicts and worth having on
record given how often representative dependence has bitten this project.

## Caveats

- Greedy smallest-first is one branch of a large tree. Every floor recorded here
  is a floor for the paths this run tried.
- The 33 multiplicity-drift stops are not evidence that those states cannot
  sparsify. They mean the natural path leaves the vertex's own fiber, and a
  different path might not.
- All residuals quoted are of the minimal-polynomial certificate and the sorted
  eigenvalue error. Neither is a proof.
