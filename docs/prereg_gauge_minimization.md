# Pre-registration: what predicts natural-orbital gauge sparsifiability

WRITTEN BEFORE THE RUN. Committed before `src/gpc_census/gauge.py` was executed
on any census state; the git history is the evidence. Scoring goes in
`docs/gauge_minimization.md` after the run, and a FAIL is a publishable finding
exactly as in `docs/prereg_bracket_3_11_3_12.md`.

## Background, stated so the prediction is not retrofitted

Natural orbitals are unique only up to unitary rotation WITHIN degenerate
occupation eigenspaces. Every census vertex has a degenerate spectrum, so every
certified state carries a gauge freedom

    G(lambda) = U(m_1) x ... x U(m_k),   m_k the eigenvalue multiplicities,

acting on the state by the N-th compound (exterior power) of the block-diagonal
unitary. Every point of that orbit has the SAME 1-RDM, exactly diag(lambda), so
every point is an equally legitimate natural-orbital representative of the same
vertex. The support recorded in the census library is the support of ONE point
of this orbit and was never minimized over it.

The quantity being computed is therefore

    s_Q^NO(v) <= min over the gauge family of the support,

still an upper bound on the census quantity s_Q^NO (other fiber components may
do better), but a strictly better one than the library column.

Two facts drive the prediction below and were known before the run:

- At the v103 support-14 certified state, all 42 off-diagonal U(4) x U(6)
  generators move the state OFF its support. The useful rotation, if there is
  one, is LARGE, not infinitesimal.
- Every fiber tangent computed in this program (silence-rank lemma, unified
  dimension formula, cascade) is SUPPORT-RESTRICTED. The gauge orbit generically
  leaves any fixed support, so support-restricted first-order data is
  structurally blind to gauge sparsification.

## The lemma that is NOT a prediction (proved, used as a screen)

Because the gauge unitary is block-diagonal, its compound preserves each
determinant's BLOCK PROFILE, the vector (|T ∩ B_1|, ..., |T ∩ B_k|). The
compound is unitary and preserves the profile decomposition, so the norm carried
by each profile class is a gauge invariant, and

    support of any gauge representative >= number of OCCUPIED profile classes.

Consequently a state whose support determinants lie in pairwise distinct profile
classes is EXACTLY gauge-minimal, with no computation. This is a lemma, not a
prediction, and it is scored as neither; it is used as an exact certificate of
minimality wherever it applies.

## Predictions, scored on the 799 certified library states

Outcome variable: DROP = (library support) - (support after gauge
minimization), an integer >= 0. A state "gauge-sparsifies" iff DROP > 0.

**P-G1 (main).** The support-restricted FIXED-RHO tangent dimension at the
certified state does not predict gauge sparsifiability. Scored PASS for this
prediction (i.e. the predictor FAILS) if at least one state with fixed-rho
tangent dimension 0 has DROP > 0. Scored FAIL for this prediction (the
predictor survives) if fixed-rho tangent 0 implies DROP = 0 across all 799.

**P-G2.** The same for the support-restricted FIXED-SPECTRUM tangent dimension,
by the same rule. Expectation: the predictor fails here too, i.e. some state
with fixed-spectrum tangent 0 gauge-sparsifies.

**P-G3.** The cascade drop count is not the gauge drop count. Expectation: the
per-state correlation is weak, and specifically at least one state has cascade
drop 0 and gauge DROP > 0. Scored on the joined table.

**P-G4 (the predictor I expect to work).** Gauge sparsifiability is governed by
the degenerate-block structure against the support incidence pattern, through
the profile-class lemma above. Concretely: define EXCESS = (library support) -
(number of occupied profile classes), the lemma's slack. Prediction: DROP <=
EXCESS for every state (this direction is forced by the lemma, so it is a
consistency check, not a discovery), and, the contentful half, EXCESS > 0
predicts DROP > 0 for a MAJORITY of the states where EXCESS > 0. Scored FAIL if
fewer than half of the EXCESS > 0 states sparsify.

**P-G5.** The cancellation-regime membership of a state (v89, v103 only) does
not predict gauge sparsifiability either way. Scored descriptively.

## Power caveat, stated before scoring

The census is not a random sample: it is one certified representative per
vertex, produced by a pipeline that had no reason to prefer any point of the
gauge orbit. That is good for P-G1/P-G2 (nothing in the construction correlates
with tangent dimension) and bad for P-G4 (the profile-class structure is a
property of the support, which the construction pipeline did shape). P-G4 is
therefore the weakest of the five and should be read as descriptive.

## What would make any of this moot

If the gauge minimization returns DROP = 0 for all 799 states, only the lemma
survives and P-G1 through P-G4 are unscoreable for want of positives. That
outcome is itself worth reporting: it would say the census library states are
already gauge-minimal, and that the v103 support-10 endpoint is reached by
genuine fiber motion rather than by gauge motion.
