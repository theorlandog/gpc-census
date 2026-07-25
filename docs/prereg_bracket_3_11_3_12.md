# PRE-REGISTRATION: out-of-sample test of the first-order fiber theory on the (3,11) and (3,12) bracket vertices (2026-07)

Status at time of writing: PREDICTIONS ONLY. Nothing below has been measured.
The measuring script (`scripts/oos_bracket.py`) is committed in the SAME state
as this file and its output does not exist yet. Score after running, before
writing any narrative. A FAIL is the more valuable outcome and must be written
up as such.

## What is being tested

The first-order fiber theory was fit on ranks 6 to 10 (the 799-vertex census):

- the unified fixed-spectrum dimension formula (`docs/unified_fiber_dimension.md`),
- the silence-rank lemma (`docs/silence_rank_lemma.md`) where the cancellation
  regime applies,
- the sparsification cascade (`gpc_census.cascade`, Task 1).

The holdout is the bracket settlement: 19 certified true vertices of (3,11)
(`docs/bracket_3_11_settlement.json`) and 19 certified true vertices of (3,12)
(`docs/bracket_3_12_true_vertices.json`). Code is run UNCHANGED, no tuning.

## Power caveat, stated before scoring

This holdout is weaker than it looks, and saying so afterwards would be
worthless. Of the 19 (3,11) true vertices, 17 carry the certificate
STATE-TRANSPORT: they are trailing-zero face embeddings of (3,10) census
vertices, and their witness states are the donor states with an unoccupied
orbital appended. All 19 (3,12) vertices are in turn face embeddings of (3,11)
vertices. So:

- genuinely new states in the holdout: TWO, both (3,11) EXPLICIT-STATE
  vertices, index 22 (frozen-core lift of an N=2 paired state, support 5) and
  index 49 (Z11 difference design on base block {0,1,3}, support 11).
- the remaining 36 test PADDING INVARIANCE of the theory and of the
  implementation at ranks the code has never run at, which is a real but much
  weaker test.

Therefore: a clean PASS here is weak confirmation. A FAIL, especially on B2 or
B3, indicates a rank-dependent implementation defect or a genuine gap in the
theory, and is the finding worth having.

## Inputs computed before predicting (support combinatorics only)

These are inputs to the formulas, not outcomes. For the two explicit states:

| vertex | support |S| | incidence kernel dim K | gauge rank g | one-hop channels |
|---|---|---|---|---|
| (3,11) index 22, [5,1^10]/5 | 5 | 0 | 5 | 0 |
| (3,11) index 49, [3^11]/11 | 11 | 0 | 11 | 0 |

## Predictions

B1. UNIFIED DIMENSION FORMULA, channel-free branch. For a one-hop-free support
the 1-RDM is exactly diagonal and the only first-order fixed-spectrum
conditions are the diagonal ones, so the weight sector contributes dim K and
the phase sector contributes |S| - g, which equals dim K because the gauge map
is the transpose of the incidence map. Predicted measured fixed-spectrum
tangent dimension = 2 dim K for every channel-free out-of-sample state, and in
particular:
    (3,11) index 22: predicted 0.
    (3,11) index 49: predicted 0.

B2. PADDING INVARIANCE. Every transported (3,11) vertex and every
face-embedded (3,12) vertex reproduces its donor's fixed-spectrum tangent
dimension, gauge dimension, and incidence kernel dimension EXACTLY. This is not
automatic: padding raises the multiplicity of the eigenvalue 0, which enlarges
the spectral projector block P_0 and therefore adds constraint rows
P_0 d(rho) P_0 = 0. The prediction is that those rows are all vacuous, because
d(rho) vanishes identically on any orbital unoccupied throughout the support.
A FAIL here is a rank-dependent implementation defect.

B3. THE TWO FIBERS AGREE ON THE HOLDOUT. For every one of the 38 out-of-sample
states, fixed_rho_tangent and fixed_spectrum_tangent report the SAME tangent
dimension. Rationale: the two differ exactly by the cross-block channel
conditions (silence-rank lemma), and B4 predicts no state here has channels at
all.

B4. SILENCE-RANK LEMMA SCOPE. The number of cancellation-regime members
(exactly diagonal 1-RDM carrying at least one one-hop channel) in the holdout
is ZERO. Consequence: the lemma is scored NOT-TESTED unless this prediction
fails. The two known members (v89, v103) are not among the 17 donors, whose
(3,10) indices are 0, 1, 2, 3, 29, 41, 64, 85, 104, 105, 106, 107, 108, 109,
110, 111, 112.

B5. CASCADE. No out-of-sample state sparsifies. Specifically support_upper
equals the certified support for both explicit states (5 and 11), and for every
transported or face-embedded vertex support_upper equals its donor's
support_upper. Rationale for the explicit states: dim K = 0 means the incidence
system has a unique solution on that support, and it is strictly positive, so
deleting any determinant makes the diagonal conditions infeasible. This is a
low-risk prediction and is scored as such.

B6. NO CRASHES OR SILENT DEGRADATION AT NEW RANK. Every one of the 38 states
returns a finished record: certificate residual below 1e-10 at the starting
point and spectrum error below 1e-8. A numerical failure at d = 11 or d = 12
that did not occur at d <= 10 counts as a FAIL of B6 even if B1 to B5 pass.

## Scoring rules

PASS / FAIL / NOT-TESTED per prediction, decided by
`results/data/oos_bracket.json` alone. Predictions whose antecedent is empty
(B4 if no cancellation states exist) are NOT-TESTED, never PASS. Any
disagreement between a prediction and a measurement is a FAIL and gets written
up with the measured numbers, not explained away.
