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

## SCORED (2026-07, after running; predictions above are unedited)

Predictions and `scripts/oos_bracket.py` were committed in b41d3e6 BEFORE the
script was run. Measurements: `results/data/oos_bracket.json`.

Score: 4 PASS, 2 FAIL. Both FAILs are the same vertex, (3,11) index 43, the
state transport of (3,10) index 108, spectrum [5,5,5,5,5,1,1,1,1,1,0]/10.

| prediction | verdict | scope |
|---|---|---|
| B1 unified dimension formula, channel-free branch | PASS | 36 states |
| B2 padding invariance | PASS | 19 embeddings |
| B3 the two fibers agree on the holdout | FAIL | 2 of 38 |
| B4 no cancellation-regime member | PASS | 0 members found |
| B5 no out-of-sample state sparsifies | FAIL | 2 of 38 |
| B6 numerics healthy at d = 11 and d = 12 | PASS | 38 states |

B1 PASS. Every channel-free state has measured fixed-spectrum tangent exactly
2 dim K, including both genuinely new states (index 22 and index 49, both 0).
Weak evidence by design, see the power caveat.

B2 PASS. All 19 (3,12) face embeddings reproduce their (3,11) inner vertex on
tangent dimension, gauge dimension, incidence kernel dimension AND cascade
outcome. The enlarged zero-eigenvalue block contributes only vacuous
constraint rows, as predicted. Notably this held through a support drop: the
inner vertex and its embedding both sparsify 8 to 7.

B3 FAIL. At (3,11) v43 and its (3,12) embedding, fixed_rho_tangent is 1 and
fixed_spectrum_tangent is 2. The prediction was wrong, and the reason is a
SCOPE ERROR in its stated rationale, not a numerical accident: B3 argued from
B4 that no holdout state has channels, but B4 only excludes the CANCELLATION
regime (diagonal 1-RDM with channels). v43 has one channel with a NONZERO
1-RDM off-diagonal, 0.2, so it sits in the MAGNITUDE regime, where the
cross-block channel conditions that separate the two fibers are live.

The finding this FAIL buys is worth more than the prediction was: the
fixed-rho / fixed-spectrum gap is NOT a peculiarity of the cancellation regime
where it was discovered (v103). It occurs in the magnitude regime too, at a
plain kdim-1 interference vertex, with gap 1 rather than v103's 6. Any text
that presents the two-fiber distinction as a cancellation-regime phenomenon is
too narrow. This is the second independent confirmation, on a holdout, that
rule 3 of the fiber program is load bearing.

B5 FAIL on its headline, PASS on its sub-clause. The headline "no out-of-sample
state sparsifies" is falsified: (3,11) v43 and its (3,12) embedding both drop
from support 8 to support 7, certificate residual 4.7e-17, spectrum error at
machine precision, dropping determinant (0,2,4). The sub-clause "transported
vertices reproduce their donor's support_upper" holds: the donor (3,10) v108
independently cascades 8 to 7 dropping the same determinant. So the FAIL is
not an out-of-sample effect at all. It is a failure of the predictor's
assumption that certified supports are hard to beat, and the census-wide run
(Task 1) confirms that assumption was wrong in-sample as well. The two
low-risk clauses of B5 that did hold: both explicit states stayed at their
certified supports, 5 and 11, as the dim K = 0 argument required.

B4 PASS, and therefore the SILENCE-RANK LEMMA is NOT-TESTED on this holdout.
Zero cancellation-regime members appeared, exactly as predicted. The lemma's
genericity question (sample size two) is untouched by this holdout, and the
count of known cancellation states is unchanged by it.

B6 PASS. All 38 states start at certificate residual below 7e-18 and spectrum
error below 7e-17 at d = 11 and d = 12; no rank-dependent numerical defect.

HONEST SUMMARY. The holdout was weak, as pre-registered, and it still produced
one genuine correction (B3) to how the two-fiber distinction is described.
Nothing here tests the silence-rank lemma, and B1's PASS should not be quoted
as out-of-sample confirmation of the unified formula in its interesting
regimes: the 36 states it scored are channel-free, which is the formula's
easiest branch.

## Scoring rules

PASS / FAIL / NOT-TESTED per prediction, decided by
`results/data/oos_bracket.json` alone. Predictions whose antecedent is empty
(B4 if no cancellation states exist) are NOT-TESTED, never PASS. Any
disagreement between a prediction and a measurement is a FAIL and gets written
up with the measured numbers, not explained away.
