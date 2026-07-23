# RANK-10 CLOSURE PLAN (authoritative checklist)

Status: 797/799 certified; v89, v103 NUMERICALLY attained (contraction attack,
residuals 2.0e-16 / 3.5e-16; real-amplitude solutions also exist numerically).
Closure = exact certified states for both, ledger regenerated, paper finalized.

## Step 1 -- certify v89 and v103 (rig; first route to succeed wins)

Route A (recommended first): REAL-manifold exactification.
  - Re-run scripts/contraction_attack.py --real at high precision (mpmath
    refinement of the float64 solution via Newton on the real quadratic system).
  - Gauge-fix (O(2) x O(8) resp. O(4) x O(6) stabilizer + det signs), then
    PSLQ/minpoly on squared amplitudes and pairwise ratios at 60+ digits.
  - Reconstruct exact state; certify by exact char-poly identity.
Route B: sparsity-seeking. L1/iterative-thresholding reoptimization on the
  solution manifold toward small support; feed the standard exactifier.
Route C: equivariant. Subgroup G with C^10 = W2 (+) W8 (resp. W4 (+) W6),
  gamma = scalar + projector by Schur; tune within the invariant space in
  closed form. Most elegant; try if A/B stall.
Gate: both states pass verify-style exact identity AND independent
  reproduction (second instance, own reconstruction), same as every crack.

## Step 2 -- score the pre-registration (BEFORE any narrative)

Mark P1-P6 PASS/FAIL/NOT-TESTED in docs/RESEARCH.md ("SKEPTICAL AUDIT +
PRE-REGISTRATION" section) against the certified states: holonomy fields (P1),
channel counts (P2), ratio law scope (P3), deformability (P4), Aut (P5),
CM-conic model (P6). A FAIL is a finding; write it as one.

## Step 3 -- regenerate the ledger

All crack artifacts (including v89/v103 records + the 17 wall states as
fiber-enrichment entries) into results/data/states.jsonl via the standard
regeneration; CI green; counts become quotable ONLY after this step.

## Step 4 -- finalize paper 1 (one pass)

- Flip the five FINAL-COUNT markers (799/799; 564/564; delete the four-vertex
  open problem; residual section becomes the self-audit narrative).
- Add the bounded v_B absorption (its fiber, the ellipse--now conic--resolution,
  the two real wall states; Remark 3's open question answered with pointer to
  paper 2). Hard cap: nothing census-wide from the fiber era.
- Soften the blog's "can't be real" line; update the blog table to 799/799.

## Step 5 -- ship

arXiv submission via the standing endorsement; the one email (paper + PhD
question + spin-fiber joint direction); Zenodo release tag; EV application
same week (momentum is its currency).

## Parallel (not blocking closure)

- SOS certificate session for cand 44 (outputs/claude_code_prompt_sos_
  certificate.md) -- settles rank 11 on the refutation branch.
- Paper 2 Stage-B closeout (10 pending walls, surfaces, quasi-fiber theorem).
