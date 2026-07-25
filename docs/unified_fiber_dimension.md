# Unified First-Order Fiber Dimension + the v103 Sparsification Cascade (2026-07)

Companion to docs/silence_rank_lemma.md. Everything below is verified against
the shipped ledger with independent code; numerical claims carry residuals.

## 1. The unified formula (T1's first-order form, both regimes)

Object: at a certified state ψ over vertex λ with support S, the non-gauge
first-order fixed-spectrum tangent dimension

    dim_1 = 2|S| − g − rank(D),   D(δc) = ⊕_i P_i δρ(δc) P_i ⊕ (norm),

with P_i the spectral projectors of ρ and g the gauge rank. The two census
regimes are combinatorial evaluations of rank(D):

- CANCELLATION regime (diagonal ρ, silent channels; the silence-rank lemma):
  dim_1 = (dim K − rank J_within^Re) + (|S| − g − rank J_within^Im).
- MAGNITUDE regime (blocky ρ, active channels): each independent loop
  contributes one weight direction (along the incidence kernel) and one free
  phase direction; touched rigid 1-term classes subtract.

## 2. Corpus verification (all 156 interference states, independent sweep)

Measured non-gauge tangent vs incidence-kernel dimension:

    kdim=0 : tangent 0 for 82 states, tangent 1 for 11 states
             (the 11 carry a free block phase with no weight kernel)
    kdim=1 : tangent 2 for 46 states (1 weight + 1 phase), including v_B,
             whose measured dim_1 = 2 independently reproduces the archive's
             "reduced spectrum-fiber is a SURFACE";
             tangent 0 for exactly ONE state: v89, the cancellation-regime
             rigid case, where the within-block silence condition eats the
             kernel (lemma corollary): the lone crossover point of the two
             regimes in the census.
    kdim=2 : tangent 4 for all 15 states (2 + 2)
    kdim=5 : tangent 6, v103 only (5−2 weights + 5−2 phases; lemma formula)

The loopy/loop-free split (63 states with kdim > 0) matches the archive's
loop census (63 loopy supports) exactly. The regime picture in one line: the
tangent dimension is (surviving weight-kernel) + (surviving free phases), and
the two published laws are the two ways channels kill those sectors:
magnitude pinning vs cancellation silence.

## 3. v103: the sparsification cascade (fiber escapes the certified support)

Continuation on the exact-spectrum variety ((ρ−aI)(ρ−bI) = 0 certificate,
gauge pinned), driving the smallest amplitude's squared modulus to zero and
iterating:

- round 0: det (2,6,9) zeroed, residual 3.3e-16 → SUPPORT-13 exact-spectrum
  state exists.
- round 1: det (1,3,6) zeroed, residual 3.9e-15 → SUPPORT-12 state exists.
- round 2: det (2,3,6) would not zero on this path (floor 1.4e-6); greedy
  cascade stops at 12. Other fiber branches untested.

Consequences:
- s_Q(v103) ≤ 12 (was: certified support 14). The route-B "isolated
  support-14 point, dropping any amplitude breaks feasibility" is a statement
  about that search's landscape/gauge slice, NOT about the fiber: the fiber
  is a 6-manifold whose boundary strata contain strictly sparser states.
- The defect hierarchy at the closures now reads
      v89 :  s_M = 4 < s_I = 7 ≤ s_Q ≤ 10
      v103:  s_M = 4 < s_I = 6 ≤ s_Q ≤ 12
  with s_M the (weak) majorization baseline, s_I the diagonal-incidence
  baseline (exact MILP + rational certificates), and only s_M ≤ s_Q
  guaranteed a priori. Ladder status for v89 at size 7: 41 supports
  enumerated (incomplete), 34 killed exactly by the projector-pattern
  argument, 6 distinct survivors all floor at 1.5e-2..6e-2 over 24
  multi-starts (strong numerical evidence, exactifiable by the same branch
  analysis that killed the v103 witness).

## 4. Exactification leads for the v103 fiber (Problem 2 input)

- INVARIANT COORDINATE: the (0,2,7) weight equals EXACTLY 5/34 (the small
  spectrum value itself) at the certified point AND at every point of the
  cascade including the support-12 endpoint (err < 3e-17). A fiber-constant
  weight is a candidate exact invariant of the reduced fiber and a
  coordinate to quotient out before computing the fiber ideal.
- The support-12 endpoint's other weights are not small-denominator
  rationals; the distinguished exact point of the support-12 stratum (if
  any) remains to be found: PSLQ over quadratic surd lattices is the next
  tool, seeded from this numerical point.
- The germ is a smooth 6-manifold at the certified point (constant tangent
  rank at continued points, multiplicities (4,6) preserved); the cascade
  shows its closure is stratified by support, directly instantiating the
  quantum-stratification program on the census's own closing vertex.

## 5. Corrections this file induces elsewhere (merge notes)

- "Cancellation rigidity exists (v89, v103)" → v89 only (fixed-spectrum);
  v103 is fixed-ρ-rigid and fixed-spectrum-deformable (lemma).
- "Isolated support-14" for v103 → landscape statement; fiber statement is
  FALSE (support 12 reached on-fiber).
- The census min-support ledger column (support of the certified state) is
  an UPPER bound on s_Q, now known strict at v103; treat it as such wherever
  quoted.
