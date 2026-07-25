# RESEARCH_NEXT.md
> **Purpose:** Candidate replacement for `RESEARCH.md`.
>
> This document is intentionally **separate** from the existing `RESEARCH.md` so it can be merged by hand.
> `RESEARCH.md` is to be ARCHIVED on merge (renamed `RESEARCH_ARCHIVE.md`), not deleted:
> it is the evidence log this document's status tags point into. See Section 9 for
> content that MUST survive the merge in some living document.
>
> House rule inherited from the archive (non-negotiable): every 🟦/❌ status tag
> carries a pointer to its certifying artifact (script, jsonl, or archive section).
> A claim with no in-repo artifact cannot carry an evidence tag.
>
> Status legend:
> - ✅ Theorem (proved in project or standard reference)
> - 🟨 Theorem candidate (appears to follow from standard mathematics but needs a written proof)
> - 🟦 Confirmed by certified project computations
> - ❌ Disproved by certified project evidence
> - 🔬 Open
> - 💡 Speculative

# 1. Mission

Develop a mathematical theory of the **inverse geometry of the fermionic one-body spectral (moment) map**, centered on fibers and sparse realization rather than isolated representative states.

# 2. Central Objects

For a determinant support S:

- Incidence equations: A_S p = λ
- Coherence equations: C_S(p, θ) = 0

Support-restricted fiber:

F_{S,λ} = {(p,θ): A_S p = λ, C_S(p,θ)=0}/gauge

**Gauge convention (mandatory precision).** "/gauge" means the quotient by
U(1)^d single-particle phases and the global phase ONLY. The further quotient by
the full degeneracy stabilizer U(m_1) × ... × U(m_k) of the target spectrum
defines the REDUCED fiber F_red, a different object. Every fiber-dimension
statement in this program must name which of the two it is about. This is not
pedantry: the fixed-ρ fiber vs fixed-spectrum fiber conflation was an actual
bug (levi.py audit, RESEARCH_ARCHIVE "BUGS FIXED" item 4), and the T1 upper
bound fails precisely because vertices are singular values of the moment map
(THEOREMS.md T1).

Working philosophy:

- Image geometry (GPC polytope) describes which spectra are possible.
- Fiber geometry describes ambiguity, rigidity, coherence, deformation, and sparse realization complexity.

# 3. Current State of the Theory

## 🟨 Theorem candidates

- Compact semialgebraic support layers Δ^(k). (Route: Δ^(k) is the image of a
  compact semialgebraic set under a polynomial map; Tarski–Seidenberg +
  compactness of continuous images. Near-free; write it.)
- Finite inverse stratification. (Route: already reduced in THEOREMS.md T7
  item 1 to Hardt trivialization + Prop 1 (semialgebraicity of the moment map,
  a polynomial map on the unit sphere). The written-proof gap is Prop 1 plus a
  hypothesis check, not new mathematics. Cross-reference T7; do not re-derive.)
- Fixed-support realization is decidable by real quantifier elimination.
  (Route: the fixed-support fiber equations are polynomial; Tarski. Near-free.)
- Local constancy of unreduced fibers (via semialgebraic triviality).
  (Route: Hardt again — local triviality over each stratum. Same dependency
  as finite stratification; prove them together.)

🔬 Universal spectral-distance → fidelity lower bound. NOT yet a theorem
candidate: no precise statement exists in the repo (direction, distance
measures, and quantifiers unfixed). Candidate route: contractivity of the
partial trace + Weyl/Lidskii eigenvalue perturbation + Fuchs–van de Graaf.
Promote to 🟨 only once the inequality is written with quantifiers.

These should not be claimed as established until formal proofs are written into the project.

## 🟦 Confirmed by project evidence

- Positive-dimensional REDUCED vertex fibers exist (v_B). [THEOREMS.md
  "v_B ANOMALY RESOLVED" (reduced spectrum-fiber is a surface);
  RESEARCH_ARCHIVE "THE FIBER IS THE 1-RDM → 2-RDM INFORMATION GAP";
  scripts/fiber_sampler.py, scripts/vb_fiber_ideal.py.]
- Cancellation rigidity exists (v89): the certified state is first-order
  fixed-SPECTRUM rigid with no touched 1-term class — its single silent
  channel is within-block and has full rank on the 1-dim incidence kernel;
  full tangent kernel = gauge (dim 9), verified independently 2026-07.
  CORRECTION (2026-07, supersedes the P4 scorecard's reading for v103): v103
  is first-order rigid only in the fixed-ρ sense. Three of its five silent
  channels are CROSS-BLOCK (9/17 vs 5/34 eigenspaces) and impose no
  first-order fixed-spectrum constraint (degenerate perturbation theory:
  cross-block δρ rotates eigenvectors, not eigenvalues). The certified v103
  state is fixed-spectrum DEFORMABLE: 6 non-gauge first-order directions,
  none of them in-span stabilizer-orbit directions (all 42 off-diagonal
  U(4)×U(6) generators leave the support), and three independent
  deformation curves continued to t = 0.10 at residual ~1e-17 with the exact
  (ρ−aI)(ρ−bI) = 0 certificate — the deformations switch on cross-block
  1-RDM coherences while holding the spectrum exact. "Isolated support-14
  point" is a statement about the route-B SEARCH LANDSCAPE and its gauge
  slice, not about the fiber. This is the fixed-ρ vs
  fixed-spectrum conflation (levi.py bug log, item 4) recurring in the P4
  scoring; the archive's P4 narrative should be annotated on merge.
  The distinction is now STRUCTURAL in the code rather than in prose:
  gpc_census.fiber exposes fixed_rho_tangent and fixed_spectrum_tangent,
  each naming its fiber in the JSON keys it writes, and a regression test
  fails on any unqualified fiber-dimension name.
  [scripts/reproduce_next_claims.py::v103_deformability;
  src/gpc_census/fiber.py; tests/test_fiber_notions.py;
  scripts/reconcile_claims.py; results/data/reconciliation.json]
- The two fibers differ in the MAGNITUDE regime too, not only in the
  cancellation regime where the gap was discovered. Found out of sample:
  at (3,11) v43 (transport of (3,10) v108, one channel with 1-RDM
  off-diagonal 0.2) the fixed-ρ tangent is 1 and the fixed-spectrum
  tangent is 2. Any text presenting the two-fiber gap as a
  cancellation-regime phenomenon is too narrow.
  [docs/prereg_bracket_3_11_3_12.md B3, scored FAIL;
  results/data/oos_bracket.json]
- The census support column is an UPPER BOUND on minimal support s_Q, and
  it is strict at 71 of the 799 certified states, all INTERFERENCE. The
  census-wide cascade removes 89 amplitudes in total; v103 goes 14 → 10
  (superseding the ≤ 12 of docs/unified_fiber_dimension.md), 55 states drop
  one amplitude, 15 drop two. NUMERICAL except where noted: certificate
  residuals ≤ 2.2e-16 with eigenvalue multiplicities gated at each round,
  not proofs.
  [docs/cascade_census.md; src/gpc_census/cascade.py;
  scripts/census_cascade.py; results/data/cascade.jsonl;
  results/data/cascade_summary.json]
- ✅ 69 of the 71 improved support bounds are EXACT, not numerical. The v103
  recipe (high-precision Gauss-Newton refine → integer relation detection on
  the squared weights → gauge-fix → exact characteristic-polynomial identity)
  certifies 69 cascade endpoints: 55 one determinant below the library
  support, 13 two below, v103 four below. The two holdouts, (5,10) v144 and
  (5,10) v256, are partly rational (four of eight weights at denominator 37)
  and partly not, the signature of a stratum whose distinguished point the
  cascade missed, exactly as v_B's real states sit at quadratic surds;
  extending recognition to quadratic surds is the open piece.
  [scripts/exactify_cascade.py; results/data/cascade_exact.json]
- ✅ s_Q(v103) ≤ 10, EXACT. The support-10 cascade endpoint has all ten
  squared weights rational (13/34, 5/34, 53/442, 195/1802, 5/68, 5/68,
  35/901, 1/34, 18/901, 84/11713), loop holonomy π so it is real up to
  gauge, and det(ρ − xI) = (x − 9/17)^4 (x − 5/34)^6 verified in exact
  rational and radical arithmetic. State denominator 46852 = 2²·13·17·53
  against the library state's 195364 = 442² = 2²·13²·17²: the prime 53 is
  new, and appears in neither the spectrum denominator nor the library
  state, which is why grid searches over multiples of the natural
  denominator could not reach it. Upper bound only; not a minimality claim.
  [scripts/v103_endpoint_precision.py; scripts/v103_support10_certificate.py;
  tests/test_v103_support10.py; results/data/v103_support10_certificate.json]
- The cancellation regime (exactly diagonal 1-RDM with silent channels) has
  exactly TWO members across all 799 certified representatives, v89 and
  v103. The census-wide cascade adds none, so the silence-rank lemma's
  genericity question stays on a sample of size two. Note the regime is a
  property of the REPRESENTATIVE, not the vertex: v103 leaves it by cascade
  round 4. [results/data/cascade_summary.json]
- Real boundary states can exist inside fibers whose generic/dense
  representatives are complex. [Both rank-10 closures are all-real with
  rational squared weights, found by reweighted-L1 continuation after the
  contraction attack returned only dense points; RESEARCH_ARCHIVE "RANK-10
  CLOSED ... the near-theorem is FALSE".]
- Kernel-dimension = silent-channel-count coincidence at both closure
  vertices (1 = 1 at v89, 5 = 5 at v103). [Same artifacts as above;
  flagged in the archive as "worth a lemma of its own" — see Problem 2b.]

- Coherence obstruction at the incidence optimum exists (the witness-level
  content of "positive realization defect"), NOW WITH IN-REPO CERTIFICATES
  (2026-07, scripts/reproduce_next_claims.py). DEFINITION (adopted): for a
  vertex spectrum λ, s_I(λ) = min |S| such that the diagonal incidence system
  A_S p = λ, p ≥ 0, Σp = 1 is feasible; s_Q(λ) = min support of a pure state
  with spec(ρ) = λ; realization defect = s_Q − s_I. CERTIFIED: s_I(v89) = 7
  and s_I(v103) = 6, exact (MILP optimum + exact rational witness weights),
  against certified quantum supports 10 and 14; AND both minimal incidence
  witnesses are PROVED quantum-infeasible in exact arithmetic — v89's by a
  projector-algebra contradiction ((P²)₀₄ = P₀₈P₈₄ with both factors forced
  nonzero by singleton channels while P₀₄ = 0), v103's by exhausting the
  eigenvalue-assignment branches (scalar-block kill; multiplicity kill; and
  the surviving rank-one branch's necessary magnitude system has zero
  nonnegative real solutions, exact solve).

🔬 STILL OPEN, do not overclaim: the GLOBAL statement s_Q > s_I (defect > 0)
  requires killing every support of size s_I..s_Q−1, not just the minimal
  witnesses; v89's s_Q = 10 currently rests on the archive's NUMERICAL
  support-9 floor. Completion path: enumerate incidence-feasible supports per
  size (MILP no-good cuts) and automate the projector-pattern +
  branch-analysis kills. CAVEAT on the baseline: Schur–Horn permits
  diag(ρ) ≺ λ, so the airtight necessary-condition baseline is
  s_M(λ) = min |S| with A_S q ≺ λ (majorization), s_M ≤ s_I; the defect
  claim against s_M is strictly stronger and untested.

## ❌ Disproved

- Fiber dimension = incidence-kernel dimension, AS A UNIVERSAL LAW. [P4, both
  closures.] SCOPE NOTE: the law was fit on magnitude-target states and fails
  in the cancellation regime; the scoped (magnitude-sector) version is the
  surviving candidate and is exactly the unobstructed sector of the T1
  Jacobian program. Do not discard the scoped version with the universal one.
- Singleton channels are the only rigidity mechanism. [P4: both closures
  rigid via multi-term cancellation constraints, no 1-term class.]
- Universal ≤2-channel interference law. [P2: v103 carries 5 exchange
  channels; independently re-verified.]
- State denominator = spectrum denominator (and its ratio-∈{1,2} refinement).
  [P3: ratios 5 (v89) and 5746 (v103), rational-weight states, no escape
  clause. LEDGER NOTE 2026-07: the regenerated ledger ships all 12
  DESIGN-REAL states at ratio 2 (archive's mined "9 of 12" is stale —
  representative-dependent, as the P2/P3/P4 caveat predicts). Re-derived in
  exact rational arithmetic, per state, by scripts/reconcile_claims.py →
  results/data/reconciliation.json.]
- Complex representative ⇒ intrinsically complex fiber. [Both closures real;
  archive "the near-theorem is FALSE".]
- Positive-dimensional fibers are literal torus orbits. [Archive retraction:
  rational curve fibered over the gauge torus, not a torus orbit.]
- Generic one-loop fibers are elliptic. [Archive retraction; census of
  kdim-1 families: 18 genus-1 vs 26 genus-0.]
- Dense numerical representatives imply no sparse exact representative.
  [The contraction attack lands on dense generic fiber points; reweighted-L1
  over the full fiber found support 10 and 14; archive "WHAT OVERTURNED THE
  NEAR-THEOREM".]

- Two spectral classes require higher-order interference — under the ADOPTED
  reading "spectra with exactly two eigenvalue blocks force complex (or
  higher-degree algebraic) interference phases". DISPROVED decisively by
  ledger census (2026-07, scripts/reproduce_next_claims.py): ALL 64 two-block
  vertices ship all-real certified states, including all four two-block
  interference vertices (v89, v103, v108, and one (5,10) vertex); trisected
  (cube-root) phases occur only at 3–5 blocks. Existence of a real state
  refutes necessity, so this direction is representative-safe.
- Spectral border support differs from exact support — under the ADOPTED
  reading "the minimal support attainable as a LIMIT of exact-spectrum states
  can be smaller than the minimal exactly-attained support". DISPROVED at the
  tested instance (v89): the certified support-10 state is the c₍₀₃₆₎ → 0
  endpoint of a support-11 one-parameter exact-spectrum family (11th
  determinant (0,3,6); continuation verified to residual ~1e-17 over
  c_x ∈ [0, 0.2]; the other five first-order openers are second-order
  obstructed), and the endpoint IS attained — border = exact here. NOTE: both
  readings above are reconstructions; if the parallel instance's original
  definitions differ, re-tag against those.

# 4. Highest-Value Open Problems

(Each problem cites its existing groundwork; new agents start there, not from zero.)

1. Jacobian formula for local fiber dimension. [= THEOREMS.md T1. Proof route
   on file: MGS local normal form in the Sjamaar–Lerman stratified setting;
   the obstruction quadric is computable per state. The cancellation regime
   (v89/v103) is the new boundary case the formula must absorb.]
2. Explicit reduced-fiber invariant for v_B. [Groundwork: v_B anomaly
   resolution, holonomy/wall data in REVISION_NOTES.md, vb_holonomy.py,
   vb_fiber_ideal.py.]
   2b. Silence-rank lemma — STATEMENT NOW FIXED BY DATA (2026-07), proof is
   ~1 page of degenerate perturbation theory, assigned: at a diagonal-1-RDM
   state over a vertex, first-order fixed-SPECTRUM weight deformations are
   ker A_S ∩ ker J_within (within-block silent channels only), while
   first-order fixed-ρ deformations are ker A_S ∩ ker J_all. Corollaries:
   v89 spectrum-rigid; v103 fixed-ρ-rigid but spectrum-deformable (≥3-dim
   verified fiber germ, tangent dim 6). The v103 kernel-dim = channel-count
   coincidence (5 = 5) is a fixed-ρ statement. Remaining open: the
   second-order structure (which of the 6 tangent directions integrate — at
   least 3 do), and whether within-block silence rank = kernel dim is forced
   for census cancellation states or accidental (sample of size 2).
3. Interior inverse walls. [Groundwork: wall_solver.py, wall_test.py,
   per-family endpoint theorem (T4).]
4. Graph-theoretic coupled-phase realizability criterion. [Groundwork:
   one-hop class combinatorics in classify.py; loop census in
   REVISION_NOTES.md.]
5. Infinite positive-defect family. [BLOCKED on the 💡 defect definition
   above; do not start before Problem 0 of Section 7.]
6. Infinite deformable and cancellation-rigid families. [Candidate seeds:
   the 44 kdim-1 deforming families; the two closure states.]
7. Levi-theoretic lower bounds. [Groundwork: levi.py post-audit, the
   rank_deficiency == stab − 1 identity (799/799), the hard/soft
   plethysm criterion.]
8. Complexity classification. [BLOCKED on choosing the complexity measure:
   support size, state-denominator, phase field degree, and clique budget
   are all in play and provably non-equivalent (den-ratio falsification).
   Define the measure first.]

# 5. Dependency Graph

Incidence geometry
    ↓
Weight polytope
    ↓
Coherence equations
    ↓
Jacobian
    ↓
Fiber dimension / rigidity
    ↓
Inverse walls
    ↓
Support complexity atlas

# 6. Dead Ends

Future agents should NOT assume:

- kernel dimension predicts fiber dimension (outside the magnitude-target scope);
- representative properties characterize an entire fiber (P2/P3/P4 and the
  projective-stabilizer demotion are all representative-dependent);
- toric-orbit descriptions;
- denominator preservation;
- two-channel universality;
- that a dense numerical fiber point says anything about sparse exact points
  (the contraction-attack artifact that nearly closed the census wrongly).

# 7. Immediate Tasks

Priority:

0. Land scripts/reproduce_next_claims.py (added with this document): exact
   s_I certificates + witness infeasibility proofs (defect), the two-block
   phase census, and the v89 support-11 family continuation. Confirm the
   adopted definitions against the parallel instance's originals; the only
   remaining unevidenced piece is the GLOBAL defect inequality (see 🔬 above).
1. Formalize theorem-candidate proofs (the two Hardt-dependent ones together).
2. Prove Jacobian dimension theorem (T1), cancellation regime included.
3. Prove the silence-rank lemma (Problem 2b).
4. Produce stabilizer-invariant for v_B.
5. Search for interior inverse walls.
6. Develop practical realization algorithms from the corrected theory
   (reweighted-L1 full-fiber continuation → PSLQ exactification is the
   validated pipeline; productionize scripts/repro_and_sparsify.py-class
   tooling out of scratch). PARTLY DONE: gpc_census.cascade is the
   productionized minimal-polynomial continuation, run census-wide in
   scripts/census_cascade.py; the PSLQ exactification half is not.
7. Certify the cascade endpoints (alpha-theory or exact recognition). The
   best target this program has produced is the v103 support-10 endpoint:
   isolated (first-order rigid in BOTH fiber senses), correct eigenvalue
   multiplicities, five of ten squared weights already exact rationals at
   denominator 68 with the remaining five summing to exactly 5/17, and the
   (0,2,7) weight still exactly 5/34. Route: high-precision re-solve seeded
   from results/data/cascade.jsonl, then integer relation detection.
   [docs/cascade_census.md]
8. Test the candidate law "positive fixed-spectrum tangent ⇔ the cascade
   sparsifies" out of sample. It holds 726/726 in the rigid direction and
   71/73 in the deformable direction across the census, with both exceptions
   flooring at 2 to 3e-6, i.e. blocked on the path rather than blocked. It is
   a post-hoc pattern in the corpus the cascade itself explored and must be
   pre-registered before it is believed. [docs/cascade_census.md]

# 8. Guiding Principle

The census should now be viewed as **experimental evidence**.

The primary research product is the mathematical theory of inverse fibers, lifting complexity, and coherence obstructions.

The falsification protocol is retained as method, not history: mined-law
mortality ran about one in two, and both rank-10 closures falsified
pre-registered laws with mechanism. Every new fiber law gets a pre-registered
scored prediction before the narrative is written, and a FAIL is a publishable
finding. (Protocol text: RESEARCH_ARCHIVE "SKEPTICAL AUDIT + PRE-REGISTRATION".)

# 9. Merge checklist — content that must survive from RESEARCH.md

Archive the file; keep these alive (here or in the named homes):

- The validation law (structural invariants gate every ship) — merge into
  AGENTS.md if not already there.
- The pre-registration protocol and final P1–P6 tally — Section 8 pointer
  plus the archive.
- The (3,11)/(3,12) campaign state: bracket settlement, partial-families
  oracle, plethysm inner hull, and the Stage-1 spec pointer — these are the
  constraint-generation frontier and are NOT superseded by the fiber program.
- The source-document index (docs/ inventory with caveats, e.g. the
  ak2008_extracted.md surd re-check warning).
- Adjacent-literature boundaries (Liebert et al. spin-adapted GPCs;
  Maciążek–Tsanov inner bounds).
- Known-stale items to fix on merge, not copy: the "785 of 799" paragraph
  (ledger is 799/799, CI-green), and the "9 of 12 DESIGN-REAL at ratio 2"
  mining figure (shipped ledger: 12 of 12).
- ANNOTATE on merge, do not silently copy: the P4 scorecard's "v103
  FIRST-ORDER RIGID" is a fixed-ρ verdict; the fixed-spectrum fiber at v103
  is positive-dimensional (see the corrected cancellation-rigidity item and
  Problem 2b). P4's v89 verdict stands in both senses.

## Reconciliation ledger (2026-07 fiber session)

All four items below are re-derived by scripts/reconcile_claims.py into
results/data/reconciliation.json, so a merge can check them rather than
trust this list.

| superseded reading | corrected reading |
|---|---|
| v103 is first-order rigid | rigid for the FIXED-RHO fiber, deformable (tangent 6) for the FIXED-SPECTRUM fiber; v89 rigid in both |
| v103 is an isolated support-14 point | true of the route-B search landscape and gauge slice, false of the fiber: support 10 is reached and CERTIFIED exactly |
| the census support column is the minimal support | it is the certified support, an upper bound on s_Q, strict at 71 of 799 vertices |
| 9 of 12 DESIGN-REAL states at denominator ratio 2 | the shipped ledger has 12 of 12; the mined figure is representative dependent |

Two further scope corrections from this session, both from evidence rather
than re-reading:

- The fixed-ρ / fixed-spectrum gap is not confined to the cancellation
  regime. It appears in the magnitude regime at (3,11) v43 (fixed-ρ 1,
  fixed-spectrum 2), found by the out-of-sample holdout.
  [docs/prereg_bracket_3_11_3_12.md B3]
- The cancellation regime is a property of the REPRESENTATIVE, not the
  vertex: v103's own cascade leaves the regime. "Exactly two cancellation
  states" counts certified representatives. [docs/cascade_census.md]
