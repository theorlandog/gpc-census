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
- Cancellation rigidity exists (v89, v103): certified states are first-order
  rigid with NO touched 1-term class; the silence conditions have full rank on
  the incidence kernel. [RESEARCH_ARCHIVE pre-registration scorecard P4, both
  closures; docs/hybrid_cracks/v89.jsonl, v103.jsonl;
  scripts/verify_hybrid_state.py. Kernel dims 1 and 5 and channel class sizes
  {2} and {4,3,2,2,2} independently re-verified 2026-07.]
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
  representative-dependent, as the P2/P3/P4 caveat predicts).]
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
   2b. Silence-rank lemma: when do the zero-target channel conditions have
   full rank on the incidence kernel (first-order rigidity), and is the
   kernel-dim = channel-count coincidence at v89/v103 forced or accidental?
   [Finite exact rank statement; most tractable new item on this list.]
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
   tooling out of scratch).

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
