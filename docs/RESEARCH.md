# Research state and continuity notes

This is the living research document. Read it before touching the science.
House rules live in AGENTS.md. The full evidence log of the census era
(pre-registration scorecards, campaign session notes, retractions with
mechanism) is frozen in `docs/RESEARCH_ARCHIVE.md`; status tags below point
into it by section name. The archive is the record; this file is the map.

House rule inherited from the archive (non-negotiable): every 🟦/❌ status
tag carries a pointer to its certifying artifact (script, jsonl, doc, or
archive section). A claim with no in-repo artifact cannot carry an evidence
tag.

Status legend:

- ✅ Theorem (proved in project or standard reference)
- 🟨 Theorem candidate (appears to follow from standard mathematics but needs a written proof)
- 🟦 Confirmed by certified project computations
- ❌ Disproved by certified project evidence
- 🔬 Open
- 💡 Speculative

# Part I. The census program (established base)

## What this project is

Complete classification program for fermionic moment polytopes (generalized
Pauli constraints, GPCs). The pure-state N-representability problem: which
one-body occupation spectra arise from N-fermion pure states in d orbitals.
Altunbulak and Klyachko (CMP 282, 287 (2008)) computed constraint systems
through rank 10 and left two vertices of the wedge^4 H_9 polytope resolved
only numerically. This program closed both (paper in results/report), then
classified every vertex of every known system, and now extends the frontier.

## The trichotomy (introduced by this program)

Every polytope vertex is exactly one of:

- DESIGN-INT: an integer weighted design exists at the natural denominator
  (nonnegative integer determinant weights, prescribed mode sums, support
  one-hop free, meaning no two support determinants share N-1 modes).
- DESIGN-REAL: a real nonnegative design exists but no integer one.
- INTERFERENCE: neither exists; complex phase cancellation is forced.

Verdicts are solver certificates. The one-hop exclusion is essential; a
degree-sums-only model is strictly weaker and misclassifies (this bug was
made and caught once; see tests pinning v_A and v_B).

## Key results encoded in results/data

- Census complete for every determinate system, ranks 6 through 10
  (799 vertices; duals covered by particle-hole verdict transport).
- Certified closed-form extremal states for ALL 799 vertices (ledger
  regenerated, CI-green): every design vertex from its witness, every
  interference vertex through the routed pipeline (per-phase recognition,
  the constructive off-diagonal-target exactifier of stage 3b, PSLQ),
  14 closed by state transport (scripts/transport_states.py), and the two
  rank-10 cancellation-regime closures v89 and v103, found off the rational
  search grid after the attainability audit self-corrected
  (RESEARCH_ARCHIVE "RANK-10 CLOSED").
- v_A (16,16,16,6,6,6,6,6,6)/21 is DESIGN-INT; v_B (20,14,14,14,14,4,4,4,4)/23
  is INTERFERENCE with cos(gamma) = 3/(4*sqrt(14)). The v_B phase is not
  special: it is the k=2 instance of the off-diagonal-target mechanism that
  fixes every interference phase.
- Interference onset: absent at rank 8 in the N=4 series, first appears at
  rank 9. Fractions stabilize across ranks 9 to 10 within each series
  (34-37 percent at N=3, 16 at N=4, 14 at N=5).
- Functorial structure: trailing-zero padding and frozen-core lifts
  (lambda -> (1, lambda)) preserve verdicts and generate most interference
  from lower-rank originals; at (4,10) they generate all of it.
- LEDGER NOTE (2026-07): all 12 DESIGN-REAL states ship at
  state-den/spectrum-den ratio 2 (the archive's mined "9 of 12" figure is a
  stale snapshot; the ratio is representative-dependent, exactly as the
  P2/P3/P4 caveat predicts). The census min-support column (support of the
  certified state) is an UPPER bound on the true minimal quantum support
  s_Q, now known strict at v103; treat it as such wherever quoted.

## The validation law (non-negotiable)

Every real error in this program was caught by structural invariants, never
by inspection: a sign-parse bug (caught by particle-hole self-duality of
(5,10)), a wrong classify port (caught by a Farkas feasibility
contradiction), a wrong solver gradient (caught by finite differences), and
the levi.py conjugation and fiber-conflation bugs (caught by the invariance
gate; RESEARCH_ARCHIVE "BUGS FIXED IN src/gpc_census/levi.py").
Consequences:

- No result ships without passing gpc_census.validate and the test suite.
- The preflight gate: any state-solving change must reconstruct v_B and
  certify its closed form end to end (scripts/solve_all.py --preflight)
  before campaigns run. The gate uses the weights-first solver
  (solve_vertex_exact_first with certify_tier_b): Tier A must attain the v_B
  spectrum exactly and Tier B must exactify to a certified closed form. It
  completes in minutes. The historical attain-based path
  (--legacy-preflight) is a Tier-A-only regression check, not the gate.
- Pin solvers. Census verdicts were produced under ortools 9.15.6755
  (a required runtime dependency); the CBC fallback labels its backend in
  every verdict.

## The falsification protocol (retained as method, not history)

Every new fiber law gets a pre-registered scored prediction before the
narrative is written, and a FAIL is a publishable finding. Mined-law
mortality ran about one in two, and both rank-10 closures falsified
pre-registered laws with mechanism. Protocol text: RESEARCH_ARCHIVE
"SKEPTICAL AUDIT + PRE-REGISTRATION". Final tally over the closure vertices
(RESEARCH_ARCHIVE "PRE-REGISTRATION SCORED: v103"): P1, P2, P5 pass; P3 and
P4 falsified with mechanism; P6 never in scope. The archive's P4 verdict
for v103 ("first-order rigid") is a fixed-rho statement and is annotated in
place there; the fixed-spectrum fiber at v103 is positive-dimensional (see
the silence-rank lemma below). P4's v89 verdict stands in both senses.

## Current campaigns: the (3,11)/(3,12) constraint-generation frontier

NOT superseded by the fiber program; this is where the polytope frontier
stands.

1. Stage 1, the constraint generator for d >= 11: docs/stage1_klyachko_spec.md
   holds the extracted algorithm (Theorem 3.2.1, cubicle extremal edges,
   Schubert coefficient test via lrcalc). Test-first: it must reproduce the
   five known N=3 systems before any new output counts. Until then
   constraints ship as a lookup (src/gpc_census/data/constraints.json).
2. The (3,11) bracket (docs/bracket_3_11.json): 46 of the 50 outer candidates
   are settled WITHOUT Stage 1, in exact arithmetic and independently
   verified (scripts/settle_bracket_3_11.py regenerates
   docs/bracket_3_11_settlement.json). 19 TRUE (17 by state transport from
   the rank-9/10 census, the uniform (3/11)^11 by an explicit Z11 difference
   design, and a frozen-core paired state); 27 REFUTED (frozen-core pinning
   to N=2 plus the even-degeneracy pairing theorem, and zero-restriction
   against known lower-rank GPCs); 4 OPEN (cands 23, 26, 34, 44, all
   genuinely rank-11 full-support; cand 23 admits no one-hop-free design, so
   INTERFERENCE if attainable). CAVEAT: the 27 refutations do NOT finish
   (3,11) as a polytope; cutting them creates new true vertices absent from
   the outer list, so Stage 1 (or a tighter bracket) is still required. The
   settled points are a mandatory acceptance oracle for any (3,11)
   generator.
3. Partial-family cross-check and the lead on the 4 OPENs: all four violate
   a published but CLAIMED level-5 series inequality (AK-RMK-4.2.1 or its
   lambda_11-extended quadruple), and none violate any PROVED family
   (docs/partial_families_3_11_3_12.json, scripts/partial_families.py). Two
   ways to close them with proof: (a) prove the level-5 inequality inside
   Stage 1 (a finite Schubert coefficient computation); or (b) attain
   idx 44, which would FALSIFY a published Klyachko inequality, itself a
   citable result. Either branch is a win; until then the four stay OPEN
   (validation law: a claim is not a certificate).
4. The (3,12) footholds: 19 certified true vertices by face embedding of the
   (3,11) settlement (docs/bracket_3_12_true_vertices.json), and the exact
   plethysm Stage-0 inner cloud (docs/bracket_3_12_stage0_inner.json,
   scripts/plethysm_inner_hull.py), whose cloud-extremality is NOT
   vertexhood until the cloud converges.

Other census-era leads (vertex-cone arithmetic, holonomy Galois structure,
state transport, the v96 campaign, core/hole isomorphisms, Levi symmetry
theorem and its audit) live in the archive under their section names; the
Levi identity rank_deficiency == stab - 1 (799/799) and the hard/soft
plethysm criterion feed Problem 7 below.

## Source documents in docs/

- RESEARCH_ARCHIVE.md: the frozen census-era evidence log (this file's
  status tags point into it).
- THEOREMS.md: numbered theorem statements (T1 Jacobian program, T4
  endpoint theorem, T7 stratification reduction, v_B anomaly resolution).
- holonomy_rationality.md: the arithmetic of loop holonomies. Theorem A
  (PROVED) settles the Holonomy Rationality conjecture for one-hop channel
  classes of size at most 2, with the adjoined radicand a weight monomial;
  Proposition B (PROVED) is the reduction step at larger classes and yields
  only a 2-tower, which is weaker than multiquadratic, so the conjecture
  stays open there on 82/82 evidence. Also carries the trisection
  dissolution and the analysis of the one non-rational-weight ledger record.
- prereg_denominator_arithmetic.md: the pre-registered denominator-arithmetic
  predictions, committed before the census run, scored FAIL/FAIL/PASS in
  results/data/denominator_arithmetic.json with mechanisms for both FAILs.
- silence_rank_lemma.md: the silence-rank lemma, PROVED (first-order
  statements, degenerate perturbation theory), with verified census
  corollaries: v89 fixed-spectrum rigid, v103 fixed-rho rigid but
  fixed-spectrum deformable.
- unified_fiber_dimension.md: the unified first-order fiber-dimension
  formula (both regimes), its corpus verification over all 156 interference
  states, the v103 sparsification cascade, and merge notes. Its s_Q(v103)
  <= 12 is superseded by <= 10, certified exactly; see cascade_census.md.
- cascade_census.md: the census-wide sparsification cascade, the exact
  certification of 69 of the 71 improved support bounds, the v103
  support-10 state, the cancellation-regime census, and the
  representative-swap check against the v0.12 real v_B record.
- prereg_bracket_3_11_3_12.md: the scored out-of-sample pre-registration on
  the (3,11) and (3,12) bracket holdout (4 PASS, 2 FAIL), including the
  power caveat stated before scoring.
- prior_art_roadmap.md: the outside-GPC prior-art map (phase retrieval,
  rigidity theory, structured inverse eigenvalue problems, numerical
  algebraic geometry, tensor rank, ED degrees, SOS) with an ordered action
  list mapped to the open problems below.
- stage1_klyachko_spec.md: the constraint-generation algorithm, verbatim
  from the Altunbulak thesis, with the implementation plan.
- ak2008_extracted.md: inequality and vertex tables for ranks 6-8 from the
  2008 paper, including extremal states (surd coefficients need re-checking
  against the PDF before use as ground truth).
- bracket_3_11.json / bracket_3_11_settlement.json: the (3,11) Stage-0
  bracket and its 46/50 settlement (19 true, 27 refuted, 4 open);
  reproducible via scripts/settle_bracket_3_11.py.
- partial_families_3_11_3_12.json: the published stable Grassmann families
  for N=3 at ranks 11 and 12, tagged PROVED / WEAK / CLAIMED. NECESSARY
  conditions, an outer bound only; use as a Farkas-style acceptance oracle
  for any Stage-1 output (scripts/partial_families.py).
- bracket_3_12_true_vertices.json: 19 certified true vertices of (3,12) by
  face embedding of the (3,11) settlement.
- bracket_3_12_stage0_inner.json: the (3,12) exact-plethysm inner cloud
  (scripts/plethysm_inner_hull.py); cloud-extremality is not vertexhood
  until convergence.

## Adjacent literature to respect

Liebert, Lemke, Altunbulak, Maciazek, Ochsenfeld, Schilling,
PRResearch 7, 023247 (2025), arXiv:2502.15464: spin-adapted GPCs for
larger systems. Different lattice of problems (spin sectors) from the
spinless census here; draw the boundary explicitly in any new paper.
Maciazek and Tsanov (J. Phys. A 50, 465304): doubly-excited inner bounds
coincide with the polytope for (3, d <= 7); the census confirms this
empirically (those systems are interference-free).
For the fiber program's outside-field siblings, see
docs/prior_art_roadmap.md.

# Part II. The inverse-fiber program (primary research product)

## Mission

Develop a mathematical theory of the **inverse geometry of the fermionic
one-body spectral (moment) map**, centered on fibers and sparse realization
rather than isolated representative states.

## Central objects

For a determinant support S:

- Incidence equations: A_S p = lambda
- Coherence equations: C_S(p, theta) = 0

Support-restricted fiber:

F_{S,lambda} = {(p,theta): A_S p = lambda, C_S(p,theta)=0}/gauge

**Gauge convention (mandatory precision).** "/gauge" means the quotient by
U(1)^d single-particle phases and the global phase ONLY. The further
quotient by the full degeneracy stabilizer U(m_1) x ... x U(m_k) of the
target spectrum defines the REDUCED fiber F_red, a different object. Every
fiber-dimension statement in this program must name which of the two it is
about, and must name FIXED-SPECTRUM vs FIXED-RHO. This is not pedantry: the
fixed-rho fiber vs fixed-spectrum fiber conflation was an actual bug
(levi.py audit, RESEARCH_ARCHIVE "BUGS FIXED" item 4), recurred in the P4
scoring (annotated in the archive), and the T1 upper bound fails precisely
because vertices are singular values of the moment map (THEOREMS.md T1).

Working philosophy:

- Image geometry (GPC polytope) describes which spectra are possible.
- Fiber geometry describes ambiguity, rigidity, coherence, deformation, and
  sparse realization complexity.

## Current state of the theory

### ✅ Proved

- The silence-rank lemma (first-order; docs/silence_rank_lemma.md, proof by
  degenerate perturbation theory). At a diagonal-1-RDM state over a vertex:
  (i) first-order fixed-SPECTRUM weight deformations are
  ker A_S ∩ ker J_within (within-block silent channels only; cross-block
  silence imposes no first-order fixed-spectrum condition);
  (ii) first-order fixed-RHO weight deformations are ker A_S ∩ ker J_all;
  (iii) the non-gauge fixed-spectrum tangent dimension is
  (dim K - rank J_within^Re) + (|S| - g - rank J_within^Im).
  Verified corollaries (scripts/reproduce_next_claims.py): v89 first-order
  fixed-spectrum rigid (formula gives 0, direct tangent kernel = gauge);
  v103 fixed-rho rigid (all-channel silence rank 5 = kernel dim 5) but
  fixed-spectrum DEFORMABLE (formula and direct tangent both give 6).

- Holonomy quadraticity at small channel classes (Theorem A,
  docs/holonomy_rationality.md). If a state has rational squared weights and
  rational channel targets, every one-hop class of size 2 gives a loop
  holonomy whose cosine lies in Q(sqrt(N)) for N a rational MONOMIAL in the
  weights. Proof by expanding the channel modulus; the v_B phase
  cos(gamma) = 3/(4*sqrt(14)) is exactly this formula.
  [docs/holonomy_rationality.md; results/data/holonomy_fields.json]

### 🟨 Theorem candidates

- Compact semialgebraic support layers Delta^(k). (Route: Delta^(k) is the
  image of a compact semialgebraic set under a polynomial map;
  Tarski-Seidenberg + compactness of continuous images. Near-free; write
  it.)
- Finite inverse stratification. (Route: already reduced in THEOREMS.md T7
  item 1 to Hardt trivialization + Prop 1 (semialgebraicity of the moment
  map, a polynomial map on the unit sphere). The written-proof gap is Prop 1
  plus a hypothesis check, not new mathematics. Cross-reference T7; do not
  re-derive.)
- Fixed-support realization is decidable by real quantifier elimination.
  (Route: the fixed-support fiber equations are polynomial; Tarski.
  Near-free.)
- Local constancy of unreduced fibers (via semialgebraic triviality).
  (Route: Hardt again, local triviality over each stratum. Same dependency
  as finite stratification; prove them together.)
- Unified first-order fiber-dimension formula, magnitude and cancellation
  regimes together: dim_1 = 2|S| - g - rank(D) with D the block-projected
  differential (docs/unified_fiber_dimension.md section 1). The
  cancellation-regime evaluation is the proved lemma above; the
  magnitude-regime evaluation (each independent loop contributes one weight
  and one phase direction; touched rigid 1-term classes subtract) is
  corpus-verified but needs a written proof. This is T1's first-order form.

🔬 Universal spectral-distance to fidelity lower bound. NOT yet a theorem
candidate: no precise statement exists in the repo (direction, distance
measures, and quantifiers unfixed). Candidate route: contractivity of the
partial trace + Weyl/Lidskii eigenvalue perturbation + Fuchs-van de Graaf.
Promote to 🟨 only once the inequality is written with quantifiers.

These should not be claimed as established until formal proofs are written
into the project.

### 🟦 Confirmed by project evidence

- Holonomy fields are 2-elementary, census-wide. All 82 gauge-invariant loop
  holonomies of all 63 loopy interference states have cosine of degree 1 (37),
  2 (32) or 4 (13) over Q, with Galois group C1, C2, V4 respectively; no C4,
  D4 or Q8 anywhere. Groups computed with sympy galois_group, not inferred
  from the polynomial having only even terms; loop-kernel bases certified
  saturated so nothing is missed; every PSLQ relation recognized at 120 digits
  and re-verified at 240 with the residual recorded (worst 1e-251).
  PRECISION PROTOCOL, now mandatory for this area: a recognized integer
  relation must be re-verified at at least double the working precision and
  the residual recorded; a relation that fails re-verification is DELETED, not
  downgraded. An exploratory scan at 60 digits produced a false counterexample
  without it. [scripts/holonomy_field_scan.py;
  results/data/holonomy_fields.json; tests/test_holonomy_fields.py]
- The theorem's hypotheses hold census-wide, checked rather than assumed: all
  84 one-hop channels of the loopy states have rational |rho_AB|^2, and for
  all 63 states the within-channel difference vectors generate the loop
  lattice, saturated over Z. The latter is a HYPOTHESIS, not a theorem: a
  support could carry a loop with no one-hop structure, and none does here.
  [results/data/holonomy_fields.json, channel_hypotheses]
- The TRISECTED phase tier is a printing artifact, not a phase-complexity
  tier. The tier is assigned syntactically ("atan" and "/3" both present in
  the printed closed form); the written phase (pi + 3*atan(t))/3 equals
  pi/3 + atan(t), and across all 13 states every representative amplitude has
  2-power degree over Q while every loop holonomy has degree 1, 2 or 4 with
  2-elementary group. No cubic subextension exists anywhere in the class.
  This is stronger than "the cube roots are a gauge artifact"; there was never
  any cubic algebra to be a gauge artifact of.
  [results/data/holonomy_fields.json, trisection_audit]
- Denominator minting is a two-stage arithmetic. The incidence solve
  A_S p = lambda, sum p = 1 and then the pinning solve by the channel
  conditions; the state denominator's new primes come from the second stage.
  At v89 the single silence condition restricted to the incidence line loses
  its quadratic terms and the surviving linear equation is literally
  130*p9 - 1 = 0, minting the prime 5. At v103 the same cancellation happens
  BETWEEN two channel equations, whose difference is 13*p12 = p9 + p13,
  minting 13; the incidence solve there is unimodular over the spectrum
  denominator 34 and mints nothing.
  [scripts/denominator_arithmetic.py; results/data/denominator_arithmetic.json]
- Exactly ONE of the 799 certified records has non-rational squared weights:
  the shipped v_B record at (4,9) index 65, whose weights are quadratic
  irrationals in Q(sqrt 35). It is a genuine non-rational-weight state, not a
  formatting inconsistency: it is the endpoint of a wall family cut out by the
  quadratic 85t^2 - 70t - 119, whose two roots are a Galois conjugate pair
  with only one physical, so the distinguished point has orbit size 2. Its own
  loop holonomy is a sign (minimal polynomial x + 1), consistent with the
  record being real up to gauge. It is excluded from denominator predictions,
  having no state denominator. [docs/holonomy_rationality.md;
  tests/test_holonomy_fields.py]

- Positive-dimensional REDUCED vertex fibers exist (v_B). [THEOREMS.md
  "v_B ANOMALY RESOLVED" (reduced spectrum-fiber is a surface);
  RESEARCH_ARCHIVE "THE FIBER IS THE 1-RDM -> 2-RDM INFORMATION GAP";
  scripts/fiber_sampler.py, scripts/vb_fiber_ideal.py.]
- Cancellation rigidity exists (v89): the certified state is first-order
  fixed-SPECTRUM rigid with no touched 1-term class; its single silent
  channel is within-block and has full rank on the 1-dim incidence kernel;
  full tangent kernel = gauge (dim 9), verified independently 2026-07.
  v103 is NOT a second instance: it is first-order rigid only in the
  fixed-rho sense; three of its five silent channels are cross-block and
  the certified state is fixed-spectrum DEFORMABLE, with 6 non-gauge
  first-order directions (none of them in-span stabilizer-orbit directions)
  and three independent deformation curves continued to t = 0.10 at
  residual ~1e-17 under the exact (rho-aI)(rho-bI) = 0 certificate.
  [docs/silence_rank_lemma.md;
  scripts/reproduce_next_claims.py::v103_deformability; archive P4 sections
  annotated accordingly.]
- Corpus-wide first-order tangent census (all 156 interference states):
  tangent dimension = (surviving weight kernel) + (surviving free phases);
  the loopy/loop-free split matches the archive loop census exactly; v89 is
  the lone crossover of the two regimes (kdim 1, tangent 0); v_B's measured
  dim 2 independently reproduces the archive's "reduced spectrum-fiber is a
  surface". [docs/unified_fiber_dimension.md section 2.]
- v103 sparsification cascade: on-fiber continuation reaches exact-spectrum
  states of support 13 and then 12 (residuals 3.3e-16, 3.9e-15), so the
  certified support-14 state is NOT isolated in the fiber. The fiber germ at
  the certified point is a smooth 6-manifold whose closure is stratified by
  support. [docs/unified_fiber_dimension.md section 3.]
  SUPERSEDED 2026-07: the floor at 12 was a property of that continuation
  path. With a geometric schedule and fallback drop targets the cascade
  continues to support 10, and that endpoint is EXACT, not numerical: all
  ten squared weights are rationals (13/34, 5/34, 53/442, 195/1802, 5/68,
  5/68, 35/901, 1/34, 18/901, 84/11713), the single loop holonomy is pi so
  the state is real up to gauge, and det(rho - xI) = (x - 9/17)^4
  (x - 5/34)^6 holds in exact rational and radical arithmetic. So
  s_Q(v103) <= 10, certified. Its state denominator 46852 = 2^2*13*17*53
  carries a prime found in neither the spectrum denominator nor the
  library state, which is why no grid search over multiples of the natural
  denominator could reach it. Upper bound only, not a minimality claim.
  [docs/cascade_census.md; scripts/v103_endpoint_precision.py;
  scripts/v103_support10_certificate.py; tests/test_v103_support10.py;
  results/data/v103_support10_certificate.json]
- Census-wide sparsification: the support column is an UPPER BOUND on the
  minimal support s_Q, and it is STRICT at 71 of the 799 certified states,
  all INTERFERENCE. The cascade removes 89 amplitudes in total: 55 states
  drop one, 15 drop two, v103 drops four. 69 of the 71 improved bounds are
  then EXACT (high-precision refine, integer relation detection on the
  squared weights, gauge-fix, exact characteristic-polynomial identity);
  the two holdouts, (5,10) v144 and v256, have partly irrational weights,
  the signature of a stratum whose distinguished point the cascade missed,
  as v_B's real states sit at quadratic surds. Quadratic-surd recognition
  is the open piece. [docs/cascade_census.md; src/gpc_census/cascade.py;
  scripts/census_cascade.py; scripts/exactify_cascade.py;
  results/data/cascade.jsonl; results/data/cascade_exact.json]
- Rigidity predicts sparsifiability across the census: every state with
  fixed-spectrum tangent 0 fails to sparsify (726 of 726) and 71 of the 73
  with positive tangent do, the two exceptions flooring at 2 to 3e-6, which
  is blocked-on-this-path rather than blocked. This is a POST-HOC pattern in
  the corpus the cascade itself explored; it needs a pre-registered
  out-of-sample test before it is believed. [docs/cascade_census.md]
- The fixed-rho / fixed-spectrum gap is NOT confined to the cancellation
  regime where it was discovered. It occurs in the magnitude regime at v_B
  itself (fixed-rho 1, fixed-spectrum 2) and, out of sample, at (3,11) v43
  (same 1 vs 2). Any text presenting the two-fiber distinction as a
  cancellation-regime phenomenon is too narrow.
  [docs/prereg_bracket_3_11_3_12.md B3, scored FAIL;
  results/data/oos_bracket.json; results/data/cascade.jsonl]
- The cancellation regime (exactly diagonal 1-RDM with silent channels) has
  exactly TWO members across all 799 certified representatives, v89 and
  v103; the census-wide cascade adds none, so the silence-rank lemma's
  genericity question stays on a sample of size two. The regime is a
  property of the REPRESENTATIVE, not the vertex: v103 leaves it by cascade
  round 4. [results/data/cascade_summary.json]
- Real boundary states can exist inside fibers whose generic/dense
  representatives are complex. [Both rank-10 closures are all-real with
  rational squared weights, found by reweighted-L1 continuation after the
  contraction attack returned only dense points; RESEARCH_ARCHIVE "RANK-10
  CLOSED ... the near-theorem is FALSE".]
- Kernel-dimension = silent-channel-count coincidence at both closure
  vertices (1 = 1 at v89, 5 = 5 at v103); the lemma shows the v103 instance
  is a fixed-rho fact. [Same artifacts as above; see Problem 2b for the
  open genericity question.]
- Coherence obstruction at the incidence optimum exists (the witness-level
  content of "positive realization defect"), with in-repo certificates
  (scripts/reproduce_next_claims.py). DEFINITION (adopted): for a vertex
  spectrum lambda, s_I(lambda) = min |S| such that the diagonal incidence
  system A_S p = lambda, p >= 0, sum p = 1 is feasible; s_Q(lambda) = min
  support of a pure state with spec(rho) = lambda; realization defect =
  s_Q - s_I. CERTIFIED: s_I(v89) = 7 and s_I(v103) = 6, exact (MILP optimum
  + exact rational witness weights); AND both minimal incidence witnesses
  are PROVED quantum-infeasible in exact arithmetic: v89's by a
  projector-algebra contradiction ((P^2)_04 = P_08 P_84 with both factors
  forced nonzero by singleton channels while P_04 = 0), v103's by
  exhausting the eigenvalue-assignment branches (scalar-block kill;
  multiplicity kill; and the surviving rank-one branch's necessary
  magnitude system has zero nonnegative real solutions, exact solve). The
  defect hierarchy at the closures now reads
      v89 :  s_M = 4 < s_I = 7 <= s_Q <= 10
      v103:  s_M = 4 < s_I = 6 <= s_Q <= 12
  with s_M the (weak) majorization baseline and only s_M <= s_Q guaranteed
  a priori. [docs/unified_fiber_dimension.md section 3.]

🔬 STILL OPEN, do not overclaim: the GLOBAL statement s_Q > s_I (defect > 0)
requires killing every support of size s_I..s_Q-1, not just the minimal
witnesses. Ladder status for v89 at size 7: 41 supports enumerated
(incomplete), 34 killed exactly by the projector-pattern argument, 6
distinct survivors all floor numerically (1.5e-2..6e-2 over 24
multi-starts), exactifiable by the same branch analysis that killed the
v103 witness; sizes 8..9 untouched. Completion path: enumerate
incidence-feasible supports per size (MILP no-good cuts) and automate the
projector-pattern + branch-analysis kills (SOS pipeline,
docs/prior_art_roadmap.md section 7). CAVEAT on the baseline: Schur-Horn
permits diag(rho) majorized by lambda, so the airtight
necessary-condition baseline is s_M(lambda) = min |S| with A_S q majorized
by lambda, s_M <= s_I; the defect claim against s_M is strictly stronger
and untested.

### ❌ Disproved

- Holonomy radicands are square roots of rational MONOMIALS in the weights.
  This was the stronger sub-clause of the Holonomy Rationality sketch; it
  fails on 15 of the 82 loops, whose radicands (42, 235, 517, 3266, 4899,
  5538) lie outside the group generated by the weights and the denominator in
  Q*/(Q*)^2. Mechanism: at a channel class of size 3 or more the elimination
  adjoins sqrt(alpha^2 + beta^2 - delta^2), a POLYNOMIAL and not a monomial in
  the weight data, and later loops see the sines of earlier ones. The
  mechanism predicted exactly where the failures must live (only at states
  with a class of size >= 3) and that prediction is confirmed 63 of 63. The
  conjecture's actual conclusion, that the field is multiquadratic, is
  untouched. [docs/holonomy_rationality.md;
  results/data/holonomy_fields.json, radicand_weight_group]
- Ratio-1 states have unimodular incidence solves (P-DEN-2, pre-registered).
  FAIL, 18 exceptions. Unimodularity of A_S is SUFFICIENT for the incidence
  solve to stay over the spectrum denominator but not NECESSARY: an invariant
  factor f > 1 obstructs only right-hand sides outside the image sublattice,
  and (lambda, 1) can lie inside it anyway. P-DEN-3 (the converse) also
  failed, as pre-registered, with its predicted shape: the single exception
  v103 is pinned. P-DEN-4, the version with the escape route named in advance,
  PASSED on 736 states with 0 exceptions.
  [docs/prereg_denominator_arithmetic.md;
  results/data/denominator_arithmetic.json]
- Fiber dimension = incidence-kernel dimension, AS A UNIVERSAL LAW. [P4,
  both closures.] SCOPE NOTE: the law was fit on magnitude-target states
  and fails in the cancellation regime; the scoped (magnitude-sector)
  version is the surviving candidate and is exactly the unobstructed sector
  of the T1 Jacobian program. Do not discard the scoped version with the
  universal one.
- Singleton channels are the only rigidity mechanism. [P4: both closures
  rigid via multi-term cancellation constraints, no 1-term class. v89's
  rigidity is fixed-spectrum; v103's is fixed-rho.]
- Universal <=2-channel interference law. [P2: v103 carries 5 exchange
  channels; independently re-verified.]
- State denominator = spectrum denominator (and its ratio-in-{1,2}
  refinement). [P3: ratios 5 (v89) and 5746 (v103), rational-weight states,
  no escape clause. LEDGER NOTE 2026-07: the regenerated ledger ships all
  12 DESIGN-REAL states at ratio 2 (the archive's mined "9 of 12" is stale;
  representative-dependent, as the P2/P3/P4 caveat predicts).]
- Complex representative implies intrinsically complex fiber. [Both
  closures real; archive "the near-theorem is FALSE".]
- Positive-dimensional fibers are literal torus orbits. [Archive
  retraction: rational curve fibered over the gauge torus, not a torus
  orbit.]
- Generic one-loop fibers are elliptic. [Archive retraction; census of
  kdim-1 families: 18 genus-1 vs 26 genus-0.]
- Dense numerical representatives imply no sparse exact representative.
  [The contraction attack lands on dense generic fiber points; reweighted-L1
  over the full fiber found support 10 and 14; archive "WHAT OVERTURNED THE
  NEAR-THEOREM". The v103 cascade sharpens this: even a certified sparse
  point can fail to be the sparsest in its fiber.]
- Two spectral classes require higher-order interference, under the ADOPTED
  reading "spectra with exactly two eigenvalue blocks force complex (or
  higher-degree algebraic) interference phases". DISPROVED decisively by
  ledger census (scripts/reproduce_next_claims.py): ALL 64 two-block
  vertices ship all-real certified states, including all four two-block
  interference vertices (v89, v103, v108, and one (5,10) vertex); trisected
  (cube-root) phases occur only at 3-5 blocks. Existence of a real state
  refutes necessity, so this direction is representative-safe.
- Spectral border support differs from exact support, under the ADOPTED
  reading "the minimal support attainable as a LIMIT of exact-spectrum
  states can be smaller than the minimal exactly-attained support".
  DISPROVED at the tested instance (v89): the certified support-10 state is
  the c_(036) -> 0 endpoint of a support-11 one-parameter exact-spectrum
  family (continuation verified to residual ~1e-17 over c_x in [0, 0.2];
  the other five first-order openers are second-order obstructed), and the
  endpoint IS attained; border = exact here. NOTE: both adopted readings
  are reconstructions; if the parallel instance's original definitions
  differ, re-tag against those. Formal definitions of support-rank and
  border-support-rank, with the tensor rank vs border rank analogy,
  are in docs/prior_art_roadmap.md section 5.

## Highest-value open problems

(Each problem cites its existing groundwork; new agents start there, not
from zero. docs/prior_art_roadmap.md maps each to its outside-field sibling
and names the tools to import.)

1. Jacobian formula for local fiber dimension. [= THEOREMS.md T1. Proof
   route on file: MGS local normal form in the Sjamaar-Lerman stratified
   setting; the obstruction quadric is computable per state. First-order
   form now on file with corpus verification
   (docs/unified_fiber_dimension.md); the cancellation regime is proved
   (silence-rank lemma); remaining: the magnitude-regime proof and the
   second-order/integrability structure.]
2. Explicit reduced-fiber invariant for v_B, and exactification of the v103
   fiber. [Groundwork: v_B anomaly resolution, holonomy/wall data,
   vb_holonomy.py, vb_fiber_ideal.py. New v103 lead: the (0,2,7) weight
   equals 5/34 (the small spectrum value itself) at the certified point and
   at every cascade point to err < 3e-17, a candidate exact invariant and a
   coordinate to quotient out before computing the fiber ideal; PSLQ over
   quadratic surd lattices, seeded from the support-12 endpoint, is the
   next tool (docs/unified_fiber_dimension.md section 4). Witness sets +
   numerical irreducible decomposition + monodromy Galois groups are the
   industrial method (prior_art_roadmap section 4).]
   2b. Silence-rank genericity. The lemma is PROVED
   (docs/silence_rank_lemma.md); what remains open: the second-order
   structure at v103 (which of the 6 tangent directions integrate; at least
   3 do), an exact constant-rank argument upgrading the numerical
   germ-smoothness, and whether within-block silence rank = kernel dim is
   forced for census cancellation states or accidental (sample of size 2).
   The rigidity-matroid formulation (coherence rigidity matroid,
   generic-rank theorem or counterexample search) is the imported frame
   (prior_art_roadmap section 2).
3. Interior inverse walls. [Groundwork: wall_solver.py, wall_test.py,
   per-family endpoint theorem (T4).]
4. Graph-theoretic coupled-phase realizability criterion. [Groundwork:
   one-hop class combinatorics in classify.py; loop census in the archive.]
5. Infinite positive-defect family. [The defect definition is now fixed and
   witnessed (s_I / s_Q / s_M above); still BLOCKED on the global defect
   inequality at one instance (the 🔬 ladder above) before hunting a
   family. Crystallographic sparse phase retrieval is the template
   (prior_art_roadmap section 1).]
6. Infinite deformable and cancellation-rigid families. [Candidate seeds:
   the 44 kdim-1 deforming families; the two closure states.]
7. Levi-theoretic lower bounds. [Groundwork: levi.py post-audit, the
   rank_deficiency == stab - 1 identity (799/799), the hard/soft plethysm
   criterion.]
8. Holonomy Rationality at channel classes of size 3 or more. [Theorem A
   settles size <= 2; Proposition B gives only a 2-tower, and a 2-tower can be
   dihedral (sqrt(1 + sqrt 2) has D4 closure). What is needed: an argument
   that the eliminated discriminants alpha^2 + beta^2 - delta^2 are always
   squares times elements already in the field, or a structural reason the
   towers cannot go dihedral. A single C4 or D4 holonomy anywhere settles it
   the other way; none occurs in 82 loops. Second open piece: is hypothesis
   (H3) forced, i.e. can a vertex support carry a loop with no one-hop
   structure? docs/holonomy_rationality.md.]
9. Complexity classification. [BLOCKED on choosing the complexity measure:
   support size, state-denominator, phase field degree, and clique budget
   are all in play and provably non-equivalent (den-ratio falsification).
   Define the measure first.]

## Dependency graph

Incidence geometry
    -> Weight polytope
    -> Coherence equations
    -> Jacobian
    -> Fiber dimension / rigidity
    -> Inverse walls
    -> Support complexity atlas

## Dead ends

Future agents should NOT assume:

- kernel dimension predicts fiber dimension (outside the magnitude-target
  scope);
- representative properties characterize an entire fiber (P2/P3/P4 and the
  projective-stabilizer demotion are all representative-dependent);
- fixed-rho rigidity implies fixed-spectrum rigidity (v103; name the fiber
  in every rigidity claim);
- the certified ledger support is the minimal support s_Q (strict at v103);
- toric-orbit descriptions;
- denominator preservation;
- two-channel universality;
- that a dense numerical fiber point says anything about sparse exact
  points (the contraction-attack artifact that nearly closed the census
  wrongly).

## Immediate tasks

Done and landed (2026-07): scripts/reproduce_next_claims.py (exact s_I
certificates, witness infeasibility proofs, two-block phase census, v89
support-11 family, v103 deformability); the silence-rank lemma proof
(docs/silence_rank_lemma.md); the corpus tangent census and v103 cascade
(docs/unified_fiber_dimension.md).

Also done and landed (2026-07, fiber session): the two fiber notions made
structural in code (src/gpc_census/fiber.py, tests/test_fiber_notions.py);
the census-wide sparsification cascade (src/gpc_census/cascade.py,
scripts/census_cascade.py, docs/cascade_census.md); exact certification of
the improved support bounds (scripts/exactify_cascade.py); the scored
out-of-sample bracket pre-registration (docs/prereg_bracket_3_11_3_12.md);
and the reconciliation artifact (scripts/reconcile_claims.py).

Also done and landed (2026-07, arithmetic session): the census-wide holonomy
field scan (scripts/holonomy_field_scan.py,
results/data/holonomy_fields.json, tests/test_holonomy_fields.py); the
Holonomy Rationality theorem note with Theorem A proved and the conjecture
scoped (docs/holonomy_rationality.md); the two-stage denominator derivation
with its pre-registered census test (scripts/denominator_arithmetic.py,
docs/prereg_denominator_arithmetic.md,
results/data/denominator_arithmetic.json); the trisection dissolution; and
the manuscript's Galois section updated from the stale 142-state tally to the
full 82-loop result.

Priority:

1. Formalize the theorem-candidate proofs (the two Hardt-dependent ones
   together).
2. Alpha-certify what exact recognition could not reach. MOSTLY OVERTAKEN:
   69 of the 71 cascade endpoints were certified outright by high-precision
   refinement plus integer relation detection plus an exact
   characteristic-polynomial check, which is stronger than an alpha
   certificate and needs no external tool. What remains for alphaCertified
   or for quadratic-surd recognition: the two endpoints with irrational
   weights ((5,10) v144 and v256) and the v89 (0,3,6) family points.
   [results/data/cascade_exact.json]
3. Witness set + numerical irreducible decomposition of the v103 fiber and
   the v_B surface; monodromy Galois groups feed the holonomy program
   (prior_art_roadmap action 2).
4. Prove the Jacobian dimension theorem (T1), unifying the two regimes
   (Problem 1); coherence rigidity matroid for Problem 2b's remainder.
5. Complete the v89 defect ladder (sizes 7..9) with the automated
   projector-pattern + branch-analysis + SOS pipeline; then the global
   defect theorem.
6. Produce the stabilizer-invariant for v_B and exactify the v103 fiber
   (Problem 2; the 5/34 invariant coordinate is the seed).
7. Search for interior inverse walls.
8. Develop practical realization algorithms from the corrected theory
   (reweighted-L1 full-fiber continuation -> PSLQ exactification is the
   validated pipeline; isospectral-flow explorers are the principled
   replacement for ad-hoc continuation; productionize
   scripts/repro_and_sparsify.py-class tooling out of scratch).

## Reconciliation ledger (2026-07 fiber session)

Each item is re-derived by scripts/reconcile_claims.py into
results/data/reconciliation.json, so a reader can check it rather than trust
this table.

| superseded reading | corrected reading |
|---|---|
| v103 is first-order rigid | rigid for the FIXED-RHO fiber, deformable (tangent 6) for the FIXED-SPECTRUM fiber; v89 rigid in both |
| v103 is an isolated support-14 point | true of the route-B search landscape and gauge slice, false of the fiber: support 10 is reached and certified exactly |
| the census support column is the minimal support | it is the certified support, an upper bound on s_Q, strict at 71 of 799 vertices |
| 9 of 12 DESIGN-REAL states at denominator ratio 2 | the shipped ledger has 12 of 12; the mined figure is representative dependent |
| the trisected states are a cube-root phase-complexity tier | a printing artifact: (pi + 3*atan(t))/3 is pi/3 + atan(t), and no cubic subextension exists in any of the 13, invariant or representative |
| holonomy radicands are square roots of weight monomials | true only for channel classes of size at most 2 (proved); false at 15 of 82 loops, all at classes of size 3 or more, exactly as the mechanism predicts |
| the manuscript's Galois tally covers 142 states / 62 holonomies | superseded by the full corpus: 63 states, 82 loops, all 2-elementary |

REPRESENTATIVE-SWAP CHECK (v0.12, the real v_B record). Shipping the real
theta=0 endpoint in place of the phase-carrying fiber point over
(4,9) index 65 leaves every fiber-program result unchanged. Re-run on the new
ledger: identical census totals (71 sparsifying states, 89 amplitudes, 33
multiplicity-drift stops, 2 cancellation-regime members), identical per-vertex
verdicts and tangent dimensions for all 799 records, and identical exact
certification (69 of 71, same weights throughout). v_B's own numbers are
stable too: incidence kernel 1, fixed-spectrum tangent 2, fixed-rho tangent 1
for BOTH representatives, and both cascade to the same support-7 point,
dropping determinant (0,1,2,4) and certifying at the same squared weights
(9/23, 4/23, 3/23, 12/115, 8/115, 1/23, 2/23), state denominator 115. The two
endpoints differ only by a gauge phase, so the real representative simply
lands on the real form directly. The one visible change is elsewhere: v_B's
phase tier in the block-count census moves from COMPLEX to REAL, so the
3-block row reads 152 REAL / 13 COMPLEX rather than 151 / 14
(scripts/reproduce_next_claims.py claim 2).

## Guiding principle

The census should now be viewed as **experimental evidence**.

The primary research product is the mathematical theory of inverse fibers,
lifting complexity, and coherence obstructions. The inverse-fiber program
imports its questions from phase retrieval (ambiguity geometry), its
rigidity from framework theory, its algorithms from structured inverse
eigenvalue problems and numerical algebraic geometry, and its pathology
expectations from tensor rank, applied to a map on which none of those
fields has yet worked (docs/prior_art_roadmap.md).
