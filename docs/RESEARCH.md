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

Feasible design verdicts carry exact witnesses. The 11 negative verdicts in
Proposition 1 and the v_B verdict in Theorem 2 now carry exact global
disjunctive Farkas/branch certificates: 192 primitive integer vectors, 5141
symmetry-orbit hitting clauses, and 12 branch DAGs verified without a solver.
The other 144 INTERFERENCE verdicts still rest on the pinned two-solver census
and structural transport checks, so their verdict-level status remains
computer-assisted. The one-hop exclusion is essential; a degree-sums-only
model is strictly weaker and misclassifies (this bug was made and caught once;
see tests pinning v_A and v_B).

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
  certified state) is an UPPER bound on s_Q^NO for the 645 records that pass
  the diagonality gate, and on s_Q^free only for the other 154; see the
  support taxonomy below before quoting it, and use
  results/data/states_natural_orbital.jsonl for an s_Q^NO bound at an
  interference vertex. It is NOT known strict at v103.

## The support taxonomy (adopted 2026-07; do not collapse these)

Four different minima, with different lower bounds. Conflating them produced a
wrong reinstatement of a retracted bound, so every support claim in this
project must name which one it is about.

- s_M(lambda): least |S| with A_S q majorized by lambda. The airtight
  a-priori baseline.
- s_I(lambda): least |S| with the diagonal incidence system A_S p = lambda,
  p >= 0, sum p = 1 feasible. A DIAGONAL notion. s_M <= s_I.
- s_Q^NO(lambda): least support of a state whose 1-RDM is EXACTLY
  diag(lambda). s_I <= s_Q^NO. CORRECTION 2026-07: this is NOT what the
  library column bounds for every record. Only the 645 records that pass the
  diagonality gate (all 643 design states plus v89 and v103) have an exactly
  diagonal 1-RDM; the other 154 INTERFERENCE records attain the spectrum in
  another basis, so their support column bounds s_Q^free instead. Check
  results/data/attainer_gate_audit.json before quoting a support against
  s_Q^NO or s_I. [scripts/attainer_gate_audit.py; docs/rigidity_rationality.md]
- s_Q^free(lambda): least support of any state with spec(rho) = lambda.
  s_M <= s_Q^free <= s_Q^NO. **s_I does NOT bound s_Q^free**, so the two must
  never appear in one inequality.

Per-STATE (not per-vertex) there is a fifth: the ORBIT-MINIMAL support of psi,
the least support over the full U(d) orbit of that one state. It is a
basis-free complexity measure of the state, bounded below by s_M, and it upper
bounds s_Q^free for the vertex psi attains.

## The validation law (non-negotiable)

Every real error in this program was caught by structural invariants, never
by inspection: a sign-parse bug (caught by particle-hole self-duality of
(5,10)), a wrong classify port (caught by a Farkas feasibility
contradiction), a wrong solver gradient (caught by finite differences), and
the levi.py conjugation and fiber-conflation bugs (caught by the invariance
gate; RESEARCH_ARCHIVE "BUGS FIXED IN src/gpc_census/levi.py").
Consequences:

- No result ships without passing gpc_census.validate and the test suite.
- THE ATTAINER GATE (added 2026-07, enforced by
  gpc_census.validate.check_vertex_attainer). A state offered as attaining a
  vertex must have (1) a diagonal 1-RDM, (2) diagonal exactly lambda, (3) unit
  norm, and, if it is claimed to be a NEW attainer, (4) a 2-RDM spectrum
  DIFFERENT from the library representative's. Conditions 1 and 2 are what make
  its support comparable with s_I at all; a state passing only
  spec(rho) = lambda attains the vertex in some other basis and bounds
  s_Q^free instead. Condition 4 is basis-free and cheap: the 2-RDM spectrum is
  invariant under every one-body unitary, so a match means the state is the
  library one in rotated orbitals, which is a gauge statement and must be
  reported as one. This gate would have caught the mislabelled v103 support-10
  bound immediately.
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
2. The (3,11) bracket (docs/bracket_3_11.json) is CLOSED: all 50 outer
   candidates are settled WITHOUT Stage 1, in exact arithmetic, at
   19 TRUE + 31 REFUTED + 0 OPEN.
   scripts/settle_bracket_3_11.py regenerates docs/bracket_3_11_settlement.json
   at 19/31/0 using five exact methods. The first four settle the same 46 they
   always did: 19 TRUE (17 by state transport from the rank-9/10 census, the
   uniform (3/11)^11 by an explicit Z11 difference design, and a frozen-core
   paired state); 27 REFUTED (frozen-core pinning to N=2 plus the
   even-degeneracy pairing theorem, and zero-restriction against known
   lower-rank GPCs). The fifth, LEVEL5-RESSAYRE, is applied last and refutes
   exactly the four that were OPEN (cands 23, 26, 34, 44), each record naming
   the certified row it violates and the excess; see item 3.
   CAVEAT unchanged: the refutations do NOT finish (3,11) as a polytope;
   cutting them creates new true vertices absent from the outer list, so
   Stage 1 (or a tighter bracket) is still required. The settled points are a
   mandatory acceptance oracle for any (3,11) generator.
3. The level-5 series inequalities are PROVED, which closed the 4 OPENs
   (docs/level5_ressayre_3_11.md, scripts/level5_ressayre_3_11.py,
   tests/test_level5_ressayre_3_11.py). The four rows that kept 23, 26, 34
   and 44 open, AK-RMK-4.2.1 and KLY09-EXT-1/2/3, all of the form
   lambda_A + lambda_11 <= 2 for a second-level quadruple A, carry exact
   Ressayre certificates from the repository's own generation package. That
   is branch (a) of the old fork, and it lands: each candidate violates at
   least one certified row, so each is outside the polytope. Branch (b),
   attaining idx 44 to falsify KLY09-EXT-1, is dead; its design route was
   already closed by exact counting in
   docs/rank11_open_candidate_designs.md.
   scripts/partial_families.py tiers these four RESSAYRE at rank 11 ONLY. The
   rank-12 rows share index sets but are wedge^3 C^12, a different
   representation, and stay CLAIMED: restriction carries a valid (3,12) row
   down to (3,11), never the reverse.
   The certificate proves VALIDITY of each row, not facetness or
   irredundancy, and the screen records both as not established. Validity is
   what refutes a candidate, so that is enough here.
   The four certificates are replayed by
   scripts/verify_level5_ressayre_3_11_standalone.py, which imports nothing
   from gpc_census, derives the fermionic sign by inversion counting rather
   than position bookkeeping, and recomputes the determinant by exact
   Fraction elimination rather than Bareiss. 64 checks, and the test suite
   asserts it fails on a tampered determinant or hyperplane. This closes the
   gap the rank-7 checkpoint records, that the Ressayre determinants sit
   outside its independent verifier's trust boundary.
   The screen is calibrated at ranks 9 and 10, the largest published fixed-N=3
   systems, by a COMPLETE test: 2,898 subset rows screened, 50 certified, and
   every certified row maximized by LP over the published polytope with
   0 invalid. Negative control: 1,857 demonstrably false rows, 0 certified.
   Caveat recorded in the artifact and the doc: the trace count did nearly all
   the rejecting, only 6 false rows reached the tangent determinant, and a
   targeted hunt over perturbed published facets found no more, so the
   determinant step is the least stressed part of the calibration.
   Superseded by this, retained as history: the normalized contraction audit
   in results/data/rank11_contraction.json, which left positive numerical
   floors on all four and identified cand 44 at 11/4320 as the sharpest
   target. It pointed the right way and was never a certificate.
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

- holonomy_two_elementary.md: the PROOF ATTACK on 2-elementarity, and the
  current state of the classes-3-and-4 gap. ✅ REDUCTION THEOREM: the
  multiquadratic conjecture follows from two conditions on one quantity per
  loop, (R1) deg_Q(cos^2 Phi) <= 2 and (R2) N(cos^2 Phi) a rational square at
  degree 4, the second being exactly the classical V4 criterion for an even
  quartic (D4 and C4 are the alternatives it excludes). ✅ THEOREM A',
  strengthening the shipped Theorem A: at one-hop class size <= 2 the cosine
  is n/sqrt(M) with n, M rational, so cos^2 is RATIONAL and the minimal
  polynomial is EVEN, which is the part that generalises. 🟦 (R1) holds 82/82
  and (R2) 13/13, and the reduction DERIVES every certified Galois group
  (C1 37, C2 32, V4 13) with no group computation. 🟦 Theorem A' has 0
  violations, and every cos^2-irrational loop sits in a state with a class of
  size 3. ❌ "The minimal polynomial is always even" is FALSE: two (5,10)
  loops have deg cos = deg cos^2 = 2; 2-elementarity survives trivially there,
  but any proof of evenness must explain them. 🔬 (R1) at class sizes 3 and 4
  is the remaining open step, and Proposition B does not supply it (squaring
  leaves a cross term). 💡 The conjugation-involution/Galois-commutant
  mechanism is stated precisely: time reversal and the fiber sign flips are
  commuting INVOLUTIONS, so the group they generate is elementary abelian and
  has no element of order 4, which is what D4 and C4 both need; the
  obstruction is that Galois permutes COMPLEX solutions while those symmetries
  act on the REAL ones, the same gap Rigidity-Rationality's descent step hits.
  [scripts/holonomy_two_elementary.py;
  results/data/holonomy_two_elementary.json;
  tests/test_holonomy_two_elementary.py]
- conditional_2rdm_program.md: steps 4-7 built on the order-three theorem.
  🟦 HIDDEN PAIR CURRENT: states with exactly the same complete 1-RDM carry
  different two-body pair currents, over the exact flux-INDEPENDENT range
  [-sqrt(15)/10, sqrt(15)/10]; the energy/current joint body is convex and is
  recovered by sweeping the support function. 🟦 The unit-edge discriminant
  gives interior branch-exchange fluxes at Phi ~ 0.0593 and ~ 3.0823, the
  genuine non-analyticities of the constrained-search functional. 🟦 The
  Slater-Condon carrier projection reproduces the full CI submatrix with max
  error 0.0, and its fermionic signs are recomputed from occupation-string
  algebra rather than assumed. 🟦 QUASIPINNING: eps_H = 2||H||(eps+eps^2) +
  ||w-w0||_1 ||y||_inf, where the second Lipschitz constant is the dual
  certificate itself; no quasipinning constant is assumed. Pinning plus the
  1-RDM narrows the achievable range by 1-2 orders of magnitude. 🟦 ATLAS: all
  799 supports tiered, 75 exactly spectrahedral (m<=3) and 724 outer-bounded,
  and only 72 have a complete coherence graph, so MOST census coherences are
  invisible to any two-body observable. ❌ Two model choices failed
  instructively and are recorded: spinful Hubbard (degenerate ground state, so
  the 1-RDM is undefined) and the symmetric ring (degenerate occupations, so
  the carrier is ambiguous). 🔬 Only 5 of 12 scan points satisfy the
  nondegeneracy preconditions, and none of them is a genuinely complex-flux
  case, so the complex engine still has NO physical benchmark.
  [gpc_census.pinned_physics; scripts/conditional_2rdm_program.py;
  results/data/conditional_2rdm_program.json;
  tests/test_conditional_2rdm_program.py]
- conditional_2rdm_body.md: the REFRAMING of the observable-range work, and
  the one to read first. The endpoints of a two-body observable's range are the
  SUPPORT FUNCTION of K(gamma) = conv{Gamma^(2)_Psi : Psi -> gamma}, so the
  solver is an exact support oracle for the conditional 2-RDM body, and the
  headline is an exact conditional N-representability statement rather than two
  solved triangle Hamiltonians. 🟦 ORDER-THREE THEOREM: the variable part of
  K(gamma) is an affine image of the scaled complex elliptope E_w, whose
  extreme points obey r^2 <= m (Li and Tam), so at m = 3 every extreme point is
  a phase matrix and the support function is an EXACT SDP, not a relaxation.
  🟦 The sextic of complex_fixed_1rdm_observable_ranges is exactly the KKT
  elimination ideal of that SDP (verified, up to primitive content), which also
  explains structurally why only |h_ij|^2 and Re(chi) enter. 🟦 Exact dual
  certificates for four carriers: optimality is witnessed POINTWISE by PSD-ness
  of Diag(y)-H from principal-minor signs, so no exhaustive-root argument is
  needed, and the zero/pi-flux cases the sextic engine excludes as degenerate
  are covered uniformly. ❌ The theorem does NOT extend past order three: an
  explicit rank-two extreme point at m = 4 gives a relaxation gap of 6.9e-2, so
  m >= 4 supports ship as outer bounds until proved otherwise. 🔬 Still open,
  in attack order: energy/pair-current joint body, Slater-Condon projection,
  quasipinning enlargement, census-wide atlas.
  [gpc_census.elliptope; scripts/conditional_2rdm_body.py;
  results/data/conditional_2rdm_body.json; tests/test_conditional_2rdm_body.py]
- complex_fixed_1rdm_observable_ranges.md: exact constrained observable ranges
  on the pinned Borland-Dennis fiber for rational Gaussian Hermitian carriers.
  🟦 The gauge-invariant cycle flux chi = h01*h12*conj(h02) is the new
  invariant: edge phases are gauge-dependent but chi is not, so two carriers
  with equal diagonals and edge magnitudes can have different ranges. Phase
  elimination gives a quartic branch curve whose critical values are a
  primitive sextic in Z[z]; endpoints are its extreme real roots by exact
  Sturm isolation, with a subresultant linear lift proving every real root has
  a real (not spuriously complex) critical preimage. ✅ Verified independently
  here: brute-force optimisation of the original two-phase objective
  reproduces both ranges to 1e-15, the fiber's 1-RDM is the stated target and
  is diagonal (the support is one-hop free), A*D-C^2 = delta^2 holds, and the
  three fermionic channel signs of the at-most-two-body realisation check out
  by hand. 🟦 The unit-edge flux family contains BOTH earlier real-symmetric
  regressions as exact u = +-1 boundary points, agreeing with
  fixed_1rdm_observable_ranges.json, which a different engine produced; pinned
  by test_flux_family_reproduces_the_real_symmetric_regressions. Scope limits
  worth respecting: rational Gaussian entries only, generic (nonzero flux,
  squarefree sextic, coprime lift) only, and one fixed three-determinant
  fiber. Complex flux is NOT generic for ordinary real-orbital molecular
  Hamiltonians; it wants magnetic fields, twisted boundaries, complex Bloch
  orbitals, or spin-orbit structure.
  [scripts/complex_fixed_1rdm_observable_ranges.py;
  scripts/verify_complex_fixed_1rdm_observable_ranges_standalone.py;
  results/data/complex_fixed_1rdm_observable_ranges.json;
  tests/test_complex_fixed_1rdm_observable_ranges.py]
- lambda3_c7_polytopes.md: the ORBIT-CLOSURE layer (entanglement polytopes),
  one polytope per SLOCC class rather than one per (N,d). 🟦 CERTIFIED EXACT
  for 10 of 13 class polytopes: all four of Lambda^3 C^6 (reproducing Walter,
  Doran, Gross, Christandl) and six of the nine of Lambda^3 C^7, including the
  decisive gate that the OPEN ORBIT reproduces the ambient (3,7) polytope
  exactly. Classes VI, VII, VIII ship as certified inner/outer brackets. Both
  bounds are exact: the inner from torus orbits of one-hop-free supports (a
  COMPLETE enumeration, because one-hop-free supports have vanishing
  coefficient moduli) plus certified degenerations; the outer from a
  Schur-Horn/Newton dominant-weight family that generalises Ky Fan. 🔬 The
  outer family is SATURATED (sweeping the weight bound 3 to 7 changes
  nothing), so the last three classes need Ressayre inequalities FOR THE ORBIT
  CLOSURE, or covariant nonvanishing. Do NOT read stage1_ressayre_3_7 as
  closing them: that manifest derives the AMBIENT cone, which is the open
  orbit's polytope (class IX), already EXACT here. The agreement is a
  cross-check on the ambient input and is pinned by a test; the bracketed
  classes are untouched by it. 🔬 Whether class V (orbit dim 26) lies in the closure of class
  VI (orbit dim 28) is open and recorded, not guessed.
  [gpc_census.orbit; scripts/lambda3_c7_polytopes.py;
  results/data/lambda3_c7_polytopes.json; tests/test_lambda3_c7_polytopes.py]
- HANDOFF.md: the note to the next agent from the 2026-07 gauge session, which
  supersedes two instructions of the handoff that preceded it (the v103
  endpoint is not diagonal, and 8 of the 71 cascade endpoints are new states).
  Read it before acting on any support claim.
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
  <= 12 is superseded by s_Q^free(v103) <= 10, certified exactly; see
  cascade_census.md.
- gauge_minimization.md: minimization of every certified state's support over
  its natural-orbital gauge family, the exact profile-class/exterior-contraction
  lower bound (670 of 799 certified gauge-minimal), the null result on the other
  129, the power calibration of the search, and the scored pre-registration.
- prereg_gauge_minimization.md: that pre-registration, committed before the run.
- cascade_census.md: the census-wide sparsification cascade, the exact
  certification of all 71 improved support bounds, the v103
  support-10 state, the cancellation-regime census, and the
  representative-swap check against the v0.12 real v_B record.
- prereg_bracket_3_11_3_12.md: the scored out-of-sample pre-registration on
  the (3,11) and (3,12) bracket holdout (4 PASS, 2 FAIL), including the
  power caveat stated before scoring, AND the 2026-07 PROVENANCE CORRECTION:
  34 of the 38 rows are transports or face embeddings of 17 rank-10 donors,
  so the effective independent sample is 2 states, not 38 rows. Read the
  corrected readings, not the original verdicts.
- prereg_template.md: the standing rules for every pre-registration in this
  project. R2 is the one added by that correction: a measurement derived by
  transport, padding, or face embedding is NOT an independent out-of-sample
  observation, and every scorecard must report independent sample size
  separately from row count.
- two_block_family.md: the four two-block interference vertices in full, and
  the census-wide result that a CROSS-BLOCK one-hop channel is NECESSARY for
  the two fibers to disagree (799/799, zero exceptions). Also kills the
  proposed cross-block counting rules and shows two-block is not the
  operative invariant.
- fiber_symmetries.md: the pair-incidence separation proposition, exact
  nonzero-minor certificates for all 799 shipped supports, the time-reversal
  involution on the exact v_B slice, and the cubic 2-RDM moment that recovers
  its family parameter modulo every one-body unitary.
- prereg_two_block_inversion.md: that study's pre-registration, committed
  before the cross-block counts were computed.
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
- prereg_rank11_candidate44.md: the frozen theorem/falsification fork for the
  normalized candidate-44 floor, exact projected spectrum, and required proof
  objects. The numerical audit is results/data/rank11_contraction.json.
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
  with both physical and one shipped, so the distinguished point has orbit
  size 2. Its own
  loop holonomy is a sign (minimal polynomial x + 1), consistent with the
  record being real up to gauge. It is excluded from denominator predictions,
  having no state denominator. [docs/holonomy_rationality.md;
  tests/test_holonomy_fields.py]

- Positive-dimensional REDUCED vertex fibers exist (v_B). [THEOREMS.md
  "v_B ANOMALY RESOLVED" (reduced spectrum-fiber is a surface);
  RESEARCH_ARCHIVE "THE FIBER IS THE 1-RDM -> 2-RDM INFORMATION GAP";
  scripts/fiber_sampler.py, scripts/vb_fiber_ideal.py.]
- Pair occupations separate every weight-sector direction on every shipped
  support. For a support S, let B_S be its pair-by-determinant incidence
  matrix. The identity sum_{j != i} (B_S w)_{ij} = (N-1)(A_S w)_i gives
  ker B_S contained in ker A_S, so full column rank of B_S makes pair
  occupations injective on all squared amplitudes and, in particular, on the
  incidence kernel. Exact nonzero maximal minors certify full column rank for
  799 of 799 supports, spanning 403 transport classes. The positive-kernel
  population is 63 records in 39 transport classes, with zero failures. This
  is a support-restricted fixed-frame statement, not a universal theorem:
  signed combinatorial 2-trades are abstract counterexamples, and pair
  occupations are not invariants of the full degeneracy-stabilizer quotient.
  [docs/fiber_symmetries.md; scripts/fiber_symmetries.py;
  results/data/fiber_symmetries.json;
  scripts/verify_fiber_symmetries_standalone.py]
- The first exact full-one-body-unitary invariant on the v_B family is now
  explicit. On its known one-parameter fixed-spectrum slice, with parameter t,
  tr(Gamma_2^2) = 1438/529 is constant but
  tr(Gamma_2^3) = 6(2899-8t)/12167. The cubic moment therefore recovers t and
  proves that distinct t-values survive as distinct points of the reduced
  fiber. It also separates the two Q(sqrt(35)) Galois wall states by
  1728*sqrt(35)/1034195, so Galois conjugacy is not one-body-unitary
  equivalence. Complex conjugation acts by theta -> -theta and fixes the two
  real walls modulo orbital-phase gauge. Whether the two interior conjugate
  sheets at fixed t are equivalent under the full U(1) x U(4) x U(4)
  stabilizer remains OPEN. At the exact interior point t=0, the full
  support-restricted fixed-spectrum constraint differential has rank 7 in 16
  real amplitude coordinates, its orbital-phase gauge has rank 7, and the
  support tangent modulo that phase gauge has dimension 2. Appending
  d tr(Gamma_2^3) and
  d tr(Gamma_2^4) raises the exact rank from 7 to 9; a nonzero 9 by 9 minor is
  299925504*sqrt(8211)/1035662780341225. Thus the cubic and quartic moments
  form a basis of the phase-gauge-quotient cotangent at t=0. Since both are
  full U(9) invariants, no additional one-body-orbit direction lies inside
  this support tangent. This is FIRST ORDER only:
  it does not prove integrability, germ smoothness, or global injectivity on
  the second surface direction.
  [docs/fiber_symmetries.md; results/data/fiber_symmetries.json;
  tests/test_fiber_symmetries.py]
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
  s_Q^free(v103) <= 10, certified. Its state denominator 46852 = 2^2*13*17*53
  carries a prime found in neither the spectrum denominator nor the
  library state, which is why no grid search over multiples of the natural
  denominator could reach it. Upper bound only, not a minimality claim.
  [docs/cascade_census.md; scripts/v103_endpoint_precision.py;
  scripts/v103_support10_certificate.py; tests/test_v103_support10.py;
  results/data/v103_support10_certificate.json]
- The natural-orbital GAUGE was never explored before 2026-07, and exploring
  it does not move the library column. Every vertex has a degenerate spectrum,
  so every certified state carries a gauge freedom U(m_1) x ... x U(m_k) whose
  every point has 1-RDM identically diag(lambda) and is an equally legitimate
  natural-orbital representative. Minimizing support over that family:
  0 of 799 states sparsify, and 670 of 799 are EXACTLY gauge-minimal, certified
  by a gauge-invariant lower bound. A block-diagonal unitary preserves each
  determinant's block profile. Within a profile class, exterior/Koszul
  contraction maps give support at least
  ceil(rank(C_q)/product_i binom(k_i,q_i)). The symbolic rank certificates
  strictly strengthen 82 records and newly close 41, and all 799 shipped bounds
  are exact (floating ranks only select the contraction map). The other 129 are
  open in both directions. The search's power is pinned by a recovery test (a provably
  rank-1 class must collapse to one determinant; it does, at residual 4e-16),
  and two earlier search methods that failed that test were removed rather than
  shipped.
  [docs/gauge_minimization.md; docs/prereg_gauge_minimization.md;
  src/gpc_census/gauge.py; results/data/gauge_min.jsonl]
- NONE of the 71 cascade endpoints is in a natural-orbital basis: 0 of 71 have
  a diagonal 1-RDM. The fixed-spectrum continuation preserves spec(rho), which
  is all it promises, and not diagonality. So every improved cascade bound is a
  bound on s_Q^free and NEVER on s_Q^NO, and none of them may be quoted beside
  s_I. At v103 the exact certificate is arithmetically sound but was labelled
  with the wrong quantity: its state carries rho_39 = -5/68 exactly, so it
  certifies s_Q^free(v103) <= 10 while s_Q^NO(v103) stays bounded by the
  library support 14. [results/data/orbit_check.json;
  scripts/v103_support10_certificate.py; tests/test_v103_support10.py]
- The cascade is mostly an ORBITAL ROTATION ENGINE, but not entirely. By the
  2-RDM spectrum (invariant under every one-body unitary), 63 of the 71
  endpoints ARE the certified library state in different orbitals (agreement
  ~1e-16), so their supports describe the STATE, not the vertex: the
  natural-orbital basis is not the sparsest basis for those states. The other
  8 differ by 3.4e-3 to 4.7e-2, which no one-body rotation can produce, so
  those are genuinely NEW extremal states over their vertices ((4,10) v93,
  (5,10) v65, v256, v144, v140, v230, v263, and (4,9) v65 = v_B). The split is
  bimodal with nothing between 1e-12 and 1e-3. At v103 the connecting unitary
  is exhibited explicitly (residual 7.5e-14) and its off-block entries reach
  0.196, so it is a U(10) rotation and not a gauge one.
  [results/data/orbit_check.json; scripts/orbit_check.py]
- Census-wide sparsification: the support column is an UPPER BOUND on
  s_Q^free, and it is STRICT at 71 of the 799 certified states,
  all INTERFERENCE. The cascade removes 89 amplitudes in total: 55 states
  drop one, 15 drop two, v103 drops four. All 71 endpoints are exactly
  recognizable, but on the SPECTRUM locus only: under the strengthened
  criterion (rho diagonal AND diag = lambda) 0 of 71 certify as attainers, 63
  are the library state in rotated orbitals and 8 are new states bounding
  s_Q^free. The final two exactifications, (5,10) v144 and v256, share
  Q(sqrt(1065)); their exact characteristic polynomials and nonzero
  off-diagonal 1-RDM entries are verified symbolically. [docs/cascade_census.md;
  src/gpc_census/cascade.py;
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
  itself (fixed-rho 1, fixed-spectrum 2) and at 73 of the 799 census states
  (41 transport classes), all INTERFERENCE. Any text presenting the two-fiber
  distinction as a cancellation-regime phenomenon is too narrow.
  PROVENANCE NOTE (2026-07): the (3,11) v43 instance once cited here as
  out-of-sample evidence is NOT independent; it is (3,10) v108 padded with an
  empty orbital, identical on every structural field. The claim stands on the
  census, not on the bracket.
  [results/data/cascade.jsonl; results/data/two_block_family.json;
  docs/prereg_bracket_3_11_3_12.md PROVENANCE CORRECTION]
- A CROSS-BLOCK one-hop channel is NECESSARY for the two fibers to disagree:
  gap > 0 implies X > 0, on all 799 states across 403 transport classes with
  ZERO exceptions. This promotes the silence-rank lemma's cross-block clause
  from a dimension formula verified at two states to a census-wide necessary
  condition. The CONVERSE IS FALSE: 82 states have a cross-block channel and
  no gap, every one of them with incidence kernel 0 and rigid in both fibers,
  so the condition exists but has no weight freedom to cut. The counting
  rules (gap == X, gap == 2X) are both FALSIFIED.
  [docs/two_block_family.md; docs/prereg_two_block_inversion.md;
  results/data/two_block_family.json; tests/test_two_block_family.py]
- The v103 inversion is not unique. 14 states are fixed-rho rigid with
  positive fixed-spectrum tangent, in 7 transport classes, patterns
  (0,1) x 11, (0,2) x 2, (0,6) x 1. v103 is the unique inversion of
  MAGNITUDE 6, not the unique inversion. Its six directions decompose as
  3 weight + 3 phase, reproducing the lemma's formula
  (5-2) + (14-9-2) term by term, with the sectors decoupling exactly as the
  proof's real-state hypothesis requires.
  [results/data/two_block_family.json, tangent_sector_split]
- Two blocks is NOT the operative invariant for fiber behaviour. The gap
  population spans block counts 2 through 6 and peaks at 4; two-block
  vertices are UNDER-represented (3 of 64 carry a cross-block channel against
  155 of 799 census-wide), and all 60 two-block DESIGN vertices have zero
  channels because design supports are one-hop free by definition. The four
  two-block interference vertices span the known behaviours because they span
  (X, dim K), not because two blocks causes anything.
  [docs/two_block_family.md]
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
  system A_S p = lambda, p >= 0, sum p = 1 is feasible; s_Q^NO(lambda) = min
  support of a pure state whose 1-RDM is EXACTLY diag(lambda); realization
  defect = s_Q^NO - s_I. (The defect is against the NATURAL-ORBITAL minimum.
  Taking it against s_Q^free instead, the minimum over all states with
  spec(rho) = lambda, is not meaningful: s_I does not bound s_Q^free. The
  original wording defined s_Q by the spectrum and then subtracted s_I from
  it, which is the conflation the support taxonomy above exists to stop.) CERTIFIED: s_I(v89) = 7 and s_I(v103) = 6, exact (MILP optimum
  + exact rational witness weights); AND both minimal incidence witnesses
  are PROVED quantum-infeasible in exact arithmetic: v89's by a
  projector-algebra contradiction ((P^2)_04 = P_08 P_84 with both factors
  forced nonzero by singleton channels while P_04 = 0), v103's by
  exhausting the eigenvalue-assignment branches (scalar-block kill;
  multiplicity kill; and the surviving rank-one branch's necessary
  magnitude system has zero nonnegative real solutions, exact solve). The
  defect hierarchy at the closures now reads
      v89 :  s_M = 4 <= s_Q^free <= 10;  s_I = 7 <= s_Q^NO <= 10
      v103:  s_M = 4 <= s_Q^free <= 10;  s_I = 6 <= s_Q^NO <= 14
  with s_M the (weak) majorization baseline. Only s_M <= s_Q^free and
  s_I <= s_Q^NO hold a priori. v89's upper bound serves both columns because
  its certified state is diagonal and did not sparsify; v103's support-10
  state serves only the free-basis column (its 1-RDM is not diagonal), and no
  gauge rotation improves its natural-orbital column.
  [docs/unified_fiber_dimension.md section 3; docs/gauge_minimization.md.]

🔬 STILL OPEN, do not overclaim: the GLOBAL statement s_Q^NO > s_I (defect > 0)
requires killing every support of size s_I..s_Q^NO-1, not just the minimal
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

0. 💡 TWO-TORSION IN THE FIBER. Relocated here from the manuscript outlook
   before the freeze: it is an interpretation resting on a single vertex, so
   it belongs in the open-questions register and not in the paper. The two
   counts it rests on ARE certified and stay in the manuscript's results
   section as data (results/data/vb_galois_orbit.json). Verbatim as it stood:

   > One observation ties the arithmetic to the geometry, and we flag it as
   > speculative because it rests on a single vertex. At v_B the amplitude
   > signs that attain the vertex are not arbitrary: for each of the two
   > Galois embeddings, exactly 64 = 2^6 of the 128 sign patterns modulo
   > global sign attain, so the attaining set is a 2-torsion coset over a
   > single weight vector rather than a scattered subset [E]. That is the
   > Kummer structure of the holonomy fields, every field in the census being
   > multiquadratic, hence generated by square roots, appearing geometrically,
   > as 2-torsion in the fiber over a rational point. Whether the fixed-rho
   > fiber of a rigid vertex is in general a finite 2-torsion set over a
   > unique weight vector is open, and we do not assert it: the one other
   > vertex where multiple isolated fiber points have been exhibited, the
   > denominator-34 vertex, has points that share a single weight vector but
   > are separated by holonomies taking values across (-1,1) rather than by
   > signs, so it does not fit the picture as stated.

   Status: the v_B leg is certified. The v103 leg no longer refutes anything,
   because v103's fiber has exactly ONE point modulo gauge, holonomy
   (0, pi, pi, 0, 0), certified by exhausting the 5-torus of gauge-invariant
   holonomies (16807 grid starts, all converging, one point). The earlier
   report of dozens of points at v103 was an artifact of a non-integer loop
   basis and is withdrawn; see correction 11 in docs/closing_ledger.md. So the
   open question is not whether v103 fits the 2-torsion picture, it is whether
   the picture says anything at all beyond v_B, where it is a single certified
   data point. Attack it as a mechanism, not by mining two vertices.
   [results/data/vb_galois_orbit.json; results/data/v103_fiber_count.json] So whatever unifies the two vertices is not
   2-torsion in the holonomy. This is arc two's first question, and it should
   be attacked as a mechanism, not mined from the two data points.
   [results/data/vb_galois_orbit.json;
   results/data/uniqueness_census.json, v103_discreteness_probe]

1. Jacobian formula for local fiber dimension. [= THEOREMS.md T1. Proof
   route on file: MGS local normal form in the Sjamaar-Lerman stratified
   setting; the obstruction quadric is computable per state. First-order
   form now on file with corpus verification
   (docs/unified_fiber_dimension.md); the cancellation regime is proved
   (silence-rank lemma); remaining: the magnitude-regime proof and the
   second-order/integrability structure.]
2. Explicit reduced-fiber invariant for v_B, and exactification of the v103
   fiber. [PROGRESS: on the known fixed-spectrum v_B slice,
   tr(Gamma_2^3) = 6(2899-8t)/12167 exactly recovers the family parameter and
   separates the two Galois walls under the full one-body group. This gives
   one exact reduced coordinate. At t=0 its differential together with
   d tr(Gamma_2^4) exactly separates the complete two-dimensional support
   tangent modulo orbital-phase gauge. Global integration and separation of
   the second surface direction,
   and the full-stabilizer status of the conjugate phase sheets, remain open.
   Groundwork: docs/fiber_symmetries.md, the v_B anomaly resolution,
   holonomy/wall data, vb_holonomy.py, vb_fiber_ideal.py. New v103 lead: the
   (0,2,7) weight
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
   witnessed (s_I / s_Q^NO / s_M above); still BLOCKED on the global defect
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
- the certified ledger support is the minimal support (it bounds s_Q^NO from
  above; the natural-orbital gauge does not improve it anywhere in the census,
  and it is NOT known strict at v103, whose support-10 witness is a
  free-basis one);
- that a support bound may be compared with s_I without naming which minimum
  it is (s_I bounds s_Q^NO only, never s_Q^free);
- that a sparser state reached by continuation is a NEW state (63 of 71
  cascade endpoints are the library state rotated; check the 2-RDM spectrum);
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

Also done and landed (2026-07, gauge session): the natural-orbital gauge
minimization and its exact lower bound (src/gpc_census/gauge.py,
scripts/gauge_census.py, docs/gauge_minimization.md); the 2-RDM acceptance gate
and the census-wide orbit check of the cascade endpoints
(scripts/orbit_check.py); the support taxonomy above; and the relabelling of
every cascade bound from s_Q to s_Q^free.

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

0c. Extend the exact global no-design certificate atlas from the 12 numbered
   theorem/proposition cases to the other 144 INTERFERENCE verdicts. The
   disjunctive Farkas method is complete and fast on all 12 published targets;
   a highly symmetric rank-10 two-block case shows that orbit clauses must be
   stored and verified symbolically rather than expanded naively.
0. Continue the gauge lower-bound campaign on the 129 open rows. The
   exterior/Koszul contraction family removed the former Lambda^a fallback,
   improved 82 records, and certified 41 more; the remaining slack needs a
   different invariant.
0b. Exactly recognize the remaining 6 genuinely new extremal states found by
   the orbit check. The v144/v256 Hodge pair is now exact over Q(sqrt(1065)).
1. Formalize the theorem-candidate proofs (the two Hardt-dependent ones
   together).
2. Alpha-certify what exact recognition could not reach. OVERTAKEN for the
   cascade: all 71 endpoints are exactly recognized on the spectrum locus by
   high-precision refinement, degree-at-most-two integer relation detection,
   and an exact characteristic-polynomial check. The v89 (0,3,6) family points
   remain candidates for alpha theory.
   [results/data/cascade_exact.json]
3. Witness set + numerical irreducible decomposition of the v103 fiber and
   the v_B surface; monodromy Galois groups feed the holonomy program
   (prior_art_roadmap action 2).
4. Prove the Jacobian dimension theorem (T1), unifying the two regimes
   (Problem 1); coherence rigidity matroid for Problem 2b's remainder.
5. Complete the v89 defect ladder (sizes 7..9) with the automated
   projector-pattern + branch-analysis + SOS pipeline; then the global
   defect theorem.
6. Integrate the exact cubic-quartic tangent coordinates over the second
   direction of the v_B reduced surface, test global injectivity, decide the
   full-stabilizer status of the conjugate sheets, and exactify the v103 fiber
   (Problem 2; the 5/34 invariant coordinate is the v103 seed).
7. Search for interior inverse walls.
8. Develop practical realization algorithms from the corrected theory
   (reweighted-L1 full-fiber continuation -> PSLQ exactification is the
   validated pipeline; isospectral-flow explorers are the principled
   replacement for ad-hoc continuation; productionize
   scripts/repro_and_sparsify.py-class tooling out of scratch).

## Rigidity and rationality (2026-07)

Full write-up: docs/rigidity_rationality.md. This is where the weight-rationality
hypothesis that docs/holonomy_rationality.md and the denominator arithmetic both
assume gets its grounding.

- ✅ DESIGN CASE PROVED. On a one-hop-free support every off-diagonal 1-RDM
  entry vanishes termwise, so rho = diag(lambda) is exactly the rational
  linear system A_S p = lambda with sum p = 1; when its solution is unique the
  weights are rational by linear algebra alone. All 643 design states have
  incidence kernel 0, so the hypothesis holds throughout and the theorem
  covers every one of them. No rigidity input and no Galois argument needed.
  The 156 interference states are not reachable this way, because phases make
  the system nonlinear in p. [docs/rigidity_rationality.md]
- 🟦 FIXED-RHO RIGID IMPLIES RATIONAL WEIGHTS, 740 of 740, zero
  counterexamples. The single irrational record (v_B, weights in Q(sqrt 35))
  is deformable. Deformability does not force irrationality: 58 of the 59
  deformable records are rational anyway.
  POWER WARNING, not optional when quoting this: all 645 natural-orbital
  representatives are rigid and all 59 deformable records are among the 154
  that fail the diagonality gate, so the discriminating power lives entirely
  in the non-natural-orbital population and there is no rigid-versus-
  deformable contrast among gate-passing records. The count is right and the
  independent content is smaller, exactly as with the bracket holdout.
  [scripts/rigidity_rationality_census.py;
  results/data/rigidity_rationality.json; tests/test_rigidity_rationality.py]
- 🎯 TARGET, not proved: if the support-restricted fixed-rho fiber modulo
  gauge is a single point (true isolation, not first-order rigidity), that
  point is defined over Q. Two open gaps, both stated in the write-up:
  (a) first-order rigidity is not isolation, needing stress-matrix or exact
  isolation certificates; (b) isolation gives finiteness, not uniqueness, so
  the honest form is "algebraic of degree at most the number of isolated
  fiber points modulo gauge, rational when unique". v_B instantiates (b):
  degree 2 over Q, orbit size 2, from the quadratic wall 85t^2 - 70t - 119.
- ❌ CORRECTION: the diagonality gate fails on 154 records, not one. The
  handoff that opened this line expected exactly one violation (v_B) and
  instructed that the expectation be verified; it is wrong. Every
  INTERFERENCE record except v89 and v103 has a non-diagonal 1-RDM. All 154
  still attain the spectrum exactly at unit norm, the library certificate
  never asserted diagonality for them (it is the characteristic-polynomial
  identity), and RESEARCH.md already recorded the same fact from the other
  side in the cancellation-regime count. So this is how the interference
  library was built, not a regression. What genuinely needed repair was the
  support-taxonomy claim, corrected above.
  [scripts/attainer_gate_audit.py; results/data/attainer_gate_audit.json]
- 🟦 REPAIRED (2026-07), additively. All 156 interference records now have a
  natural-orbital representative in results/data/states_natural_orbital.jsonl,
  each with 1-RDM exactly diag(lambda) (worst off-diagonal 2.8e-16) and each
  passing the gate, re-verified from the shipped amplitudes. The 2-RDM
  spectrum is unchanged to 1.6e-15, so the rotation moved the basis and not
  the state. Convention note: the compound is of U^T, NOT U^H; they agree for
  real rotations, so a wrong convention fails only on the complex records.
  THE COST: support grows on 142 of 156 (mean 6.8 -> 12.1, max 14 -> 79) and
  the exact closed forms do NOT survive the rotation, so the new ledger is
  NUMERICAL. states.jsonl is therefore NOT replaced. The two ledgers bound
  different quantities: the census column bounds s_Q^free, the new ledger
  bounds s_Q^NO, and any s_Q^NO claim must quote the latter.
  [scripts/natural_orbital_ledger.py;
  results/data/states_natural_orbital.jsonl;
  tests/test_natural_orbital_ledger.py]
- 🟦 CLOSED, PREDICTION CONFIRMED: the v_B Galois embedding test. Both roots
  of the wall polynomial 85t^2 - 70t - 119 give positive weight vectors that
  attain the vertex EXACTLY, each by 64 of 128 sign patterns modulo global
  sign. Orbit size 2, weight field degree 2: degree equals orbit size, which
  instantiates gap (b) of the Rigidity-Rationality target. The test resolves
  the Galois action on the wall parameter t rather than on nested surds,
  which is precisely what the naive attempt got wrong.
  [scripts/vb_galois_orbit.py; results/data/vb_galois_orbit.json]
- (superseded, kept for the record) PREREG-PENDING: the v_B embedding test. The naive substitution
  sqrt(35) -> -sqrt(35) fails attainment (spectrum error 0.068), and that
  failure is DIAGNOSTIC, not refuting: the amplitudes carry nested surds like
  sqrt(120 - 18 sqrt(35)) and independent principal-branch substitutions are
  not a field embedding. Since the fiber equations are polynomial over Q, a
  genuine embedding carries exact solutions to exact solutions, so a
  substitution that breaks attainment proves only that it was not an
  embedding. The proper test (build the field, take its real embeddings,
  apply each consistently, test attainment, orbit size pre-registered) is
  specified and NOT YET RUN.

## Reconciliation ledger (2026-07 fiber session)

Each item is re-derived by scripts/reconcile_claims.py into
results/data/reconciliation.json, so a reader can check it rather than trust
this table.

| superseded reading | corrected reading |
|---|---|
| v103 is first-order rigid | rigid for the FIXED-RHO fiber, deformable (tangent 6) for the FIXED-SPECTRUM fiber; v89 rigid in both |
| v103 is an isolated support-14 point | true of the route-B search landscape and gauge slice, false of the fiber: support 10 is reached and certified exactly |
| the census support column is the minimal support | it is the certified support, an upper bound on s_Q^NO; the cascade's 71 improvements are improvements of s_Q^free, a different quantity |
| the certified support column was never minimized over the natural-orbital gauge, so it is loose | it was never minimized, and minimizing it changes nothing: 0 of 799 sparsify, 670 of 799 are certified gauge-minimal |
| s_Q^NO(v103) <= 10 (the support-10 endpoint read as a natural-orbital state) | FALSE. That endpoint's 1-RDM carries rho_39 = -5/68 exactly, so it certifies s_Q^free(v103) <= 10 only; s_Q^NO(v103) stays bounded by the library support 14 |
| the cascade endpoints are all the library state in another basis | 63 of 71 are; the other 8 have 2-RDM spectra no one-body rotation can reach, so they are genuinely new extremal states |
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
certification (all 71 after quadratic recognition, same numerical endpoints
throughout). v_B's own numbers are
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

## Exact RDM observability atlas (2026-07)

The census now has a second role: a certified sparse correlator-reconstruction
atlas with a phase-synchronization interpretation.
For 771 of 799 shipped determinant supports, including all 156 INTERFERENCE
supports, the 2-RDM rationally reconstructs every pure full-support state on
the fixed carrier and therefore every higher RDM. The proof uses two exact
support gates: full column rank of pair incidence and a connected spanning
tree of collision-free 2-RDM coherence channels. Pair incidence has full
column rank on all 799 supports.

This is a fixed-frame, carrier-restricted, pure-state statement. It does not
extend to arbitrary mixed states, infer a carrier, prove non-injectivity for
the 28 unresolved supports, or supply fixed-rho dynamics. Its immediate use
is exact higher-RDM reconstruction and carrier-preserving reduced-dynamics
benchmarks. The high-value dynamic gate, a correlation-changing two-body
trajectory that remains in a certified carrier and is not a one-body orbit,
is met below for \(v_B\).
The remaining gates are an infinite-family and stability theorem.

[docs/rdm_observability_census.md;
scripts/rdm_observability_census.py;
scripts/verify_rdm_observability_standalone.py;
results/data/rdm_observability_census.json;
tests/test_rdm_observability_census.py]

## Exact scheduled v_B carrier dynamics (2026-07)

The main dynamic gate above is now closed for one exact finite model. A
Hermitian pair kernel with a seven-entry symbolic support drives the known
\(v_B\) conic, preserves its complete eight-determinant carrier in the full
126-dimensional four-particle sector, and satisfies the Schrödinger equation as a vector
identity. A diagonal orbital lift gives a second at-most-two-body schedule
whose reference trajectory has one constant complete 1-RDM, not merely a
fixed diagonal or spectrum.

The basis-free invariant
\(\operatorname{Tr}(\Gamma_2^3)=6(2899-8t)/12167\) varies strictly with
\(t\). Therefore the loop is not a one-body orbital orbit, including an orbit
inside the stabilizer of the fixed 1-RDM. Since the disconnected two-body
contribution is fixed by that complete 1-RDM, the changing part is the
two-body cumulant.

The carrier has a unimodular pair-incidence minor and a connected seven-edge
2-RDM coherence tree. Therefore the scheduled flow remains on a chart with an
exact rational \(\Gamma_3(\Gamma_2)\) map and gives exact closed nonautonomous
TD-2RDM evolution. This is a reference-schedule result, not an autonomous
Hamiltonian, a mixed-state closure, or a scalable speedup. The package proves
closure but does not yet implement an amplitude-free reduced propagator. The
next gates are that executable reconstruction, conditioning, additional
nontransport carriers, an infinite family, and a beyond-full-wavefunction
benchmark.

[docs/vb_carrier_hamiltonian.md;
scripts/vb_carrier_hamiltonian.py;
scripts/verify_vb_carrier_hamiltonian_standalone.py;
results/data/vb_carrier_hamiltonian.json;
tests/test_vb_carrier_hamiltonian.py]
