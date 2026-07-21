# Research state and continuity notes

This file encodes the working understanding of the research program so any
agent or collaborator can continue from the command line. Read this before
touching the science. House rules live in AGENTS.md.

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
- Certified closed-form extremal states for 785 of the 799 vertices: every
  design vertex (integer and real) plus every interference vertex that has been
  reduced to a closed form, including all corners solved by the constructive
  off-diagonal-target exactifier (stage 3b), v_B itself, and 14 closed by state
  transport (scripts/transport_states.py: 10 by padding/frozen-core lift, 4 by
  particle-hole duality; see the leads section). The 14 uncertified vertices are
  all TIMEOUT / SOLVE-FAIL at Tier A (state-finding), a compute frontier, not an
  exactify or open-math frontier. A wider clique sweep
  (max_clique=4) closed many that the k=3 run left as SOLVE-FAIL.
- v_A (16,16,16,6,6,6,6,6,6)/21 is DESIGN-INT; v_B (20,14,14,14,14,4,4,4,4)/23
  is INTERFERENCE with cos(gamma) = 3/(4*sqrt(14)). The v_B phase is not
  special: it is the k=2 instance of the off-diagonal-target mechanism that
  fixes every interference phase (stage 3b).
- Interference onset: absent at rank 8 in the N=4 series, first appears at
  rank 9. Fractions stabilize across ranks 9 to 10 within each series
  (34-37 percent at N=3, 16 at N=4, 14 at N=5).
- Functorial structure: trailing-zero padding and frozen-core lifts
  (lambda -> (1, lambda)) preserve verdicts and generate most interference
  from lower-rank originals; at (4,10) they generate all of it.

## The validation law (non-negotiable)

Every real error in this program was caught by structural invariants, never
by inspection: a sign-parse bug (caught by particle-hole self-duality of
(5,10)), a wrong classify port (caught by a Farkas feasibility
contradiction), and a wrong solver gradient (caught by finite
differences). Consequences:
- No result ships without passing gpc_census.validate and the test suite.
- The preflight gate: any state-solving change must reconstruct v_B and
  certify its closed form end to end (scripts/solve_all.py --preflight)
  before campaigns run. The gate uses the weights-first solver
  (solve_vertex_exact_first with certify_tier_b): Tier A must attain the v_B
  spectrum exactly and Tier B must exactify to a certified closed form (v_B
  has realizations with a pi/8 interference phase). It completes in minutes.
  The historical attain-based path (--legacy-preflight) is a Tier-A-only
  regression check: it reaches the spectrum numerically but lands on an
  arbitrary orbital gauge, so exactify certifies design vertices (trivial
  phases) but not interference vertices, and it is minutes-to-hours slow;
  it is not the gate.
- Pin solvers. Census verdicts were produced under ortools 9.15.6755
  (a required runtime dependency); the CBC fallback labels its backend in
  every verdict.

## The state-construction algorithm (stitched pipeline)

State construction is a routed pipeline, not a single solver. scripts/solve_all.py
drives it; each stage below is a named function, and the stages compose so that
cheap, certain work is never redone by an expensive stage.

0. Route by verdict (already computed in the census). A DESIGN-INT vertex is
   built directly from its classification witness: solve_design_vertex takes the
   integer design weights k and returns psi = sum sqrt(k_t/den) |t>. The design's
   support is one-hop free, so the 1-RDM is diagonal and equals the sorted
   spectrum by inspection, exact on the natural grid, no iterative solve. Only
   DESIGN-REAL and INTERFERENCE reach the solver.

1. Block-budget preflight (min_block_count). The combinatorial half of the
   solver run as pure feasibility, no phase solve: the smallest number of 2x2
   natural-orbital blocks for which some block ansatz has a degree-feasible
   support that is one-hop free off the blocks. 0 = design, k >= 1 = the block
   budget interference needs, None = outside the block-ansatz family (the
   extended-ansatz frontier). Validated: v_A -> 0, v_B -> 1, and the 4_9
   interference vertices the solver cannot yet reach -> None. On None the solver
   fails fast (seconds) with an explicit frontier reason instead of sweeping.

2. Weights-first solve (solve_vertex_exact_first). Sweeps the block ansatz
   family at the predicted budget: block_ansatze enumerates 2x2 blocks that mix
   two spectrum values e1 > e2 by an integer split (a, b) with a + b = e1 + e2
   and off-diagonal sqrt(a*b - e1*e2)/den, so the block eigenvalues are exact by
   construction. Support is confined to the degenerate-signature closure (a
   sound superset of the true support) with off-block hops forbidden. Moduli are
   fixed at sqrt(k/den); only the phases are optimized, as the smooth quartic
   |rho - T|_F^2 with analytic gradients (no eigendecomposition, immune to the
   degeneracy flatness). Verdicts of the census are the input; the census-time
   attain solver is retained only as a Tier-A cascade fallback for frontier
   vertices.

2b. Block-size generalization (k x k), clique_ansatze / phase_solve_clique,
   enabled by max_clique. Some vertices need a natural-orbital rotation that
   mixes k >= 3 distinct eigenvalue classes at once, not just 2x2 blocks; the
   2x2 sweep flags them min_blocks = None. The general block is a k-mode clique
   whose canonical diagonal (integer occupations) is any vector majorized by the
   block's k eigenvalues. The Schur-Horn theorem (I. Schur, Sitzungsber. Berl.
   Math. Ges. 22, 9 (1923): the diagonal of a Hermitian matrix is majorized by
   its eigenvalues; A. Horn, Amer. J. Math. 76, 620 (1954): the converse, every
   majorized vector is a diagonal) characterizes exactly the realizable
   diagonals, and _schur_horn_diagonals enumerates them (the 2x2 split is the
   k=2 case; permutation diagonals are kept, since for k >= 3 they admit
   non-diagonal realizations). A clique mixes distinct eigenvalues, so its block
   is non-degenerate and phase_solve_clique matches it by its
   characteristic-polynomial (elementary-symmetric) coefficients via power sums
   and Newton's identities, with analytic gradients (L-BFGS), again without
   eigendecomposition. A clique-size preflight (min_clique_count, pure
   feasibility from k = 3 up) sets the starting size and fails fast off the
   single-clique family. Validated end to end and UNAIDED: idx 24 of (4,9),
   (9:6:5:5:5:2:2:1:1)/9, which no 2x2 configuration reaches, reconstructs from
   the spectrum alone (residual 1e-31, ~160 s) with a single 3x3 clique. The
   pipeline finds its own realization (clique on modes {1,4,8}, eigenvalues
   (6,5,1)), not the one an exploratory prototype had been handed (modes {0,1,4},
   (9,6,5)), which is the evidence that the search is unbiased. It certifies end
   to end: idx 24's closed form is REAL, psi = (|0125>+|0134>+|0237>+|0245>)/3 +
   sqrt(2)/3 |0268> + sqrt(3)/3 |0348>, Tier B EXACT by the existing gauge-fixed
   exactify with no extension. This is a structural point about Tier B, not a
   one-off: a k-clique has a (k-1)(k-2)/2-dimensional Schur-Horn fiber, so for
   k >= 3 a REAL realization is often available (idx 24 is one), and a real state
   gauge-fixes to phases 0/pi and certifies trivially. So "interference" at k >= 3
   frequently means only that the support cannot be made one-hop free, not that
   complex phases are forced (contrast v_B's genuine pi/8 at k = 2, where the
   tight 2x2 fiber leaves no real solution). Consequences: Tier B already reaches
   most k >= 3 vertices for free through real-realization detection; the
   genuinely hard cases needing higher-degree algebraic phases (no real, no
   p*sqrt(q)/r) are a smaller residual, now also handled by PSLQ (see stage 3).
   Multi-clique ansatze (multi_clique_ansatze, several disjoint cliques on
   disjoint orbitals, one-hop free off the union) are implemented behind
   max_cliques, with a per-clique diagonal cap and the escalation 1 -> 2 -> ...
   cliques. Status: the machinery runs, but does not yet crack the (4,9) FAILs
   (idx 15, 40, 42). The obstruction is that degree feasibility of a multi-clique
   structure does NOT imply phase solvability: for idx 15 a probe found a
   degree-feasible two-3-clique structure, yet 425 diagonal combinations x
   several skeletons each all failed to phase solve. So finding the RIGHT
   multi-clique structure (which classes, which placement) is a search problem
   for which feasibility is too weak a guide; that search, and clique placement
   beyond one representative orbital per class (WLOG for a single clique by
   degenerate-class symmetry, not for several sharing a class), are open.

3. Exactify (Tier B, gpc_census.exactify). Squared amplitudes snap to k/den.
   The state is defined only up to the single-particle U(1)^d phase gauge
   (c_t -> c_t prod_{m in t} exp(i phi_m), a diagonal-unitary conjugation of the
   1-RDM), so exactify gauge-fixes first (projects the phases orthogonal to the
   gauge orbit) and then recognizes the residual interference phases (rational
   multiples of pi, cosines on the p*sqrt(q)/r lattice, then general low-degree
   algebraic cosine and sine via PSLQ integer-relation detection on their powers,
   recognize_algebraic). The PSLQ branch is a proposal only, gated by
   verify_exact, so a spurious hit from float precision cannot certify a wrong
   state. verify_exact is gauge-invariant (a characteristic-polynomial identity),
   so the certificate
   stands in the fixed frame. Failures fall through labeled TIER-C for hand
   analysis. v_B certifies as a single 14/4 block with a pi/8 interference
   phase; solve_vertex_exact_first(certify_tier_b=True) searches realizations
   until one exactifies, since a vertex has many realizations differing only in
   interference phase and some carry a clean lattice.

3b. Constructive interference exactify (exactify_interference). Per-phase
   recognition (stage 3) fails on a residual of interference corners whose
   phases are coupled polygon angles: after gauge-fixing, each determinant's
   absolute phase is a path sum of several relative phases and so of high
   algebraic degree, while a generic numeric solve lands on an arbitrary point
   of the gauge orbit, so recognize_phase and PSLQ see only noise. The fix is to
   stop treating the phases as the unknowns. The moduli are already exactly
   rational (|c_t|^2 = k_t/den), so the diagonal occupations of the 1-RDM are
   known; wherever two orbitals p, q carry an equal occupation that the spectrum
   must split, the 2x2 block [[occ_p, x],[conj(x), occ_q]] is forced to
   eigenvalues lo, hi with lo + hi = occ_p + occ_q, which pins the off-diagonal
   magnitude |x| = sqrt(occ_p occ_q - lo hi) to an exact algebraic number
   (rational or a low-degree surd). This is the interference analogue of a
   design's diagonal 1-RDM (_active_offdiagonals computes these targets). Each
   active off-diagonal is then a closed polygon whose sides are the fixed term
   moduli, so its relative phase is an exact arccos of that target; phases
   propagate across edges by constraint propagation (an edge with one unassigned
   determinant fixes it as arg(P) + arccos(rhs/|P|) of its already-fixed partial
   sum P), and edges that share every determinant are closed by a small bounded
   joint search over the same algebraic angles. Every candidate is gated by
   verify_exact, so the construction cannot certify a wrong state. This is the
   same mechanism as v_B's cos(gamma) = 3/(4 sqrt(14)); v_B was the first visible
   instance, not a special case. It closed the entire NO-EXACT residual: all 17
   interference corners that had a numeric state but no recognized closed form
   (single 2x2 edge with phases 2pi/3, arccos(+-3/4), arccos(-7/8),
   arccos(sqrt(6)/6), arccos(7 sqrt(30)/60), and coupled two-edge corners
   combining an arccos with a pi/3 triangle closure) now certify in about two
   seconds total. The remaining uncertified vertices are all TIMEOUT / SOLVE-FAIL
   at the state-finding stage (Tier A), not exactify: they are compute-bound
   (bigger clique-timeout, wider block search on stronger hardware), not open
   research.

   verify_exact robustness (needed for any genuinely complex closed form): the
   char poly is taken as the Berkowitz determinant of (x I - rho), because
   sympy's charpoly() block-factorization tries to order complex exp() factors
   and raises; and rho is expanded to rectangular a + b i form first, and each
   coefficient of the difference reduced with expand_complex before the zero
   test, because sp.simplify alone leaves true zeros like 1 + exp(2 i pi/3)
   unreduced and had been silently rejecting valid interference closed forms.

Historical note: Altunbulak and Klyachko (2008) resolved the (4,9) vertices and,
per Sec 6.2.2, determined extremal states for all but v_A and v_B, which were
left numerical-only. Those rank-9 state data are not in the paper text (dead zip
/ paywalled ESM), so this pipeline reconstructs them independently, and the
design/interference trichotomy is this program's framing, not theirs: a vertex
can be interference (phase cancellation forced, as in the surd-coefficient
rank-7 states they did publish) yet still have a known closed form. v_A and v_B
were the two they could not close at all.

## Current campaigns and open items

1. State construction (Tier A) and exactification (Tier B) of all vertices:
   scripts/solve_all.py (--all writes the full census to
   results/data/states.jsonl; the interference-only default writes
   states_interference.jsonl); scripts/build_states.py is the checkpointed,
   resumable campaign driver. Design vertices certify from their witness. Tier B
   now certifies every interference vertex whose support Tier A finds: the
   per-phase recognition layer plus the constructive off-diagonal-target solver
   (stage 3b) leave no NO-EXACT residual. What remains uncertified (14 of 799)
   is Tier A: TIMEOUT and SOLVE-FAIL vertices whose sparse support the block /
   clique search has not yet found, cleared by a longer clique-timeout and wider
   block search (--max-cliques 0) on stronger hardware, or by transporting a
   certified sibling's state (scripts/transport_states.py), not by new
   mathematics.
   Recognition layer is unit tested on the v_B cosine and the pi/8 realization;
   the constructive solver is exercised by the shipped closed forms it certifies.
2. Stage 1, the constraint generator for d >= 11: docs/stage1_klyachko_spec.md
   holds the extracted algorithm (Theorem 3.2.1, cubicle extremal edges,
   Schubert coefficient test via lrcalc). Test-first: it must reproduce the
   five known N=3 systems before any new output counts. Until then
   constraints ship as a lookup (src/gpc_census/data/constraints.json).
3. The (3,11) bracket (docs/bracket_3_11.json): 46 of the 50 outer candidates
   are now settled WITHOUT Stage 1, in exact arithmetic and independently
   verified here (scripts/settle_bracket_3_11.py regenerates
   docs/bracket_3_11_settlement.json). 19 TRUE (17 by state transport from the
   rank-9/10 census, all verify_exact-checked; the uniform (3/11)^11 by an
   explicit Z11 difference design on base block {0,1,3}, one-hop-free; and a
   frozen-core paired state for (1,(1/5)^10)); 27 REFUTED (18 by frozen-core
   pinning to N=2 plus the even-degeneracy pairing theorem, 9 by
   zero-restriction, each restricted point violating the known lower-rank GPCs,
   checked against gpc_census.constraints); 4 OPEN (cands 23, 26, 34, 44, all
   genuinely rank-11 full-support; cand 23 admits no one-hop-free design by
   weight counting on its two 6/7 modes, so INTERFERENCE if attainable). CAVEAT:
   the 27 refutations do NOT finish (3,11) as a polytope: cutting them creates
   new true vertices absent from the outer list, so Stage 1 (or a tighter
   bracket) is still required. The refuted and certified points are a mandatory
   acceptance oracle for any (3,11) constraint generator: it must exclude every
   refuted point and retain every certified one.
4. Real-attainability certification per vertex is open research; the
   historical scripts/interference8_1.py is the v_B-specific instance
   (exhaustive two-block real stratum sweep). enumerate_designs in
   classify.py is the generalizable exhaustive layer.
5. Vertex-cone arithmetic (positive-geometry lead): the exact lattice
   invariants of simple vertex cones are computed for all systems; see the
   section below. Open: prove den(v) = normal-cone index for simple
   interference vertices; extend edge-cone SNF to non-simple vertices via
   pulling triangulation; compute the full canonical form (adjoint) of
   Pi_{3,7} as the smallest non-simplex.
6. Holonomy Galois structure: every certified interference holonomy in the
   analysis snapshot is 2-elementary abelian; see the section below. Open:
   prove the norm-square law from the stage-3b polygon-closure mechanism;
   certify the remaining Tier-A stragglers (14 after the k=4 sweep and state
   transport, down from 73 at analysis time) and re-test (the only place a
   non-abelian D4 holonomy could hide); Groebner elimination for the exact
   antiunitary overlap of psi_B.
7. State transport (VERIFIED and APPLIED, scripts/transport_states.py):
   padding and frozen-core lift transport states, not just verdicts; this
   closed 14 Tier-A stragglers for free. See the dedicated section below.

## Vertex-cone arithmetic and holonomy Galois structure (July 2026 leads)

Two exact-computation campaigns, run against the shipped census data
(reproduction: leads_analysis.py; all statements exact rational/symbolic
unless marked numerical). These are LEADS, not shipped results: the
arithmetic layer (norm identities, minimal polynomials, prime
factorizations) is hand-verified in this repo, but the group-theory and
lattice conclusions (Galois groups via PARI polgalois, Smith normal forms)
are as computed by the analysis script and are not independently re-run
here. Snapshot note: Lead B's per-vertex counts are as of a 770-certified
snapshot (127 certified interference vertices); the certified set has since
grown (128 interference, 14 total residual and shrinking under the k=4
sweep). Lead A's counts are classification/geometry, independent of how
many states are certified, so they are stable. Any use of these in the
paper (results/report/main.md) is DEFERRED and must first refresh the
certification-dependent counts.

### Lead A: positive-geometry / lattice structure of vertex cones

Full H-representation = published GPC inequalities + ordering walls
(lambda_i >= lambda_{i+1}) + lambda_d >= 0 in the slice sum = N. Exact
vertex-facet incidence and edge graphs for all nine systems give:

- Simpliciality census: 179 full-dimensional simple vertices across the
  systems (e.g. 25 of 103 at (4,9)). v_A and v_B are both simple.
  Simpliciality does NOT correlate with the trichotomy, so local
  canonical-form degeneracy (adjoint at a non-simple vertex) is not an
  interference detector.
- v_B tangent cone: edge-ray lattice index 23^7, Smith normal form
  diag(1, 23 x7), i.e. the affine toric variety is the quotient
  singularity C^8/(Z/23)^7. Dual (normal) cone: index exactly 23 = the
  natural denominator, SNF (1,...,1,23), cyclic. Every nonzero
  facet-ray pairing <W_j, r_k> equals exactly 23 ("isotypic purity").
  The local canonical form in lattice-normalized facet coordinates is
  Omega = 23 d^8y / prod_j <W_j, y> (numerator 23^7 in the dual
  vertex-expansion frame); no adjoint numerator, all structure is
  arithmetic at the prime 23.
- LAW (exact, all 179 simple vertices checked): den(v) = normal-cone
  index holds for EVERY simple INTERFERENCE vertex (54/54) and fails
  only at design vertices (68 of 125 simple design vertices violate it,
  typically nidx = 2 den or 3 den). Theorem candidate; no proof yet.
- Isotypic-pure vertices (all pairings = den, edge quotient
  (Z/den)^{dim-1}): 23 across all systems; every one has PRIME
  denominator (7, 11, 13, 17, 23, 37). All pure interference vertices
  have prime den >= 11; all pure design vertices have den = 7 (the
  (3,3,3,3,3,3,3)/7 vertex and its pads/lifts). Converse fails:
  prime den does not imply purity ((4,9) v91, den 13, is not pure).
- rad(den) divides the edge-cone index at every simple vertex checked.
- Pi_{3,6} is a 3-simplex (adjoint degree 0), consistent with the
  paper's positive-geometry table; Pi_{3,7} (dim 6, 10 vertices,
  10 facets, adjoint degree 3) is the smallest target for a full exact
  canonical form with nontrivial adjoint.

### Lead B: Galois structure of interference holonomies

Amplitude phases are diagonal-gauge dependent; the invariant is the loop
holonomy: left-kernel vectors of the support incidence matrix (dets x
orbitals) applied to the phase vector, in the natural-orbital basis the
census fixes.

- v_B: kernel dim 1, kernel vector (1,-1,0,0,-1,1,0,0) on the 8 certified
  support dets; holonomy = atan(sqrt(119)/3), cos = 3 sqrt(2)/16 (so the
  certified-form phase acos(-3 sqrt(2)/16)+pi is exactly the invariant
  content). Minimal polynomial of e^{i hol}: 32x^4+55x^2+32 (hand-verified);
  splitting field Q(sqrt(2), sqrt(-119)); Galois group V4 (Klein four);
  119 = 7*17. The block angle cos(gamma) = 3/(4 sqrt(14)) is a different
  parametrization (field Q(sqrt(14), sqrt(-215)), also V4; 215 = 5*43).
  Note the primes decouple from the geometric prime 23 of Lead A.
- Census-wide (127 certified interference vertices in the snapshot): 80
  supports are loop-free (kernel dim 0; signs suffice), 47 carry 59
  nontrivial holonomies with 24 distinct minimal polynomials. Every
  holonomy field is 2-elementary abelian: Z/2 (degrees 1-2), V4 (all 11
  quartics), and (Z/2)^3 for ALL FIVE degree-8 polynomials (PARI polgalois,
  package pari-galdata). Degree-4 => V4 is forced by unit-circle
  reciprocity (splitting field Q(sqrt(a(2a-b)), sqrt(-a(2a+b))) for
  ax^4+bx^2+a); the degree-8 result is NOT forced and is the contentful
  finding: no Z/4, no D4, no Q8 anywhere in the certified census snapshot.
- Norm-square law (exact, all three occurrences, hand-verified): every
  nested radical sqrt(a + b sqrt(d)) in a holonomy cosine satisfies
  a^2 - b^2 d in (Q^x)^2: 11^2-21 = 10^2; 817^2-1633 = 816^2;
  2651^2-48^2*2769 = 805^2. This is exactly the criterion separating
  abelian V4 quartic subfields from generic D4. Conjecture: forced by
  stage 3b (the conjugate radical is itself a physical off-diagonal
  magnitude; the product is a rational norm of the exact moduli). If
  proved: certified interference holonomies are always abelian of exponent
  2. A D4 counterexample can only hide in the still-uncertified Tier-A
  vertices (14 now, 73 at analysis time), which is a shrinking target.

### Correction found: the antiunitary overlap of psi_B (numerical)

The paper's Remark on non-realifiability reports max overlap ~ 0.933184
over the reduced group U(1) x U(4) x U(4). Re-running the optimization
(compound-matrix objective |<psi, Lambda^4 V conj(psi)>|, multi-start
L-BFGS, independent parametrizations):
- Within the SAME reduced group: 0.936497506 from 6+ independent seeds
  (also under the symmetric/involutive restriction). The published
  0.933184 appears to be a premature convergence.
- Over full U(9): 0.991732701478395695713606553348762758... (43 digits,
  50-dps Newton on the exact symbolic reduced objective). The maximizer
  is a single symmetric 2x2 block MIXING modes {4,8} across the two
  degenerate eigenspaces, plus diagonal phases.
- REVIEW NOTE (rigor): the non-realifiability CONCLUSION rests on the
  REDUCED-group max being < 1, because any real form would give overlap
  exactly 1 and its V must map eigenspaces to eigenspaces (hence lie in
  the reduced group). So 0.936498 < 1 is the sufficient, complete, and
  safer evidence. The full-U(9) value 0.991733 is the true orbit invariant
  ("closest any rotation gets") but is (a) not what the conclusion depends
  on and (b) uncomfortably close to 1; and the value has crept upward under
  successive search (0.9332 -> 0.9365 -> 0.9917), the signature of a hard
  optimization whose supremum could be higher still. The paper should lead
  with the reduced-group 0.936498 as load-bearing and present 0.991733 as a
  secondary, explicitly numerical remark. Both are < 1, so the paper's
  CONCLUSION (no real form under any orbital rotation) is unaffected; only
  the number and its labeling need updating. Paper edit DEFERRED to end of
  session.
- The value is algebraic in principle (critical value of a polynomial
  system) but resists PSLQ at degree <= 8 / height <= 1e10 with 43
  digits; exact identification needs Groebner elimination on the
  stationarity system (open).

### State transport (VERIFIED and APPLIED): clears 14 Tier-A stragglers

Unlike the two leads above, this is verified in this repo and applied to the
shipped data. The verdict-transport maps also transport STATES, not just
verdicts. Trailing-zero padding leaves the state untouched (one more empty
orbital, eigenvalue 0); a frozen-core lift lambda -> (1, lambda) prepends an
orbital occupied in every determinant (1-RDM becomes 1 (+) rho, eigenvalue 1).
So two vertices with the same CORE (the eigenvalues strictly in (0,1)) are
transport-equivalent: strip a certified donor to its core state (drop the
frozen and empty orbitals) and re-add the target's frozen/empty orbitals.

Verification is the same exact characteristic-polynomial identity Tier B uses
(verify_exact against the target spectrum), so a transport is accepted only if
it certifies exactly. Applied over the residual, this cleared 14 vertices (10 by padding/frozen-core lift, 4 by particle-hole duality in the self-dual (5,10) system),
several from a HIGHER-rank donor than the failing root (e.g. the k=4 sweep
certified (3,10) v47, whose core state transports down to certify (3,8) v15,
which had SOLVE-FAILed at rank 8): the rank-10 run succeeded where the rank-8
run did not, and the state moves down for free. Certified transports (tierB
EXACT, provenance in a `transport` field), donor -> target:
  (3,10)v47 -> (3,8)v15, (4,9)v15, (5,10)v15
  (3,10)v63 -> (3,8)v20, (4,9)v20, (5,10)v20
  (3,10)v30 -> (3,9)v11, (4,10)v11
  (4,10)v58 -> (4,9)v38, (5,10)v38
(3,9) v22 was already certified directly by the k=4 sweep. After transport the
residual is 14, and (3,8) and (3,9) are fully certified. Implementation:
scripts/transport_states.py (dry run by default, --apply to merge). Future:
fold this in as a stage of solve_all.py that runs before Tier A on any failing
vertex (both directions), and add a nonzero-eigenvalue count check to validate
alongside the power sums.

### Geometry of the remaining 14 (lead)

The 14 survivors skew strongly non-simple (tangent-cone complexity is a genuine
difficulty predictor: the uncertified set is far more non-simple than the
certified interference set). The simple survivors carry full lattice data and a
speculative sweep prior from Lead A: their non-den facet-ray pairings (a prime
sub-pairing such as 13 | 26, 17 | 34, 13 | 39) suggest trying block/clique
moduli denominators at the pairing values, not only at den. A second empirical
prior from Lead B: most certified interference states are loop-free (signed-real
suffices), so order the skeleton sweep signed-real-first before phase solves.

### Why the 14 fail: two strata (diagnosis, verified against the shipped data)

The 14 SOLVE-FAILs in states.jsonl (785/799 certified) split into two distinct
failure modes, both checked against the current data and the solver source. This
sharpens the standing "cleared by a longer clique-timeout, not new mathematics"
claim: it holds for one stratum and not the other.

- Off-clique cancellation gap (the non-simple majority). These vertices carry no
  feasible 2x2 block ansatz (min_block_count = None) but a single 3-clique is
  degree-feasible (min_clique_count = 3), matching the >=2-channel synthesis
  prediction. But _solve_via_cliques (src/gpc_census/states.py, the one-hop cut
  at the AddBoolOr over hop_pairs) forbids EVERY support one-hop pair that lies
  off the clique union, whereas 80 of the 127 certified interference states
  settle their off-diagonals by SIGNED CANCELLATION on such pairs, not by
  avoiding them. So the k-clique family as implemented cannot express those
  states at all: for these the failure is an ANSATZ-EXPRESSIVENESS gap, not a
  walltime budget, and more hardware alone will not close them. LEAD (proposed,
  not implemented, not verified): add an ansatz axis that ADMITS off-clique
  one-hop pairs required to cancel by signs. With integer weights on the natural
  grid the sqrt(k_a k_b) products cancel iff the off-clique pairs partition into
  equal +/- halves, which is CP-SAT-encodable per off-clique mode pair. This
  would also explain why min_clique_count passes while the sweep fails: the
  preflight tests only the degree system, without the one-hop-off-clique cut.
- max_card cap (the two simple survivors). (3,10) v89 (den 26) and v103 (den 34)
  are the two SIMPLE, isotypically impure survivors from the Lead-A lattice
  census (pairings {2,13,26} and {17,34}); both have min_block_count = 1 yet
  SOLVE-FAIL. Their denominators exceed the solver's default max_card = 24, so a
  legitimate support may need more determinants than the cap admits. This IS the
  genuine longer-walltime / raise-the-cap case: rerun with max_card >= den + 10
  and the Lead-A prior (prioritize block splits aligned with the 13- and
  17-structure). Not yet run to completion here.

### Synthesis: at most two exchange channels (VERIFIED, census-wide)

The cleanest cross-cutting pattern, re-checked here at current counts: every
certified interference state uses at most TWO one-hop exchange channels (active
off-diagonal 1-RDM pairs). Distribution over the 142 certified interference
states: 99 with one channel, 43 with two, ZERO with three or more (max = 2). If
a constant channel bound survives the remaining Tier-A stragglers and higher
ranks, the statement is strong: extremal N-representability requires only O(1)
quantum interference, ever. This bears directly on discussion question (ii)
(the sign-problem / interference gradation), and is a candidate for the paper's
Discussion. Falsifiable prediction for the residual: since the certified solver
sweeps small channel budgets first, the survivors should need >= 3 channels
(equivalently, the <=2 law is either a theorem or a portrait of solver reach).

More speculative synthesis leads (NOT verified here, recorded for follow-up):
the elementary-abelian pattern across three unrelated computations (toric
quotients (Z/p)^k at simple vertices, holonomy Galois groups (Z/2)^k with
k<=3, channel budget <=2) suggests a toric/binomial structure on the extremal
fiber after the Schur-Horn substitution, which would tie the multiquadratic
phase fields, the norm-square law, and the elementary quotients together; the
only measured invariant that is NOT elementary is the antiunitary overlap of
psi_B (no minpoly of degree <=8, height <=1e10). First concrete test: write the
polygon-closure equations of a kernel-dim-1 vertex and check binomiality.

## The selection rule is basis-relative (degeneracy lemma)

At a vertex, every active facet forces the extremal state onto
determinants achieving its bound (the superselection rule; cf. Liebert
et al., PRR 7, 023247). Crucial subtlety discovered empirically on v_B:
the rule holds in the state's OWN natural orbital basis, and a sparse
extremal state generally does not have its natural orbitals aligned with
the canonical modes. Diagnostic history: the strict canonical filter
leaves 9 of 126 determinants for v_B, the resulting skeletons plateau at
residual 1e-4, and none carries the historical weight multiset
[7,4,3,2,2,2,2,1]. A first attempted fix (closure of the admissible set
under degenerate-block signatures) was insufficient, and the failure is
provable: the degree system pins diag(rho) to the spectrum, and if the
eigenvalues also match, Ky Fan equality forces rho = diag(lambda)
exactly, so every canonical-degree skeleton must cancel all off-diagonal
1-RDM entries. Exhaustive CP-SAT enumeration shows no cardinality-8
canonical-degree skeleton with the historical multiset passes even the
polygon necessary condition. The sparse v_B state has DIFFERENT mode
sums: its 1-RDM is block diagonal with 2x2 blocks mixing a 14-mode and a
4-mode, exactly the structure of the historical interference8_1.py
ansatz.

The implementation (states.block_ansatze, states.solve_vertex_exact_first):
sweep a discrete family of exact 1-RDM targets. Each target is the
canonical diagonal or a block diagonal with 2x2 blocks mixing two
spectrum values e1 > e2 via an integer split (a, b): a + b = e1 + e2,
off-diagonal sqrt(a*b - e1*e2)/den, which has exact eigenvalues e1/den,
e2/den by construction. For each target the degree system uses the split
mode sums; support is confined to the degenerate-signature closure (a
sound superset of the true support) with off-block one-hop pairs
forbidden so the 1-RDM off-diagonal stays on the blocks (the historical
interference8 ansatz); skeleton enumeration is deduplicated to one
representative per orbit of the degree-preserving mode permutations, a
polygon inequality per mode pair prunes skeletons that cannot cancel, and
the phase problem is solved as the smooth quartic |rho - T|_F^2 with
analytic phase gradients (no eigendecomposition, hence no degeneracy
flatness). This reconstructs v_B end to end: a SINGLE 2x2 block mixing a
14-mode and a 4-mode (split (6,12), off-diagonal 4/23), support size 8,
weights [8,4,3,3,2,1,1,1]/23, phase residual 1e-32.

Ruled out: the continuous OPERATOR form of the selection rule (put the
state in the b-eigenspace of dGamma(U diag(a) U^T), U the target's
natural orbitals, and take the joint kernel's support). It is ill posed
at degenerate vertices, because the facet coefficients differ across a
degenerate lambda-block while the natural orbitals inside the block are
free, so U diag(a) U^T is not determined by <a, lambda>; checked
directly, the exact v_B state fails dGamma(A_nat) psi = b psi by O(1) on
both active facets. A correct continuous form would project each facet
coefficient vector into the blocks and characterise the invariant
subspaces; that is open. This is a lemma about the method and belongs in
the states paper.

## Source documents in docs/

- stage1_klyachko_spec.md: the constraint-generation algorithm, verbatim
  from the Altunbulak thesis, with the implementation plan.
- ak2008_extracted.md: inequality and vertex tables for ranks 6-8 from the
  2008 paper, including extremal states (surd coefficients need re-checking
  against the PDF before use as ground truth).
- bracket_3_11.json: the Stage-0 inner/outer bracket for (3,11).
- bracket_3_11_settlement.json: 46/50 candidates settled (19 true, 27
  refuted, 4 open) by state transport + frozen-core/N=2 pairing +
  zero-restriction; see campaign item 3. Reproducible via
  scripts/settle_bracket_3_11.py; all 46 witnesses verified in this repo.
- bracket_3_12_true_vertices.json: 19 CERTIFIED true vertices of (3,12),
  each a face embedding of a settled (3,11) true vertex (the {lambda_12=0}
  face of Delta(3,12) is Delta(3,11); a vertex of a face is a vertex of the
  polytope). Rigorous and reproducible from the (3,11) settlement; the
  first vertices of a rank-12 system. LEAD (not shipped, not re-derived
  here): a separate exact-plethysm approach (multiplicity of s_lambda in
  h_m[e_3] via Murnaghan-Nakayama) generates a (3,12) Stage-0 INNER cloud
  of attainable spectra without Stage 1; its extreme points are vertex
  CANDIDATES only (cloud-extremality does not certify true vertexhood, cf.
  (3,7) needing M >= 7), so that cloud needs an in-repo generator and the
  usual acceptance-oracle checks before it ships. The plethysm generator
  complements Stage 1 (inner side), it does not replace it.

## Adjacent literature to respect

Liebert, Lemke, Altunbulak, Maciazek, Ochsenfeld, Schilling,
PRResearch 7, 023247 (2025), arXiv:2502.15464: spin-adapted GPCs for
larger systems. Different lattice of problems (spin sectors) from the
spinless census here; draw the boundary explicitly in any new paper.
Maciazek and Tsanov (J. Phys. A 50, 465304): doubly-excited inner bounds
coincide with the polytope for (3, d <= 7); the census confirms this
empirically (those systems are interference-free).
