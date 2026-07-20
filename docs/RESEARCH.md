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
- Certified closed-form extremal states for 726 of the 799 vertices: every
  design vertex (integer and real) plus every interference vertex that has been
  reduced to a closed form, including all corners solved by the constructive
  off-diagonal-target exactifier (stage 3b) and v_B itself. The 73 uncertified
  vertices are all TIMEOUT / SOLVE-FAIL at Tier A (state-finding), a compute
  frontier, not an exactify or open-math frontier.
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
   (stage 3b) leave no NO-EXACT residual. What remains uncertified (73 of 799)
   is Tier A: TIMEOUT and SOLVE-FAIL vertices whose sparse support the block /
   clique search has not yet found, cleared by a longer clique-timeout and wider
   block search (--max-cliques 0) on stronger hardware, not by new mathematics.
   Recognition layer is unit tested on the v_B cosine and the pi/8 realization;
   the constructive solver is exercised by the shipped closed forms it certifies.
2. Stage 1, the constraint generator for d >= 11: docs/stage1_klyachko_spec.md
   holds the extracted algorithm (Theorem 3.2.1, cubicle extremal edges,
   Schubert coefficient test via lrcalc). Test-first: it must reproduce the
   five known N=3 systems before any new output counts. Until then
   constraints ship as a lookup (src/gpc_census/data/constraints.json).
3. The (3,11) bracket (docs/bracket_3_11.json): 17 vertices already
   certified true by inner-hull membership; 33 gap vertices await either
   restriction tests against rank-10 or Stage 1.
4. Real-attainability certification per vertex is open research; the
   historical scripts/interference8_1.py is the v_B-specific instance
   (exhaustive two-block real stratum sweep). enumerate_designs in
   classify.py is the generalizable exhaustive layer.

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

## Adjacent literature to respect

Liebert, Lemke, Altunbulak, Maciazek, Ochsenfeld, Schilling,
PRResearch 7, 023247 (2025), arXiv:2502.15464: spin-adapted GPCs for
larger systems. Different lattice of problems (spin sectors) from the
spinless census here; draw the boundary explicitly in any new paper.
Maciazek and Tsanov (J. Phys. A 50, 465304): doubly-excited inner bounds
coincide with the polytope for (3, d <= 7); the census confirms this
empirically (those systems are interference-free).
