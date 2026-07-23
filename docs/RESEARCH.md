# Research state and continuity notes

This file encodes the working understanding of the research program so any
agent or collaborator can continue from the command line. Read this before
touching the science. House rules live in AGENTS.md.

## FIBER EXPLORATION III: anholonomy, the correlation dial, and the
frontier vertices carry fibers too (2026-07)

1. BERRY PHASE around the v_B fiber loop (natural gauge, phase carried
   on the loop determinant): gamma_B = -2.1555 rad = -0.6861 pi, i.e.
   +1.3139 pi mod 2pi -- a NONTRIVIAL anholonomy (neither 0 nor pi).
   Adiabatically cycling the fiber once returns the state with a ~0.69 pi
   geometric phase. Caveats: defined mod one-body gauge windings (the
   mod-2pi value is the invariant); "adiabatic cycling" presupposes a
   Hamiltonian path tracking the fiber -- the selection-transition
   coupling circle is the natural candidate driver, unverified.
2. THE CORRELATION DIAL: mutual information I(2:4) between modes varies
   by a factor ~2000 along the fiber -- 6e-5 bits at the t=0 rational
   point (modes essentially uncorrelated; n24 within 0.2% of n2*n4) to
   0.107 bits at the right wall -- while n2, n4 and ALL single-mode
   entropies are frozen. The fiber coordinate is a correlation dial
   invisible to occupations: correlations redistribute between mode
   pairs with no one-body trace. (Curious datum: the shipped rational
   point sits almost exactly at the uncorrelated point of its family.)
3. FRONTIER FIBERS: the achievable width of <n_2 n_3> at v89's exact
   1-RDM is 0.1917, range [0.0000, 0.1917] -- the newest vertex carries
   the same O(1) underdetermination the day it was attained, including
   states where the two degenerate modes NEVER co-occur (p11 = 0
   exactly attainable). Width includes stabilizer-orbit directions
   (orbit-vs-intrinsic split not yet separated); v103 pending (bigger
   space, optimizer budget).


## SOS CERTIFICATE FOR CAND 44: the Coleman-hull obstruction and the
symmetry-reduced route (session spec shipped)

Lemma (proved by formulation, saves a wasted level): LINEAR duality
cannot certify any rank-11 candidate's non-attainability. Constant
multipliers separate gamma_0 only from conv{gamma(pure)} = Coleman's
set {0 <= gamma <= 1, tr = N}, which contains every candidate (they all
satisfy Pauli). The nonconvexity of the pure-state 1-RDM image is the
entire content, so the certificate requires genuine degree-4
Positivstellensatz (level-2 Lasserre on 330 real variables). Raw size
(~55k moment basis) is intractable; the S5 x S6 spectrum-stabilizer
symmetry (plus conjugation and gauge) block-diagonalizes it to
expected-hundreds-size blocks -- the symmetry-adapted-SOS regime
(SumOfSquares.jl + SymbolicWedderburn.jl / TSSOS). Session prompt with
staged gates shipped (outputs/claude_code_prompt_sos_certificate.md):
G1 pipeline validation on a known-exterior (3,6) point WITH an
attainable-vertex negative control; G2 symmetry-machinery equivalence;
G3 cand 44 with exact-arithmetic rationalization of the dual
certificate (SDP proposes, exact algebra disposes -- census
discipline); G4 the remaining candidates. A rational certificate with
delta > 0 candidate-settles (3,11) at 19 vertices and makes the
claimed level-5 inequality theorem-adjacent without Schubert calculus.


## EXPERT SCRUTINY PASS: claim repairs before the papers (2026-07)

Each recent claim attacked as a referee would; repairs binding.

1. VALLEY CLAIM -- three repairs. (a) WORDING ERROR, fixed: the
   measurement fixes the full 1-RDM gamma = diag(lambda), not merely the
   spectrum. The correct (and STRONGER, and exactly RDMFT-relevant)
   statement: THE 1-RDM ITSELF underdetermines pair correlations by
   O(1) -- the constrained-search domain at fixed gamma contains the
   measured width. Fixed-spectrum-only phrasing would be trivially true
   by basis rotation and must not appear. (b) The widths are certified
   LOWER bounds (inner approximation by local optimization); the true
   ranges may be larger; endpoints are achieved values, and interval-
   ness of the achievable set is UNPROVEN (fiber connectivity unknown)
   -- sample intermediate targets or prove connectivity before writing
   "range". (c) "Maximal pinning locally minimizes the ambiguity" rests
   on one ray and three epsilons: downgraded to observed-along-this-ray
   + conjecture; the eps=0.10 failure is unattributed (optimizer vs
   feasibility) and must not be glossed.
2. SELECTION PRINCIPLE -- gap found and closed with data. The endpoint
   theorem is PER-FAMILY; the physical minimization runs over the FULL
   fiber, which the valley measurement showed is larger than any
   family. New test (this session): global fiber extremes of <n2 n4>
   at v_B are attained by REAL states (complex vs real-restricted
   extremization agree to 8e-7 / 1e-9). Paper language: theorem
   (per-family) + global numerical verification at v_B (one observable,
   one vertex) + conjecture (general), with a mechanism sketch: c ->
   conj(c) preserves gamma and all pair occupations, so extreme points
   plausibly admit time-reversal-symmetric representatives. Prove or
   census this before claiming generality.
3. 2-RDM EMBEDDING -- upgrade the conjugate-sheet separation from a
   spot-check to a one-line proof: the T-odd coherences are
   proportional to sin Theta, nonzero on the open interior, vanishing
   exactly at the walls where the sheets glue. Injectivity of the
   circle then follows from monotone t-resolution + interior sheet
   separation.
4. RANK-11 FLOORS -- strengthen from evidence toward certificate: the
   floors are unconstrained local minima; a Lasserre/SOS relaxation of
   min ||gamma(psi) - gamma_0||^2 yields a CERTIFIED positive lower
   bound (dual certificate) = rigorous non-attainability, independent
   of and faster than the Schubert route. Recommend running the
   degree-2 relaxation on cand 44 (smallest floor) first: if even the
   nearest candidate certifies infeasible, the level-5 claim is
   theorem-adjacent without new mathematics.
5. UNCHANGED after attack: v89/v103 numerical attainment (discipline
   already in place), wall certificates, endpoint theorem (scope
   stated), pre-registration protocol.


## QUASI-FIBER STABILITY DEMONSTRATED: the valley persists and WIDENS
off-vertex (2026-07)

Measured the achievable range of <n_2 n_4> over states with EXACT
spectrum lambda (contraction attack + kappa-continuation extremization,
residual < 1e-13), for lambda walking from v_B into the polytope
interior along (1-eps) v_B + eps uniform:

  eps = 0.00: width 0.2522   (range [0.217, 0.470])
  eps = 0.02: width 0.2627
  eps = 0.05: width 0.2899   (eps = 0.10: optimizer budget, not physics)

CONTENT: the occupation-blindness is not a knife-edge vertex
phenomenon. States sharing an identical one-body spectrum differ in
pair correlations by O(1) across a whole neighborhood of the vertex --
and the width GROWS off-vertex, so maximal pinning is where the
ambiguity is SMALLEST, and even there it is ~0.25. This is the
data-backed form of the stability theorem and the deflation-proof
version of the physics claim: in the quasipinning regime real systems
occupy, occupations underdetermine two-body physics by a bounded-below
amount. (Note the full-fiber width 0.25 at eps=0 exceeds the single
v_B family's range -- the complete fiber, stabilizer orbits and all
components included, is bigger than any one family; the family gave
the lower bound, the attack measures the whole thing.)

THEOREM TO WRITE (now with its numbers): valley width is continuous in
lambda, locally minimized at the vertex, bounded below near it by the
vertex width. Proof route: incidence surjectivity (linear algebra) +
open polygon conditions + the measured monotonicity. PHYSICS QUEUE:
(1) Berry phase around the fiber circle (closed-form afternoon);
(2) exact RDMFT functional table F[gamma](H) along the valley for a
    family of two-body H -- first exact functional data at a GPC
    vertex, addressed to the ensemble-DFT/exchange-force program;
(3) preparation-and-measurement sketch for a programmable simulator
    (6-8 determinant states, shallow circuits; observables: pair
    occupations for the radius, T-odd pair coherences for the
    orientation, the discontinuous selection jump as the signature);
(4) valley-width census over other deforming vertices (same script).


## v89/v103 NUMERICALLY ATTAINED (verified here) -- but NOT certified; census
unchanged at 797/799 (2026-07)

VERIFIED (ran scripts/contraction_attack.py independently). The ansatz-free
contraction-map attack -- minimize ||gamma(Psi) - gamma_0||^2 over the FULL 120-dim
Lambda^3 C^10 tensor with an analytic Wirtinger gradient, no support/sparsity/
rationality restriction -- attains BOTH holdouts to machine precision from random
starts: v89 residual 4.3e-17 (complex) and 7.6e-17 (REAL), v103 residual 2.2e-16
(complex) and 5.0e-16 (REAL), eigenvalues exactly (15,15,6^8) and (18^4,5^6). So
BOTH v89 and v103 numerically admit REAL states -- which CORRECTS an earlier hedge
here (my real-diagonal attain floored at ~5e-6; that was a weak numeric gradient,
not a real obstruction: the analytic-gradient version nails it). The analytic
gradient is the whole difference from every prior stalled search.

WHAT THIS IS, PRECISELY. Numerical ATTAINABILITY was never the open question --
these are vertices, so states exist and my earlier attain already hit ~1e-14. The
open problem is a CERTIFIED CLOSED FORM (an exact state passing the char-poly gate),
and the contraction solutions do NOT provide one: they are DENSE generic fiber
points (support 100-120 of 120), exactly the "dense, no exploitable symmetry" state
my analysis predicted (design/split/cyclic/neighbor/Pauli all fail; see below). So
this REINFORCES the near-theorem rather than overturning it: the state exists and is
findable numerically, but it is generic in the degenerate subspace and yields no
closed form. Census stays 797/799; do NOT read "numerically attained" as "solved".

CERTIFICATION ROUTES (all still open; a distinguished, recognizable representative
is required): (a) sparsity-seeking reoptimization on the solution manifold toward a
small support, then the standard exactifier -- but the manifold may carry no sparse
point (the numeric support is full); (b) the real-amplitude manifold (now known
nonempty) with PSLQ/minpoly recovery of algebraic weights -- viable only if the
weights are low-degree algebraic, unknown; (c) the equivariant route (a subgroup
forcing gamma = scalar+projector by Schur) -- ALREADY SHOWN TOO RIGID for v89 (the
order-128 Pauli construction gives only 2-dim covariant spaces, residual 0.34), so
(c) is effectively closed. Net: attainability + reality settled numerically;
certification remains the genuine open problem, and the honest odds of a clean
closed form are low given the density.

## RANK-11 FORK: contraction attack floors on all four open (3,11) candidates
(verified here) -- strong evidence for the refutation branch, NOT proof (2026-07)

VERIFIED (ran the same scripts/contraction_attack.py on the (3,11) settlement
candidates). The attack that reaches ~1e-16 on true vertices FLOORS 11-14 orders
higher on every open candidate. Controls (TRUE-VERTEX): idx0 7.1e-17, idx7 4.1e-16,
idx14 2.8e-16, idx20 1.6e-16 -- machine precision, the solver has no trouble with
attainable (3,11) spectra. Open candidates (docs/bracket_3_11_settlement.json,
verdict OPEN): idx23 (den 7, spec 6,6,1^9) residual 9.2e-2; idx34 (den 17, 11^4,1^7)
2.5e-2; idx26 (den 27, 21,21,11^3,1^6) 9.9e-3; idx44 (den 12, 6^5,1^6) 2.2e-3. My
numbers reproduce the parallel run's to two figures. Targeting diag(spectrum) is
WLOG (orbital rotation), so this tests attainability of the SPECTRUM directly; the
floor is the squared distance from the target to the nearest pure-N-representable
1-RDM, so a positive floor means the spectrum is OUTSIDE the (3,11) polytope.

READING (evidence, discipline preserved): all four candidate spectra are very likely
outside the true (3,11) polytope, hence the claimed level-5 second-kind inequality
of Refs. AK2008/Klyachko2009 is very likely TRUE, and the bracket settles at 19
true vertices + 27 refuted + 4 (now evidentially-refuted) = the falsification branch
of the paper's "decisive fork" is heavily DISFAVORED. The floor magnitudes rank the
candidates by distance-to-attainable (idx23 farthest, idx44 nearest), useful for a
Stage-1 facet hunt. CAVEATS (do not overstate): nonconvergence is EVIDENCE, not
proof -- a thin missed optimization basin is disfavored by the true-vertex control
but not excluded; rigorous closure still requires the Schubert-coefficient proof of
the level-5 series (which doubles as the Stage-1 constraint-generator verifier). So
in the paper this shifts the fork's evidential weight toward refutation but must NOT
be written as "resolved" -- it is a strong numerical indicator, and the exact proof
is the open task. (Recorded in RESEARCH.md; not injected into main.md, which
correctly keeps the fork open pending the Schubert computation.)

## SKEPTICAL AUDIT + PRE-REGISTRATION for the rank-10 closure (2026-07)

Selection-bias audit of the fiber era + anti-overfitting protocol (adopted from the
parallel instance; sound practice, recorded verbatim in spirit):
1. Every fiber law was mined from the 44 kdim-1 single-touched-class families of the
   SOLVED corpus -- states our ansatz family finds. "44/44 conics" is a THEOREM for
   single-2-term-class, +-1-loop families PLUS an empirical fact about the solver's
   output distribution, not a law of nature. A v89/v103 state with multi-class or
   3-term structure sits outside every fiber law's proven scope while violating
   nothing. The endpoint theorem and wall-interiority are proved only in the
   single-class scope.
2. Mined-law mortality to date is about one in three (class-count, toric-orbit,
   elliptic-fiber-as-the-rule -- all retracted). Presume survivors carry similar
   out-of-sample risk until tested.
3. Weapon vs law: the ENGINE features (symmetry reduction, kernel quotienting,
   rigid-class arithmetic, continuous/wall-first search, and now the contraction
   map) are sample-independent and bear on the holdouts; the LAWS describe solved
   states post hoc and crack nothing. The x1.32/denominator scaling forecast is a
   6-point fit -- an estimate, not a promise.

PRE-REGISTERED PREDICTIONS (committed BEFORE any CERTIFIED v89/v103 state exists;
score PASS/FAIL/NOT-TESTED against the certified states only, before writing any
narrative; a FAIL is itself a publishable finding):
- P1 every gauge-invariant holonomy of the v89/v103 states generates an exponent-2
  abelian extension of Q (Conjecture 2).
- P2 each state uses at most 2 exchange channels (Conjecture 3).
- P3 rational-weight states have state-den = spectrum-den (26, 34); irrational
  (continuous-search) states invoke the scope clause and do NOT test P3.
- P4 loopy supports with no touched 1-term class deform (fiber dim = kernel dim).
- P5 the states break their spectrum symmetry (Aut != G).
- P6 if the family is single-2-term-class +-1-loop, its CM curve is a conic with
  both walls real-type.
(Caveat consistent with the above: the contraction solutions are DENSE, so a
certified state -- if one is ever extracted -- may itself be dense/multi-class and
land in several predictions' scope-exclusion clauses rather than testing them.)

## THE FIBER IS THE 1-RDM -> 2-RDM INFORMATION GAP (verified at v_B, 2026-07)

Framing (parallel instance's, adopted; candidate central claim for a second paper,
NOT injected into the current main.md): even with the one-body spectrum EXACTLY
pinned at a vertex -- the sharpest boundary the polytope offers -- the many-body
state is not fixed. Natural occupation numbers are incomplete coordinates at maximal
pinning; the fiber is precisely the information the 1-RDM cannot see but the 2-RDM
can. This is STRUCTURAL and independent of the empirical quasipinning-deflation
debate: pinning fixes the FORM of the state, the fiber is what lives inside the form,
generically a continuum.

VERIFIED HERE (v_B, exact 2-RDM built by second quantization, scratch vb_2rdm):
- the RADIAL fiber coordinate t is read by PAIR OCCUPATIONS (diagonal of Gamma^2) --
  already established (8 pair occupations vary linearly along the family; in the
  paper).
- the two conjugate sheets +Theta and -Theta (the state and its time-reversed image)
  are INDISTINGUISHABLE in every T-even two-body datum: the diagonal-Gamma^2
  (pair-occupation) difference is EXACTLY 0.0 and the real-part difference of the
  full 2-RDM is EXACTLY 0.0.
- they differ ENTIRELY in the T-ODD (imaginary) part of the 2-RDM: the whole
  conjugate-sheet 2-RDM distance ||Gamma^2(+) - Gamma^2(-)||_F = 1.40 is imaginary,
  max |Im Gamma^2| ~ 0.12. (The parallel run quotes 1.30 and ~0.24 for these; the
  gap is a norm/normalization convention, the STRUCTURE -- all-difference-is-T-odd --
  is exact and reproduced.)
So the fiber orientation is a purely time-reversal-odd, two-body observable: the
1-RDM is blind to the whole fiber; pair occupations resolve where on the arc; T-odd
pair coherences resolve which sheet. Clean, correct, and paper-2 material.

MISSING THEOREM (named honestly, unproven): quasi-fiber stability -- does the moduli
continuum persist as a near-degenerate valley for spectra NEAR the vertex (the
quasipinning regime, where real systems live)? Plausible by perturbing the exact
families, but NOT proved. Until then the physicality claim is "exactly at the
boundary"; with it, it reaches real systems. Do not state it as established.

QUASI-FIBER STABILITY NUMERICALLY DEMONSTRATED (verified here; the theorem itself
is still unproven). Extremized the pair occupation <n_2 n_4> over states with EXACT
1-RDM diag(lambda), for lambda walking from v_B into the interior along
(1-eps) v_B + eps*uniform (constraint-weighted contraction attack, residual < 1e-4):
- eps=0.00 (at v_B): <n_2 n_4> in [0.216, 0.470], WIDTH 0.254
- eps=0.05:          <n_2 n_4> in [0.200, 0.494], WIDTH 0.294
My numbers reproduce the parallel run's (0.252 -> 0.290) to ~1%, and the lower end
is exactly the Frechet bound occ_2 + occ_4 - 1 = 2*(14/23) - 1 = 0.2174 (analytic,
matches). CONTENT confirmed: the occupation-blindness is NOT a knife-edge -- states
with an identical FULL 1-RDM (gamma = diag(lambda), the constrained-search domain of
RDMFT, NOT merely the spectrum) differ in a pair occupation by ~0.25 AT the vertex,
and the width GROWS off-vertex, so along this ray maximal pinning is where the
ambiguity is SMALLEST and even there it is bounded below by ~0.25. This is the
data-backed, deflation-proof form of the physics claim: at fixed 1-RDM, two-body
physics is underdetermined by a bounded-below amount.

SCRUTINY REPAIRS (expert pass, binding -- corrects wording and scope of the above):
- WORDING (important): the measurement fixes the full 1-RDM gamma = diag(lambda),
  NOT the spectrum alone. Fixing only the spectrum is trivially large by basis
  rotation and must never be the stated claim; the correct (stronger, RDMFT-relevant)
  statement is "the 1-RDM underdetermines pair correlations". My computation targeted
  diag(lambda) throughout, so the NUMBERS are for the fixed-1-RDM reading and stand;
  only the "spectrum" phrasing was loose and is fixed here.
- The widths are certified LOWER bounds (inner approximation by local optimization);
  the true achievable ranges may be LARGER, and interval-ness of the achievable set
  is UNPROVEN (fiber connectivity unknown) -- so "range [a,b]" means "endpoints a,b
  are achieved values", not "every value between is achieved" (sample interior
  targets or prove connectivity before asserting a full interval).
- "maximal pinning locally minimizes the ambiguity" rests on ONE ray + three
  epsilons: it is observed-along-this-ray + conjecture, not established; the eps=0.10
  point failed for an unattributed reason (optimizer budget vs feasibility) and is
  excluded, not glossed.
- DISCIPLINE: numerically DEMONSTRATED, not a proved theorem -- the stability theorem
  (valley width continuous in lambda, locally minimized at the vertex, bounded below)
  still needs its proof (route: incidence surjectivity + open-polygon conditions +
  the measured monotonicity).
(Caveat on magnitude: the full-1-RDM-fiber width 0.25 exceeds the single ledger
family's <n_2 n_4> range -- the whole fiber, all components and stabilizer orbits, is
larger than any one deforming family; the family is a lower bound, the ansatz-free
attack measures the whole thing.)

SELECTION-PRINCIPLE GAP (patch-6, partially checked here): the endpoint theorem is
PER-FAMILY, but the physical minimization runs over the FULL fiber (larger than any
family), so "density-density interactions select real wall states" needs the fiber
EXTREMES to admit real representatives. The parallel run reports complex and
real-restricted extremization of <n2 n4> at v_B agree to ~1e-8. NOT INDEPENDENTLY
CONFIRMED HERE: my quick penalty-extremizer is optimizer-limited (real vs complex
differ by ~1e-2, and it even returned real-max > complex-max, impossible at the true
extreme -- so it is not resolving the extremes to the claimed 1e-8). The mechanism
is plausible (c -> conj(c) preserves gamma and all pair occupations, so extreme
points plausibly have time-reversal-symmetric representatives), but treat "extremes
are real" as parallel-reported + conjecture until a reliable constrained solver or a
proof confirms it; scope any paper claim to "per-family theorem + v_B numerics".

2-RDM EMBEDDING, proof upgrade (patch-6, consistent with my check): the
conjugate-sheet T-odd coherences are proportional to sin Theta -- nonzero on the
open interior, vanishing exactly at the |cos|=1 walls where the sheets glue. This
matches my earlier finding that the +/- sheets differ ENTIRELY in the imaginary
2-RDM part (real-part difference exactly 0). Injectivity of the fiber circle then
follows from monotone t-resolution (pair occupations) + interior sheet separation
(T-odd coherences). Upgrades the spot-check toward a one-line proof.

RANK-11 CERTIFICATE PATH (patch-6 recommendation; NOT run -- no cvxpy/SDP solver in
this environment): the contraction floors are UNCONSTRAINED local minima (evidence,
not proof). A Lasserre/SOS (moment-relaxation) dual of min ||gamma(psi)-gamma_0||^2
yields a CERTIFIED positive lower bound = rigorous non-attainability, independent of
and faster than the Schubert route. Run the degree-2 relaxation on cand 44 (smallest
floor, nearest to attainable) first: if even the nearest candidate certifies
infeasible, the level-5 claim becomes theorem-adjacent without new mathematics. This
is the recommended next concrete step for the rank-11 closure and needs an SDP-capable
environment.

## RANK-11 CERTIFICATE PATH: SDP RUN, method validated, scaling wall found (2026-07)

Executed the recommended SDP path (Clarabel/cvxpy installed in the venv;
scripts/sos_nonattain.py). Three concrete results, discipline preserved.

1. CONE FIX (the naive relaxation is wrong). f = sum_ij |gamma_ij - g0 delta|^2 is
   U(1)-gauge invariant, so it is a sum of squared magnitudes of the SESQUILINEAR
   forms gamma_ij - g0 delta, NOT a holomorphic SOS. A holomorphic {1,c,c^2}-block
   SOS cannot even represent it (fails to certify the trivial delta=0). The correct
   cone squares the CHARGE-0 space B = {1} U {Re(c_a cbar_b), Im(c_a cbar_b)} (size
   1 + nd^2, nd = C(d,N)): f - delta = B^T Q B + mu(||c||^2 - 1), Q >> 0, mu a free
   charge-0 degree-<=2 multiplier. This is the level-2 relaxation on the sphere.

2. VALIDATED on (2,4), whose pure spectra are the degenerate pairs (a,a,1-a,1-a).
   EXTER (0.8,0.6,0.4,0.2) certifies delta = 0.04, EXACTLY the squared distance to
   the nearest attainable (0.7,0.7,0.3,0.3); ATTAIN (0.7,0.7,0.3,0.3) gives
   delta ~ 0. Two independent solvers (cvxpy and a hand-assembled native Clarabel
   program) agree; the certificate reconstructs f to ~1e-14. So a positive delta IS
   a rigorous non-attainability certificate and the pipeline is sound on a case with
   a known answer.

3. G1 PASSES un-reduced at (3,6) with the right solver. The single PSD block has
   side 1 + nd^2: 401 for (3,6), 27226 for (3,11). Clarabel's dense-IPM
   Nesterov-Todd scaling is O(npsd^2) = 80601^2 * 8 B ~ 52 GB and OOMs at (3,6);
   SCS (ADMM, eigendecomposition PSD projection) is O(nz^2) and clears it. With
   SCS the G1 gates both pass: EXTER diag(0.9,0.8,0.7,0.35,0.2,0.05) certifies
   delta = 0.00416 > 0 (non-attainable, violates the (3,6) pairing), the attainable
   CONTROL diag(1,1,1,0,0,0) returns delta ~ 1.5e-6 ~ 0 (feasible); certificate
   residual ~1e-8. So the pipeline is validated on the actual (3,6) target, not
   only the (2,4) toy.

4. (3,11) STILL NEEDS SYMMETRY. SCS is O(nz^2) memory but O(nz^3) per PSD
   projection: side 27226 is ~2e13 flops/iteration -- time-infeasible. cand 44
   g0 = diag(6^5,1^6)/12 has spectrum stabilizer S5 x S6; block-diagonalizing Q
   under the fermionic determinant rep collapses side 27226 into irrep blocks.
   Burnside counts (computed): dim End_G(V) = sum m_lambda^2 = 18 (tiny
   multiplicities in the 165-determinant rep), and dim End_G(Herm V) = sum
   mult_nu(Herm V)^2 = 28884, so the largest reduced PSD block is <= ~170 -- easily
   in SCS range. That symmetry-adapted assemble-and-solve (scripts/sos_symmetry.py)
   is the next step; the validated un-reduced assembler and G1 self-test ship in
   scripts/sos_nonattain.py.

TARGET DATUM (for the reduced solve to reproduce). By orbital-rotation invariance of
the pure-1-RDM set, min_Psi ||gamma(Psi) - diag(g0)||^2 = squared Euclidean distance
from g0's spectrum to the (N,d) moment polytope (Hoffman-Wielandt: the nearest pure
1-RDM aligns eigenbases). So a positive delta certifies the spectrum is OUTSIDE the
polytope, and delta equals dist^2 to it. For cand 44 the contraction floor is
~2.155e-3 (~1/464 to 7 figures, real and complex agree); that is the number the
symmetry-reduced certificate must lower-bound. NOTE this reduction also means any
single violated GPC facet inequality gives delta >= ((a.g0 - b)/||a||)^2 by exact
rational arithmetic -- but the (3,11) level-5 facet is precisely the open Schubert
object, which is why the ansatz-free SOS route (no facet list needed) is the target.

## RANK-11 CERTIFICATE PATH: cand 44 symmetry-reduction analysis (2026-07)

Groundwork for the S5 x S6 reduced solve (G3). All numbers below are computed
(Burnside character sums + numerical isotypic split of the 165-determinant rep),
not estimated.

DECOMPOSITION of V = signed permutation rep on the 165 = C(11,3) determinants under
G = S5 (modes 0-4) x S6 (modes 5-10). V splits into 9 distinct S5 x S6 irreps: three
with multiplicity 2 and six with multiplicity 1 (so sum m_lambda = 12,
sum m_lambda^2 = dim End_G(V) = 18). By (mult, dim): (2,4), (2,6), (2,20), (1,5),
(1,10) x three distinct dim-10 irreps, (1,30), (1,40); check sum m*d = 165.

The SOS Gram Q is indexed by the charge-0 basis {1} U {c_a cbar_b}, transforming as
trivial (+) V (x) Vbar. Its G-invariant part block-diagonalizes with one PSD block
per irrep of size = mult_nu(Herm V); dim End_G(Herm V) = sum mult_nu^2 = 28884, so
the largest block is <= ~170 -- trivially in SCS range (vs the un-reduced side
27226 that is time-infeasible).

KEY COMPUTATIONAL UNLOCK (makes the matrix-free reduction cheap): because G is a
DIRECT PRODUCT and rho_V(sigma,tau) = rho_V(sigma) rho_V(tau) with the character
factorizing chi_nu(sigma,tau) = chi_a(sigma) chi_b(tau), the isotypic projector
factorizes P_nu = P_a^{S5} . P_b^{S6}, each a sum over only 120 (resp. 720) group
elements, not |G| = 86400. Applying rho_V(g) (x) conj(rho_V(g)) to a 165x165 matrix
X is X with rows/cols signed-permuted -- O(nd^2). So building every isotypic
projector costs ~840 * nd^2, seconds not hours. This is the lever that turns cand 44
from intractable into a bounded computation.

G2 DONE -- the reduction is BUILT and VALIDATED (scripts/sos_symmetry.py). The
numerical-Wedderburn reduction (character-free: eigendecompose a random G-invariant
operator for the isotypic split, link equal-irrep copies via a second invariant,
parametrize the invariant Gram by signed index-pair orbits with y_o = y_{o^T} for
symmetry) block-diagonalizes the charge-0 Gram and reproduces the un-reduced answer.
Gate: on the S2 x S4-symmetric EXTERIOR (3,6) point diag(0.9,0.9,0.3,0.3,0.3,0.3)
the reduced delta = 0.0882338 matches the un-reduced delta = 0.0882353 to 1.5e-6,
with 10 blocks of sizes up to 26 (vs the un-reduced side-401 block). The
group-action rho_B on the SOS basis is checked to be a genuine orthogonal rep
before use.

G3 SCALE-UP: primitives BUILT and VALIDATED, one bottleneck remains
(scripts/sos_symmetry_scaled.py). The scalable swaps needed for cand 44 are done and
each is checked:
  * FACTORIZED REYNOLDS. For G = A x B on DISJOINT modes (A,B commute),
    Reynolds_G = Reynolds_A o Reynolds_B -- |A|+|B| = 840 conjugations, not 86400.
    Verified to reproduce the dense full-|G| average to 9e-16, so the
    block-diagonalization scales.
  * GENERATOR-BFS signed orbits (over adjacent transpositions), and exact S_n
    characters (scripts/sn_characters.py, Murnaghan-Nakayama, verified by
    dim-sum-of-squares and orthogonality).
  * SUPPORT (term-sparsity) parametrized by hop distance N - |a & b|. FINDING: the
    certificate needs TWO-HOP support. On the (3,6)/S2xS4 exterior point the scalable
    pipeline reproduces the un-reduced delta EXACTLY at hop<=2 (0.08822 vs 0.08823)
    but COLLAPSES to ~0 at hop<=1 -- one-hop pairs alone cannot carry the
    Positivstellensatz. So the whole scalable pipeline (factorized Reynolds +
    generator-BFS + hop-restricted reduced SDP) is validated end-to-end at (3,6).
  * WHY cand 44 is not yet solved: at hop<=2 the (3,11) basis is nz = 17986 (vs full
    27226). The current assembly enumerates all nz^2 index pairs (BFS + coefficient
    matching) and forms dense constraint matrices (nk x norb) -- ~3.2e8 pairs and a
    multi-GB constraint matrix, which does not fit. The remaining piece is a
    GROUP-SUM assembly: compute each orbit's coefficient contribution as the Reynolds
    of the single monomial product B_{P0} B_{Q0} (840 elements) instead of
    enumerating orbit members, and keep the constraint map sparse. That replaces the
    O(nz^2) enumeration with O(norb * 840 * small) and is the one component between
    here and a numerical cand-44 delta.

Then RIGORIZE per repo discipline: round the dual to rational and verify
f - delta = B^T Q B + mu(||c||^2 - 1) as an EXACT polynomial identity with each
reduced block exactly PSD (block-diagonal => exact LDL is cheap). That exact identity
is the theorem and a false-positive-proof gate independent of any reduction-code bug.
Method proven (G1 + G2), primitives validated; the group-sum assembly + rational
endgame remain.

## SOLVER UPGRADES: moduli/symmetry-informed search (2026-07)

Four of the five corpus-mined solver upgrades implemented. All are exact
reformulations or theorem-grade constraints, never support filters; each ships
with tests, and each that encodes a law ships with a corpus-containment proof
that no certified state is excluded.

- F1 SYMMETRY-REDUCED SKELETON SEARCH (states._skeleton_model): lex-leader orbit
  reduction within equal-nv mode classes, with block modes colored as singletons
  so the hop constraints stay invariant. Keeps the lex-minimum representative of
  every orbit, so no orbit is emptied (sound; asymmetric support_filters skip the
  transposition they break). Wired behind states._SYMMETRY_BREAK (default on).
  Full suite runs faster with it on (about 27s vs 34s); v96 canary
  min_block_count = 1 in 0.4s, production solve OK in 1.8s.
- F3 RIGID ONE-HOP-CLASS CONSTRAINT (states.one_hop_classes,
  rigid_class_products): a class with a single hop pair forces k_A k_B = x2 (the
  block eigenvalue identity), a sound search prune. test_rigid_class proves the
  identity over the whole corpus, so the prune never drops a certified state.
- F4 DENOMINATOR-GRID TIERING (states.MAX_DENOMINATOR_TIER): the ratio
  state_den/spectrum_den is 1 or 2 corpus-wide and 2 only for DESIGN-REAL, so the
  interference/design-int paths stay m=1 and DESIGN-REAL is the sole m=2 consumer.
  The design-real record carries denominator_tier; test_denominator_tier proves
  the ratio law over the corpus.
- F5 EXPONENT-2 LIVE TEST (exactify.is_exponent_two, conjecture2_scan): every
  certified state's gauge-invariant loop holonomies are classified by the Galois
  group of minpoly(2 cos theta); exactify records holonomy_exponent_two and raises
  CONJECTURE2-CANDIDATE-COUNTEREXAMPLE on a definitive non-2-elementary holonomy
  without ever rejecting the verify_exact-certified state. The check is on the
  gauge invariant, not amplitude phases (trisections are expected). v_B, v113,
  v17 confirmed exponent-2; a C3 angle (cos 2pi/7) is correctly flagged.
- F2 KERNEL QUOTIENTING (states._KERNEL_QUOTIENT, _incidence_kernel_dim): on a
  fixed support the integer weight solutions form a coset of the incidence
  kernel; free-kernel translates are phase-equivalent (fiber-dimension law), so
  solve_vertex_exact_first phase-solves one representative per (support,
  rigid-class products) coset. Keying on the rigid products preserves the F3
  coupling (rigidly-pinned distinctions, which cut the kernel, keep separate
  keys), so nothing sound is skipped; verify_exact gates. Records
  fiber_kernel_dim. Validated: v40 certifies with quotient on/off.
- DEGENERATE-BLOCK ANSATZ (gpc_census.degenerate): the 2x2/clique families mix
  DISTINCT eigenvalue classes and miss repeated-eigenvalue blocks (v89's
  eightfold 6). The new generator emits degenerate k-block degree vectors (a
  size-k eigen sub-multiset with a repeat, Schur-Horn-majorized integer
  diagonal) and solves each with the shared skeleton enumeration (F1/F2 apply)
  plus attain (eigenvalue matching is basis-free, so a block-diagonal 1-RDM is
  fine) gated by exactify. crack_vertex_exact falls back to it when the
  block/clique family exhausts. Soundness rests on the exact gate.

Acceptance: G1 full suite green (75 tests). G3 v96 canary passes. G2 regression
interference vertices still certify with symmetry on, suite faster.

G4 OUTCOME (upgraded sweep, max_clique = 4, max_cliques = 0, max_card above the
denominator, 7200 s budget each): both open vertices EXHAUSTED the block/clique
ansatz family with no phase-solvable skeleton, and both finished well under the
budget, so this is an ansatz-family EXHAUST, not a timeout:
- (3,10) v89 = (15,15,6,6,6,6,6,6,6,6)/26: FAIL in 705 s (exhausted).
- (3,10) v103 = (18,18,18,18,5,5,5,5,5,5)/34: FAIL in 1278 s (exhausted).
The symmetry reduction (F1) is what let the full block search COMPLETE within
budget where it previously timed out, so the negative is now informative: the
2x2-block + clique family (through max_clique 4) provably does not reach v89 or
v103. This is NOT a proof of non-attainability. v89's eightfold-6 degeneracy
points to the degenerate-kxk ansatz (scripts/kxk_degen.py), which is outside
solve_vertex_exact_first's family; v103 (den 34, two classes) likewise needs a
richer or degenerate block. The census stays at 797/799; the ledger is untouched.

FOLLOW-UP (both next levers now implemented): F2 kernel quotienting shipped, and
the degenerate-block ansatz was folded into the family and wired into
crack_vertex_exact as a fallback.

DEGENERATE-FALLBACK RUN OUTCOME (v89 3600 s + v103 3600 s degenerate budget): no
solve, and crucially the degenerate pass TIMED OUT for both (v89 3608 s, v103
3746 s), it did NOT exhaust. So this is inconclusive, not a family-does-not-reach
result: brute degenerate enumeration (per-config support enumeration times a
per-support attain) is too slow to finish even with F1 symmetry and F2 kernel
quotienting. Census stays 797/799; ledger untouched.

SYMMETRY REDUCTION OF v89 -- a provable obstruction and the exact next ansatz
(2026-07; scripts/v89_symmetry.py). v89 = (15,15,6^8)/26 has symmetry S_2 x S_8
(head = the two 15-orbitals, octet = the eight 6-orbitals). The 3-subset
determinants form THREE orbits ({3 octet} 56, {1 head 2 octet} 56, {2 head 1 octet}
8), so a fully symmetric state has just THREE amplitudes (x,y,z) and its 1-RDM
block-diagonalizes: antisym-head eigenvalue a-b, octet standard rep c-d (mult 7),
and a 2x2 [[a+b,4e],[4e,c+7d]]. The five entries are exact quadratic forms:
a=28y^2+8z^2, b=28y^2, c=21x^2+14y^2+z^2, d=15x^2+12y^2+z^2, e=7y(3x-z).
RESULT: the fully-symmetric ALL-POSITIVE (trivial-rep = "design") ansatz PROVABLY
cannot attain v89. The octet off-diagonal d is large and same-sign as c, so the
mult-7 octet eigenvalue c-d = 6x^2+2y^2 is STRICTLY below the required (6/26) n2
(gap = (90x^2+142y^2+24z^2)/13, positive-definite), yet it must equal 6/26 (part of
the eightfold 6). So no symmetric design attains v89: the vertex demands
INTERFERENCE that REDUCES the octet off-diagonal, i.e. intra-orbit sign/phase
structure transforming under a NONTRIVIAL S_8 representation (a per-orbit phase or
an all-positive signing cannot do it -- the reduction must live inside the octet
3-subset orbit). This is a constructive re-derivation of v89's interference
character (sharper than the MILP) and pinpoints the search: find the S_8-rep-
theoretic signing of the 56 octet 3-subsets that lowers rho_{octet} to the target.
That is a small structured combinatorial-design problem (56 signs modulo an
S_8-orbit action), NOT a 120-determinant brute enumeration -- the first attack that
uses v89's own symmetry. Still OPEN; this is progress and a precise next ansatz,
not a solve.

FULL-STEAM FOLLOW-UP (two more ansatze eliminated; the target structure pinned).
The exact 1-RDM target is now explicit: eight eigenvalues equal 6/26, so
rho - (6/26) I has rank EXACTLY 2, i.e. rho = (6/26) I + (9/26) P_2 with P_2 the
rank-2 projector onto the two 15/26 natural orbitals. This is a maximally-mixed
(scalar) 1-RDM on an 8-dim subspace plus a rank-2 bump -- a very rigid target.
(a) A magnitude-S_8-symmetric octet block is gauge-equivalent to the all-positive
design (rho_{ij} = lambda eps_i eps_j with eps a sign vector is the gauge
a_i -> eps_i a_i), so it inherits the design obstruction: reducing the octet
off-diagonal REQUIRES non-uniform octet off-diagonal magnitudes, i.e. the state
BREAKS the octet magnitude-symmetry. (b) The natural symmetry-broken ansatz,
splitting the octet into two S_4 blocks (S_2 x S_4 x S_4, nine determinant orbits,
signs free), does NOT attain v89 either: an exact Powell search over the nine
orbit amplitudes floors at residual 0.14 (eigenvalues spread 16.7,14.7,10.4,...,0),
nowhere near the eightfold 6. So neither the fully symmetric nor the 4+4-split
symmetric sparse ansatz reaches the rank-2 target. CONCLUSION: v89's extremal state
is not captured by a small symmetric-orbit ansatz -- the rank-2 = (6/26)I + bump
structure forces the state to be generic in the 8-dim degenerate subspace, which is
exactly what defeats sparse symmetric constructions and the rational-grid search.
This is honest negative progress: it rules out the two natural structured ansatze
and identifies the precise obstruction (rank-2 target vs sparse symmetry). A crack
likely needs either the cluster-mutation route (a certified neighbor plus an exact
mutation, see below) or an explicit rank-2 1-RDM inverse construction, not another
symmetric ansatz. v103 = (18^4,5^6)/34 is the analogous rank-2 case (rho -
(5/34)I has rank 4: four 18/34 orbitals) and is expected to be harder still.

NEIGHBOR-MUTATION ROUTE ALSO CLOSES (checked against shipped vertex data). v89's
nearest CERTIFIED neighbor sharing its (15,15,6...) architecture is v81 =
(15,15,6^5,4^3)/24 (INTERFERENCE, 9-determinant support). Its orbital occupations
are [15,15,6,6,6,6,5,4,4,5]/24: it makes FOUR of its 6's DIAGONALLY (orbitals 2-5
are occupation-6 outright) and only ONE by an interference 2x2 block (orbitals
(6,9): (5,5) -> (6,4)). That sparse "mostly-diagonal-plus-one-block" mechanism is
incompatible with v89's requirement that ALL EIGHT generic orbitals be degenerate
at 6/26 via the maximally-mixed rank-2 block -- so mutating/extending v81's support
cannot reach v89 (different construction in kind, not degree). Every concrete route
tried -- symmetric design (provably out), 4+4-split symmetric (residual 0.14),
Z_8 cyclic complex (search did not converge), and neighbor mutation from the
nearest certified vertex (mechanism-incompatible) -- converges on the same reading:
v89's extremal state is dense/generic in its 8-dim degenerate subspace and is not
reachable by any sparse, symmetric, or neighbor-derived construction available here.
v89/v103 stay OPEN; what is now firmly established is WHY (rank-2/rank-4
maximally-mixed target), which is itself a paper-worthy characterization of why
exactly two rank-10 vertices resist closed forms. A genuine crack needs dedicated
(non-interactive) machinery: an exact rank-2 pure-state N-representability inverse
for rho = (6/26)I + (9/26)P_2, or a real cluster-algebra mutation engine -- not more
ansatz search.

RANK-2 INVERSE ATTEMPT via the group-theoretically CORRECT construction (a sharp
negative result). To force the generic 8-block scalar one needs a NONABELIAN group
with an IRREDUCIBLE 8-dim rep (abelian groups have only 1-dim irreps -- which is why
the Z_8 cyclic attempt was doomed). The 3-qubit Pauli/Heisenberg group (order 128)
is exactly that: generic orbitals ~ F_2^3, X_i/Z_i generators, irreducible on C^8.
A 3-fermion state covariant under a linear character chi of this group has, by
Schur, rho_generic proportional to I AND rho_{head,generic} = 0 (an intertwiner
trivial -> irreducible vanishes), so rho is forced block-diagonal with a scalar
octet and the two 15/26's in the head 2-block -- exactly the rank-2 target with
P_2 = head. Computed: 28 of the 64 characters carry a nonzero covariant 3-fermion
space, each exactly 2-DIMENSIONAL. RESULT: none reaches v89 -- scanning the
one-parameter family in each, the best residual over all 28 is 0.34. The covariant
states have the right STRUCTURE (scalar octet by Schur) but their rigid 2-dim
freedom cannot tune the octet scalar and head block to the required (6/26, 15/26).
CONSEQUENCE (a real characterization): v89's extremal state is NOT of the natural
"two pure head orbitals + symmetric scalar octet" form; its rank-2 plane P_2 must be
a MIXED generic 2-plane spanning head and generic orbitals together, with a
NON-scalar octet block. Since the smallest group that could force a scalar octet is
already the full order-128 Pauli group and it is too rigid, NO symmetry-aligned
construction produces v89: the state is genuinely generic in its degenerate
subspace, with no aligning symmetry to exploit. This closes the natural
rank-2-inverse route and upgrades the open status to a near-theorem -- v89/v103
resist closed forms because their extremal states carry NO exploitable symmetry.
Any future crack must be a non-symmetric exact construction (a direct real-algebraic
solve of the rank-2 N-representability system over a FREE mixed P_2), a dedicated
computation rather than an ansatz.

The right next move is structural, not more brute force. Under the moment-map
lens the extremal supports are candidate Grassmannian cluster seeds and the
one-hop graph is the hypersimplex 1-skeleton, so v89/v103 should be reachable by
CLUSTER MUTATION from a certified neighbor rather than by enumerating supports.
Cheap falsifiable pre-checks (both reuse shipped data): (1) are the 154
interference supports weakly separated collections (Leclerc-Zelevinsky /
Oh-Postnikov-Speyer)? if so the fibers are cluster varieties and mutation search
applies; (2) does the Frobenius-Schur indicator of each state under its Aut group
track the design/signed/complex trichotomy? if so reality is a Schur-index
computation, not a numerical overlap. These are the recommended next probes.

WEAKLY-SEPARATED PROBE RESULT (run over the corpus): the literal "supports are
cluster seeds" hypothesis is FALSIFIED. 0 of 154 interference supports and only
50 of 643 design supports are weakly separated collections (Leclerc-Zelevinsky).
But a clean one-directional law is verified: weakly-separated (a Gr(N,d) cluster
seed) implies DESIGN-attainable, and every interference support contains a
crossing pair (154/154). So interference lives exactly in the non-cluster
(crossing) regime. Consequence: single-seed mutation does not directly apply to
v89/v103 (their supports are crossing). The refinement "crossing count = channel
count" is ALSO falsified: crossings range 3 to 24 across the 154 supports,
uncorrelated with the <=2 channel bound. So the only surviving verified fact from
this probe is the one-directional law weakly-separated implies design; the
cluster reading of interference supports (single seed, or crossing = channel) is
refuted. Recorded as a falsified lead so it is not re-proposed.

FROBENIUS-SCHUR PROBE RESULT (run over the complex corpus, implementation
validated on design states): also does NOT unlock the reality question. Every
complex certified state, v_B included, has NO finite antiunitary permutation
symmetry (49/49 tested; real states trivially have one via complex conjugation
K with pi = identity). So the real/complex split does coincide with
antiunitary-symmetric/not, but there is NO Frobenius-Schur = -1 (quaternionic /
Kramers) case among finite symmetries, so no exact Kramers proof of complexity
is available and the reality obstruction for v_B stays CONTINUOUS (the open
U(4)xU(4) block-rotation question, same as the numerical overlap faces). Net of
both structural probes: the census's empirical rigidity does not reduce to
Grassmannian-cluster or Frobenius-Schur structure in a form that decides reality
or reaches v89/v103. The moment-map / GIT-symplectic framing remains the correct
organizing lens (it is established: Klyachko, Berenstein-Sjamaar, Walter et al),
but the mechanism behind the observed laws (fiber continua, denominator ratio,
2-elementary holonomy) is still open. Both probes are logged as tested leads.

NORMAL-FAN PROBE RESULT (this one PAID OFF, exact): computed the normal-cone
lattice index at every simple vertex of all nine systems (Smith/|det| of the
incident primitive facet normals, projected to the polytope lattice
Z^d / Z(1..1); validated on v_B, |det| = 23 = den). Findings, all exact integer
arithmetic:
- den = normal-cone index at ALL 51 simple INTERFERENCE vertices across every
  system (the paper conjectured this and had verified only v_B; it now holds
  51/51). Interference vertices are exactly the minimal, cyclic-order-den toric
  quotient singularities.
- DESIGN simple vertices have index = m * den with m in {1,2,3,4,5,6}
  (m=1: 60, m=2: 45, m>=3: 23). So the law sharpens to a one-directional
  characterization: normal-cone index > denominator implies DESIGN (never
  interference). Interference forces m = 1.
This is the toric-quotient-singularity structure the moment-map lens predicts,
and unlike the cluster and Frobenius-Schur probes it is CONFIRMED and exact. It
strengthens the paper's simple-vertex denominator conjecture from one vertex to
all 51 simple interference vertices, and adds the design index-multiple law.
Full Ehrhart quasi-polynomials (which would explain the counting sequences
directly) still need a lattice-point counter (LattE/Barvinok), not vendored.

SNF REFINEMENT (the honest limit of the normal-fan result): among the saturating
(index = den > 1) simple vertices, the Smith normal form is CYCLIC for BOTH
classes (51 interference, 60 design, all one non-unit invariant factor = den). So
the normal-cone arithmetic does NOT separate interference from design: those 60
designs and 51 interference vertices have identical normal cones. Consequence for
the epistemics: den = index is a general saturation fact (den | index always, by
Cramer on the incident facet equations; equality is the tight case), and it is
too COARSE to carry the interference/design distinction. That distinction is a
function of the full determinant combinatorics of the vertex (classify.py decides
it by integer programming, no state solve), not of the normal-fan arithmetic. The
genuinely new content of the census is therefore NOT the base-polytope arithmetic
(a function of the known Klyachko polytope) but the FIBER data over the vertices:
the explicit extremal states, their fiber dimensions, holonomies, and the
design/interference/complex trichotomy of the fibers. Klyachko proved the base
polytope exists; this work computes the fibers, which the base does not encode at
the level of its normal-fan arithmetic. The normal-cone result is a clean exact
statement about the KNOWN object; the new information lives one level up, in the
states.

FIBER PROBES (pushing into the new-information side). Two computations over the
state data.

(A) BASE-FIBER BRIDGE, combinatorial not arithmetic (this works): the vertex's
facet-incidence combinatorics partially predicts the fiber dimension, where the
normal-cone arithmetic could not. Over the certified corpus (fiber = incidence
kernel dim of the support; excess = facets on the vertex minus dim):
- fiber dim 2 (surfaces) occurs ONLY at simple or excess-1 vertices (15/15), a
  clean necessary condition: big fibers live at near-simple vertices.
- inverse tendency: non-simple vertices skew rigid (fiber 0: 70/105), simple
  vertices carry the positive-dimensional fibers (26/49 have fiber > 0).
- designs are always rigid (643/643).
So facet-incidence combinatorics (not normal-cone arithmetic) is the base
quantity that bridges to fiber dimension. Partial (a correlation plus the hard
dim-2 constraint), not a determination.

(B) v_B FIBER as a variety: it is a genus-1 ELLIPTIC CURVE (now derived, not
sampled; scripts/vb_fiber_ideal.py). Displace the eight support weights along the
kernel vector, w(t) = w0 + t (1,-1,0,0,-1,1,0,0). This keeps every orbital
occupation fixed (the kernel is exactly the cycle that leaves the diagonal 1-RDM
invariant: occ = (20,14,14,14,5,4,4,4,13)/23 for all t), so the whole 1-RDM stays
block-diagonal with a SINGLE nontrivial 2x2 block on modes (4,8), diagonal
(5,13)/23, which must split into eigenvalues (14,4)/23. That forces the one
off-diagonal magnitude to the fixed value 3/23. Exactly two determinant pairs
drive that entry, (D0,D1) and (D4,D5), and they carry OPPOSITE fermionic signs,
so the modulus condition is
  |sqrt(w0 w1) e^{i theta} - sqrt(w4 w5)|^2 = 9,
i.e. the defining polynomial (c = cos theta, the gauge-invariant holonomy cosine)
  F(t,c) = 4 c^2 (1+t)(8-t)(4-t^2) - (2 t^2 - 7 t - 3)^2 = 0,
irreducible over Q, bidegree (deg_t, deg_c) = (4, 2). The physical branch is
  c(t) = (3 + 7 t - 2 t^2) / (2 sqrt((1+t)(8-t)(4-t^2))),
so at t=0 cos = +3 sqrt2 / 16 (the earlier -3sqrt2/16 was a sign error: the two
pairs interfere with a relative minus) and at t=1/2 cos = 4 sqrt3 / 15, both now
carrying an EXACT verify_exact certificate, and numerically the re-solved state
matches the SPEC eigenvalues to 1e-15 across the entire real domain
t in (-1, 2) (positivity of the weights and the radicand). Genus: substituting
s = sqrt((1+t)(8-t)(4-t^2)) makes the curve birational to s^2 = (1+t)(8-t)(2-t)(2+t),
a quartic with FOUR DISTINCT roots {-2,-1,2,8}, hence a nonsingular genus-1
elliptic curve (not the rational curve a naive law-of-cosines would have given).
Its quartic invariants are I=1116, J=-66528, discriminant 4I^3 - J^2 =
1133740800 != 0, and j-invariant = 1906624/225 = 2^6 * 31^3 / (3^2 * 5^2). So the
extremal fiber over a single interference vertex carries genuine modular
structure (for THIS slice). (The earlier retraction stands as to the wrong
one-phase reduction; this derivation replaces it and pins the magnitude condition
to the correct opposite-sign form.)

CAVEAT (genus is slice-dependent, not a vertex invariant -- corrects an earlier
overclaim). "The v_B fiber is an elliptic curve" is too strong: the genus is a
property of the chosen one-parameter kernel family, not of the vertex. The
family above lives on the LEDGER support (weights (1,8,4,3,2,2,1,2)/23, block
(5,13)->(14,4), x^2 = 9) and is genus 1. The PAPER's Theorem-3 support is a
DIFFERENT 8-determinant support attaining the same vertex (weights
(4,2,7,3,2,2,1,2)/23, block on orbitals (8,9) = (7,11)->(14,4), x^2 = 21,
cos gamma = 3/(4 sqrt14)); its kernel cycle is (0,1,-1,0,-1,1,0,0) and its
one-parameter family is s^2 = (2+t)^2 (7-t)(2-t) -- a DOUBLE root, hence
genus 0 (RATIONAL). So two natural slices through v_B have different genus; the
elliptic structure is real but is a statement about the ledger slice, not "the
fiber." (This also resolves the parallel instance's "(7,11), x^2 = 21" remark: it
is not a phantom, it is exactly the paper's own support -- I was wrong to file it
as unconfirmed.)

ARC ENDPOINTS ARE REAL STATES -- v_B admits a REAL extremal state, so the paper's
open question is ANSWERED (YES). Where the physical branch reaches cos theta = +-1
the holonomy is 0 or pi and every amplitude is real. Solving F(t, +-1) = 0 gives
85 t^2 - 70 t - 119 = 0, so t = 7/17 -+ 18 sqrt35 / 85 (= (35 -+ 18 sqrt35)/85),
BOTH interior to the positivity range (-1, 2) (t ~ -0.84105 and +1.66458). At each
wall the real-sign state certifies EXACTLY over Q(sqrt35) (verify_exact True,
all-real amplitudes; scripts/vb_fiber_ideal.py -> real_wall_states):
  t = 7/17 - 18 sqrt35/85 (cos theta = -1), signs (+,-,+,+,+,+,+,+), weights
    (24/17 - 18sqrt35/85, 129/17 + 18sqrt35/85, 4, 3,
     27/17 + 18sqrt35/85, 41/17 - 18sqrt35/85, 1, 2)/23;
  t = 7/17 + 18 sqrt35/85 (cos theta = +1), signs all +.
So v_B (the vertex/spectrum) IS realized by a real extremal state; what fails is
realifying the HISTORICAL ledger state psi_B(t=0), which sits at cos = 3sqrt2/16
!= +-1 and is genuinely complex up to gauge. Corrections this forces on the
narrative (do BEFORE any paper/blog line):
1. "v_B has no real extremal state" / "the first vertex that can't be real" is
   FALSE as a vertex-level statement and must be softened. The correct residual
   theorems are class-relative: no real state with RATIONAL squared amplitudes in
   the searched classes (every exhaustive sweep ran on rational grids; the real
   states live at quadratic-irrational Q(sqrt35) weights, which is exactly why the
   grid searches missed them), and psi_B(t=0)'s U(9) orbit is non-realifiable.
   "The phase is forced" is CLASS-RELATIVE, not a fact about the vertex.
2. The interference CLASSIFICATION is untouched: the wall states are real but NOT
   designs (an active exchange pair and irrational weights), so the design
   infeasibility certificates stand.
3. Clean fiber picture: one arc, complex in its interior, real exactly at its two
   endpoints, containing the previously-known certified points. Conjecture
   (testable on every other deforming interference state): interference fibers
   generically interpolate between real endpoint states through a forced-complex
   interior.

REAL STATE ON THE PAPER'S OWN SUPPORT (the cleanest witness, use this in the
paper). Displace psi_B's squared weights along its support's incidence cycle
(0,1,-1,0,-1,1,0,0): k(t) = (4, 2+t, 7-t, 3, 2-t, 2+t, 1, 2). The single
off-diagonal magnitude |rho_{89}|^2 = 21/23^2 is fixed; both exchange-pair
fermionic signs are +1, so |rho_{89}|^2 = k1 k2 + k4 k5 + 2 sqrt(k1 k2 k4 k5)
cos theta. Real means cos theta = +-1: solving with cos theta = +1 gives
t = 55/109 -+ 42 sqrt15 / 109, both interior to the positivity range and BOTH
certified real over Q(sqrt15) (verify_exact True, all-real amplitudes;
scratchpad support_A_certify). At the canonical t=0 the same relation gives
cos theta = 3/(4 sqrt14) = cos gamma, the Theorem-3 phase. So the narrative is
exact and self-contained: deform psi_B's weights along its one free support
direction and the phase gamma relaxes to 0 -- the endpoint is a REAL extremal
state for v_B with weights in Q(sqrt15). This is what should go in the paper
(same support as Theorem 3), with the Q(sqrt35) ledger-support state as an
independent second witness.

CORPUS-WIDE (real fiber-wall states are GENERIC, verified 2026-07;
scripts/wall_test.py). The v_B mechanism is not special. Ran the exact wall test
over single-block, kernel-dim-1 interference vertices whose SHIPPED states are
complex: displace the weights along the incidence cycle, solve the active
two-term block's law-of-cosines wall polynomial exactly, verify_exact the
real-sign state at each real root inside the positivity range. SEVENTEEN certify
(v_B plus 16), each verify_exact True with all-real amplitudes at quadratic-
irrational weights:
- (4,9) v_B: t = 7/17 + 18 sqrt35/85
- (3,8) v28 (Q sqrt14), v31 (Q sqrt30)
- (3,9) v34 (Q sqrt3), v37 (Q sqrt2)
- (4,10) v34 (Q sqrt3), v37 (Q sqrt2)  [transported siblings of the (3,9) pair:
  wall t-values coincide, v34 -> -sqrt3/2 and v37 -> -8/49 - 40 sqrt2/49, exactly
  as functorial padding transport of the whole fiber predicts]
- (3,10) v40 v49 v57 v73 v81 v92 v96 v99 v101 v108 (fields sqrt21, sqrt5, sqrt6,
  sqrt15, sqrt6, sqrt5, sqrt5, sqrt10, sqrt6, sqrt5)
So an interference vertex being REALIZED BY A COMPLEX shipped state does NOT mean
it forces complexity: real extremal states at the fiber walls are generic.

SCOPE / WHAT IS NOT CLAIMED (discipline). Independently verified here: the 17
positives above (exact t, verify_exact, all-real). NOT re-derived here: a parallel
run reported testing 26 of 44 kdim-1 families and leaving 10 with no first-wall
real state, but flagged that only the FIRST common root was sign-searched (second
walls untried) and multi-class fermionic sign-compatibility was not analyzed --
so those 10 are NOT negative results and are not recorded as such. Also untested:
kdim-2 vertices (fiber is a surface, walls are curves not points) and multi-active-
block interference. The clean, defensible statement is the one above: for every
single-block kdim-1 interference vertex tested, a certified real fiber-wall state
exists.

## FIBER-CURVE ARITHMETIC: honest genus census (a "CM surprise" REFUTED, 2026-07)

Both this instance and a parallel run computed the genus of the fiber curve
w^2 = Q(t) (Q = product of the four weight factors of the single active block) for
every single-block kdim-1 interference family, and BOTH first over-reached before
correcting. Recording the corrected, agreed result and the DEAD claims so neither
is re-applied from a pasted patch.

REFUTED (do NOT re-add): "40 of 44 fiber curves are genus 1", "j = 1728 nineteen
times (lemniscatic CM by Z[i])", and the "harmonic wall" cross-ratio verification
(published twice with self-corrections). ALL ARTIFACTS of a wrong reduction of
w^2 = Q when Q has a DOUBLE root: a repeated linear factor is an even-order zero,
so it is NOT a branch point -- divide it out entirely (w^2 = (t-r)^2 g gives a
conic v^2 = g, genus 0). The bad pass kept one power of the repeated factor,
inventing a phantom cubic with a spurious branch point at infinity whose 4-point
cross-ratio came out harmonic, hence the fake j = 1728. Riemann-Hurwitz settles
each case: branch points = odd-multiplicity roots of Q (+ infinity iff their
count is odd); genus 1 iff that count is 3 or 4.

HONEST CENSUS (scripts/fiber_curves.py, committable source of truth): of the 44
single-block kdim-1 families, GENUS 1 = 18, GENUS 0 (rational) = 26. The elliptic
fiber is a structured MINORITY. The 26 rational ones share a weight between the
two exchange terms (a repeated factor -> conic); one direction of the earlier
"guessed genus law" is that trivial fact, the converse fails on coincidental root
collisions. (Reconciliation note: a looser filter that keyed only on active 2-term
off-diagonals counted 20/9-distinct-j by wrongly admitting (5,10) v140 and v263,
which have TWO touched exchange classes -- multi-block, not a single w^2 = Q curve;
the strict "exactly one touched class" filter gives 18.)

WHAT SURVIVES (independently verified here): 8 distinct rational j-invariants over
the 18 genuine genus-1 families, and transport towers preserve the curve exactly
(the face-map functors are isomorphisms of fiber curves):
- j = 5088448/441 x7: (3,8) v31 -> (3,9) v48 -> (3,10) v95 -> (4,9) v31 ->
  (4,10) v48 -> (5,10) v31 (+ pad)
- j = 1906624/225 x2: v_B (4,9) v65 and its (5,10) v65 Hodge partner
- j = 19248832/3969 x2: (4,9) v42, (4,10) v62;  964430272/480249 x2;
  3631696/2025 x2;  three singletons (3,10) v40, v73, v101
- v_B's j = 1906624/225 = 2^6 * 31^3 / (3^2 * 5^2) reconfirmed by both the
  quartic-invariant and cross-ratio routes.

MORDELL-WEIL RANKS (validated in PARI ellrank by the parallel run; not re-run in
THIS environment, which has no PARI/Sage/gp -- so exact numbers below are on that
authority, not a local recomputation): over the 8 genuine isomorphism classes,
v73 rank 0 (torsion order 8), five classes rank 1 (including v_B), two rank 2. An
independent corroboration of the corrected genus split falls out of the same run:
ellrank REFUSED exactly 4 curves, and those four are precisely the fake
double-root "cubics" -- PARI flagged the genus artifact on its own. Locally
checkable and consistent: the branch points are rational (roots -w0/u of rational
weight factors), so the degree-4 models carry full rational 2-torsion
(Z/2 x Z/2, order 4); "order 8" and the ranks are the parts that needed the
descent. Point-locus dictionary: 2-torsion (w = 0) = the support-degeneration
walls; the |cos| = 1 locus = the real wall states (generally irrational t);
rational (t, w) points = rational-holonomy-cosine states, whose supply the rank
governs. The 26 genus-0 fibers' (simpler) point arithmetic is uncatalogued.

## REALITY IS A POSITIVITY BOUNDARY: the Cayley-Menger bridge (verified 2026-07)

The reality question is not a separate kind of constraint from positivity -- it is
the BOUNDARY of a positivity region, one level up from the moment polytope. The
Schur-Horn condition fixes each interference off-diagonal MAGNITUDE |x_e| = sqrt(tau_e).
The amplitudes on the one-hop-connected determinants are vectors that must sum to
x_e: a closed polygon (for a 2-term block, a triangle with sides A, B, sqrt(tau),
apex angle = the holonomy). Then:
- COMPLEX state exists  <=>  the polygon CLOSES  <=>  triangle inequality
  |A - B| <= sqrt(tau) <= A + B  <=>  the amplitude Gram matrix is PSD. This is a
  convex/positivity condition, the exact analogue of "the 1-RDM is PSD" that cuts
  out the moment polytope, but applied to the amplitude vectors.
- REAL state  <=>  the polygon DEGENERATES (collinear, area 0)  <=>  the
  CAYLEY-MENGER DETERMINANT VANISHES. Reality is the boundary of the positivity
  region, not its interior.

EXACT IDENTITY (verified, scripts reality_positivity): for the 2-term block with
side-squares (p1, p2, tau), 16 Area^2 = -[(tau - p1 - p2)^2 - 4 p1 p2] = -(wall
polynomial). So the Cayley-Menger determinant IS the wall polynomial: for v_B it is
85 t^2 - 70 t - 119, whose roots are exactly the walls 7/17 -+ 18 sqrt35/85, and
16 Area^2 > 0 on the physical (complex-allowed) interior, < 0 outside (no triangle,
|cos| > 1). General for every single 2-term-block family, not just v_B.

WHAT THIS BUYS (the transferable machinery):
- A basis-free UNIVERSAL reality certificate: reality <=> a Gram/Cayley-Menger
  determinant vanishes. Replaces the per-family wall-polynomial solve with one
  determinantal object, and it is the SAME object.
- The multi-edge case (the "10 pending" sign-compatibility gap) becomes a clean
  real-algebraic question: do the several edge Cayley-Menger hypersurfaces have a
  COMMON real zero inside the positivity region? That is a positivity/semialgebraic
  problem with real tools (real Nullstellensatz, SOS certificates), not an ad hoc
  sign search.
- The k-term block generalizes: reality = the k amplitude vectors are COLLINEAR =
  the amplitude Gram matrix drops to rank 1 = all 2x2 minors vanish = a
  determinantal variety; positivity (state exists) = Gram PSD. Reality is the
  deepest boundary stratum of that spectrahedron. Positivity of weights is
  polyhedral (LP, the polytope); reality is spectrahedral (PSD-boundary): the
  transfer is polytope -> spectrahedron.
- Positive geometry: the physical region is {Cayley-Menger <= 0} and reality is its
  boundary component, so the reality locus is a natural boundary of a positive
  geometry -- a concrete handle for the Section-on-positive-geometry connection the
  paper flags (cos gamma as "new data").

This is already half-present in the paper's exactifier language ("each off-diagonal
is a closed polygon whose relative phase is an exact arccosine"): the polygon is the
positivity object, its degeneration is reality. Now IN the paper's positive-geometry
section (reality = a spectrahedral positivity boundary on the amplitudes, the
Cayley-Menger determinant = the wall polynomial).

THE BRIDGE IS ALSO A BETTER SOLVER, AND THE CENSUS IS NOW COMPLETE
(scripts/wall_solver.py): taking candidate t as the common real zeros of the edges'
Cayley-Menger determinants and letting verify_exact adjudicate the global sign
assignment certifies ALL 44 genuine single-parameter deforming kdim-1 interference
families with an exact real extremal state (up from 17 with the single-block,
single-wall, single-sign scripts/wall_test.py). This closes the "10 pending /
multi-class sign-compatibility" gap wholesale: framing reality as simultaneous
Cayley-Menger degeneration turned an ad hoc sign search into a common-zero
computation.

THE TWO REMAINING ARE NOT A GAP -- they are structurally different. The
multi-touched-class (5,10) families v140 and v263 each carry a 1-TERM active edge
(the off-diagonal rho_{2,5}, target 7/1521). A one-term off-diagonal has no phase
to adjust, so its magnitude sqrt(w_i w_j)/den is RIGID; along the incidence kernel
that magnitude equals its Schur-Horn target only at a single point (t = 0, the
shipped complex state). Hence there is NO 1-parameter spectrum-preserving fiber to
carry real endpoints: the deformation is pinned, and the wall/reality method is
inapplicable by construction, not by solver deficiency. Their reality is genuinely
OPEN (a numeric real-restricted attainment floors at ~1e-11, indistinguishable from
the v_B control that provably has a real state, so it settles nothing). They are the
census's candidate GENUINELY-COMPLEX interference vertices -- the analogue of what
v_B was once thought to be -- and deciding them needs a construction outside the
single-kernel deformation (a different support, or an exact real-infeasibility
certificate à la the design MILP). So the clean statement: EVERY genuine
single-parameter interference family (44/44) admits a certified real extremal state,
the certificate being a positivity (Cayley-Menger) condition; exactly two nominal
kdim-1 families are rigid and their reality is the remaining open problem.

## PHYSICS OF THE FIBER: the fiber is a 2-body observable, and a REALITY
SELECTION PRINCIPLE (verified here for v_B, 2026-07)

VERIFIED HERE (scripts/vb_fiber_ideal.py -> fiber_physics_vb): the fiber coordinate
t is invisible to 1-body data but is a TWO-BODY observable. Along the v_B ledger
family exactly EIGHT pair occupations <n_i n_j>(t) = sum_{T contains i,j} w_T(t)/23
vary, all linear in t and all on the exchange-channel modes, e.g.
<n_2 n_4> = (3 + t)/23, <n_1 n_8> = (9 - t)/23. Two states at the same maximally
pinned vertex (identical natural occupations) therefore have measurably different
density-density correlations: occupation spectroscopy is blind to the fiber, pair
spectroscopy reads it.

REALITY SELECTION PRINCIPLE (verified for v_B, general derivation, scope stated).
The PHYSICAL interval of t is where |cos theta(t)| <= 1, and for v_B its endpoints
are exactly the |cos| = 1 REAL walls t = 7/17 -+ 18 sqrt35/85, which sit strictly
INSIDE the positivity interval (-1, 2) (checked: |cos| < 1 at t=0, |cos| > 1 at
t = -0.9 and t = 1.8, both still weight-positive -- so reality, not positivity, is
the binding constraint). Since every <n_i n_j>(t) is linear, any purely
density-density Hamiltonian H = sum J_ij n_i n_j has <H>(t) LINEAR in t, hence its
minimum over the family is at an endpoint of the physical interval -- a REAL wall
state. So density-density interactions pin a pinned state to REALITY; the
forced-complex interior is selected only by Hamiltonians with genuine exchange
(nonlinear-in-t phase) terms. This gives the wall states a physical role
(interaction-selected fiber representatives) and is the natural RDMFT statement:
the exact functional's boundary degeneracy along the fiber is lifted LINEARLY by
density-density terms. SCOPE (per-family, not yet census-wide): the "endpoints are
real walls" step needs the |cos| = 1 walls to bind before positivity, verified for
v_B and consistent with the 17 corpus vertices (all had interior real walls) but
NOT checked for every family; linearity of <H> is general.

ENDPOINT THEOREM (verified here, 44/44 -- the scope caveat above is now DISCHARGED
for the "walls interior" half). For a single active two-term block,
cos^2 theta(t) = (tau - p1 - p2)^2 / (4 p1 p2) with p1, p2 the two weight products
(linear in t) and tau the fixed exchange target. As t approaches any positivity
boundary some weight -> 0, so a p_i -> 0 while the numerator -> (tau - p_j)^2 != 0,
hence |cos theta| -> infinity: the |cos| = 1 walls are ALWAYS crossed strictly
before positivity fails. Checked across every single-block kdim-1 family: 44/44
(46/46 incl. the two multi-touch (5,10) cases) have the physical interval strictly
interior to positivity, bounded by |cos| = 1 walls on both sides (|cos| > 1 at both
positivity boundaries, physical region nonempty). So density-density minima land on
|cos| = 1 boundary states (real: holonomy 0 or pi) for EVERY deforming single-block
family, not just v_B. (What remains per-family is only the exact symbolic
verify_exact of each wall's sign realizability: 17 done, the rest structural.)

EXCHANGE SELECTS COMPLEX (mechanism verified here; exact v_B numbers parallel-
reported). Adding a correlated-exchange term makes <H>(t) nonlinear: for v_B the
exchange expectation <n_2 (a_4^dag a_8 + h.c.)>(t) reduces exactly to a multiple of
(3 + 7t - 2t^2)/(2 sqrt(4 - t^2)), a nonlinear function of t. Minimizing
[linear density-density] + lambda [exchange] over the physical interval: at
lambda = 0 the min is the hi wall (real); at intermediate lambda the min moves
INTERIOR (a genuinely complex extremal state is selected); at larger lambda it
snaps to the lo wall, the jump between basins being discontinuous -- a first-order-
like, two-basin selection crossover. So the "only genuine exchange reaches the
complex interior" clause is backed by an explicit two-body mechanism. Parallel run
reports the exact instance: exchange n_2(a_4^dag a_8 + h.c.) at mixing ratio 0.10
puts the global min at t = -0.465 (interior/complex), with the transition near
ratio 0.08-0.10 (their specific couplings; I reproduced the mechanism, not the
exact threshold).

PARALLEL-REPORTED, NOT INDEPENDENTLY CHECKED (a null result, low stakes): the
quantum-metric / elliptic-differential bridge is COLD -- neither g(t) Q(t) nor
g(t) sqrt(Q(t)) is constant along the v_B family (g the Fubini-Study metric), so
no naive identity between the physical metric and the curve's invariant
differential. Parked without decoration.

CANDIDATE FOR THE PAPER (flag, do not inject unilaterally): the reality selection
principle + the 2-body-observable framing is the most paper-ready physics to come
out of the fiber work -- it connects the reality result to RDMFT and pinning
(Section on applications) with a falsifiable, elementary statement. Hold for an
explicit decision on scope.

## FIBER-DIMENSION CENSUS: the loopy-state question ANSWERED (Stage A full, Stage B sampled)

Method (two stages). STAGE A, exact linear algebra per state: kernel of
the support incidence = the weight-deformation space (mode sums frozen);
classify each one-hop class by (term count, target); a 1-TERM class is
RIGID (its off-diagonal magnitude sqrt(k_i k_j)/den is pinned by the
target spectrum, so any kernel direction with v_i + v_j motion changing
the product is cut -- first-order test: v_i k_j + v_j k_i = 0); classes
with >= 2 terms are OPEN (phases/signs absorb motion in the polygon
interior). Predicted fiber dim = kernel dim minus the rank of the rigid
cuts. STAGE B, certification: displace to a rational t, re-solve the
phase, verify the exact char-poly identity (or exhibit the rigid-class
obstruction for blocked states).

RESULTS (797 states):
- NEW LAW (exact): ALL 643 design states are LOOP-FREE -- design
  supports never carry incidence kernels, so design fibers are rigid
  points (up to symmetry orbits). Continuous moduli are exclusive to
  interference.
- Interference: 93 of 154 loop-free (rigid); of the 61 loopy states,
  46 have kernel dim 1 and 15 have kernel dim 2, and Stage A predicts
  FULL kernel-dimension families for 59 of 61 -- i.e. 44 certified-
  candidate CURVES (the 46 kernel-dim-1 states minus the two blocked)
  and 15 SURFACES of extremal states. Exactly TWO states are blocked:
  (5,10) v140 and v263, each killed by a rigid 1-term class.
- Stage B certifications: v_B (t=1/2, cos = 4 sqrt 3/15) and (3,10)
  v40 (t=1/2, cos = sqrt(165)/33) displaced states verified by exact
  char-poly identity; v140's blockage PROVED exactly (its 1-term class
  (2,5) has product k_0 k_1 = 7t + 7 along the kernel, non-constant,
  so any t != 0 moves the block off-diagonal with the diagonal fixed
  and shifts the eigenvalues -- QED).
- FIBER-DIMENSION LAW (empirical, 59/61 + 2 proved exceptions):
  fiber dim = kernel dim, except where a rigid 1-term class cuts it.
  Remaining Stage B work: certify the other 57 displaced points (same
  script), certify one 2-parameter surface point, and prove v263's
  blockage (same argument shape as v140).

Consequence for the paper (do not add until Stage B completes): the
extremal fibers of these polytopes are generically positive-dimensional
at loopy interference vertices -- the census's one-state-per-vertex
library is sampling curves and surfaces, and the multiplicity structure
(symmetry orbits x continuous families) is now computable exactly.

VERIFICATION (this session, against the 797 ledger): the design loop-free
law reproduced exactly (643/643 kernel dim 0); the interference kernel-dim
distribution is {0: 93, 1: 46, 2: 15} (61 loopy) -- NOTE the source patch
wrote "44 have kernel dim 1", which is the CURVE count (46 dim-1 minus the
2 blocked), not the dim-1 count; 44 + 15 = 59 predicted families is
consistent, 46 + 15 = 61 loopy is the correct partition. v140 and v263
confirmed kernel-dim-1 interference; the v_B (cos 4 sqrt 3/15) and v40
(cos sqrt(165)/33) t=1/2 displaced states both pass verify_exact. The
v140 first-order blockage arithmetic (7t+7) was not re-derived here.

## Toric-fiber conjecture: both tests RUN; a continuum of extremal states at v_B

Ran the two named tests. Results are a genuine mixed verdict plus one
discovery that supersedes the question as posed.

- DISCOVERY (exact, certified): v_B's extremal fiber is POSITIVE-
  DIMENSIONAL. The degree system on its 8-det support has a 1-dim
  kernel (the vector (1,-1,0,0,-1,1,0,0)), and moving the weights along
  it PRESERVES extremality because the block phase re-solves: at
  t = 1/2 the displaced state with weights (3/2,15/2,4,3,3/2,5/2,1,2)/23
  and cos(theta) = 4 sqrt(3)/15 satisfies the exact char-poly identity
  (verified symbolically). The family lives on t in (-1, 2); at both
  boundaries a weight vanishes, the exchange class drops to one term,
  and NO real sign choice certifies -- so the family is genuinely
  complex throughout and the paper's open question (does v_B admit any
  real extremal state) STANDS, now over a strictly larger known search
  space (the whole family, each point a distinct U(9) orbit; the
  non-realifiability optimization was run only at t = 0).
  FRAMING CAUTION (own goal avoided): weight-deformation directions and
  holonomy loop vectors are ker of the SAME incidence map -- their
  coincidence is definitional. The nontrivial content is PHASE
  ABSORPTION: the char-poly conditions do not cut the kernel direction;
  the phase absorbs it. First known continuum of extremal states at a
  GPC vertex; the certified psi_B is one point of a curve.
- LITERAL TORUS-ORBIT FORM: REFUTED. Along the family the loop monomial
  prod k^v is NOT constant (1/8 at t=0, 1/3 at t=1/2), so the fiber
  curve is not a naive torus orbit in amplitude coordinates.
- BINOMIAL-SHAPE FORM: MOSTLY SUPPORTED, exceptions identified. Across
  all 154 certified interference states, one-hop classes have 1 term
  (123), 2 terms (61), or 3 terms (15); the >= 3-term (non-binomial)
  states are exactly the 15 trisection states -- whose gauge-invariant
  content nonetheless collapses to the abelian catalog. So the phase
  ideals are binomial for ~92% of classes, trinomial precisely where
  cube roots appear, with invariants lattice-clean even there.
- TEST 2 AS DESIGNED: INCONCLUSIVE (the 24 v96 fiber states occupy 24
  distinct supports; no shared-support pairs to test loop-monomial
  constancy). Redesign: scan the v_B family instead -- its defining
  curve in (k, z) coordinates is rational, cut by
  m1^2 + m2^2 + 2 m1 m2 cos(theta) = x2 with m_i^2 linear in t.
REVISED CONJECTURE DIRECTION (speculative, flagged): the fiber is not a
torus orbit but a rational curve (or higher family) fibered over the
kernel lattice, with holonomy varying along it -- the lattice controls
deformations and phases jointly, cluster-algebra-like rather than
naively toric. Census-wide question now well-posed: for every loopy
certified state, does the kernel direction always survive phase
absorption (fiber dim = kernel dim), or do 1-term rigid classes block
it? That is a finite exact computation over the corpus.

VERIFICATION (this session): the t=1/2 displaced state passes verify_exact
exactly; the loop monomial is 1/8 at t=0 and 1/3 at t=1/2; 154 certified
interference states and 15 trisection states confirmed against the ledger.
Not recomputed here: the full 123/61/15 one-hop term-count histogram.

## Symmetry census at 797: spontaneous symmetry breaking at the corners

Computed the FULL automorphism group of every certified state (mode
permutations within spectrum-degeneracy classes preserving the weighted
support; 1-WL color refinement + in-cell brute force; 788/797 computed,
9 skipped on cell size). New invariant per vertex: the SYMMETRY-BREAKING
INDEX [G : Aut(psi)], G = the spectrum's class-permutation group.

- SLATER RIGIDITY LAW (exact): exactly NINE states realize the full
  symmetry of their spectrum (Aut = G) -- and they are precisely the
  nine Slater corners (index 0 of each system, single determinant).
  Every correlated extremal state in the census breaks its spectrum's
  symmetry. Theorem-shaped; worth a short proof attempt (a G-invariant
  correlated state on a vertex may be impossible for elementary
  reasons).
- STRATIFICATION BY CLASS: rigidity (Aut = 1) hits 63% of interference
  states vs 40% of integer designs -- interference correlates with
  symmetry breaking. ALL 12 DESIGN-REAL states have nontrivial Aut:
  combined with the ratio-2 law, real designs are doubly special
  (symmetric AND on the doubled grid).
- LARGE EXTREMAL MULTIPLETS: the breaking index concentrates at
  16-144, so each such vertex carries an orbit of that many distinct
  extremal states -- massive extremal degeneracy, previously invisible.
- HODGE BREAKING: of the 14 Hodge-self-dual spectra at half filling
  with certified states, only TWO ((4,8) v21 and (5,10) v291, both
  two-determinant states) have complement-closed (self-dual) supports;
  the other 12 are attained by non-self-dual states, whose Hodge images
  are DIFFERENT extremal states in the same fiber. (Shipped-state
  statement only: whether every self-dual spectrum also admits a
  self-dual representative is a clean open question -- the two
  existing ones show it is sometimes possible.)

Speculative framing, flagged: this is spontaneous symmetry breaking at
polytope corners -- symmetric spectra, asymmetric extremal states,
degenerate multiplets of order [G:Aut]. A toric extremal fiber
(the standing conjecture) would organize exactly these orbits.

VERIFICATION (this session): the nine single-determinant Slater corners
(all at index 0) and the two Hodge self-dual supports ((4,8) v21,
(5,10) v291, both two-determinant, complement-closed) confirmed against
the ledger. Not recomputed here: the [G:Aut] automorphism statistics
(the 63%/40% stratification, the 16-144 multiplet range, all-12-DESIGN-
REAL-nontrivial), which rely on the automorphism tooling not shipped in
the repo.

## Corpus-wide algebraic mining at 797/799 (2026-07)

Three exact findings and one honest negative from mining every certified
state for algebraic relations and generating structure:

1. DENOMINATOR-RATIO LAW (exact, 797/797): state-denominator over
   spectrum-denominator is ALWAYS 1 or 2. All 785 DESIGN-INT and
   INTERFERENCE states sit at exactly 1; the ratio-2 states are
   precisely 9 of the 12 DESIGN-REAL states (12/6, 26/13, 28/14; the
   other three DESIGN-REAL sit at 1). So the m=2 amplitude grid IS
   realized in nature -- exclusively by real designs. Refines the old
   "ratio-1 law" and gives DESIGN-REAL a sharp arithmetic signature.
2. GAUGE TRISECTION, INVARIANT ABELIANNESS: 16 certified interference
   states carry TRISECTED amplitude phases (atan-form arguments like
   (pi + 3 atan sqrt(7))/3 -- cube roots of complex numbers, from the
   cubic characteristic polynomials of 3-clique blocks). These look
   non-abelian at the amplitude level, but the GAUGE-INVARIANT loop
   holonomies collapse every time into the exponent-2 catalog (checked
   at (3,10) v17: loops give minpolys 2x^2-3x+2, field Q(sqrt(-7)), and
   the known octic 4x^8-6x^6+5x^4-6x^2+4). Conjecture 2 is doing real
   work: parametrizations may trisect; invariants may not.
3. FIRST NON-EVEN RECIPROCAL HOLONOMY + NORM-SQUARE LAW EXTENSION: the
   fresh v113 crack carries holonomy minpoly 12x^4-2x^3-19x^2-2x+12 --
   reciprocal but NOT even, the first holonomy outside the forced-V4
   lemma's hypothesis in the census. Abelianness survives via the
   norm-square mechanism one level up: with y = x + 1/x satisfying
   12y^2-2y-43, the norm (y1^2-4)(y2^2-4) = 1/16 is a perfect square,
   forcing V4 (confirmed by direct Galois computation). Conjecture 2's
   evidence now includes its first non-even instance, and the
   norm-square law generalizes: for reciprocal holonomy quartics, the
   product of (y_i^2 - 4) over the y-conjugates is a rational square.
   (v113's other loop is the even quartic 32x^4-17x^2+32, forced-V4,
   field Q(sqrt 2, sqrt -94).)
4. NEGATIVE, recorded honestly: the census counting sequences (vertices
   4,10,38,58,113; designs; interference; facets 1,4,31,52,93 along
   (3,d)) admit NO low-order linear recurrence on five terms -- second
   differences are irregular (22,-8,35) -- so no generating function is
   proposed. One qualitative regularity: growth alternates by parity of
   d (jumps into even d are consistently larger: +28 into d=8, +55 into
   d=10, vs +6/+20 into odd d), suggesting half-filling-adjacent
   enrichment; parked as an observation, not a law.

VERIFICATION (this session, against the shipped 797/799 ledger): law 1
reproduced exactly (788 states at ratio 1, 9 at ratio 2, all DESIGN-REAL,
spectrum-denominators 6/13/14); the v17 minpolys (2x^2-3x+2, octic
4x^8-6x^6+5x^4-6x^2+4) and the v113 minpolys (12x^4-2x^3-19x^2-2x+12 and
32x^4-17x^2+32) reproduced exactly, and the norm-square (y1^2-4)(y2^2-4)=1/16
confirmed; the vertex sequence 4,10,38,58,113 and its second differences
(22,-8,35) reproduced. Not independently confirmed here: the exact "16"
trisected-state count (a phase-pattern scan found 15) and the (3,d) facet
sequence 1,4,31,52,93 (no per-(3,d) facet-count artifact ships for d<=8).

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
   PARTIAL-FAMILY CROSS-CHECK and a lead on the 4 OPENs (verified here in exact
   arithmetic, 2026-07). The only rank-11/12 constraint data in print are the
   stable Grassmann families of Altunbulak-Klyachko (CMP 2008, Thms 4.2.1 and
   4.3.1) and Klyachko (arXiv:0904.2009): the first-kind pairs
   lambda_{k+1}+lambda_{r-k} <= 1, the four second-kind quadruples (2345, 1346,
   1256, 1247 <= 2), and a level-5 series extending each quadruple by lambda_11
   (indices 1,2,4,7,11 etc.). Consistency checks all PASS: all 19 certified true
   (3,11) vertices satisfy every family; the proved families restrict validly to
   the {lambda_11=0} face (all 113 certified (3,10) vertices satisfy them). LEAD:
   all four OPEN candidates violate a published level-5 series inequality (idx 23
   and 26 and 34 violate AK-RMK-4.2.1 (1,2,4,7,11 <= 2); idx 44 violates only the
   extended quadruple (2,3,4,5,11 <= 2)); none violate any PROVED family. If the
   series inequalities hold, (3,11) settles completely, 19 TRUE + 31 REFUTED.
   CAVEAT (validation law): the level-5 series is CLAIMED not proved: AK-RMK-4.2.1
   is a remark whose proof CMP 2008 explicitly defers, and the lambda_11-extended
   quadruples are an unrefereed letter claim. Neither is a certificate here, so
   the four stay OPEN. Two ways to close them with proof: (a) prove the level-5
   second-kind inequality inside Stage 1 (a finite Schubert/Grassmann coefficient
   computation of the kind docs/stage1_klyachko_spec.md targets; a positive c_gamma
   for the (1,2,4,7,11) diagram and its three variants converts all four to
   REFUTED); or (b) attain idx 44 (it rests on the weakest, letter-only claim) --
   an exact attainment would FALSIFY a published Klyachko inequality, itself a
   citable result. Either branch is a win.
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
(reproduction: scripts/leads_analysis.py, now in-repo and runnable over all
nine systems; all statements exact rational/symbolic unless marked numerical).
These are LEADS, not shipped results: the arithmetic layer (norm identities,
minimal polynomials, prime factorizations) is hand-verified in this repo, but
the group-theory conclusions (Galois groups via PARI polgalois) are as computed
by the analysis script and not re-run here.

REFRESHED COUNTS (2026-07, regenerated by scripts/leads_analysis.py against the
current 785-certified data; the prose below still carries the older snapshot
numbers where noted). Lead A: 183 simple vertices census-wide (the earlier 179
dropped the Borland-Dennis system (3,6), the sole system with equality
constraints, which full_hrep now includes; (3,6) contributes 4 simple design
vertices, all facet-degenerate so excluded from the normal-index tally). The
den = normal-cone-index law holds for all 51 simple interference vertices
(51/51; the earlier 54 predates a 3-vertex classification refresh, all of which
satisfy the law either way) and fails at 68 of the 128 simple design vertices
with a computable normal index; 24 isotypically pure vertices (23 without
(3,6)). Lead B: over the current 142 certified interference vertices, kernel
dims {0: 92, 1: 37, 2: 13} give 92 loop-free supports and 50 carrying 62
nontrivial holonomies with 24 distinct minimal polynomials (11 quartics, 5
degree-8), superseding the 770-snapshot 127/80/47/59. Any use of these in the
paper (results/report/main.md) must cite these refreshed numbers.

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
  MEASURED CAVEAT (live runs): the stratum is not uniform. (3,10) v96 EXHAUSTS
  its sweep (returns FAIL, not TIMEOUT) at max_clique=4, max_cliques=0
  (capacity 5), max_card=36 -- a true expressiveness gap; reproduced here in
  ~73 s (an earlier handoff note reported ~10 s at a smaller cap; 73 s is the
  measured value at these settings). But (4,9) v42 at the same settings is
  reported COMPUTE-BOUND past 200 s (no exhaustion), not reproduced here. So
  before routing a root to the extension, classify it by exhaust-vs-timeout at
  a short fixed budget: fast-exhaust roots need the new ansatz axis; slow roots
  may still fall to walltime. Preflight table (min 2x2 blocks / min single
  k-clique, max 5/5; reproduced exactly by scripts/leads-adjacent preflights):
  (3,10) v40 None/3, v49 None/3, v57 None/3, v60 None/3, v73 None/3,
  v89 1/None, v96 None/3, v103 1/None; (4,9) v40 None/3, v42 None/3;
  (5,10) v113 None/3, v261 None/3. For the max_card stratum: v89 at
  max_card=44 is reported to run past 280 s (real search space opens),
  consistent with the raise-the-cap-and-wait plan.
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

## Plethysm probe: cross-validation and a saturation warning (2026-07)

Ran the shipped plethysm engine (scripts/plethysm_inner_hull.py) against the
partial-families oracle and the four OPEN (3,11) candidates. Exact arithmetic
throughout.

- CROSS-VALIDATION (PASS): the certified-attainable clouds (217 points at
  (3,11) M<=7 and 88 points at (3,12) M<=6) satisfy every published family in
  partial_families_3_11_3_12.json with ZERO violations, across all three tiers
  (PROVED, WEAK, CLAIMED). Mutual validation of the plethysm engine and the
  transcribed families, and weak positive evidence for the CLAIMED level-5
  series (attainable data keeps respecting it).
- ATTAINABILITY PROBES of the OPENs at m = den: cand 23 (m=7) mult 0; cand 23
  (m=14) mult 0; cand 44 (m=12) mult 0. INCONCLUSIVE (see the calibration below
  before reading these as refutation).
- SATURATION WARNING (the key finding): m = den is NOT sufficient even for TRUE
  vertices. Control probes: cand 22 (TRUE, N=2-lift certificate) has mult 1 at
  m=5, but cand 49, the uniform (3/11)^11 vertex with an explicit Z11
  difference-design certificate, has mult 0 at m=11. A point can lie in the
  moment polytope (mult(k*lambda) > 0 for some stretch k) while the primitive
  multiplicity vanishes: plethysm coefficients are famously non-saturated, and
  cand 49 is a concrete instance in the wild (attainable, mult 0 at the
  primitive lattice point). Consequences: (a) negative probes at small m carry
  almost no evidential weight, so the falsification branch on cand 44 needs m at
  2-4x den (degree 72+, rig-scale, parallel over the mu-support); (b) the
  convergence cost of the whole Stage-0 inner approach is set by the SATURATION
  STRETCH of each vertex, not its denominator, so charting mult vs m for the
  known rank<=10 vertices would calibrate that constant; (c) if cand 49 survives
  scrutiny it may be independently citable (a physically meaningful
  non-saturation example in the fermionic plethysm cone).
- VERIFICATION STATUS (this repo): the saturation control probes (cand 22 mult 1
  at m=5, cand 49 mult 0 at m=11) and the (3,12) M<=6 zero-violation check are
  reproduced here. The (3,11) M<=7 count (217) and the OPEN-candidate probes are
  from the handoff engine and not re-run at this budget (the (6,6,1^9)-type
  characters are expensive).

## Structure scan for the next campaigns (2026-07)

Three empirical structures from re-mining states/vertices/facets, including the
rank-11/12 data. These are HANDOFF LEADS: spot-checked (the flagship signature
pair below is confirmed) but not fully reproduced in this repo; treat the
per-count figures as the scanning engine's, pending an in-repo generator.

1. VERTEX/FACET FUNCTORIAL ASYMMETRY. Vertices are mostly functorial (transport
   dominates); facets are only HALF functorial, even after adding a second facet
   map, SERIES EXTENSION (append a new index with coefficient +-1 or +-2, rhs
   unchanged: the shape of the published level-5 series): pad+extend explains
   48/93 facets at (3,10) and 53/125 at (4,10); appends (-1,0) and (-2,0) alone
   hit 19 and 15 of the (4,10) facets from (4,9). The residual genuinely-new
   facets per rank are characterized: dense support (6-8 nonzero coefficients at
   d=10) and small height (max |coeff| <= 4, growing ~+1/rank). Generator
   consequence: candidate facets = functorial closure (free, ~half) + a bounded
   box (height <= prev+1, >=6 nonzero, decreasing-compatible) for the Schubert
   test. The tower's irreducible novelty lives on the facet side.
2. SIGNATURE GRAMMAR (multiplicity patterns of integer forms). TRUE vertices are
   signature-hereditary: 17/19 rank-11 true vertices reuse a rank-10 signature
   (novel: the boundary shapes (1,10) and (11,)). Signature-NOVELTY predicts
   hardness: all 4 OPEN rank-11 candidates and 6 of the 14 open rank-10 vertices
   ((3,10) v49/v57/v89/v96/v103, (5,10) v261) have signatures never certified
   before. NEW SOLVER TECHNIQUE: ansatz transfer by signature -- warm-start the
   block/clique sweep from certified same-signature interference siblings (same
   mixed classes, split ratios rescaled to the new denominator). 8 of the 14
   have donors; flagship pair (signature (1,4,2,2), confirmed here): (4,9) v42
   (14,9,9,9,9,3,3,2,2)/15 with certified sibling (4,9) v30
   (16,9,9,9,9,4,4,2,2)/16, near-identical architecture at the neighboring
   denominator (the v_B-family phenomenon, now algorithmic). Implement as a
   warm-start stage before the general sweep and before the signed-cancellation
   extension.
3. SUPPORT RIGIDITY. Only 267/785 certified supports admit even one
   equal-occupation-class transposition symmetry as stored (caveat: stored
   supports are orbit representatives, so this underestimates).
   Equivariant/difference-design search is a useful prior (it built the (3,11)
   uniform vertex), not a universal law.

Hardness ordering of the 14 implied: 8 warm-startable by signature transfer; 6
signature-novel (route to the signed-cancellation extension and, failing that,
treat as the genuinely new mathematics of rank 10).

## Generator blueprint: the Edge Law and rhs Law (2026-07, exact)

Reverse-engineered the cubicle/Schubert layer of Stage 1 from the 542 known GPC
facet inequalities (all nine systems, exact rational arithmetic). Two laws:

- RHS LAW (VERIFIED HERE, 542/542): the right-hand side of every known GPC facet
  equals an N-subset sum of its own sorted coefficient vector (the pure-state
  specialization b = a_T). Reproduced in this repo across all nine systems
  (per-system 1/4/31/15/52/60/93/125/161). Its DEPTH (position among the
  distinct subset sums, 1 = maximum = trivial) is reported bounded and slow
  growing: max depth 2 through rank 7, 5 at ranks 8-9, 7 at rank 10 (the single
  depth-1 facet per half-filled system is the Pauli facet lambda_1 <= 1). Depth
  is what the Schubert coefficient certifies, so the generator only needs depths
  2..~8 at rank 11. (The depth analysis is the handoff engine's, not re-run
  here; the b = a_T law itself is reproduced.)
- EDGE LAW (VERIFIED IN-REPO, 542/542): every facet's sorted coefficient vector
  is an EXTREMAL RAY of the arrangement cut on the ordered cone
  {a_1 >= ... >= a_d} by the subset-sum tie hyperplanes a_T = a_T' (tight rank
  exactly d-2, a ray modulo the constant shift, under which the inequality
  system is invariant). Sharper than Klyachko's remark that edges SUFFICE:
  empirically facets ARE edges. Reproduce both laws with
  `python scripts/facet_laws.py --verify` (exact rational arithmetic; also
  prints the per-system depth tables). Worth a proof attempt from Thm 3.2.1.
- Facet grammar corollary (handoff): facet vectors have few distinct values
  (typically 4-6 at d=10, heights <= 7), so facets carry their own
  multiplicity-signature grammar, mirroring the vertex grammar.

STAGE-1 BLUEPRINT this yields (handoff proposal; parameters reported measured):
 1. Enumerate integer extremal rays of the ordered subset-sum arrangement with
    height <= H (H = previous-rank max + 1; H ~ 8 at rank 11): a
    value-multiplicity composition (<= ~6 blocks), decreasing integer values
    <= H, closed under subset ties to tight rank d-2; lrs-style ray enumeration.
 2. Prescreen against functorial closure (pad + series-extension supplies
    ~half the candidates with known provenance).
 3. Per ray, per depth k = 2..K (K = previous max + 2): Schubert test c_vw(a)
    via lrcalc for rhs = the depth-k subset sum.
 4. Assemble; verify inner = outer against the plethysm cloud, the vertex
    attainability decision procedure, and the (3,11)/(3,12) oracles (19+19
    certified true retained, 27 refuted excluded).
Acceptance gates: reproduce 1/4/31/52/93 (the (3,d) counts) at ranks 6-10 first;
the RHS law is a free invariant (any generated facet violating b = a_T is a
bug), and the Edge law becomes one once reproduced.

## Core/hole isomorphisms and decomposability of the 14 (2026-07)

Fresh sweep of functorial reductions over the residual, two results:

- DUAL PAIR (VERIFIED IN-REPO): (5,10) v113 (17,17,16,8,8,8,4,4,4,4)/18 and
  (5,10) v261 (14,14,14,14,10,10,10,2,1,1)/18 are HODGE DUALS of each other
  (lambda -> 1 - lambda at half filling), previously carried as separate roots.
  Checked against states.jsonl: v261's complement (18 - lambda, sorted) is
  exactly v113. They are ONE problem -- an exact state for either yields the
  other by the Hodge star immediately, so the residual is 11 independent
  problems, not 12: the 14 SOLVE-FAIL vertices minus the two (4,10) entries
  (v60, v62), which strip their trailing zero to become (4,9) v40, v42 exactly
  (same denominator 15, verified), minus the merged dual pair. (5,10) is the
  only self-dual system in the census (d - N = N), so this is the only
  within-census dual: the (3,10), (4,9), (4,10) roots dualize to (7,10),
  (5,9), (6,10), none of which is a solved system.
- NO DECOMPOSITIONS (HANDOFF sweep, not re-run in-repo): the campaign reports
  that none of the 14 splits as an exact wedge product of lower states on
  disjoint mode sets, checked over all particle/mode splits with the full
  reduction toolkit (frozen-core peeling, trailing-zero stripping,
  particle-hole reduction, N=2 pairing, certified-library lookup) on each
  part -- zero valid splits. Independently confirmed here only the necessary
  conditions this rests on: no root (the 12, excluding the two (4,10) pads
  where the trailing zeros live) has lambda_1 = 1 or a trailing zero at the
  root level. The exhaustive negative itself is the handoff engine's, pending
  an in-repo reproduction. If it holds, every residual vertex is "connected":
  its extremal state entangles across the full mode set, so no product
  construction cracks any of them.

## v96 campaign: the signed-design family (NO-basis formulation, 2026-07)

New ansatz family for the residual, formulated in the state's OWN natural-orbital
basis, where rho is diagonal and the problem is exact: integer weights k_T >= 0
with mode sums EXACTLY the integer spectrum, plus signs (or phases) making EVERY
one-hop pair class cancel. This is a third family, strictly between the two the
pipeline searches: designs (one-hop-free, rational weights) and block/clique
states (one-hop pairs confined to blocks, algebraic weights). SIGNED DESIGNS
(one-hop pairs present but cancelling, rational weights in the NO basis) were
never searched by any census stage.

v96 RESULT (VERIFIED IN-REPO, exhaustive; supersedes the time-capped CP-SAT runs
of 9,056 and 6,395 skeletons). The structured enumeration
(scripts/signed_design_v96_full.py: tail-cover decomposition over the five
incidence-1 modes with a memoized head solve) covers the ENTIRE m=1 family. The
enumeration is exhaustive because the total weight is fixed at 9 (mode sums sum
to 27, three per determinant), so every support has <= 9 determinants; the tail
(the five incidence-1 modes, each in exactly one weight-1 det) plus the head
(weighted dets inside {0..4}) enumerate every integer m=1 skeleton. REPRODUCED
on this machine: 116,916 tail covers, 350,980 weight skeletons, ZERO hits, 44 s
(handoff run 41 s). Epistemic status per rung:
 - signs (rung 1): EXACT and EXHAUSTIVE, theorem-grade, VERIFIED IN-REPO -- v96
   admits NO signed design at m=1 (no extremal state with rational weight squares
   in its own NO basis and sign-only cancellation).
 - phase cancellation (rung 2-3): polygon filter exact-necessary; surviving
   skeletons excluded by multi-start numeric phase solve (12 starts, residual
   1e-10). High confidence, NOT certificate-grade. An exact algebraic exclusion
   would upgrade to: v96's extremal states have IRRATIONAL NO-basis weight
   squares, a structural first. Note the elegant self-test here: a class of three
   equal-magnitude terms cancels at cube-root-of-unity phases, which carry a Z/3
   holonomy -- so a phase-design hit would FALSIFY the exponent-2 holonomy law
   (Conjecture 2), and its absence is weak positive evidence for it.
Every searched family is now empty for v96 at m=1: 2x2 blocks (preflight None),
k-cliques with the one-hop cut (solver exhaust), single 3-clique + signed
off-clique cancellation (prototype slice), and the full signed/phased-design
family (this run). Remaining rungs: (a) the HYBRID family -- rational weights in a
block basis where designated block mode-pairs carry NONZERO off-diagonal targets
(Schur-Horn magnitudes) and all off-block classes cancel, the capped prototype
completed by the same tail-cover technique. The lone-pair funnel gives this
search a new exact skeleton gate: a lone one-hop class cannot cancel (this is
exactly the P1 prune in scripts/signed_design_fast.py), so in the hybrid it must
sit on a block pair, its value becoming that block's off-diagonal; hence BLOCKS
MUST COVER THE LONE-PAIR SET, and a skeleton whose lone pairs cannot be covered by
admissible block pairs (distinct-eigenvalue mixing, few blocks) dies before any
solve. The solver core is a POLYGON-TARGET SOLVER (BUILT: scripts/polygon_target.py):
the coupled one-hop classes with mixed targets (prescribed magnitudes in-block,
zero off-block), solved by triangle propagation (3-term classes pin relative
phases by the law of cosines, up to reflection branches), 2-term propagation, and
a bounded exact branch search over the remaining coupling; the stage-3b
off-diagonal exactifier is its single-class special case. (b) m=2 rational (mode
sums doubled; no interference precedent -- all 142 certified interference states
have state-den = spectrum-den, VERIFIED in-repo, the 9 den-doubling states being
DESIGN-REAL); (c) exactify rung 2-3. The same pipeline applies to the other five
signature-novel roots.

Polygon-target solver, verification (scripts/polygon_target.py, 2026-07). The
solver only PROPOSES exact symbolic phases; a state is returned only if
gpc_census.exactify.verify_exact certifies it (gauge-invariant char-poly
identity), so an incomplete propagation can only miss a solution, never certify a
wrong one -- soundness is architectural, the tests measure completeness. Two
constraint geometries, handled distinctly: a MAGNITUDE target |P|=X is one real
equation (the off-diagonal's own phase is free inside the block), reduced to a
single-unknown law-of-cosines step (the stage-3b exactifier generalized); a
CANCELLATION target P=0 is two real equations, so a 1-term class is infeasible
(the lone-pair funnel, identical to signed_design_fast.py's P1 prune), a 2-term
class pins its phase difference exactly, and a 3-term class is a rigid triangle
(law of cosines, two reflection branches). Empirically the census never needs
cancellation: across all 142 certified interference states every one-hop class
with terms carries a NONZERO Schur-Horn target (active-class arity 1/2/3 =
122/50/13, zero off-block classes), so the corpus tests the magnitude path only.
VERIFIED IN-REPO: `python scripts/polygon_target.py --recertify-all` re-solves all
142 from weights + spectrum alone (phases discarded), 142/142. The cancellation
path is tested on a constructed positive (the 4-cycle signed design
(|01>+|02>+|13>-|23>)/2, spectrum (1/2)^4, two 2-term classes cancelling) and the
funnel negative (tests/test_polygon_target.py). Arity <= 3 (all that occurs in
the census) is complete; arity >= 4 polygons carry internal freedom and fall to
the bounded branch search, not proven complete. Next: feed hybrid block/target
specifications for the v96 siblings (Task 1) into solve(..., targets=...).

## THE FILTER WAS UNSOUND: block path now filter-free, residual collapsing (2026-07)

The real root cause, three bugs deep. The block-ansatz selection-rule support
filter (admissible_support signature closure) is UNSOUND for a block target: it
is valid only when the 1-RDM is diagonal in the canonical basis (the 0-block
design probe), but a block target's 1-RDM is non-diagonal, and the closure
provably DROPS true-solution determinants. Verified directly: the census's own
exact v96 solution has 2 of its 7 support determinants OUTSIDE even the
class-merged closure. Merging classes across block pairs (the earlier fix) was
necessary but insufficient; the only sound choice is to enumerate ALL
determinants for block ansatze. FIX applied to both min_block_count and the
solve_vertex_exact_first block stage (adm = None / all determinants; the strict
filter remains only for the 0-block design probe, where it is valid). The clique
path already ran filter-free for exactly this reason.

RESULT (this SUPERSEDES the earlier "gate fix does not change the census output"
note): the census's OWN production pipeline now cracks the false negatives, no
external tooling. Each solve returns an exact certificate that is ALSO checked by
an independent from-scratch 1-RDM (scripts/verify_hybrid_state.py):
 - (3,10) v96 (den 9): min_block_count 1, SOLVED 7 s, indep-verified.
 - (3,10) v60 (den 12): SOLVED 3 s, indep-verified.
 - (3,10) v49 (den 13): SOLVED 42 s, indep-verified.
 - (3,10) v73 (den 14): SOLVED 103 s, indep-verified.
 - (3,10) v40 (den 18): SOLVED 159 s, indep-verified.
 - (3,10) v57 (den 28): SOLVED 1042 s, indep-verified -- a DEDICATED run past the
   300 s sweep cap; proves the high-denominator "timeouts" are compute limits, not
   FAILs, and that even den 28 is a false negative.
 - the test's own "off-family frontier" (9,6,5,5,5,2,2,1,1)/9: SOLVED 185 s (a
   false negative too; the min_block_count budget test was corrected accordingly).
Full filter-free sweep of all 14 (per-vertex 300 s cap; every OK independently
verified): SIX of the (3,10) roots SOLVED (v40 v49 v57 v60 v73 v96); v89 (den 26)
and v103 (den 34) hit the 300 s cap but v57 shows those are compute-bound, not
FAIL; the (4,10)/(4,9)/(5,10) roots were still running. Even (4,9) v42 -- the
supposed mixed-ansatz frontier -- flips to min_block_count 1 filter-free, so the
"v42 needs a mixed clique+block ansatz" conclusion is retracted pending its solve
verdict: likely another filter false negative, not a genuine family gap.

CAVEATS kept honest: (a) filter-free costs CP-SAT model size, not correctness --
verify_exact still gates every certificate, so no false positive is possible, but
high-denominator vertices are much slower and may hit the cap (inconclusive, not
FAIL). (b) min_block_count is now SOUND but a weaker gate: it rarely returns None
(all determinants make degree-feasibility easy), so its fail-fast value is
reduced. (c) the 785 previously-SOLVED states are unaffected (independently
certified); this only turns former FAILs into SOLVEs. The paper's residual count
must be revised DOWN to whatever survives the full sweep -- provisionally most of
the 14 are false negatives, not a genuine residual.

## CLASS-COUNT LAW: RETRACTED (two-class interference IS 2x2-solvable)

RETRACTION (found by a QA check the same day it was posted). The claimed law --
"the two holdouts v89/v103 are two-class, and two classes structurally require a
degenerate k x k block the 2x2 family cannot express" -- is FALSE. Counterexamples
in the shipped census: (3,10) v108 = (5^5,1^5)/10 and (3,8) v32 = (15^4,6^4) are
BOTH two-class INTERFERENCE states already SOLVED by the ordinary 2x2 family (v108
with an 8-determinant support and a pi/8 phase). A two-class vertex has one (top,
bottom) block ptype, which is a perfectly good 2x2 block (v96's (5,1) block is
exactly this) -- two classes do NOT preclude 2x2 interference. The only TRUE part
is the correlation within the residual: v89/v103 happen to be the two-class roots.
The MECHANISM ("needs k x k") does not follow and is refuted.

REVISED READING of v89/v103. v89 = (15^2,6^8)/26 has a valid single (15,6) block
ansatz (diagonal e.g. (10,11)/26, off-diagonal sqrt(20)/26, eigenvalues 15,6),
structurally identical to v108's (5,1) block one class-size up and one
denominator up. v108 (den 10) solved; v89 (den 26) did not. That is a DENOMINATOR
/ ENUMERATION-SCALING difference, not a family gap: at den 26 the census's
enumerate_weight_vectors (a 60 s enumerate_all_solutions CP-SAT pass) cannot emit
even one support that min_support_cardinality proves feasible, so the ansatz
contributes no support and the vertex reads as FAIL. So v89/v103 are almost
certainly COMPUTE ARTIFACTS (a third false-negative source, after the support
filter and the max_card cap), not genuine residuals. The fix is to make support
enumeration scale (incremental find-and-cut, or a targeted search seeded from the
known lower-denominator analog v108), NOT a new ansatz family.

kxk_degen.py (the degenerate 3x3 solver) stands as a correct, math-validated tool
for genuinely degenerate cases, but it is NOT the keystone for v89: v89 is
2x2-solvable in principle, so the degenerate-block detour is unnecessary for it.

Two supporting laws from the shipped cracks (verified): (1) ONE-CHANNEL -- every
shipped residual state (v96, v60, (4,9) v40/v42) uses exactly one exchange
channel on 6-7 determinants; the residual was minimal interference hidden behind
the filter bug, supporting Conjecture 3 where its falsification was predicted.
(2) CONJECTURE 2 clean -- all new phases are exponent-2 abelian: v96 cos -1/4 ->
2x^2+x+2 (Q(sqrt-15)); (4,9) v42 cos -sqrt(10)/20 -> 10x^4+19x^2+10; (4,9) v40
cos sqrt(3)/6 -> 3x^4+5x^2+3 (both biquadratic, V4/abelian); v60 loop-free. Zero
non-abelian holonomies. (A conjectured "scaling law", solve time ~1.32^den, is a
handoff extrapolation, not reproduced here.)

KEYSTONE TOOL: scripts/kxk_degen.py -- the degenerate 3x3 solver. Enumerates
eigenvalue sub-multisets E (size 3, >= 2 distinct, repeats allowed), Schur-Horn
integer diagonals majorized by E, mode-first DFS with one-hop pairs confined to
the three block pairs, exact integer-surd off-diagonal sums with core-only sign
enumeration, and the two exact 3x3 char-poly conditions
(sum s^2 = e2(diag) - e2(E); cos theta from e3). Both identities VALIDATED
in-repo (3000/3000 random symmetric 3x3). Every config reports EXHAUSTED or
TIMEOUT (the exhaust-vs-timeout discipline). It runs at any rank, so it also
attacks the two-class (3,11) open candidates (cand 23 = (6,6,1^9)/7, cand 44 =
(6^5,1^6)/12), where a hit would falsify a claimed level-5 inequality -- the
handoff reports cand 23 largely exhausted with two survivors under longer budget;
those verdicts are the handoff engine's, pending a dedicated in-repo run. The one
open decision for the whole rank-10/11 frontier is now this single extension:
solve a two-class vertex with a degenerate k x k block, or exhaust it there.

## RESIDUAL SWEEP RESULT: 8 of 11 roots solved; v42 retraction; v89 the frontier

Dedicated filter-free runs (per root: full max_card = total weight, max_blocks 2,
20 min budget, every SOLVE independently re-verified by the from-scratch 1-RDM).
The 11 independent roots (14 SOLVE-FAIL vertices minus the two (4,10) trailing-
zero pads of the (4,9) roots and the v113/v261 Hodge-dual merge):
- SOLVED and VERIFIED (8): (3,10) v40 v49 v57 v60 v73 v96; (4,9) v40 v42. State
  artifacts in docs/hybrid_cracks/ for v96, v60, 4_9_v40, 4_9_v42 (each passes
  both the independent 1-RDM check and shipped verify_exact). The (4,10) v60/v62
  pads and v261 dual follow from these for free.
- COMPUTE-BOUND, NOT FAIL (2): (5,10) v113 (den 18, N=5, 252-det model exceeds
  20 min) and (3,10) v103 (den 34, largest support space). v57 (den 28, solved at
  1042 s) is the precedent that these need only more walltime.
- GENUINE FRONTIER CANDIDATE (1): (3,10) v89 = (15,15,6,6,6,6,6,6,6,6)/26. It
  returns an EXHAUSTIVE FAIL (867 s, full max_card, filter-free), not a timeout,
  and it has only TWO eigenvalue classes {15,6}, so the block family (<= 2 blocks,
  no 3-clique possible) is COMPLETE for it. Status: UNRESOLVED, not confirmed
  genuine. The census phase solve is numeric (L-BFGS) and can miss a solution the
  exact solver would find; an exact-solver check (hybrid_search, 1119 s) was
  INCONCLUSIVE -- 0 solve calls, because its skeleton enumeration never surfaced a
  block-hop-bearing support in the den-26 space (the pre-filter rejected all
  106,400 skeletons it reached). So neither run settles v89: the census searched
  supports and found no numeric phase solution, but the exact solver has not yet
  been run on those same feasible supports. The decisive test is to extract the
  census's block-feasible-but-unsolved supports and exact-phase-solve each; until
  then v89 is the single open vertex, genuine-vs-numeric-miss undetermined.

RETRACTION: the "(4,9) v42 needs a MIXED clique+block ansatz -- a genuine family
gap" conclusion (the v42-audit and no-clique-precedent sections below) is WRONG.
v42 solves via the plain filter-free block path as an ordinary interference state
(phase exp(i(acos(-sqrt(10)/20)+pi))). The mixed-ansatz-gap narrative was entirely
downstream of the support-filter bug; strike it. The clique machinery was never
needed for any residual vertex reached so far.

NET FOR THE PAPER: the residual of 14 (11 independent) is NOT a genuine residual.
At least 8 of 11 roots are false negatives (verified), 2 more are compute-bound
(expected to fall with walltime), and only v89 is a live candidate for a true
residual pending its exact-solver verdict. Do not report 14/11; report the
verified solved set and v89 as the single open frontier.

## v57 dive: why den 28 resists, and the core+completion algorithm (handoff)

v57 = (19,19,10,10,6,6,6,6,1,1)/28 is the first residual vertex where
filter-free enumeration genuinely saturates: the production sweep hit its 300 s
cap (inconclusive, verified here), and a dedicated 20-minute filter-free run was
also launched. The handoff reports two bespoke DFS implementations (det-ordered
with in-search off-block pruning, and mode-first branching over an off-block
conflict graph) both exhausting 9 s/ansatz budgets across all one-block ansatze;
the obstruction is the weight budget (Sigma k = 28), not the ansatz count.
Structure map, reproduced in-repo: 27 one-block ansatze (VERIFIED), and exactly
two incidence-1 modes 8, 9 (VERIFIED, the two 1/28 modes) that survive in every
ansatz not blocking the 1-class, forcing one weight-1 det each (exploitable as in
the tail-cover enumerations).

CORE+COMPLETION (handoff algorithm proposal, not yet implemented): factor the
search into (i) an EXCHANGE CORE -- the block-pair hop dets and their weights,
constrained by number theory (the signed surd sum over cores must equal
sqrt(x2) exactly, so core products live in the squarefree class of x2; cores are
tiny), then (ii) a DESIGN COMPLETION -- the residual degree vector met by a
one-hop-free support avoiding off-block hops with the core, which is the design
MIP the census already solves fast, plus adjacency exclusions. Enumerate the
(number-theoretically bounded) cores, MIP-complete each: (small) x (fast) in
place of one intractable joint search. In parallel: 2-block ansatze and cliques
under real walltime, and hybrid_search (which handled a 540 s pass at v96). (The
handoff's v57_solver.py prototype is not in this repo; core+completion is a
design spec here, not shipped code.)

## v60 exact state shipped (all-positive form); conjecture updates from the cracks

docs/hybrid_cracks/v60.jsonl: the (3,10) v60 = (8,8,5,5,5,1,1,1,1,1)/12 state
from the fixed production sweep, gauge-reduced to its minimal form and
independently certified (both the from-scratch 1-RDM and the exact char-poly
identity pass, scripts/verify_hybrid_state.py). The support is LOOP-FREE, so all
sweep phases were pure gauge: the state is ALL-POSITIVE, psi = (sqrt2 |012> +
sqrt3 |029> + |035> + |048> + |067> + 2 |349>)/sqrt(12), with a single exchange
channel (1,9) realizing the (8,1)-block off-diagonal sqrt(6)/12.

Conjecture ledger as the residual collapses (update per new crack):
- Conjecture 2 (exponent-2 holonomies): v96's single loop carries cos = -1/4,
  minimal polynomial 2z^2 + z + 2, discriminant -15, splitting field
  Q(sqrt(-15)), abelian Z/2 -- PASS (verified in-repo). v60 contributes no
  holonomy (loop-free) -- vacuous PASS. Each further crack: compute the loop
  kernel and minimal polynomial before the record lands.
- Conjecture 3 (<= 2 channels): v96 one channel, v60 one channel -- both PASS,
  strengthening the bound precisely where its falsification was predicted.
- The "survivors need >= 3 channels" prediction is now FALSIFIED (it was
  downstream of the filter bug); strike it wherever quoted.

## v96 SOLVED: a census false negative (2026-07)

(3,10) v96 = (5,5,5,5,2,1,1,1,1,1)/9, recorded SOLVE-FAIL, IS SOLVABLE. Exact
closed-form extremal states exist; 24 were found in a single 540 s hybrid pass
and all 24 pass BOTH an independent from-scratch 1-RDM spectrum check and the
shipped verify_exact gate (scripts/verify_hybrid_state.py; stored in
docs/hybrid_cracks/v96.jsonl; regression guard tests/test_v96_crack.py). A
representative state, weights (1,1,1,1,1,2,2)/9 on dets (1,2,6)(1,3,7)(1,8,9)
(0,1,4)(1,4,5)(0,2,3)(2,3,5), carries a single interference phase
exp(-i acos(-1/4)) and reproduces the spectrum exactly.

They live on the DEGENERATE (5,1) BLOCK ansatz: modes 0 and 5 both carry
occupation 3/9, and their 2x2 block [[3/9, 2/9],[2/9, 3/9]] has eigenvalues
5/9, 1/9. The 1-RDM is block diagonal (the supports are off-block-hop free), so
this is NOT the signed-cancellation extension -- it is the plain block ansatz
gpc_census.states.block_ansatze already generates (ptype (5,1), split (3,3),
x2 = 3*3 - 5*1 = 4, target |off-diagonal|^2 = 4/81, exactly the value the census
(3,10) interference blocks carry).

ROOT CAUSE, traced and FIXED at source (2026-07; this supersedes two earlier
readings in this file's history -- the "min_block_count preflight" note and the
"selection-rule collapse, grouping irrelevant" note, both wrong). The admission
gate min_block_count took the selection-rule signature closure over the
lambda-DEGENERATE (equal-value) classes only. A block that mixes two DISTINCT
classes (the (5,1) block: a 5-mode with a 1-mode) rotates the natural basis
across the class boundary, so the rotated support carries occupancy signatures
the within-class closure drops; the gate found the CP-SAT model infeasible and
returned None, and solve_vertex_exact_first bailed. FIX (src/gpc_census/states.py,
_block_merge_groups, wired into both min_block_count and solve_vertex_exact_first):
close over groups that MERGE the classes each block's mode pair joins. v96 is
highly saturated (11 of 93 facets active, strict set size 1 vs 16 for a solved
vertex like v17), and the earlier note wrongly concluded the size-1 strict set
was a hard wall -- but the signature closure over the MERGED groups expands from
it to a feasible admissible set (size 84 for the (5,1) ansatz). VERIFIED: the
shipped test suite passes (64) and the fixed gate re-audit flips 12 of 14 to
feasible (below).

IMPORTANT SCOPE: the gate fix is necessary but NOT sufficient. Gate feasibility
is a lower bound; the census's own SOLVE path (per-ansatz min_support_cardinality
plus the numeric L-BFGS phase solve and exactify) still returns FAIL on v96 even
with the gate open -- a further CP-SAT inconsistency at the single-block ansatz
plus the numeric solver's limits. So fixing the gate does NOT by itself change
the census SOLVED/FAIL output. The actual exact states come from the
direct-enumeration polygon-target solver (scripts/hybrid_search.py), which uses
no selection-rule prune and phase-solves exactly; that is what cracked v96.

FIXED-GATE RE-AUDIT of all 14 (max_blocks=2, time_cap 8 s; every number below
reproduced in-repo against the patched gate): TWELVE now report a feasible block
ansatz -- (3,10) v40:1 v49:2 v57:1 v60:1 v73:1 v89:1 v96:2 v103:1; (4,10) v60:1;
(4,9) v40:1; (5,10) v113:2 v261:2. Only (4,9) v42 and its (4,10) v62 pad remain
None. So the gate was a real false-negative FILTER on most of the residual.

IMPLICATION, stated carefully: the residual of 14 (11 independent) is an UPPER
BOUND and v96 is a confirmed crack, but the residual does NOT collapse to 2, and
gate feasibility must not be read as solvability. A direct-enumeration sweep
(hybrid_search, the real solve) finds v40 and v49 RESIST even a deep max_blocks=2
slice (669k and 2.1M skeletons, no hit) despite both being gate-feasible -- so
"gate artifact" explains the mislabeling but not the hardness; some residual
vertices look genuinely hard. High-denominator vertices (v57/28, v89/26, v103/34,
v113/v261/18) are not meaningfully searched yet (weight = den pushes supports to
26-34 determinants). v42/v62 is the surviving gate-frontier candidate for a true
family gap. This does not touch the 785 SOLVED states (each independently
verify_exact-certified) -- only the FAIL labels are suspect. Do NOT cite 14/11 as
final: v96 is SOLVED; the rest stay open pending a proper direct-enumeration
sweep plus phase solves, not the gate re-audit alone.

v42 audit (the surviving frontier vertex; handoff analysis, key claims
reproduced in-repo). Unlike v96, (4,9) v42's persistence is NOT a gate artifact:
its block path is exactly infeasible through 2 blocks (min_block_count = None
with BOTH fixes in, a 3 s verdict not a timeout), while its clique path is
gate-feasible (min_clique_count = 3: a single 3-clique ansatz is degree
feasible). So v42 passes the clique gate and fails only in the phase solve --
structural or compute bound, not gating. STRUCTURAL GAP located and verified: the
solver enumerates block-only ansatze (block_ansatze, 2x2 blocks) OR clique-only
ansatze (_solve_via_cliques over multi_clique_ansatze, disjoint cliques of size
>= 3), and NEVER a MIXED configuration (a 2x2 block together with a 3-clique).
v42 may need exactly such a mixed ansatz -- a genuine family gap, not a bug. The
polygon-target solver can test mixed block+clique targets directly once
hybrid_search is extended to emit them; that is the recommended next attack on
v42. (The handoff's specific "35 feasible / 5 infeasible" single-clique tally was
NOT reproduced here: clique_ansatze(sizes=(3,)) yields 374 ansatze, so that count
is from _solve_via_cliques's internal filtered enumeration, not independently
confirmed; the gate-clean conclusion stands on min_clique_count = 3.)

NO GENUINE CLIQUE PRECEDENT (verified 2026-07, and it gates the whole clique
attack). Across all 785 solved states the MAXIMUM number of active off-diagonal
edges is 2, and ZERO states have three mutually-active edges (a genuine
triangle / 3-clique): every solved interference state is 0, 1, or 2 disjoint 2x2
pairs. So the census's clique machinery, though present in code, has never
produced a solved state, and there is NO positive control for a genuine k>=3
clique solver. Consequence for v42: a clique/mixed solver can be built (the k x k
realization primitive, a Chan-Li Givens sweep giving one symmetric matrix with
the target diagonal and eigenvalues, was validated on 2836 random majorized
cases), but without a genuine-clique control a NULL result on v42 is
uninterpretable (no-solution vs incomplete solver vs too-slow), while a HIT stays
independently verifiable via the from-scratch 1-RDM check. A first-pass clique
driver was prototyped and deliberately NOT committed: it fails its own self-test
(the chosen control (3,8) v20 turned out to be a single (3,7) pair, not a
3-clique, and surd-magnitude recognition is incomplete), and high-denominator
vertices like v42 hit the same enumeration-scaling wall as v57/v89/v103. v42 is
genuinely novel territory: a real clique+mixed exact solver (surd-exact magnitude
targets, realization-family search, scalable enumeration) is a substantial
dedicated build, not a quick extension.

Generalized enumerator (scripts/signed_design_generic.py): the same three-rung
search for an ARBITRARY (N,d) integer spectrum, exhaustive DFS over determinants
ordered so the tightest (smallest-budget) modes bind first, with feasibility
pruning and an optional walltime. CROSS-CHECK VERIFIED IN-REPO: on the v96
spectrum (--ints 5,5,5,5,2,1,1,1,1,1 -N 3 --no-phases) it reproduces the
specialized run EXACTLY -- 350,980 skeletons, 0 hits (79 s here; slower than the
v96-specialized tail-cover/head-memo path, same answer). The five sibling roots
are addressable directly by their integer forms (docstring --ints lines, each
verified to match its vertex): v49 (9,9,5,5,5,2,1,1,1,1) den 13; v57
(19,19,10,10,6,6,6,6,1,1) den 28; v89 (15,15,6,6,6,6,6,6,6,6) den 26; v103
(18,18,18,18,5,5,5,5,5,5) den 34; (5,10) v261 (14,14,14,14,10,10,10,2,1,1) N=5
den 18. Roots with no incidence-1 modes (v89, v103) have a much larger skeleton
space; run under --max-seconds and read the frontier line (a capped run is a
documented partial slice, not a silent cap).

Fast variant (scripts/signed_design_fast.py): same enumeration and same hit
definitions as the generic version, with three EXACT per-skeleton prunes -- (P1)
a one-hop class with a single term can never cancel; (P2) a two-term class
cancels only if the two surd magnitudes match (m1==m2 and q1==q2); (P3) each
surviving two-term class is a linear GF(2) constraint on the sign bits, so
Gaussian elimination replaces the 2^(n-1) sign sweep, and the remaining >=3-term
classes are checked on the reduced solution set. All three are exact necessary
conditions, so the answer is unchanged. VERIFIED IN-REPO two ways:
tests/test_signed_design_fast.py is a differential guard asserting fast == generic
(skeleton count AND hit set, signs and phases) across six small spectra chosen so
the design/signed-design cases actually carry hits -- the guard that matters,
since v96's zero-hit result cannot exercise a wrongly-dropped hit; and the v96
full run reproduces 350,980 skeletons / 0 hits (60 s here vs the generic 79 s).
Use it for the larger sibling roots (v89, v103) under --max-seconds.

## Depth triage score: back-testing the facet laws on the open vertices

Define the TOTAL ACTIVE DEPTH of a vertex: the sum of rhs depths over its
saturated GPC facets (`python scripts/facet_laws.py --triage N D`; exact, no
solving; compare only within a system and within the interference class, since
design corners carry huge degenerate facet counts and are exact by construction
anyway).

- STRATIFICATION (distribution shift, NOT clean separation; reproduced here):
  (3,10) certified interference median 12 vs the six clique-family opens at
  17-34; the top scores INTERLEAVE (certified v90 = 41 SOLVED, open v96 = 34
  with eleven active GPCs, certified v70 = 33, v86/v77 = 31, all EXACT). So
  depth is a distribution shift, not a threshold: v96 is not even the deepest
  interference vertex, only the deepest OPEN one, and the whole system reaches
  depth 304 (a degenerate design vertex, v0). (4,9): opens v40/v42 at 11/21 vs
  certified v_B/v30 at 6/6. v89/v103 sit at 3-4 (one active GPC each), the
  opposite extreme. The opens are the depth-extremes of the OPEN set in both
  directions, not of the polytope.
- CALIBRATION AGAINST SOLVE COST (handoff, NOT reproducible from shipped data):
  the handoff reports Spearman(total depth, Tier-A secs) = -0.625 over certified
  (3,10) interference vertices, median 8 s at depth > 12 vs 147 s at depth <= 12
  -- high constraint load makes the family search DECISIVE (fast solve when a
  family solution exists, e.g. v90 at depth 41 solved; fast exhaust when not,
  e.g. v96 at 34), near-zero load is under-constrained blowup (v89/v103).
  states.jsonl stores no Tier-A solve times, so this correlation is the
  campaign's and is NOT re-run here; ship the solve-time logs alongside the
  states if it is to reach the paper.
- TRIAGE RULE (proposed): compute total active depth (one exact pass, no
  solving) at any rank. Near-zero => search-bound; raise max_card / walltime.
  High relative to the system's certified band => decisiveness expected; if the
  standard sweep exhausts quickly, route to the signed-cancellation extension
  without more hardware. LIMITS: no cross-system threshold (scales differ);
  exhaust-vs-compute-bound within the high band is not depth-determined on
  current evidence (v96 depth 34 exhausts, v42 depth 21 grinds).
- Quantified negative kept for the record: the naive position-space
  selection-rule filter (support inside the active facets' canonical realizer
  sets) fails 0/34 on certified (3,10) interference states; stored supports
  carry arbitrary mode labels, so this conflates label permutation with genuine
  rotation, and the degeneracy lemma's v_B diagnostic remains the clean witness.

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
- partial_families_3_11_3_12.json: the published stable Grassmann families for
  N=3 at ranks 11 and 12 (the only rank-11/12 constraint data in print), tagged
  PROVED (CMP 2008 Thms 4.2.1/4.3.1), WEAK (thesis odd-rank head), and CLAIMED
  (the level-5 series, proof deferred / letter only). NECESSARY conditions, an
  outer bound only. Regenerated and validated by scripts/partial_families.py:
  all certified (3,11) and (3,12) true vertices satisfy every family, the proved
  families restrict validly to {lambda_11=0} on all 113 (3,10) vertices, and the
  four OPEN candidates each violate a CLAIMED (never a PROVED) inequality. Use as
  a Farkas-style acceptance oracle for any future Stage-1 (3,11)/(3,12) output.
- bracket_3_12_true_vertices.json: 19 CERTIFIED true vertices of (3,12),
  each a face embedding of a settled (3,11) true vertex (the {lambda_12=0}
  face of Delta(3,12) is Delta(3,11); a vertex of a face is a vertex of the
  polytope). Rigorous and reproducible from the (3,11) settlement; the
  first vertices of a rank-12 system.
- bracket_3_12_stage0_inner.json + scripts/plethysm_inner_hull.py: the (3,12)
  Stage-0 INNER cloud, generated WITHOUT Stage 1 by exact plethysm (multiplicity
  of s_lambda in h_m[e_3] via Murnaghan-Nakayama, exact Fractions; every emitted
  point is a certified attainable spectrum). VERIFIED (2026-07): the generator's
  own --validate gate reproduces the (3,6) vertex set exactly at M<=6 and finds
  all ten (3,7) vertices by M<=10 (tests/test_plethysm.py guards the fast (3,6)
  case); an independent cross-check confirms all 88 M<=6 points satisfy the
  proved partial Grassmann families (partial_families oracle, zero violations)
  and that the 9 face-embedded (3,12) vertices reachable at M<=6 all appear (the
  other 10 have denominator >= 7 and need larger M, the documented convergence
  behavior). CAVEAT: cloud-extremality (35 exact-certificate extreme points, 2
  numerically-unverified candidates at M<=6) is NOT true-polytope vertexhood
  until the cloud converges (M >= max vertex denominator); the rigorous vertices
  are the face-embedded ones above. The plethysm generator complements Stage 1
  (inner side), it does not replace it.

## Adjacent literature to respect

Liebert, Lemke, Altunbulak, Maciazek, Ochsenfeld, Schilling,
PRResearch 7, 023247 (2025), arXiv:2502.15464: spin-adapted GPCs for
larger systems. Different lattice of problems (spin sectors) from the
spinless census here; draw the boundary explicitly in any new paper.
Maciazek and Tsanov (J. Phys. A 50, 465304): doubly-excited inner bounds
coincide with the polytope for (3, d <= 7); the census confirms this
empirically (those systems are interference-free).
