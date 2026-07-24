# THEOREM PROGRAM (response to external expert reviews, 2026-07)

## T0 -- The frame (upgraded by the second review, with corrections)
FACT, not conjecture (the reviewer's "even stronger conjecture" is the
literature's front door): gamma IS the moment map for the U(r) action on
P(Lambda^N C^r); the GPC polytope IS the Kirwan polytope; Klyachko's
construction runs on exactly this. Consequence: Atiyah-Guillemin-Sternberg,
Kirwan convexity, and Sjamaar-Lerman singular reduction apply verbatim --
paper 2's opening section cites rather than proves.

GENERIC FIBER THEOREM (reviewer's conjecture; adopt, route = standard
transversality): over regular values (interior, nondegenerate spectra) the
sorted-spectrum map is a submersion and fibers are smooth of constant
dimension dim P(H) - (r - 1). VERIFIED numerically at a generic (3,10)
point: rank D(spec) is full; fiber dimension 229. Consistent with the
measured valleys.

THE CORRECTION (our Slater/v60 computation, decisive): the reviewer's
Step-3 formula "dim F = nullity of D Gamma" holds ONLY at regular values.
The census's entire subject matter -- vertices -- consists of SINGULAR
values, where the formula fails in the direction opposite the reviewer's
Step-4 intuition: the nullity is LARGE (ker >= 2 nd - d^2, e.g. 148 at
v60) while the actual fiber is SMALL (a circle at Slater vertices) -- the
excess kernel is a second-order OBSTRUCTED CONE, and dimensions at deep
strata DROP, they do not jump. Designs are not "singular fibers with
jumping dimension"; they are the deepest strata where the REDUCED fiber
collapses to a point. The correct vertex-fiber statement is stratified
(Sjamaar-Lerman), and the constructive families are its unobstructed
sector -- see T1.


Goal: invert the 70/30 numerical/theoretical ratio. Each target theorem below
carries an honest status and its proof route. Numerical facts cited are exact
or machine-precision and live in RESEARCH.md.

## T1 -- Fiber Dimension Theorem (the centerpiece)
Statement (target): for an extremal spectrum v and a certified state psi over
it, the REDUCED fiber (states with gamma(psi') = gamma, modulo the stabilizer
U(m1) x ... x U(mk) of gamma and global phase) has dimension
   dim F_red(psi) = dim ker(incidence map of supp psi) - #(touched rigid
   1-term classes),
except at strata where higher-order obstructions activate.
STATUS:
- LOWER BOUND: PROVED CONSTRUCTIVELY for the census corpus -- the kernel
  deformation families ARE curves/surfaces in the fiber (44 + 15 explicit
  families; phase re-solve = open polygon condition).
- UPPER BOUND: OPEN, and provably NOT a first-order rank computation:
  at a Slater vertex the fiber is a circle while dim ker Dgamma >= 2*nd - d^2
  (grossly larger) -- vertices are SINGULAR VALUES of the moment map and
  kernel directions are generically second-order obstructed (adding epsilon
  on an off-support determinant moves gamma only at O(epsilon^2)).
- PROOF ROUTE: Marle-Guillemin-Sternberg local normal form at the orbit,
  in the Sjamaar-Lerman stratified-reduction setting; the obstruction
  quadric is computable (the second fundamental form of the gamma-map),
  and on the corpus it can be VERIFIED exactly per state (finite check:
  which kernel-of-Dgamma directions survive the quadratic obstruction).
  Executed lesson (this session): the naive full-space rank test at v60
  gives ker = 148, orbit = 30 -- the 118 excess is the obstructed cone,
  not moduli. The support-restricted kernel calculus we have used all
  along is precisely the unobstructed sector.

## T2 -- Holonomy formulation (referee-proof Berry phase)
Statement (target): the reduced fiber of a single-2-term-class +-1-loop
family is a smooth circle S in projective Hilbert space; the tautological
U(1) bundle restricted to S carries the Berry connection; its holonomy is
   gamma_B = -oint k_L(t) / den  dTheta
(k_L = the loop-determinant weight), parametrization-independent and
gauge-invariant mod 2pi once the lift convention (phase carried on the loop
determinant, support basis fixed) is stated. v_B value: 1.3139 pi mod 2pi.
STATUS: formulation complete above; smoothness of S = the conic model
(proved); the integral in closed form on the conic = elementary (OPEN,
one afternoon); contractibility in CP^H: S bounds a disk (CP^H simply
connected), so gamma_B = integral of Berry curvature over any bounding
disk -- making it an honest geometric invariant of the embedded circle.

## T3 -- Cumulant Separation Theorem
Statement (target): every positive-dimensional reduced fiber contains states
with distinct two-body cumulants; concretely, along any kernel deformation
some pair occupation varies linearly with nonzero slope.
STATUS: reduces to "no deforming family kernel is annihilated by the
pair-incidence map" -- an exact rational rank check per family, queued
corpus-wide (verified at v_B: 8 varying pairs). If a counterexample family
exists, the theorem gains an interesting exceptional class; either outcome
is a result.

## T4 -- Reality Selection (per-family: PROVED; global: conjecture)
Per-family endpoint theorem proved (tau > 0 forces |cos|=1 walls interior).
Global version (fiber-wide extremes of pair-linear functionals attained by
real states): verified numerically at v_B (gaps < 1e-6); MECHANISM ROUTE:
extreme points admit time-reversal-symmetric representatives (conjugation
preserves gamma and all pair occupations) -- candidate two-line proof via
averaging; assigned.

## T5 -- Quasi-fiber stability
Valley persistence off-vertex: measured (widths 0.25 -> 0.29). Theorem via
incidence surjectivity + open phase conditions; assigned (Session 3F).

## T1-NUMERICAL: the reviewer's Jacobian test, executed exactly (2026-07)

Proposal (second review): compute Jacobian nullity at certified states,
compare with observed fiber dimension; designs should have smaller kernels.
Executed with two refinements the vertex geometry forces, each taught to us
by a failed prediction:

STEP 1 (matrix map, support-restricted): nullity - gauge gave MODULI = 0 at
ALL of design/v_B/v89/v96 -- including states with CERTIFIED families. Not a
contradiction: the families preserve the SPECTRUM, not the matrix (weights
move, the degenerate-block frame rotates along the coadjoint orbit). The
matrix-fiber and spectrum-fiber are different objects; certified families
live in the latter.

STEP 2 (spectrum map = commutant-projected Jacobian; at simple spectra this
IS the reviewer's D Phi, so the two regimes unify): exact results --
  DESIGN (3,10) v0:  moduli 0   (predicted 0)  PASS
  v89 silent-channel: moduli 0  (predicted 0)  PASS -- the zero-target
     class violation appears at FIRST ORDER in the projected map: a
     differential-geometric proof-sketch of the state's rigidity, i.e. the
     mechanism behind pre-registration P4's failure, now theorem-shaped.
  v_B:  moduli 2  (predicted 1)
  v96:  moduli 2  (predicted 1)
The +1 excess at deforming states is identified: a within-support coadjoint
direction (first-order coherence between DISTINCT eigenvalue blocks,
representable inside the support) -- orbit motion the spectrum-fiber
contains but the REDUCED fiber quotients. The invariant hierarchy is
therefore three nested objects, all now computable exactly per state:
  matrix-fiber  (fixed gamma)         -- moduli 0 at all tested vertices;
  spectrum-fiber (fixed lambda)       -- gauge + orbit + families;
  REDUCED fiber (mod stabilizer/orbit) -- the physical moduli = the
     certified families (final quotient pass: subtract within-support
     orbit rank; queued for the corpus-wide run).
Reviewer's headline prediction CONFIRMED in the corrected frame: designs
(and the new silent-channel class) have vanishing moduli at first order;
interference deforming states have positive moduli. The design/interference
distinction is now differential geometry, not just classification.

## T1-NUMERICAL II: the linear hierarchy REFUTED both ways; a v_B anomaly opened

The proposed final quotient (reduced moduli = spectrum-fiber moduli minus
within-support orbit rank) FAILS in both directions, exactly (all ranks
rational/algebraic):
  DESIGN: 0 - 0 = 0 vs families 0  (consistent)
  v89:    0 - 0 = 0 vs families 0  (consistent)
  v_B:    2 - 0 = 2 vs families 1  (OVERCOUNT: neither extra direction is
          within-support orbit motion -- orbit excess is exactly 0)
  v96:    2 - 2 = 0 vs families 1  (UNDERCOUNT: the certified family's
          first-order tangent lies ALONG orbit directions; the family
          escapes the orbit only at second order)
LESSON (now demonstrated three independent ways): at singular strata no
sequence of linear quotients computes the reduced moduli -- curvature mixes
orbit and moduli directions, and first-order kernels both over- and
under-count. The MGS local normal form is not a nicety; it is the only
correct tool. The Generic Fiber Theorem (regular values) is untouched.

OPEN ANOMALY (the sharpest new thread): v_B carries a SECOND in-support,
non-orbit, spectrum-preserving first-order direction beyond the conic
family -- confirmed by exact rank AND float SVD (nullity 9 = gauge 7 + 2),
and a predictor-corrector integration along it reaches exact-spectrum
states (residual ~1e-25) drifting slowly off the known family's weight
profile (0.0007 -> 0.0031 over six steps). Candidate explanations, ordered:
(a) an UNDISCOVERED SECOND FAMILY through psi_B (the fiber richer than the
conic model -- would be a discovery); (b) a first-order direction that is
second-order obstructed, with the corrector tracking a curved shadow of
the known family (the small drift is suggestive but not conclusive);
(c) a subtle error in the commutant projection at the degenerate 4+4
blocks. Deciding requires the careful session: exact second-order
obstruction analysis along direction 2, and/or high-precision arclength
continuation with orbit-invariant distance measures. Assigned to the rig
queue with priority just below v103 scoring.

## v_B ANOMALY RESOLVED: the reduced spectrum-fiber is a SURFACE; the conic
families are its matrix-fiber slices (2026-07)

The decisive tests (constrained-corrector scaling + 2-RDM-spectrum
invariant, this session):
1. With the corrector forbidden from the gauge/family-tangent/direction-2
   span, exact-spectrum states persist at machine residual (1e-21..1e-24)
   at every displacement -- NOT the eps^4 scaling of a second-order
   obstruction. Direction 2 integrates.
2. The endpoint state's weights move the two loop pairs INDEPENDENTLY
   (deltas -+0.523 and -+0.368) -- impossible on the 1-dim incidence
   kernel -- so its 1-RDM is NON-DIAGONAL: diagonal and off-diagonal
   co-vary preserving eigenvalues. A single-phase rebuild at the new
   weights certifies the spectrum at 1e-21.
3. The unitary-invariant discriminator: its 2-RDM spectrum matches NO
   family point (min distance 6.2e-3, ~60x the weight-rounding noise) --
   a genuinely new REDUCED point, continuously connected to psi_B.

CONCLUSION (framing matters): not a second component -- a SECOND
DIMENSION. The reduced spectrum-fiber through psi_B is locally (at
least) a surface; the censused conic families are exactly its
FIXED-DIAGONAL slices (the matrix-fiber moduli). The fiber-dimension
census, the wall census, and the conic model all remain correct AS
matrix-fiber statements, and every physics claim already used fixed
full gamma (the scrutiny repair), so nothing retracts. What opens is a
new layer: spectrum-fiber directions where the natural-orbital frame
and the occupations co-rotate. Expect every deforming vertex to carry
spectrum-fiber dimension >= matrix-fiber dimension + 1 (conjecture;
same test, corpus-wide, queued). The reviewer's Local Structure Theorem
should be stated for the spectrum fiber, with the matrix fiber as the
commutant-diagonal stratum -- and since gamma is exactly quadratic, the
MGS normal form TERMINATES: the obstruction analysis is exact, with the
obstruction quadric equal to gamma(delta) plus degenerate-PT feedback.
This is the theorem's worked example, ready.

## T6 — THE GENERALIZATION QUESTION: fiber geometry vs face combinatorics
(response to the reviewer's "inverse-geometry principle", 2026-07)

Their bold statement -- "the combinatorics of a moment polytope are
encoded in the geometry of its moment-map fibers" -- splits cleanly on
the abelian/nonabelian line, and the census decides the nonabelian half:

1. ABELIAN (toric) CASE: the principle is a CLASSICAL THEOREM (Delzant).
   Fibers are tori; dimension drop equals face codimension; face lattice
   and fiber degeneration determine each other exactly. Gelfand-Tsetlin
   is (via its toric degeneration) in this regime -- so testing there
   would confirm trivially, not discriminate. GT is the NEGATIVE CONTROL,
   not the test.
2. NONABELIAN (our case): the naive equivalence is FALSE BY EXHIBIT.
   The census contains 183 simple vertices -- combinatorially the same
   local face structure -- splitting as 128+ simple DESIGN vertices
   (rigid fibers) and 51 simple INTERFERENCE vertices (including the
   deforming families). Same local combinatorics, radically different
   fiber geometry. Arithmetic agrees: the den = normal-cone-index law
   holds 51/51 at simple interference and fails at 68 of the simple
   designs -- quantities split along the STATE class, not the face class.
   Conclusion: the face lattice does NOT determine the fibers.
3. THE CORRECT (evidence-backed) STATEMENT is the refinement, running
   the other way: fiber invariants define a stratification of the
   Kirwan polytope STRICTLY FINER than the face lattice -- a "quantum
   stratification" (design/interference, rigidity, silent channels,
   fiber dimension) invisible to convex geometry -- with partial forward
   control (fiber dim 2 forces simple-or-excess-1, 15/15). Fibers know
   the faces AND MORE; faces cannot recover the fibers.
4. THE DISCRIMINATING EXPERIMENT (replacing their GT proposal): port the
   fiber census to a genuinely nonabelian, independently-understood
   moment polytope -- the natural candidate is a small ENTANGLEMENT
   POLYTOPE (Walter-Doran-Gross-Christandl), where vertices classify
   entanglement types. If the refinement phenomenon recurs (same face,
   different fiber classes), the quantum stratification is a general
   feature of nonabelian moment maps and the method is a general
   instrument. GT alongside as control, where the refinement must
   collapse to the toric equivalence.
5. Adopted claim (their phrasing, endorsed): the complete GPC census is
   a uniquely well-characterized testbed for fiber-based methods on
   moment polytopes. No stronger claim until the entanglement-polytope
   experiment reports.

## T6 PROOF PROGRAM (how the quantum stratification becomes theorems)

THEOREM A (strictly-finer; exhibit + definitions): design/interference
are vertex-level certified invariants; simple vertices have pairwise
combinatorially isomorphic vertex figures (simplicial normal cones); the
183-simple census supplies same-system design/interference pairs. Hence
the classification does not factor through the local face lattice.
Status: WRITING TASK (days). Requires: pin the definitions, name one
certified pair, cite both certificates.

THEOREM B (design rigidity) -- FINITE PART PROVED THIS SESSION:
Reduction: a design support has no one-hop pairs, so all phases are
gauge and the only candidate deformations are incidence-kernel weight
shifts, which preserve every mode sum and hence gamma itself -- giving a
continuum of designs iff the kernel is nonzero. VERIFIED EXACTLY: all
643 certified design supports have kernel dimension ZERO (rational rank,
corpus-wide). Therefore every certified design state is rigid within its
support, by proof. Remaining open sliver: whether design supports are
kernel-free NECESSARILY (a design vertex with a loopy design support
would carry a design continuum -- not excluded abstractly; in-corpus it
never occurs). Either resolution is publishable.

THEOREM C (fiber-dim 2 => excess <= 1; 15/15 observed): corollary of the
capstone -- the MGS slice representation's weight count bounds fiber
dimension against active walls. Priced into T-MAIN, not separate.

T-MAIN (the capstone, with the reviewer): the Local Structure Theorem in
algorithmic form. Because gamma is EXACTLY quadratic, mu(psi + delta) =
mu(psi) + d mu(delta) + mu(delta) with no remainder, so the reduced
fiber germ at any certified state is the solution set of an explicit
quadratic system computable from (support, weights, classes) -- the
theorem is an algorithm plus a correctness proof in the MGS frame. The
v_B surface, the silent-channel rigidity, and the conic families are its
first three worked corollaries, already computed.

TRANSFER PREREQUISITE: the portable (category-free) definitions --
e.g. "design" as: the fiber meets the toric stratum (states whose
moment image already lies in the Cartan). Needed before the
entanglement-polytope experiment can be stated cleanly.

Division of labor: A + B-sliver + algorithm implementation = us/rig;
MGS formalism = the reviewer's frame; finite verifications = CI.

## T6 FINALIZED: the Incompleteness Principle, the rigorous definition,
and the transfer experiment's pilot (2026-07)

ADOPTED (reviewer's formulation, verbatim):
  CONJECTURE (Combinatorial Incompleteness of Nonabelian Moment
  Polytopes). For generic nonabelian moment polytopes, the face lattice
  is not a complete invariant of the geometry of the associated reduced
  fibers. Consequently there exists a strictly finer stratification of
  the polytope by reduced-fiber type.
Sharpened headline sentence (adopted): THE FACE LATTICE IS NOT A
COMPLETE INVARIANT OF THE UNDERLYING QUANTUM GEOMETRY. The GPC census is
the first substantial evidence body (183 simple vertices, two fiber
classes at one combinatorial type).

RIGOROUS DEFINITION (reviewer's, adopted): a QUANTUM STRATIFICATION of a
Kirwan polytope is a decomposition into strata whose points possess
equivalent reduced-fiber geometry. Design/interference, rigidity, silent
channels, and fiber dimension are REALIZATIONS of the invariant, not its
definition. Theorem to record: the design/interference dichotomy refines
the face stratification of the GPC polytopes (Theorem A's content).

TRANSFER PILOT EXECUTED -- THREE QUBITS IS A SECOND NEGATIVE CONTROL:
sampled the fiber over the entanglement polytope's top vertex
(1/2,1/2,1/2); every convergent sample has 3-tangle = 1.000000 (GHZ
LU-orbit; consistent with the known classification of locally maximally
mixed 3-qubit pure states). All 3-qubit vertex fibers reduce to points:
the quantum stratification COLLAPSES there -- the smallest entanglement
polytope behaves torically at its vertices despite the nonabelian group.
(Sample size small, backed by the known theorem; noted.)
CONSEQUENCE FOR THE EXPERIMENT: negative controls are now GT (abelian)
and 3 qubits (nonabelian-but-collapsed). The first LIVE test is the
2x2x2x2 (four-qubit) or 2x3x3 entanglement polytope -- big enough for
fiber classes to differ, small enough for the contraction-attack
machinery to census. Rig session: port the local-spectra map to the
multipartite setting (per-party marginal Gram, same Wirtinger gradient),
enumerate the published vertex list, run the fiber-moduli test at every
vertex, and score the Incompleteness Principle on independent ground.

## T7 — THE INVERSE-IMAGE PROGRAM: foundations, the Refinement Theorem
split into provable pieces, and the gamma-2 branch (2026-07)

1. THE EXISTENCE THEOREM IS ALREADY IN THE LITERATURE (key catch): the
   reviewer's aspirational "Refinement Theorem" splits into three parts
   with very different costs, and the first is FREE. Existence: the map
   Phi is semialgebraic (Prop 1), so by HARDT'S TRIVIALIZATION THEOREM
   the Kirwan polytope admits a FINITE semialgebraic partition over
   which Phi is locally trivial -- i.e. the quantum stratification
   EXISTS and is FINITE, as a theorem of real algebraic geometry, not a
   conjecture. Strictness: proved at the vertex level by the 183-simple
   exhibit (Theorem A). Computation: our program -- computing Hardt
   strata explicitly for a nontrivial nonabelian moment map would
   itself be new. Remaining care: "strictly finer than the FACE
   stratification" globally requires fiber-type behavior on face
   interiors (constancy on the polytope interior follows from the
   Generic Fiber Theorem; boundary faces are where the content lives).
2. THE (P, S) INVARIANT -- in-corpus evidence it behaves functorially:
   the pad/lift maps embed P(N,d-1) as a face of P(N,d), and the census
   shows the stratification is CARRIED along these embeddings (classes,
   fibers, walls, and curve isomorphism classes all transport). Hodge
   gives isomorphic (P, S) pairs, as it must (induced by an isomorphism
   of the quantum problems), so it is not the discriminating test; a
   true "same P, different S" pair must come from accidental polytope
   isomorphisms or the entanglement side. Open, correctly framed.
3. THE gamma-2 BRANCH (reviewer's surprise entry) -- calibrated against
   known theory: unlike the 1-RDM case, gamma-2 fibers are GENERICALLY
   POINTS (generic pure states are determined among pure states by
   their 2-RDMs -- the parts-and-whole literature), so the fiber
   program there is the classification of the EXCEPTIONAL LOCUS where
   determination fails. Census connection already in hand: on every
   fiber we probed, gamma-2 is injective (the T-odd coherences separate
   conjugate sheets), consistent with generic determination; our fiber
   states are thus NOT in the exceptional locus, and the machinery
   (contraction attack with Gamma-2 targets) ports directly to hunting
   states that ARE. Physics-facing: the exceptional locus is exactly
   where 2-RDM methods can be blind, the most chemistry-adjacent
   question this program owns.
4. TARGET RANKING (adopted): entanglement polytopes and quantum
   marginal polytopes first (mature images, near-empty fiber theory);
   Horn/LR next; GT = negative control. Home mathematics confirmed:
   stratified symplectic geometry (Sjamaar-Lerman), now with Hardt as
   the semialgebraic backbone.

## T7-PILOT: the gamma-2 fiber program's first data point (2026-07)

Ran the ansatz-free attack with a Gamma-2 MATRIX target at psi_B over
the full 126-dimensional space: 10/10 random starts converge to residual
< 1e-16 with overlap |<psi_B|found>| = 1.00000000 every time. VERDICT:
psi_B has a POINT 2-RDM fiber -- fully determined among pure states by
its two-body marginals, hence NOT in the exceptional locus, despite
sitting on a 1-RDM fiber that is a measured surface. The hierarchy in
one sentence: the 1-RDM forgets a surface; the 2-RDM forgets nothing --
at this state. The exceptional-locus hunt (states the 2-RDM fails to
determine) now has validated machinery and its first negative control;
candidate targets for the hunt: maximal-ODLRO / eta-pairing-type states
and high-symmetry states, where determination is classically suspected
to fail. This branch is the program's chemistry-facing frontier: the
exceptional locus is exactly where 2-RDM-based methods are blind, and
mapping it is a fiber question, not an image question.

HIGH-Tc CALIBRATION (completing the reviewer's truncated paragraph, at
the honesty level this repo requires): no path from fibers to materials.
The defensible indirect chain: fiber diagnostics (uniqueness, flat
directions, hidden degeneracies) can harden v2RDM-type solvers; those
solvers are used on pairing models (Hubbard-class) relevant to
superconductivity THEORY; better instruments for model studies is the
entire claim, several removes from any material. Anything stronger is
marketing.
