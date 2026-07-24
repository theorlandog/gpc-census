# PAPER: "Inverse Geometry of Moment Maps" (retitled per final review
# round) -- framework skeleton (2026-07)
# Architecture: definitions and canonical constructions first; then
# three worked examples as APPLICATIONS -- I. generalized Pauli
# (complete instance), II. entanglement polytopes, III. Horn's problem.
# Later papers are Application papers of this framework.

INVERSION DIRECTIVE (final review round, adopted): the paper must be
readable almost without knowing what a generalized Pauli constraint is.
Structure: "here are four natural objects associated with inverse images
of moment maps" -- reduced fibers, quantum stratification, equivariant
exposedness, symmetry-resolved gap functions -- developed abstractly
with the classical existence machinery (Hardt, Tarski-Seidenberg) doing
the foundation work; THEN "the generalized Pauli setting provides the
motivating example and first complete computation." Mathematics primary,
census illustrative. Sections 1-5 below to be redrafted in that order:
objects and definitions first (currently Sections 2 and 4), setting and
worked example second (currently Sections 1 and 3), interaction agenda
last. Layer IV items are presented as THE RESEARCH AGENDA, never as
expected truths, until proved or demonstrated in an independent family.

OFFICIAL PAPER ROADMAP (adopted from review, in order):
  1. Rank-10 census (finished reference work; in final assembly).
  2. Definitions/framework paper (this skeleton).
  3. Local fiber geometry paper (PROVED results only).
  4. Transfer paper (one genuinely different moment polytope --
     the four-qubit entanglement polytope).
2-RDM exceptional-locus and Hamiltonian-design expansions come only
after these four.

PURPOSE (per external review): define the objects independently of the
GPC census; state precisely what is proved, by what kind of proof, and
what is conjectured. Minimal computation -- the census appears only as
the existence witness and the worked example. Target length: short
(15-20 pp). This is the framework paper; the census (paper 1) and the
differential geometry (paper 2) cite it or absorb it -- decide at
drafting time whether this is standalone paper 3 or paper 2's Part I.

## Section 1 -- Setting
Hamiltonian U(r) action on P(Lambda^N C^r); gamma as moment map; Kirwan
polytope. All standard; cited (Kirwan, Klyachko, Sjamaar-Lerman).
Semialgebraicity of the spectrum map (elementary proposition).

## Section 2 -- The quantum stratification (Definition 1)
Decomposition of the Kirwan polytope into strata of equivalent
reduced-fiber geometry.
- EXISTENCE: Hardt trivialization (cited; classical). [TYPE: structural
  consequence]
- STRICTNESS over the face lattice: Theorem A -- the 183-simple exhibit;
  design/interference as certified vertex invariants. [TYPE: corollary
  through the classification + certified computation; both certificates
  re-runnable]
- The Incompleteness Principle stated as the conjecture for general
  nonabelian moment polytopes; GPC = evidence body; GT and 3-qubit
  entanglement polytope = negative controls (computed); 4-qubit = the
  designated live test. [TYPE: conjecture]

## Section 3 -- Fiber-local structure
The three-object hierarchy (matrix-fiber, spectrum-fiber, reduced
fiber); exact quadraticity of the moment map; the terminating MGS
normal form; obstruction quadric = gamma(delta) + degenerate-PT
feedback.
- Design rigidity: reduction lemma (gauge/kernel argument) [TYPE:
  symbolic proof] + kernel-freeness of all 643 design supports [TYPE:
  computational theorem, exact rational rank, re-runnable].
- The v_B surface, silent-channel rigidity, conic families as worked
  corollaries. [TYPE: certified computation]
- T-MAIN (Local Structure Theorem, algorithmic form): stated; proof
  program with the reviewer's frame. [TYPE: theorem-in-progress]

## Section 4 -- Equivariant exposedness (Definitions 2-3)
H_G = G-symmetric 2-body Hamiltonians; g_G(psi) = sup gap; E_G = {g_G
> 0}. Canonicity discussion (SDP-independence).
- Semialgebraicity of E_G and of g_G: Tarski-Seidenberg. [TYPE:
  structural consequence]
- The spinless T-theorem (real H => nondegenerate eigenstates real).
  [TYPE: classical, cited]
- Corollary: certified-complex interior fiber states are never uniquely
  realizable by real 2-body parents. [TYPE: corollary through the
  classification]
- Measured values (wall 1.147 real; interior 0.582 complex): evidence,
  clearly labeled pilot-grade pending exact certification. [TYPE:
  computational evidence, not theorem]

## Section 5 -- Layer IV: the interaction questions (all open, all
precise)
Which strata admit G-parents; does E_G refine the quantum
stratification; Whitney regularity of E_G; Hardt triviality of g_G;
does symmetry detect hidden fiber structure. The symmetry matrix
(walls x symmetry classes) as the standing experiment with the one
proved column (g_real = 0 off the real locus).

## Section 6 -- Provenance table
Every numbered statement tagged: [classical/cited], [symbolic proof],
[computational theorem -- exact, re-runnable], [computational evidence
-- numerical], [conjecture]. The referee-requested typing, made a
structural feature of the paper rather than a concession.

## Process notes
- Per review: pause new computations for the drafting window, with two
  standing exceptions already in flight on the rig: the cand-44 SOS
  certificate (settles a published question; nearly done) and the
  pre-registered symmetry matrix (feeds Section 5's evidence and is
  harness-protected).
- The historical framing (knot tables -> invariants; enumerations ->
  theories) is the introduction's opening move: the census as
  laboratory, the objects as the contribution.

## Section 7 (ADDED per review) -- Methodological Portability
Explicit table of which constructions are polytope-independent:
  Construction              Depends on GPC?   Demonstrated in
  Projected Jacobian        No                GPC, 4-qubit
  Orbit quotient            No                GPC, 4-qubit
  Stabilizer computation    No                GPC, 4-qubit (exact)
  Fiber sampling            No                GPC, 4-qubit, Horn
  Obstruction analysis      No                GPC (Horn/4q queued)
  Hardt/T-S framework       No                all (cited machinery)
  Pre-registration harness  No                all three families
The claim the table supports: these are canonical constructions for
inverse images of moment maps, with the GPC census as first complete
instance -- not GPC-specific techniques.

## Section 8 (ADDED) -- Limits of applicability: the designed failure
Per the reviewer's falsification challenge, a deliberate failure search.
STRUCTURAL FAILURE EXHIBITED: the cut polytope (CUT(K4) computed as the
minimal exhibit). Its defining map has a discrete source, no group
action, finite constant-size fibers ({S, S^c}, size exactly 2 over
every vertex), and no operator class -- every instrument in the toolkit
is inapplicable or vacuous, for identifiable structural reasons. PRESENTATION CAUTION (final review, adopted): the cut polytope is an
EXAMPLE OF AN OBJECT LACKING THE INGREDIENTS the framework uses -- it
illustrates the intended scope. It is NOT a proof that those
ingredients are necessary, and must not be presented as characterizing
the exact boundary of applicability. Scope statement: the framework is
designed for inverse images of semialgebraic maps with geometric
structure -- the toolkit's natural domain is HAMILTONIAN
settings -- semialgebraic maps from smooth (projective/orbit) domains
carrying a group action and a dual operator class. Outside that class
(cut, matching, and other purely combinatorial polytopes) the framework
claims nothing.
IN-CLASS FAILURE: none found to date; the search is OPEN and logged as
a standing task (candidate stress tests: maps with wild non-proper
fibers; infinite-dimensional limits; non-compact groups). Per the
reviewer: continued survival of serious falsification attempts is the
naturalness evidence that positive examples alone cannot provide.

## Language directive (final round, adopted): "the pilot computations
exhibit ANALOGOUS rigid-versus-deformable behavior in three distinct
inverse-image problems" replaces any "same phenomenon in three
families" phrasing, until a common theorem unifies the observations.
The evolved thesis (recorded as the program's falsifiable statement):
inverse geometry of semialgebraic moment maps admits natural canonical
structures -- stratifications, symmetry-constrained realizability,
computational diagnostics -- that can be studied independently of the
particular application.

## Section ordering (final): per the closing review -- 1. Motivation;
2. Canonical definitions; 3. Classical existence (Hardt, T-S);
4. New constructions (quantum stratification, equivariant exposedness,
symmetry gap functions); 5. Applications AS EVIDENCE (GPC,
entanglement, Horn); 6. Open structural questions. The applications
are evidence for the framework, not the framework itself.

## External assessment of record (closing review round, 2026-07)
Confidence levels, as stated by the external reviewer and adopted as
the program's public posture:
- WELL ESTABLISHED (conditional on exact verification): the census;
  nontrivial GPC fibers; machinery transfer to 4-qubit entanglement
  polytopes; applicability to Horn's problem.
- PLAUSIBLE, NEEDS BROADER VALIDATION: quantum stratification as the
  canonical stratification; equivariant exposedness as a useful
  invariant; structural theorems for gap functions beyond
  semialgebraicity.
- GENUINELY OPEN (the frontier): a common theorem behind the analogous
  rigid/deformable behavior; extension beyond moment-map-type
  problems; independent adoption.
The program's remaining decisive test is INDEPENDENT ENGAGEMENT --
not computable, only maximizable through clean, intrinsic,
application-independent definitions. Identity statement of record:
the census is the experimental platform; the definitions paper is the
mathematical core; the later papers are applications and tests of
that core.
