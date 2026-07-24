# PAPER: "Quantum Stratifications and Equivariant Exposedness of
# Moment-Map Fibers" -- axiomatization skeleton (2026-07)

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
