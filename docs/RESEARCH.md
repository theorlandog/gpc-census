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
- v_A (16,16,16,6,6,6,6,6,6)/21 is DESIGN-INT; v_B (20,14,14,14,14,4,4,4,4)/23
  is INTERFERENCE with cos(gamma) = 3/(4*sqrt(14)).
- Interference onset: absent at rank 8 in the N=4 series, first appears at
  rank 9. Fractions stabilize across ranks 9 to 10 within each series
  (34-37 percent at N=3, 16 at N=4, 14 at N=5).
- Functorial structure: trailing-zero padding and frozen-core lifts
  (lambda -> (1, lambda)) preserve verdicts and generate most interference
  from lower-rank originals; at (4,10) they generate all of it.
- The published (3,9) table is missing one inequality (page-break loss);
  the repaired system (52 inequalities, 58 vertices) is what ships here.
  See results/data/PROVENANCE.md.

## The validation law (non-negotiable)

Every real error in this program was caught by structural invariants, never
by inspection: a sign-parse bug (caught by particle-hole self-duality of
(5,10)), the missing (3,9) facet (caught by embedding coherence and the
face identity P(N,d) cap {lambda_d = 0} = P(N,d-1)), a wrong classify port
(caught by a Farkas feasibility contradiction), and a wrong solver gradient
(caught by finite differences). Consequences:
- No result ships without passing gpc_census.validate and the test suite.
- The preflight gate: any state-solving change must reconstruct v_B end to
  end (scripts/solve_all.py --preflight) before campaigns run.
- Pin solvers. Census verdicts were produced under ortools 9.15.6755
  (the cpsat extra); the CBC fallback labels its backend in every verdict.

## Current campaigns and open items

1. State construction (Tier A) and exactification (Tier B) of all
   interference vertices: scripts/solve_all.py, output
   results/data/states_interference.jsonl. Tier B snaps squared amplitudes
   to k/den, recognizes phases on the p*sqrt(q)/r lattice, and verifies the
   symbolic state by exact characteristic polynomial identity; failures
   fall through labeled TIER-C for hand analysis (likely the interesting
   ones). Recognition layer is unit tested on the v_B cosine.
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
the rule holds in the state's OWN natural orbital basis. Vertex spectra
are maximally degenerate, and within a degenerate lambda-block the
natural orbitals rotate freely, so the canonical-basis admissible set is
too strict. Diagnostic that proved it: the strict filter leaves 9 of 126
determinants for v_B, the resulting skeletons plateau at residual 1e-4,
and none carries the historical weight multiset [7,4,3,2,2,2,2,1].

Fix (implemented in states.admissible_support): block-signature closure.
Group orbitals by degenerate blocks of the spectrum; a determinant is
admissible iff for every active facet the bound is achievable by some
choice of that many coefficients within each block (exact subset-sum per
block). This is a lemma about the method and belongs in the states paper.

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
