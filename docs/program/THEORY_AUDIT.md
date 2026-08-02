# Observable-manifold theory audit

## Purpose

This document separates theorem-level results, computational evidence, conjectures,
and prior-art boundaries for the proposed observable-manifold framework.

## Definitions

Let \(S=\{D_0,\ldots,D_{m-1}\}\) be a fixed determinant carrier and let
\[
|\psi\rangle=\sum_{a=0}^{m-1}c_a|D_a\rangle,\qquad
X_{ab}=\overline{c_a}c_b.
\]

A **certified order-two observable chart** consists of:

1. an exact full-column-rank pair-incidence minor recovering all weights
   \(w_a=X_{aa}\);
2. a connected set of collision-free 2-RDM coherence channels recovering a
   spanning tree of \(X_{ab}\);
3. an explicit rational reconstruction map for all remaining \(X_{ab}\).

The full-support chart domain is \(w_a>0\) for every \(a\).

## Theorem set that is ready to formalize

### T1. Rational pure-state reconstruction on a certified carrier

If pair incidence has full column rank and the collision-free 2-RDM coherence
graph is connected, then every pure full-support state on the fixed carrier is
determined, up to global phase, by the selected 2-RDM entries. The complete
carrier matrix \(X\), and hence every higher RDM, is a rational function of
those entries.

**Status:** theorem; proof is elementary exact linear algebra plus rank-one path
products.

### T2. Interior observable-manifold chart

For a carrier with \(m\) determinants, the pure full-support projective state
space is diffeomorphic to
\[
\Delta_{m-1}^{\circ}\times\mathbb T^{m-1}.
\]
A certified chart realizes this space in selected 2-RDM coordinates.

**Status:** theorem after carefully specifying the chart map and its inverse.
The topology itself is standard; novelty can only lie in the explicit
low-order-observable realization and certificate.

### T3. Pullback symplectic form

With
\[
c_0=\sqrt{w_0},\qquad c_a=\sqrt{w_a}e^{i\phi_a},
\]
the Fubini--Study symplectic form restricts, up to the chosen normalization, to
\[
\omega=\sum_{a=1}^{m-1}dw_a\wedge d\phi_a.
\]

**Status:** standard geometric quantum mechanics specialized to the certified
observable chart. Do not claim the Fubini--Study/action-angle fact as new.
Claim the exact representation in certified 2-RDM coordinates.

### T4. Boundary graph theorem

For an active support \(A=\{a:w_a>0\}\), let \(G[A]\) be the surviving
observable-coherence graph. If it has \(c(A)\) connected components, the
selected chart leaves a residual relative-phase ambiguity
\[
\mathbb T^{c(A)-1}.
\]

**Status:** theorem for the selected channel family. A proof must distinguish
the rank defect of one chart from non-identifiability by the complete 2-RDM.

### T5. Exact reduced Hamiltonian dynamics

For a known carrier-preserving Hamiltonian and a certified chart, the
Schrödinger flow descends to an exact Hamiltonian flow in observable
coordinates, and the BBGKY hierarchy closes through the exact rational map
\(\Gamma^{(3)}=F_{S,3}(\Gamma^{(2)})\).

**Status:** theorem on the pure full-support locus. Extension across boundary
strata requires the multi-chart Cartesian construction.

### T6. Multi-chart continuation

A bank of certified coherence graphs, supplemented by Cartesian coherence
coordinates, permits continuation through amplitude-zero events whenever the
active measured graph remains phase-synchronizing on each side and the
Cartesian derivatives determine the outgoing coherence direction.

**Status:** algorithmic theorem needs a precise regularity statement. Current
numerics are evidence, not a general proof.

### T7. Cyclic-window family

For cyclic \(N\)-window supports on \(d\) orbitals, singleton incidence has
rank
\[
d-\gcd(d,N)+1.
\]
Under stated coprimality and labeled-channel hypotheses, pair incidence is
full rank and adjacent-window channels yield a connected observable graph.

**Status:** theorem candidate. Channel collision-freedom and labels must be
proved with exact fermionic conventions, not inferred from examples.

### T8. Gluing theorem

For disjoint-orbital products, the cross-pair incidence block factors as a
Kronecker product of singleton-incidence matrices. The glued observable graph
contains a spanning categorical-product graph. Connectivity follows when the
factor graphs are connected and at least one is non-bipartite.

**Status:** theorem candidate. Use “contains a spanning subgraph,” not “is
exactly,” unless all extra channels are classified.

## Results that remain evidence or special constructions

- The 771/799 census count is an exact computational theorem only to the extent
  that the independent verifier and input provenance are correct.
- The 29-state invariant transport carrier is an exact finite construction.
- The multi-chart zero-crossing experiment is a numerical validation of the
  algorithm, not yet a global existence theorem. It is weaker evidence for T6
  than it looks, and in one direction stronger than it looks. The recorded
  trajectory in `results/data/multichart_symplectic_propagator.json` never
  invokes the synchronization fallback: at both amplitude zeros the vanishing
  determinant is a leaf of the selected tree, leaves are never reconstruction
  denominators, and deleting a leaf leaves `T[A]` connected. Those crossings
  are therefore a corollary of T4 with `c(A)=1`, not evidence for T6. The
  fallback path is exercised only by a forced unit test. The open case, and
  the one T6 actually needs, is a zero at an articulation point of the
  coherence graph.
- Exponential benchmark separation for products of two-level cells is valid
  only for the specified invariant product manifold and algorithm classes.
  Sharper: the Hamiltonian in that family has no inter-cell coupling and the
  benchmark state is a product state, so the Schmidt rank stays 1 at every
  tested size. The separation is against a dense `2^q` amplitude vector that
  nobody would choose for the system. Against the like-for-like factorized
  amplitude update, which is also `O(q)`, the margin is a constant near 1.3.
- “Transport rigidity” is proved for each exact finite leakage-rank audit, but
  the broad cyclic-family rigidity statement remains a conjecture.
- Determinacy closure on the unresolved 28 supports has been generated:
  `results/data/rdm_determinacy_closure.json` records `newly_resolved: 0` and
  `remaining_unresolved: 28`. The count did not move, and the artifact's own
  nonclaim field says so correctly: failure of the certificate is not
  non-injectivity, since nonlinear rank-one equations may still resolve
  further coherences. Cite the negative result, not a pending task.

## Prior-art boundaries

`docs/prior_art_roadmap.md` already maps this territory outside GPC with
primary citations: phase retrieval, quantum marginals, matrix completion and
the rest, ordered by expected payoff. Use it as the source for the paper's
closest-result table rather than starting a second list. The bare topics below
are an index into it, not a replacement.

The following ideas are established and must be cited rather than presented as
new:

- projective Hilbert space as a Kähler/symplectic manifold;
- Hamiltonian formulation of Schrödinger dynamics;
- pure states determined by collections of reduced density matrices;
- phase recovery from connected measurement graphs;
- BBGKY and approximate TD-2RDM closure;
- complete-active-space and product-state factorizations;
- matrix completion of RDM entries under additional hypotheses.

The defensible contribution is the combination of:

1. support-combinatorial exact certificates;
2. an explicit rational inverse from selected fixed-frame 2-RDM entries;
3. a census-scale exact atlas;
4. constructive infinite certified carrier families;
5. graph-controlled boundary strata and multi-chart continuation;
6. exact carrier-preserving reduced Hamiltonian dynamics;
7. independently verifiable software artifacts and scaling benchmarks.

## Claims that must not appear

- “The 2-RDM generically determines a fermionic pure state.”
- “The framework solves the \(N\)-representability problem.”
- “All unresolved census carriers are non-injective.”
- “Observable coordinates universally beat wavefunction methods.”
- “The symplectic geometry of quantum mechanics is new.”
- “The complete 2-RDM is always smaller than a carrier wavefunction.”
- “An exponential speedup holds outside the stated invariant family.”
- “Solver infeasibility evidence is an exact no-design proof.”

Three more, each of which a benchmark in this repository invited before being
corrected. They are recorded because the mistake is easy to repeat:

- Reporting a speedup over a baseline that calls the method being measured.
  The dense 2-RDM propagator extracts the same certified observables and calls
  the same `reconstruct_X`, so it is the sparse right-hand side plus a
  wrapper and cannot be faster by construction. Compare variable counts there,
  not clocks.
- Reporting a speedup of an analytic solve over numerical integration of the
  same trajectory as if it were a per-evaluation advantage. The exact
  observable propagator is one call; the ODE baselines take 314 right-hand-side
  evaluations across 26 adaptive steps.
- Quoting a benchmark factor without the host. Across two machines the
  compiled-kernel margin moved from 1.44 to 1.17 and the exact propagator's
  margin over carrier-wavefunction propagation from 2.29 to 1.76. Method
  ordering was stable; the factors were not. Quote the ordering, the variable
  counts, and the accuracy, which do reproduce.

## Adversarial audit outcome

`docs/adversarial_observable_audit.md` runs items 2 and 3 of the list below
and reports verdicts on T1 to T8. Two of its findings change what may be
written about T1, and both are recorded in the ledger.

Full pair-incidence rank does not imply order-two reconstruction. The carrier
`{{0,1,2},{3,4,5}}` in `wedge^3 C^6` is pair-rigid and its complete 2-RDM is
independent of the relative phase, because a two-body operator moves at most
two orbitals and these determinants differ in three. Weight recovery is not
half of the theorem; the coherence graph is doing separate work.

A connected collision-free graph does not imply linear identifiability from
all 2-RDM channels. At `d=6`, `N=3`, `m=3` all 1140 pair-rigid carriers have a
connected graph but only 960 have a full-rank real-linear channel matrix. The
tree certificate, linear identifiability, and nonlinear pure-state injectivity
are three separate statements and none implies the next.

The fermionic sign convention is now checked by three independent
implementations across 8,186,968 cases with zero mismatches. That is software
corroboration, not the written proof item 2 asks for.

## Required scrutiny before submission

1. Independent proof review of T1--T8.
2. Symbolic sign audit using a second implementation.
3. Exhaustive small-system counterexample search.
4. Full literature comparison with explicit closest-result table.
5. Reproducible artifact bundle with hashes and no hidden generator imports.
6. Referee-style attempt to falsify every headline claim.
7. Separate theorem, numerical-result, and conjecture labels throughout.
