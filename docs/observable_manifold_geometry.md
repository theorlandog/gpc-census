# Universal geometry of certified observable carriers

## Result

Let a determinant carrier contain `m` determinants. Suppose its pair-incidence
matrix has full column rank and choose a collision-free 2-RDM spanning tree
`T`. On the pure full-support locus, the selected 2-RDM data are not merely an
injective list of measurements. They form a global action-angle chart on the
projective carrier state space.

Write

\[
|\psi\rangle=\sum_{a=0}^{m-1}c_a|D_a\rangle,
\qquad
w_a=|c_a|^2,
\qquad
X_{ab}=\overline c_a c_b.
\]

An invertible pair-incidence minor recovers every `w_a` from diagonal 2-RDM
entries. For each oriented tree edge `e=(a,b)`, its collision-free 2-RDM
channel recovers `z_e=X_ab`, including the fermionic sign.

## Observable-chart theorem

On the normalized full-support pure-state stratum, the map

\[
[\psi]\longmapsto (w_0,\ldots,w_{m-1};z_e:e\in T)
\]

is a real-analytic diffeomorphism onto the semialgebraic set

\[
\sum_a w_a=1,
\qquad w_a>0,
\qquad |z_{ab}|^2=w_aw_b\quad ((a,b)\in T).
\]

Equivalently, the observable manifold is

\[
\Delta_{m-1}^{\circ}\times\mathbb T^{m-1}.
\]

Its real dimension is `2m-2`, exactly the dimension of
`CP^(m-1)`. The inverse is the path-product reconstruction already certified
in the RDM atlas.

### Proof

Choose a root and fix its phase. The positive weights determine the coefficient
moduli. The phase of each nonzero tree coherence gives the relative phase
across that edge. Since a tree has a unique root-to-vertex path, these edge
phases recover every coefficient phase uniquely. Conversely, coefficients
produce the displayed relations. Both directions are analytic while all
weights are positive. The incidence minor only applies an invertible linear
change to the weight coordinates. QED.

## Symplectic consequence

Writing

\[
c_a=\sqrt{w_a}\,e^{i\theta_a}
\]

and removing one global phase, the pullback of the Fubini-Study symplectic
form is, up to the conventional overall normalization,

\[
\omega=\sum_{a=1}^{m-1}dw_a\wedge d(\theta_a-\theta_0).
\]

Thus the reconstructed weights and relative phases are canonical action-angle
coordinates. For every carrier-preserving Hamiltonian, its expectation value

\[
h(w,\theta)=\langle\psi(w,\theta)|H|\psi(w,\theta)\rangle
\]

generates the exact Schrödinger flow by Hamilton's equations directly on the
observable manifold. Exact 2-RDM closure is therefore not only an algebraic
substitution for the 3-RDM: it is a finite-dimensional Hamiltonian system in
low-order observable coordinates.

This is the conceptual upgrade. The carrier wavefunction is a lift of a
symplectic observable trajectory, not an additional dynamical object.

## Boundary stratification and singular locus

Let

\[
A=\{a:w_a>0\}
\]

be the active determinant support and let `T[A]` be the induced forest. If it
has `c(A)` connected components, the selected tree channels determine phases
inside each component but do not determine the relative phases between
components. The certificate chart therefore has a residual fiber

\[
\mathbb T^{c(A)-1}
\]

and loses `c(A)-1` relative-phase directions.

Consequently:

* there is no branch locus or monodromy on the full-support interior;
* the apparent poles in path-product formulas are coordinate-boundary poles;
* the selected chart is regular on a boundary stratum exactly when `T[A]` is
  connected;
* chart singularities are governed by graph cuts, not mysterious analytic
  branching.

This statement concerns the selected certificate chart. Additional collided
2-RDM channels can remove some boundary ambiguities, so a disconnected
`T[A]` is not by itself a non-injectivity theorem for the complete 2-RDM.

## Atlas theorem

A collection of collision-free coherence trees gives an observable atlas.
For an active support `A`, any chart whose induced tree remains connected is
regular. This turns numerical chart design into a combinatorial covering
problem:

> choose a small set of certified spanning trees so that physically relevant
> active supports remain connected in at least one tree.

The exact artifact `results/data/observable_manifold_geometry.json` audits
path, star, and balanced trees through ten vertices. It demonstrates why one
fixed tree is unsuitable globally. For example, deleting the center of an
`m`-vertex star leaves `m-1` components and an `(m-2)`-torus ambiguity, whereas
other trees fail on different support boundaries.

## What is universal and what is carrier specific

Universal:

* the full-support topology `simplex interior x torus`;
* the `2m-2` dimension;
* action-angle symplectic coordinates;
* graph-controlled boundary rank loss;
* rational reconstruction of all higher RDMs.

Carrier specific:

* which 2-RDM entries realize the weight minor and tree edges;
* fermionic channel signs;
* conditioning of the incidence inverse and path products;
* which Hamiltonians preserve the carrier;
* whether extra collided channels repair boundary disconnections.

## Significance

The result reframes the census as a library of explicit symplectic reductions.
Each certified support supplies a many-fermion manifold on which a finite set
of 2-RDM observables is a complete canonical state description. The infinite
families and the 29-state transport model then become constructive examples of
Hamiltonian dynamics on these observable manifolds.

A strong Paper 2 can now be organized around three theorem layers:

1. **Observable-chart theorem:** exact 2-RDM coordinates and universal
   action-angle geometry.
2. **Construction theorem:** infinite families and gluing operations that
   preserve observability.
3. **Dynamics theorem:** entangling and hopping Hamiltonians whose exact flow
   remains on an observable manifold.

The remaining high-impact calculation is no longer another support example.
It is a multi-chart propagator that switches trees before a boundary pole and
proves stable global evolution across amplitude zeros.
