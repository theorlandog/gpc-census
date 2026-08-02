# Exact RDM observability atlas

## Result

The census now certifies when a 2-RDM is an exact coordinate system for a
pure many-body state on its shipped determinant carrier.  The certificate
holds for 771 of 799 supports, including every one of the 156 INTERFERENCE
supports.

This is not a generic solution of the 2-RDM problem.  It is a sparse,
fixed-frame, pure-state reconstruction theorem with an exact certificate for
every census record.  Its value is that the inverse is rational, constructive,
and immediately supplies exact higher-RDM closure on every certified carrier.

The generated artifact is
`results/data/rdm_observability_census.json`.  It is reproduced by
`scripts/rdm_observability_census.py` and independently checked by
`scripts/verify_rdm_observability_standalone.py`.

## The support theorem

Fix a determinant carrier

\[
S=\{D_0,\ldots,D_{m-1}\},\qquad
|\psi\rangle=\sum_{a=0}^{m-1}c_a|D_a\rangle,
\]

and write

\[
X_{ab}=\overline{c_a}c_b,\qquad w_a=X_{aa}.
\]

Let \(B_S\) be the pair-incidence matrix,

\[
(B_S)_{P,a}=\mathbf 1[P\subset D_a],
\]

where \(P\) runs over orbital pairs.  The diagonal pair occupations obey

\[
\Gamma^{(2)}_{P,P}=\sum_a(B_S)_{P,a}w_a.
\]

The convention is \(P=(p_1<p_2)\),
\(A_P=a_{p_2}a_{p_1}\), and
\(\Gamma^{(2)}_{P,Q}=\langle A_P^\dagger A_Q\rangle\). Orbital labels in
the artifact are zero-based.

For sorted orbital pairs \(P,Q\), every off-diagonal entry has the exact form

\[
\Gamma^{(2)}_{P,Q}
=\sum_R \epsilon(P,R)\epsilon(Q,R)
X_{\operatorname{idx}(P\cup R),\operatorname{idx}(Q\cup R)},
\]

where the sum is over sorted \(R\) of size \(N-2\) such that
\(R\cap(P\cup Q)=\varnothing\), \(P\cup R,Q\cup R\in S\), and
\(\operatorname{idx}\) is the carrier index. Also,
\(\epsilon(K,R)=(-1)^{\operatorname{inv}(K\mathbin\Vert R)}\).

Call a channel collision-free when this sum contains exactly one nonzero
coefficient coherence.  Make a graph whose vertices are the determinants in
\(S\) and whose edges are the collision-free channels.

**Theorem.** Suppose:

1. \(B_S\) has full column rank.
2. The collision-free 2-RDM graph is connected.
3. Every \(c_a\) is nonzero.
4. The state is pure and the orbital frame and carrier are fixed.

Then \(\Gamma^{(2)}\) rationally determines \(X\), hence the pure state up to
global phase and every higher RDM.

**Proof.** An invertible \(m\) by \(m\) minor of \(B_S\) recovers all weights
from diagonal entries of \(\Gamma^{(2)}\).  Each collision-free tree edge
recovers one \(X_{ab}\), including its fermionic sign.  Along a tree path
\(a=v_0,\ldots,v_r=b\),

\[
X_{ab}=
\frac{\prod_{j=0}^{r-1}X_{v_jv_{j+1}}}
{\prod_{j=1}^{r-1}w_{v_j}}.
\]

This reconstructs the complete rank-one carrier matrix.  Every \(k\)-RDM is
a signed linear contraction of that matrix.  All operations are rational in
the selected 2-RDM entries.  \(\square\)

The full-support hypothesis is the domain of this particular rational chart.
A different tree can often remove a coordinate singularity at a boundary.

## Exact census result

| Classification | Records | Reconstructed from 2-RDM | Unresolved by this certificate |
|---|---:|---:|---:|
| DESIGN-INT | 631 | 603 | 28 |
| DESIGN-REAL | 12 | 12 | 0 |
| INTERFERENCE | 156 | 156 | 0 |
| **Total** | **799** | **771** | **28** |

Additional exact facts:

- Pair incidence has full column rank for all 799 supports.
- The 799 record supports contain 305 unordered one-hop determinant pairs in
  total, and every corresponding 2-RDM channel is included in the scan.
- The 63 positive singleton-incidence kernels are all INTERFERENCE records,
  with dimensions \(1^{\times47},2^{\times15},5^{\times1}\).
- For the 28 supports not certified at order two, the collision-free graph
  first becomes connected at order three for 22, order four for 5, and order
  five for 1. Combined with the pair-diagonal weight certificate, these are
  sufficient higher-order reconstruction trees, not a proof that lower-order
  reconstruction is impossible.
- The \(v_B\) carrier has eight vertices and a seven-edge 2-RDM reconstruction
  tree.

The last caveat is essential.  Collided channels can contain useful joint
information.  A disconnected collision-free graph means only that this simple
certificate is silent.

## What this unlocks

### Exact carrier-conditional reconstruction

On every certified carrier, researchers can replace the wavefunction by a
small set of diagonal pair occupations and \(m-1\) complex tree coherences.
The artifact records the exact incidence minor, fermionic signs, and tree for
each state.  No nonlinear fit or phase search is required.

### Exact higher-RDM closure

For every certified support,

\[
\Gamma^{(k)}=F_{S,k}(\Gamma^{(2)}),\qquad 2\le k\le N,
\]

on the pure full-support stratum.  The map is rational and generated directly
from the certificate.  In particular, \(\Gamma^{(3)}\), the missing term in
the usual time-dependent 2-RDM hierarchy, is exactly reconstructible there.

### Closed reduced dynamics when a Hamiltonian preserves the carrier

Let a known two-body Hamiltonian preserve a certified carrier. Along a pure
full-support trajectory in that carrier, BBGKY becomes

\[
\dot\Gamma^{(2)}
=\mathcal L_H\left(\Gamma^{(2)},F_{S,3}(\Gamma^{(2)})\right).
\]

For a time-independent \(H\), this is an exact autonomous reduced equation on
the image of the carrier. For a scheduled \(H(t)\), it is exact and closed but
nonautonomous. Diagonal density-density Hamiltonians preserve every determinant
carrier and therefore provide immediate exact test cases. They do not generally preserve
the complete 1-RDM, and on incidence-rigid DESIGN supports their phase motion
can be only orbital gauge.  A scientifically stronger driver must be checked
against both carrier preservation and the intended fixed-1RDM condition.

### Ground truth for algorithms

The atlas supplies hundreds of exact, sparse correlator-reconstruction
instances with
known answers and varied channel collisions.  It can test:

- pure-state 2-RDM tomography and noise stability;
- \(\Gamma^{(3)}\) reconstruction methods;
- time-dependent 2-RDM closures;
- Hamiltonian learning on known carriers;
- variational 2-RDM relaxations when combined with exact parent Hamiltonians.

## What is not unlocked

- Mixed states on a whole carrier are not reconstructed.  Purity is used in
  the path-product identity.
- The theorem does not solve generic \(N\)-representability.
- It does not infer an unknown carrier or orbital frame.
- It is not invariant under an unresolved degeneracy-stabilizer quotient.
- It does not prove that the 28 unresolved records are non-injective.
- It does not by itself produce fixed-1RDM dynamics.

These boundaries separate an exact reusable result from the next research
gate.

## Relation to prior work

The general ideas that lower RDMs can determine pure states, that unique
two-body ground states are determined by their 2-RDMs, and that phases can be
propagated along a connected measurement graph are established. Close
reference points include [Mazziotti's complete RDM
reconstruction](https://doi.org/10.1016/S0009-2614(00)00773-9),
[Chen et al. on unique preimages of extreme
RDMs](https://arxiv.org/abs/1205.3682), and [Alexeev, Bandeira, Fickus, and
Mixon on propagating relative phases along a connected polarization
graph](https://arxiv.org/abs/1210.7752). Explicit
counterexamples to global uniqueness also exist in [two-body-blind
Hamiltonians and quantum codes](https://arxiv.org/abs/1010.2717).

The closest current completion result is [Massaccesi et al. on unique matrix
completion of RDMs](https://arxiv.org/abs/2603.13087), which is based on a
nondegenerate two-body ground state and a Hamiltonian-dependent set of known
entries. The present result is different in scope: it is Hamiltonian-free and
applies to every pure full-support state on each certified carrier, but it is
conditional on that known carrier and fixed frame. Its contribution is the
constructive support certificate, explicit rational inverse, and exact
census-scale atlas. The spanning-tree identity or the abstract existence of a
higher-RDM functional should not be claimed as new by itself.

For time-dependent 2-RDM work, current methods reconstruct the 3-RDM
approximately, as in [Lackner et al.](https://arxiv.org/abs/1411.0495), and
recent tests show that time-local closure remains regime-dependent, as in
[Egenlauf et al.](https://arxiv.org/abs/2512.13913). The exact scheduled
\(v_B\) driver in `docs/vb_carrier_hamiltonian.md` supplies a
carrier-preserving, fixed-complete-1RDM trajectory in one finite model. The
remaining physics gate is to establish stability and scalable breadth.

## Research program

The census should be used as a theorem generator, benchmark generator, and
counterexample generator in parallel.

### 1. Exact nontrivial closed 2-RDM trajectory: completed for v_B

The sparse scheduled Hamiltonian in `docs/vb_carrier_hamiltonian.md` satisfies

\[
(1-P_S)H(t)P_S=0
\]

and, with no projective remainder,

\[
H(\phi)|\psi(\phi)\rangle
=i|\dot\psi(\phi)\rangle.
\]

It uses a wall-regular conic clock, has a seven-entry symbolic pair-kernel
support, cancels the sole carrier-leakage channel, and has an aligned lift with
one constant complete 1-RDM. The carrier's unimodular pair-incidence minor and
seven-edge coherence tree then give exact nonautonomous
\(\Gamma^{(3)}(\Gamma^{(2)})\) closure for the driven loop. The certificate is
in `results/data/vb_carrier_hamiltonian.json`.

### 2. Turn the inverse into a stable algorithm

Derive conditioning bounds for the selected incidence minors and tree path
products.  Optimize the tree for noise amplification, emit alternative charts
near zero amplitudes, and compare reconstruction error with generic 3-RDM
closures.  Exact invertibility becomes computationally important only when it
is numerically stable.

First emit an executable reconstruction DAG for \(F_{S,3}\), verify the
reconstructed 3-RDM against direct contraction around the \(v_B\) loop, and
integrate the closed 2-RDM equation without reading wavefunction amplitudes.

### 3. Prove an infinite-family theorem

Identify operations on supports that preserve pair-incidence rank and
collision-free connectivity, such as controlled padding, frozen cores, or
combinatorial gluing.  A scalable family changes the result from a finite
atlas into a reusable design principle.

### 4. Build a parent and fidelity atlas

Screen the certified supports for strict two-body carrier exposers and
pair-kernel parents that isolate the target ray. Every candidate must then be
upgraded to integer minors, rational dual certificates, certified many-body
gaps, and explicit 2-RDM fidelity witnesses in a stated particle-number
sector. This is a proposed campaign, not a current census result.

### 5. Benchmark beyond full configuration interaction

Generate larger certified carriers from the infinite-family theorem, evolve
them with carrier-preserving two-body drivers, and publish full
\(\Gamma^{(2)}(t)\), reconstructed \(\Gamma^{(3)}(t)\), exact errors, and
runtimes.  The decisive computational result is an exact reduced simulation at
a size where the full wavefunction is no longer practical.

## Significance threshold

The atlas and the exact \(v_B\) driver support a focused mathematical-physics
or quantum-information paper now. The first two items below are achieved in
one finite model. A broader result should combine all five at scale:

1. an exact correlation-changing, carrier-preserving two-body trajectory that
   is not a one-body orbit, achieved for \(v_B\);
2. an exact closed 2-RDM theorem, achieved nonautonomously for that trajectory,
   followed by an executable amplitude-free propagator;
3. three nontransport primitive classes, including a higher-dimensional
   singleton kernel;
4. a stability theorem or a measured numerical advantage;
5. either an infinite support family or a benchmark beyond full wavefunction
   propagation.

That package would turn the census from evidence about exceptional fermionic
states into a constructive theory of reduced observability, control, and
certification.
