# Four-gate fiber dynamics campaign

## Result summary

This campaign attacks the four remaining gates identified by the exact 2-RDM
observability atlas:

1. a nonseparable carrier-preserving two-body Hamiltonian;
2. an amplitude-free reduced propagator with an explicit conditioning bound;
3. a benchmark beyond practical full configuration interaction;
4. a stronger exact certificate for the 28 supports left unresolved by the
   collision-free-tree test.

The first gate is closed by theorem. The second is closed for diagonal
carrier Hamiltonians, including an entangling cross-block interaction. The
third is demonstrated as a state-dimension and propagation benchmark. The
fourth is implemented as an exact census scan and comes back negative: it
resolves none of the 28, so that gate stays open.

## 1. Entangling carrier-preserving two-body dynamics

Let `S_A` and `S_B` be determinant carriers on disjoint orbital blocks. Choose
one orbital `p` in block A and one orbital `q` in block B and define

\[
H_{AB}=J n_p n_q.
\]

This is a two-body density-density interaction. It is diagonal in the
occupation basis, so it preserves every determinant carrier, in particular
`S_A \boxtimes S_B`.

Write a full-support product state as

\[
|\psi(0)\rangle=\sum_{a,b}\alpha_a\beta_b|D_aE_b\rangle.
\]

After time `t`, its coefficient matrix is

\[
C_{ab}(t)=\alpha_a\beta_b
\exp[-iJt\,u_a v_b],
\]

where `u_a=1[p in D_a]` and `v_b=1[q in E_b]`. If both indicator vectors take
both values zero and one, select rows `a_0,a_1` and columns `b_0,b_1` with
`u_{a_r}=r`, `v_{b_s}=s`. The corresponding 2 by 2 determinant is

\[
\alpha_{a_0}\alpha_{a_1}\beta_{b_0}\beta_{b_1}
(e^{-iJt}-1).
\]

It is nonzero whenever `Jt` is not in `2*pi*Z`. Therefore a generic
full-support product state acquires Schmidt rank at least two. This gives an
explicit entangling, cross-block, carrier-preserving two-body Hamiltonian.
It changes inter-block correlation while retaining the exact carrier 2-RDM
coordinate chart.

This does not transport particles between blocks. It is nevertheless a true
interaction, not a sum of factor Hamiltonians, and the evolved state is not a
product state.

## 2. Amplitude-free exact propagation and stability

A pure state on an `m`-determinant certified carrier can be represented without
amplitudes by

\[
\{w_a=X_{aa}\}_{a=0}^{m-1},\qquad
\{X_{0a}\}_{a=1}^{m-1},
\]

where `X_ab=conjugate(c_a)c_b`. That is `2m-1` stored real coordinates: `m`
weights and `m-1` complex root coherences whose moduli the rank-one
constraints fix. The underlying manifold has `2m-2` real dimensions, since
normalization removes one and the global phase another. Every carrier matrix
entry follows from

\[
X_{ab}=\frac{\overline{X_{0a}}X_{0b}}{w_0}.
\]

For any diagonal carrier Hamiltonian with determinant energies `E_a`,

\[
w_a(t)=w_a(0),\qquad
X_{0a}(t)=e^{-i(E_a-E_0)t}X_{0a}(0).
\]

Thus the cross-block entangling Hamiltonian above has an exact amplitude-free
propagator using `O(m)` state storage. Higher RDMs are generated from the
carrier matrix only when requested.

For a tree-path reconstruction with `L` measured coherence edges,

\[
X_{ab}=\frac{\prod_{j=0}^{L-1}X_{v_jv_{j+1}}}
{\prod_{j=1}^{L-1}w_{v_j}}.
\]

If every input has relative perturbation at most `epsilon`, the first-order
relative condition coefficient is at most

\[
(2L-1)\epsilon+O(\epsilon^2),
\]

provided every denominator weight is nonzero. This exposes the two stability
controls directly: choose short reconstruction paths and avoid charts whose
root or internal weights approach zero. For cyclic factors, shortest-cycle
paths have length at most `floor(d/2)`; product constructions should use
low-diameter spanning trees or switch charts near small weights.

The current implementation is exact for diagonal carrier Hamiltonians. A
general carrier-preserving non-diagonal Hamiltonian requires evaluating the
closed rational vector field and numerically integrating it; that remains a
separate engineering and stability step.

## 3. Beyond-FCI dimensional benchmark

For `q` copies of `S(5,2)`,

\[
m_{carrier}=5^q,
\quad
m_{block}=\binom52^q=10^q,
\quad
m_{FCI}=\binom{5q}{2q}.
\]

The benchmark propagates the `2m-1` reduced coordinates under a cross-block
`n_0 n_0` interaction. At `q=8`:

- carrier dimension: 390,625;
- reduced coordinate count: 781,249;
- separately fixed block-number sector: 100,000,000;
- full fixed-particle FCI sector: 62,852,101,650;
- full-FCI/carrier ratio: about 160,901;
- measured local propagation time: about 1.8 seconds;
- measured peak Python allocation: about 41 MB.

This is a reproducible dimensional and propagation benchmark. It does not yet
compare against an optimized FCI code, and it should not be advertised as a
runtime speedup until the same observables, precision, hardware, and evolution
problem are benchmarked on both sides.

## 4. Determinacy closure for the 28 unresolved carriers

The collision-free-tree certificate discards every channel containing more
than one coherence. The new scan retains the complete signed two-body channel
system.

For each unresolved support it:

1. creates one variable for every unordered carrier coherence;
2. builds the exact integer channel matrix from every oriented 2-RDM entry;
3. computes exact rational row reduction;
4. identifies every coordinate `X_ab` whose unit vector lies in the row space,
   hence is individually determined by the measured channel system;
5. forms the graph of those determined coherences;
6. applies the pure-state path identity when this graph is connected.

This certificate strictly dominates the collision-free graph: every singleton
channel is linearly identifiable, while collided channels can add independent
rows. If it still fails, that is not a non-injectivity result. Nonlinear
rank-one equations may resolve additional variables, and a full algebraic
identifiability calculation would be the next escalation.

Run:

```sh
python scripts/fiber_dynamics_campaign.py closure \
  --out results/data/rdm_determinacy_closure.json
```

### Result: a null

The scan is shipped as `results/data/rdm_determinacy_closure.json`, and its
answer is **0 of 28 newly resolved**. Keeping the collided channels adds real
equations but never enough: 17 of the 28 have no linearly identifiable
coherence at all, and the largest identifiable component reaches only 5. Most
of these carriers are small (support 2 to 4) and simply have too few two-body
channels, which is consistent with the atlas finding that they first connect
at order three, four, or five.

So gate 4 is open, not closed. As stated above this is not a non-injectivity
result; it says the linear-plus-path certificate is silent here, and the next
escalation is the nonlinear rank-one system.

## Physics interpretation

The combined result is an infinite collection of finite-dimensional
interacting fermionic manifolds for which:

- a two-body interaction creates genuine inter-block entanglement;
- the determinant carrier remains invariant;
- the pure state is exactly coordinatized by its 2-RDM;
- the entire RDM hierarchy closes at order two;
- an exact reduced propagator exists without wavefunction amplitudes;
- examples extend far beyond practical storage of the containing FCI vector.

The physically new mechanism is **interaction without loss of reduced
observability**. The 2-RDM is not merely an approximate summary on these
manifolds; it is a complete dynamical state variable. The remaining high-value
step is a non-diagonal interaction that transports occupation or excitation
between factors while preserving a certified glued carrier.
