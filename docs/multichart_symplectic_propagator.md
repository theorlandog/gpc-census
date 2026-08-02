# Multi-chart symplectic 2-RDM propagator

## Result

This implementation turns the observable-manifold geometry theorem into a
global reduced propagator on the 29-determinant invariant transport carrier.
The propagated state is represented only by quantities read from the fixed-frame
2-RDM:

- the 29 determinant weights recovered from a full-rank pair-incidence minor;
- 69 Cartesian collision-free 2-RDM coherences forming an overlapping chart bank.

No wavefunction amplitudes are read by the reduced right-hand side. Amplitudes
appear only when constructing a test initial condition and in the independent
direct-evolution validation.

## Atlas and chart selection

The collision-free coherence graph has 69 usable edges. The implementation
maintains 108 distinct certified spanning trees obtained from different roots,
orders, and BFS/DFS traversals. At each right-hand-side evaluation it scores
candidate trees by the smallest weight used as a reconstruction denominator.
The best regular chart is selected automatically.

A tree chart reconstructs the carrier rank-one matrix from

\[
X_{ab}=\overline c_a c_b,
\qquad
X_{rb}=\frac{X_{ra}X_{ab}}{w_a}.
\]

A chart is rejected before an internal denominator reaches the configured
threshold.

## Boundary chart and zero crossing

Relative phase is undefined when a determinant amplitude vanishes, but the
Cartesian 2-RDM coherences remain regular. If no precomputed tree is regular,
the propagator switches to a boundary phase-synchronization chart using the
full certified edge bank. Vertices with zero weight are omitted from the active
coherence graph, their rows and columns of `X` are set to zero, and the phases
of all remaining active vertices are synchronized from measured Cartesian
coherences.

The test trajectory was designed to pass through exact amplitude-zero events at

\[
t=\pi/4,\qquad t=3\pi/4.
\]

The minimum direct weight was

\[
7.70\times10^{-34},
\]

and the reduced trajectory continued across both events without amplitude
input. Nine chart switches occurred over the interval `0 <= t <= 3`.

The synchronization fallback was in fact never invoked on this trajectory. Its
charts are tagged `-2` in the artifact and `-1` for the dynamic tree, and only
the precomputed tree charts `30, 35, 48, 67, 77` appear. This is the
observable-manifold prediction rather than a shortfall: a vanishing weight
costs the chart nothing when the vanishing determinant is a leaf, because
leaves are never reconstruction denominators, and deleting a leaf leaves the
induced forest `T[A]` connected. At `t=\pi/4` determinant 6 vanishes and chart
35 holds it as a leaf; at `t=3\pi/4` determinant 15 vanishes and chart 77 holds
it as a leaf. The fallback is exercised directly by
`tests/test_multichart_symplectic_propagator.py`, which forces it at
`t=\pi/4`, so the code path is covered even though the selector never needs it
here. A support that vanishes at an articulation point of the coherence graph
would need it.

## Reduced equation

For the carrier matrix convention `X_ab = conjugate(c_a)c_b` and the real
symmetric carrier Hamiltonian,

\[
\dot X=i(HX-XH).
\]

The code reconstructs `X` from the current 2-RDM chart, evaluates this equation,
and returns derivatives only for the chart weights and certified Cartesian
coherences. It never advances a carrier wavefunction.

## Validation

Across 123 sample times, including the two exact zeros:

| quantity | maximum error |
|---|---:|
| carrier matrix `X` | `1.15e-13` |
| `Gamma_2` | `2.17e-13` |
| reconstructed `Gamma_3` | `1.35e-13` |

The initial chart extracted from the 2-RDM agrees with the test state to
`9.71e-17`.

Two symplectic checks are included:

1. exact-unitary propagation preserves the projective Fubini-Study form to
   `7.11e-15`;
2. finite-difference tangent propagation through the reduced integrator
   preserves the same form to `2.56e-10`.

Only the second is a statement about the reduced flow. The first compares
`<u,v>` with `<Uu,Uv>` for a unitary `U`, which is the same inner product on
both sides; what it actually measures is the unitarity of `expm` at this step,
and it is recorded as a consistency floor, not as evidence about the
propagator.

The second number includes ODE and finite-difference error. DOP853 is not itself
a symplectic time-stepper, so the implementation demonstrates preservation to
numerical tolerance rather than exact discrete symplecticity. Replacing the
adaptive integrator with an implicit midpoint or variational reduced integrator
is a clear next numerical improvement.

## Reproducibility

Reruns on one machine reproduce the artifact byte for byte. Across different
BLAS builds the reported maxima are stable but individual sample rows can move
in the last ulp, so the shipped file records one run rather than a
bit-portable computation. The structural fields, the chart sequence, the switch
count and every headline maximum are reproducible.

## Interpretation

The result establishes a functioning global reduced-dynamics method on this
carrier:

- multiple certified observable charts are maintained;
- conditioning is monitored dynamically;
- charts switch before rational denominators become singular;
- exact amplitude-zero events are crossed in Cartesian observable coordinates;
- the reduced trajectory is advanced without wavefunction amplitudes;
- `Gamma_3` is reconstructed exactly from the propagated 2-RDM chart;
- projective symplectic structure is preserved within numerical tolerance.

The representation is deliberately redundant: 167 real chart variables versus
56 physical projective degrees of freedom. The redundancy supplies overlapping
boundary coverage. A production implementation can prune the 69-edge bank to a
smaller fault-tolerant graph or use adaptive edge activation.
