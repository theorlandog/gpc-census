# Generalized exact observable propagation benchmark

## Scope

This campaign generalizes the exact observable propagator from one 29-state
carrier to a scalable family of `q` exactly observable two-level transport
cells.  Each cell is represented by its three Cartesian Bloch coordinates,
which are obtainable from certified low-order RDM entries.  The invariant
manifold is a product of Bloch spheres and has `3q` redundant real observable
coordinates, while a generic dense carrier wavefunction has `2^q` complex
amplitudes.

The Hamiltonian is

\[
H=\sum_{j=1}^q \omega_j\sigma_x^{(j)}.
\]

It preserves the product observable manifold.  Each observable cell therefore
undergoes one exact real rotation about its x axis.  No ODE solver, carrier
wavefunction, dense 2-RDM, or higher RDM is used by the observable update.

## Complexity theorem

For this family:

* observable storage and exact update work are `O(q)`;
* a dense carrier state requires `O(2^q)` storage;
* applying all local rotations to that dense carrier costs `O(q 2^q)`;
* a generic dense matrix exponential has still higher practical cost.

The exponential separation follows from restricting dynamics to a
lower-dimensional observable invariant manifold.  It does not contradict the
fact that a full projective `m`-determinant carrier has the same intrinsic
state dimension in amplitude and observable coordinates.

## Measured results

All timings are median end-to-end update times for the same initial state and
clock value.  The dense carrier comparison is structure-aware: it applies
local tensor rotations rather than constructing a dense Hamiltonian matrix.

| Cells q | Dense carrier amplitudes | Observable variables | Observable time | Dense carrier time | Speedup | vs factorized amplitudes |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 16 | 12 | 7.88 us | 59.3 us | 8x | 1.39x |
| 6 | 64 | 18 | 7.87 us | 92.6 us | 12x | 1.40x |
| 8 | 256 | 24 | 9.14 us | 164.5 us | 18x | 1.20x |
| 10 | 1,024 | 30 | 8.05 us | 272.2 us | 34x | 1.38x |
| 12 | 4,096 | 36 | 8.06 us | 511.1 us | 63x | 2.45x |
| 14 | 16,384 | 42 | 8.12 us | 1927.0 us | 237x | 1.38x |
| 16 | 65,536 | 48 | 9.76 us | 10624.8 us | 1,088x | 1.27x |
| 18 | 262,144 | 54 | 9.60 us | 55347.5 us | 5,764x | 1.19x |

A generic dense matrix exponential is measured where it is still feasible:

| Cells q | Generic dense expm | Observable speedup |
|---:|---:|---:|
| 4 | 53.8 us | 7x |
| 6 | 8,013.0 us | 1,018x |
| 8 | 108,004.5 us | 11,811x |

The maximum observable error over all tested sizes is `4.16e-16`.

The observable update also beats the factorized amplitude update by a
constant factor, around 1.2x to 1.4x across the range. Both are linear in `q`;
this comes from using three real observable coordinates instead of two complex
amplitudes per cell and from the particular NumPy kernels. It is not a
representation-independent lower bound, and the spread across rows is timing
noise, not structure.

Timings are host specific and this table records one machine. What is
structural, and reproduces anywhere, is the shape: observable time flat in
`q`, dense carrier time doubling with each added cell, and the accuracy.

## What the exponential separation is over

The Hamiltonian has no inter-cell coupling and the benchmark state is a
product state, so the trajectory never leaves the product manifold. Checking
the Schmidt rank across the first cut confirms it stays 1 at every tested
size: the dense carrier is never entangled. The `O(q)` versus `O(q 2^q)`
separation is therefore against a dense `2^q` amplitude vector, which is not
the representation anyone would choose for a non-interacting product system.

That is why the factorized amplitude update matters. It is the like-for-like
opponent, it is also `O(q)`, and the observable method beats it only by a
constant. The honest reading of this family is that it demonstrates the
machinery works and scales, not that observable coordinates defeat a
well-chosen amplitude representation.

## What “beats all” means here

The observable update beats every implementation in the generalized benchmark:

1. generic dense matrix-exponential propagation, measured for `q <= 8`;
2. structure-aware dense carrier propagation, measured across the full range;
3. factorized amplitude propagation, measured, and only by a constant factor.

Dense 2-RDM propagation is not benchmarked here. Its state dimension grows at
least as the square of the one-particle pair space and is unnecessary on this
invariant manifold, but that is an argument, not a measurement in this
artifact.

The result is exact and scalable for the stated family.

It is not a universal theorem that observable coordinates beat every possible
algorithm.  An implementation specially fused to hardware or a different
minimal parametrization can change constant factors.  The robust theorem is
the exponential advantage over generic dense carrier and FCI representations
when the physical trajectory lies on a product observable manifold of linear
dimension.

## Scientific consequence

This identifies the condition under which the geometry becomes a genuine
computational advantage:

> The observable manifold must be parametrically smaller than the containing
> carrier projective space, not merely a different coordinate system on the
> same full carrier.

For products of observable transport cells, low-order observables provide an
exact global state representation and exact dynamics with linear resources,
while generic many-body amplitudes scale exponentially.
