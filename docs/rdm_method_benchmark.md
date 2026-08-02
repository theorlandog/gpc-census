# Four-way propagation benchmark

## Scope

This benchmark compares four representations on the same exact transport
problem: the 29-determinant invariant carrier in four fermions and ten
orbitals, driven by the Hermitian cross-block hopping Hamiltonian used by the
multi-chart propagator.

Every method uses the same initial state, interval `0 <= t <= 1`, DOP853
solver, tolerances, and maximum step. Each timed run occurs in a fresh Python
process and is repeated three times.

The methods are:

1. **Full FCI:** 210 complex amplitudes in `wedge^4 C^10`.
2. **Carrier wavefunction:** 29 complex amplitudes in the known invariant carrier.
3. **Dense 2-RDM:** the complete `45 x 45` complex matrix, closed exactly by
   reconstructing the carrier matrix from the 2-RDM at each RHS evaluation.
4. **Sparse observable manifold:** 29 weights and 69 complex certified
   collision-free coherences, with multi-chart reconstruction.

The dense and sparse reduced methods do not evolve wavefunction amplitudes.

## Results

| Method | Real ODE variables | Median runtime (s) | Median peak RSS (MiB) | Maximum reference error |
|---|---:|---:|---:|---:|
| Full FCI | 420 | 0.00697 | 178.15 | 5.55e-17 |
| Carrier wavefunction | 58 | 0.00439 | 177.87 | 5.55e-17 |
| Dense 2-RDM | 4050 | 0.68326 | 180.82 | 2.14e-12 |
| Sparse observable manifold | 167 | 0.62023 | 178.13 | 1.40e-12 |

All methods took 26 accepted steps and 314 RHS evaluations, so the runtime
comparison primarily measures RHS cost rather than different adaptive-step
behavior.

The runtimes are machine specific and the table records one host. What
reproduces exactly elsewhere is the structure: the same 420, 58, 4050 and 167
variable counts, the same 314 evaluations and 26 accepted steps for all four
methods, and the reference errors to seven significant digits. On a slower
host every method scales up together and the dense-to-sparse runtime ratio
stayed above one (1.07 rather than 1.10), while the variable-count ratios are
exact integers and do not move.

## What became easier

The sparse observable method uses 24.25 times fewer propagated real variables
than the dense 2-RDM method:

```
4050 / 167 = 24.25.
```

It is also 1.10 times faster in this unoptimized implementation, but that
number needs a caveat and should not be read as a comparison of two
independently designed integrators. The dense method closes the 2-RDM by
extracting the same certified observables from `Gamma_2` and calling the same
`reconstruct_X`, then applying the sparse `vec(X) -> vec(Gamma_2)` map and
carrying 4050 variables through the stepper. It is the sparse right-hand side
plus a wrapper, so it cannot be faster by construction, and the measured
per-evaluation ratio of 1.12 is that wrapper. The defensible reading is the
variable count, not the clock.

The compression claim stands on its own: the observable atlas turns the dense
reduced problem into a much smaller ODE system without sacrificing exact
closure.

The sparse method also stores fewer propagated variables than full FCI:

```
420 / 167 = 2.51.
```

That variable-count advantage is modest in this small model but can grow in
families where the ambient FCI sector grows combinatorially.

## What did not become easier yet

The current sparse implementation is not faster than direct wavefunction
propagation on this small system. Full FCI is approximately 89 times faster,
and carrier-wavefunction propagation is approximately 141 times faster.

This is not a contradiction. The FCI Hamiltonian in this benchmark is an
extremely sparse 210-dimensional matrix, while every sparse-2RDM RHS call
currently performs:

- chart conditioning checks;
- carrier-matrix reconstruction;
- dense `29 x 29` commutators;
- Python-level graph traversal and coherence extraction.

The benchmark therefore establishes representational compression and exact
closure, not a finished high-performance solver.

Peak RSS is not yet a sensitive metric at this size because Python, NumPy, and
SciPy import overhead dominates all four processes. The variable counts are a
more meaningful memory indicator here.

## Consequence for the research claim

The defensible computational claim is:

> A certified sparse subset of 2-RDM coordinates reduces the propagated
> reduced system by a factor of 24 relative to the dense 2-RDM while retaining
> exact many-body closure and boundary-crossing capability.

The following claim is not supported yet:

> Sparse 2-RDM propagation is faster than direct wavefunction propagation.

## Optimization target

The timing profile identifies the next engineering theorem-to-algorithm step:
compile each chart into sparse linear and rational kernels so that the RHS
never constructs the complete carrier matrix. For a two-body Hamiltonian, it
should reconstruct only the coherences required by the selected observable
derivatives and required 3-RDM contractions.

A useful second benchmark must then scale across a family of carriers and
report the crossover against ambient FCI, rather than relying on this one
small sector.
