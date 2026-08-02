# Software architecture: Observable Dynamics Toolkit

## Product boundary

The research repository should remain the exact reference implementation and
artifact source. A separate package should provide a stable API, compiled
kernels, integrations, and visualization.

Suggested package name: `obdyn` (working name only).

## Modules

### `obdyn.carriers`

- determinant-carrier data model;
- frozen-core and active-space constructors;
- cyclic-window and product constructors;
- closure enlargement under supplied Hamiltonians;
- canonical hashing and provenance.

### `obdyn.certify`

- pair-incidence exact rank and minor certificates;
- signed RDM channel enumeration;
- collision-free and determinacy-closure graphs;
- chart bank generation;
- independent certificate verifier.

### `obdyn.geometry`

- observable chart maps;
- active-support boundary strata;
- symplectic forms and Jacobians;
- chart overlap maps;
- conditioning estimates;
- atlas optimization.

### `obdyn.compile`

- sparse reconstruction DAG;
- needed-coherence analysis from Hamiltonian sparsity;
- generated NumPy/Numba/JAX kernels;
- exact and floating backends;
- cache keyed by carrier and Hamiltonian hashes.

### `obdyn.propagate`

- regular minimal-chart integrator;
- redundant Cartesian boundary mode;
- event detection and chart switching;
- implicit midpoint or variational symplectic integrator;
- scheduled and autonomous Hamiltonians;
- reconstruction-on-demand for higher RDMs.

### `obdyn.integrations`

Initial adapters:

- PySCF: determinants, CI vectors, one- and two-electron integrals;
- OpenFermion: `FermionOperator` conversion;
- FCIDUMP import/export;
- NumPy/SciPy reference backend.

Later adapters:

- Psi4;
- Q-Chem data interchange;
- Julia/ITensors or Python tensor-network bridges;
- quantum-hardware RDM measurement inputs.

### `obdyn.visualize`

- carrier graph;
- incidence rank and minor diagnostics;
- coherence graph and chart bank;
- weight simplex trajectories;
- phase-torus and boundary events;
- benchmark dashboards;
- provenance and certificate explorer.

### `obdyn.verify`

- direct wavefunction cross-check;
- RDM contraction cross-check;
- symplectic residual;
- conservation laws;
- exact certificate replay;
- randomized property testing.

## API sketch

```python
carrier = obdyn.carriers.from_determinants(determinants, particle_number=n)
certificate = obdyn.certify.order_two(carrier)
atlas = obdyn.geometry.build_atlas(certificate, redundancy=3)
kernel = obdyn.compile.hamiltonian(atlas, hamiltonian, backend="numba")
result = obdyn.propagate.solve(
    atlas,
    kernel,
    initial_rdm2=gamma2,
    t_span=(0.0, 10.0),
    symplectic=True,
)
```

## Engineering standards

- Python reference implementation plus compiled backend.
- Exact-integer certificate generation.
- Standalone verifier with no imports from generator modules.
- Deterministic JSON schema with SHA-256 provenance.
- Semantic versioning.
- Typed APIs and documented numerical tolerances.
- Property-based tests for fermionic signs and chart overlaps.
- Benchmark results stored as artifacts, never hard-coded claims.
- No proprietary lock-in for certificate verification.

## Open-source versus commercial split

Open source:

- theorem reference implementation;
- certificate formats and verifiers;
- basic constructors;
- small-system propagation;
- paper reproduction.

Potential commercial value:

- optimized compiled backends;
- large-scale atlas optimization;
- enterprise integrations;
- GPU/distributed kernels;
- workflow management;
- validated support packages;
- proprietary carrier discovery heuristics.

Keeping the verifier and certificate schema open is strategically valuable:
users can trust results without receiving every optimization or integration.
