# Exact reduced evolution on the 29-determinant transport carrier

The four-determinant hopping closure of `S(5,2) x S(5,2)` is invariant under

\[
H=(a_1^\dagger a_0)(a_7^\dagger a_5)+(a_0^\dagger a_1)(a_5^\dagger a_7).
\]

It has 29 determinants, full-column-rank pair incidence, and a connected
collision-free 2-RDM coherence graph. The script
`scripts/transport_29_reduced_evolution.py` evolves a full-support state under
this non-diagonal two-body Hamiltonian and, at every sample time:

1. contracts the exact carrier rank-one matrix to the 2-RDM;
2. recovers all weights from 29 selected diagonal pair occupations;
3. recovers 28 tree coherences from collision-free 2-RDM channels;
4. reconstructs the complete carrier rank-one matrix by path products;
5. reconstructs the 3-RDM and compares it with direct contraction.

For 17 times on `0 <= t <= 4`, the largest observed errors were

- carrier matrix: `3.78e-14`;
- reconstructed 2-RDM: `4.68e-14`;
- reconstructed 3-RDM: `3.78e-14`.

The minimum carrier weight reached `5.71e-5`.  Regenerating the artifact
reproduces every headline number exactly; the per-sample errors move in the
last bits (order `1e-16`) because `expm` is not bit-deterministic, so the file
is a frozen measurement rather than a byte-reproducible artifact. The largest errors occur near
those small-weight chart regions, as predicted by the denominators in the
path-product inverse. This is numerical roundoff, not closure error.

## Did this make 2-RDM calculations easier?

There are three different comparisons.

### Against the ambient full-CI problem: yes

The model lives in `wedge^4 C^10`, whose fixed-particle sector has 210 complex
wavefunction amplitudes. The invariant carrier requires only 29 complex
amplitudes, and the certified reduced chart uses 29 real weights plus 28
complex tree coherences. The exact 3-RDM is obtained without storing or
evolving the 120 by 120 3-RDM.

### Against propagating the known 29-component carrier wavefunction: not yet

A normalized pure state on the carrier has 56 real physical degrees of
freedom. The simple certificate chart stores 85 real numbers, so it is
redundant. A full 45 by 45 2-RDM is much larger still. For this small known
carrier, direct 29-component Schrödinger propagation is cheaper.

### For reduced-density-matrix methods: substantially

The BBGKY equation for the 2-RDM normally requires the 3-RDM. Here the 3-RDM
is an exact algebraic function of selected 2-RDM coordinates. This removes the
closure approximation and turns the model into exact ground truth for TD-2RDM
algorithms. The computational win becomes plausible for scalable carriers
whose dimension is far below their ambient FCI sector, provided the reduced
chart is implemented without expanding the full 2-RDM or full carrier matrix.

The current script validates exact closure but still reconstructs dense
matrices for verification. The next optimization is a sparse reconstruction
DAG that evaluates only the 3-RDM contractions requested by the BBGKY right
hand side.
