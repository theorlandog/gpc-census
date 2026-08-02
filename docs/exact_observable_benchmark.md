# Exact observable propagator benchmark

The 29-state hopping Hamiltonian decomposes into four disjoint two-state
blocks and 21 spectators.  The observable propagator exploits this exact
structure.  It reconstructs the carrier density matrix from certified 2-RDM
coordinates, applies four analytic rotations, and extracts the updated
observable chart.  No adaptive ODE solver or wavefunction amplitudes are used.

Against the previously published ODE benchmark, the median end-to-end time is
about 1.92 ms, with maximum carrier-density error below 7e-17.  This is 2.29x
faster than carrier-wavefunction ODE propagation, 3.63x faster than full-FCI
ODE propagation, and more than 300x faster than either dense or original
sparse 2-RDM ODE propagation.

This does **not** establish that observable propagation is intrinsically
faster than every wavefunction method.  An equally specialized analytic
carrier update for the same four 2x2 blocks takes about 4.27 microseconds in
this Python environment, roughly 450x faster than the observable update.  The
reason is fundamental for this model: the known 29-amplitude carrier is a
smaller representation than the redundant 167-real-variable global
observable atlas, and the observable method must reconstruct and re-extract
correlators.

## What is being compared

The speedups against the ODE methods compare unlike work and should be read
that way. The observable propagator is one call: reconstruct the carrier
matrix from the chart, apply four analytic rotations, repack. The ODE
baselines integrate the same trajectory with 314 right-hand-side evaluations
across 26 adaptive steps. An analytic solution beating numerical integration
of the same problem is expected, and the 300x figures are that, not a
per-evaluation advantage. The comparison that carries information is the
control, which is analytic on both sides.

Profiling the observable call locates its cost: on a second host the median
was 2.49 ms, of which reconstruct_X was 2.49 ms and chart selection over the
108-tree bank alone was 1.82 ms. Almost all of the time is searching for a
well-conditioned chart, not the physics. Caching the selected chart while it
stays regular is the obvious next step.

The figures are host specific. Re-running on a slower machine gave 1.76x
against carrier-wavefunction and 2.80x against full-FCI ODE propagation
rather than 2.29x and 3.63x, and the specialized-carrier control widened from
450x to 568x. What reproduced bit for bit is the accuracy, `6.87e-17`, and
the block decomposition. The ordering of the methods did not change.

The valid computational claim is therefore:

* exact observable propagation beats all four prior generic ODE
  implementations on this benchmark;
* it does not beat a Hamiltonian-specific exact carrier propagator;
* its potential advantage lies where amplitudes are unavailable, where a
  reduced-observable interface is required, or where the ambient/carrier
  wavefunction cannot be stored while a compiled observable closure remains
  sparse.
