# Compiled sparse 2-RDM benchmark

## Question

Does compiling the observable-manifold right-hand side into Hamiltonian-specific sparse kernels turn the reduced variable-count advantage into a runtime advantage?

The original multi-chart implementation reconstructed the complete 29 by 29 carrier matrix at every right-hand-side evaluation. The transport Hamiltonian has only eight nonzero carrier matrix entries, so most of that reconstruction was unnecessary.

## Compiled kernel

For each propagated weight and certified coherence, the compiler records only the carrier coherences that occur in

\[
\dot X_{ab}=i\sum_k(H_{ak}X_{kb}-X_{ak}H_{kb}).
\]

A selected chart first reconstructs one root-coherence vector in linear time. Any needed coherence is then evaluated on demand as

\[
X_{ab}=\frac{\overline{X_{ra}}X_{rb}}{w_r}.
\]

The compiled right-hand side does not assemble the dense carrier matrix and does not traverse the full channel graph after chart selection.

## Results

Seven repeated integrations used the same 167 real observable variables, DOP853 tolerances, 26 accepted steps, and 314 right-hand-side evaluations.

| Method | Median runtime | Maximum carrier-matrix error |
|---|---:|---:|
| Original sparse 2-RDM | 0.601 s | 1.40e-12 |
| Compiled sparse 2-RDM | 0.418 s | 1.40e-12 |
| Dense 2-RDM, previous benchmark | 0.683 s | 2.1e-12 |
| Full FCI, previous benchmark | 0.00697 s | 5.6e-17 |
| Carrier wavefunction, previous benchmark | 0.00439 s | 5.6e-17 |

The compiler gives:

- 1.44-fold speedup over the original sparse observable implementation;
- 1.63-fold speedup over dense 2-RDM propagation;
- unchanged numerical accuracy;
- no synchronization fallbacks during this benchmark trajectory.

Two qualifications belong with those speedups. First, they are host specific
and the margin is not stable: on a second, slower machine the same script
measured 1.17, and an interleaved A/B run alternating the two right-hand sides
measured 1.23 there. The direction is consistent, the factor is not. Second,
the 1.63 figure divides a fresh-process dense timing from the previous
benchmark by an in-process compiled timing from this one, so it compares
across two harnesses; re-running the dense method on the same host as the
compiled one gave 1.34.

What does hold exactly, and is host independent, is equivalence: the compiled
kernel and the reference right-hand side agree to `4e-18` over the whole
trajectory including both amplitude zeros, and both methods take 314
evaluations and 26 accepted steps and end at a bit-identical error. The
speedup is therefore a genuine reduction in work per evaluation rather than a
different amount of work. The compiler exploits that the carrier Hamiltonian
has 8 nonzero entries out of 841.

## Interpretation

Compilation converts the geometric result into a better reduced-density-matrix algorithm, but not yet into a wavefunction-beating algorithm on this small system. The remaining approximately 60-fold gap to full FCI is dominated by Python-level chart reconstruction and propagation of a redundant 167-variable atlas, while the ambient FCI Hamiltonian is exceptionally sparse and only 210 dimensional.

The benchmark supports the limited claim that the compiled observable method is now both smaller and faster than dense 2-RDM propagation. It does not support a general speed advantage over direct wavefunction propagation.

The next optimization should remove redundancy by propagating one 56-dimensional regular chart and activating backup Cartesian coherences only near a chart boundary. A compiled extension should then be benchmarked across a scalable family where the ambient FCI dimension grows combinatorially.
