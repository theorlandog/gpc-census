# Hybrid rank-generator prototype and 2-RDM bridge

## Scope and result

The first exact regression target is the Borland--Dennis moment polytope
\(\wedge^3\mathbb C^6\).  The prototype takes its certified exact vertex
oracle, recovers the irredundant facets by rational arithmetic, identifies
Weyl-chamber versus nontrivial facets, and then analyzes the determinant
carrier selected by saturation of the recovered GPC.

It recovers the nontrivial inequality

\[
\lambda_4\le\lambda_5+\lambda_6
\]

together with the three chamber facets in reduced coordinates.  No
floating-point hull routine is used.

This is a successful **V-to-H regression harness**, not yet a complete
weight-to-facet Ressayre generator.  It validates exact canonicalization,
facet orientation, redundancy handling, constraint labeling, and the
downstream carrier audit needed before rank 11.

## Why this is the correct first milestone

The repository already has exact certified vertex oracles through rank ten.
A new rank generator must first reproduce their H-descriptions exactly.
Separating this regression layer from candidate generation catches:

- sign and chamber convention errors;
- trace-equation shifts;
- duplicate Weyl images;
- inherited versus genuinely new facets;
- facet-normal normalization errors.

The next module replaces the vertex input by admissible hyperplanes from
exterior-power weights, then applies the Ressayre trace and determinant tests.
The present output is the acceptance oracle for that module.

## Concrete 2-RDM result

Saturation of the Borland--Dennis facet

\[
D(\lambda)=2-(\lambda_1+\lambda_2+\lambda_4)=0
\]

has second-quantized selection operator

\[
\widehat D=2-(n_1+n_2+n_4).
\]

Its zero-eigenspace in the determinant basis is the nine-determinant carrier

\[
S_D=\{D:|D\cap\{1,2,4\}|=2\}.
\]

The exact audit finds:

- carrier size: \(9\);
- singleton-incidence rank: \(5\);
- non-toric phase complexity:
  \[
  \kappa(S_D)=9-5=4;
  \]
- pair-incidence rank: \(9\);
- collision-free 2-RDM graph: connected, with 27 unique edges.

Therefore every pure full-support state on this pinned carrier is certified
reconstructible from selected fixed-frame 2-RDM entries.

This is the first direct bridge between the GPC generator and the observable
manifold program:

> A nontrivial GPC facet can select a carrier with intrinsic interference
> complexity while still admitting exact order-two reconstruction.

## What this says about the 2-RDM question

Schubert/Bruhat/Ressayre machinery does not reconstruct a 2-RDM by itself.
It helps indirectly and systematically:

1. A facet generator produces an integral normal \(D\).
2. Pinning gives a selection operator \(\widehat D\).
3. Its zero-weight determinant support \(S_D\) is a canonical carrier.
4. Toric nullity \(\kappa(S_D)\) measures internal phase complexity.
5. Pair incidence and signed channel graphs test whether the 2-RDM resolves
   those phases.
6. The resulting carrier can support exact closed reduced dynamics.

This turns facet generation into an automatic source of physically meaningful
2-RDM test manifolds.

## New research hypothesis

For a GPC facet \(D\), define its **RDM observability profile**

\[
\mathcal O(D)=
\bigl(
|S_D|,
\operatorname{rank}A_1,
\kappa(S_D),
r_2(S_D),
q_{\min}(S_D)
\bigr),
\]

where \(r_2\) is pair-incidence rank and \(q_{\min}\) is the lowest RDM order
with a connected certified coherence graph.

The profile may distinguish facets that are spectrally similar but impose
very different wavefunction and dynamical reductions.

The Borland--Dennis facet has the provisional profile

\[
(9,5,4,9,2).
\]

Computing this profile for every known rank-\(\le10\) facet is now a concrete,
finite campaign.  It can reveal whether higher-level GPC generators
systematically produce easier exact 2-RDM problems.

## Remaining generator work

The complete hybrid generator still needs:

1. admissible hyperplane enumeration from exterior-power weights;
2. Weyl/Bruhat canonicalization and pruning;
3. Ressayre trace filtering;
4. exact sparse determinant nonvanishing witnesses;
5. relative-Schubert cross-certificates;
6. exact comparison against the rank-\(\le10\) H-oracles;
7. only then, rank-11 candidate generation.

No claim of a complete rank generator is made by this prototype.
