# Theorem candidates extracted from the Schubert/Bruhat audit

## Theorem A: carrier support-polytope theorem

Let \(S\subseteq\binom{[d]}N\) and let
\([\psi]=[\sum_{D\in S}c_D|D\rangle]\) have \(c_D\ne0\).
Under the diagonal complex orbital torus, the normalized orbit closure is the
projective toric variety associated with
\(P_S=\operatorname{conv}\{\chi_D:D\in S\}\).  The compact-torus moment image
is \(P_S\).

## Theorem B: non-toric phase-complexity formula

If \(A_1\) is singleton incidence and \(m=|S|\), then the quotient of the full
relative-phase torus by the effective diagonal orbital torus has dimension

\[
\kappa(S)=m-\operatorname{rank}A_1.
\]

Equivalently, \(\kappa(S)\) is the affine-circuit dimension of the support
configuration.

## Corollary C: cyclic windows

For the cyclic-window carrier \(C_{d,N}\),

\[
\kappa(C_{d,N})=\gcd(d,N)-1.
\]

Hence coprime cyclic carriers are support-simplex carriers.

## Proposition D: Grassmannian viability criterion

A carrier support can equal the nonzero Plücker support of a decomposable state
only if it is the basis set of a realizable matroid.  Failure of basis exchange
is therefore a certificate that the carrier has no full-support point in the
Grassmannian.

## Research conjecture E: RDM order versus internal phase complexity

For fixed particle number, the minimal RDM order required for generic
pure-state reconstruction is governed more strongly by the interaction of
\(\kappa(S)\) with the signed channel hypergraph than by carrier size alone.
This is a testable conjecture, not yet a theorem.
