# Stage 1 checkpoint: first-principles reproduction of (3,6)

## Status and scope

This checkpoint is **proved for `(N,d)=(3,6)`**. It does not claim that the
rank-7 through rank-10 systems have been regenerated, and it does not yet
produce rank-11 constraints.

The implementation is in `src/gpc_census/generation/`. The deterministic
entry point is `scripts/stage1_generator.py`, and the complete replayable
artifact is `results/data/stage1_ressayre_3_6.json`. Generation reads neither
the stored constraint table nor the stored vertex ledger. Those data are used
only by a post-generation regression test.

The checkpoint combines three logically separate routes:

1. exhaustive admissible-hyperplane enumeration and exact Ressayre
   determinants for valid outer inequalities;
2. four flag-Schubert divided-difference certificates as an independent
   Klyachko cross-check;
3. exact plethysm inner points and exact H-to-V enumeration for a completeness
   sandwich.

The Ressayre convention follows Definition 11 of Bürgisser, Christandl,
Mulmuley, and Walter: a certificate `(H,z)` proves

\[
H\mathbin{\cdot}\lambda\ge z.
\]

See [Membership in moment polytopes is in NP and coNP](https://arxiv.org/pdf/1511.03675)
and [Inequalities for Moment Cones of Finite-Dimensional Representations](https://arxiv.org/pdf/1410.8144).

## Exterior representation data

The weights of \(\wedge^N\mathbb C^d\) are the multiplicity-free determinant
weights

\[
\Phi=\left\{\omega_I=\sum_{i\in I}e_i:|I|=N\right\}.
\]

Negative type-A roots are represented in code by zero-based pairs `(i,j)`,
\(i<j\), meaning \(\alpha=e_j-e_i\). The root operator is
\(E_{ji}=a_j^\dagger a_i\). Its action on every determinant is computed with
the exact exterior sign \(\pm1\).

The central ambiguity in an affine pair is removed by requiring
\(\sum_i H_i=0\). Adding a scalar to every component of \(H\), with the
corresponding shift of \(z\), leaves all weight levels unchanged.

## Complete admissible-hyperplane enumeration at rank 6

The exterior weights lie in the trace hyperplane of affine dimension
\(r=d-1\). An admissible affine hyperplane has dimension \(r-1=d-2\), so it
contains \(d-1\) affinely independent weights. The generator enumerates every
\((d-1)\)-subset of weights. For each subset it computes the normal from the
signed maximal minors of the matrix containing \(d-2\) weight differences
and the row \((1,\ldots,1)\).

This enumeration is complete: every admissible hyperplane has an affine
basis of exactly this size, that basis occurs in the enumeration, and its
one-dimensional normal space produces the same primitive pair `(H,z)`.

For \(\wedge^3\mathbb C^6\), the exact counts are:

| Stage | Count |
|---|---:|
| Unoriented admissible hyperplanes | 362 |
| Oriented candidates | 724 |
| Trace-condition survivors | 64 |
| Nonzero determinant certificates | 21 |

## Ressayre certificate and replay

For each oriented pair `(H,z)`, define

\[
\Phi_{=}=\{\omega:\omega(H)=z\},\qquad
\Phi_{<}=\{\omega:\omega(H)<z\},
\]

and

\[
R_< = \{e_j-e_i:i<j,\ H_j-H_i<0\}.
\]

The trace filter requires \(|\Phi_<|=|R_<|\). The tangent matrix has rows
indexed by \(\Phi_<\), columns indexed by \(R_<\), and entries equal to the
signed variables of on-plane determinants that the corresponding root
operator maps into the row weight.

The determinant is expanded exactly with SymPy. If it is nonzero, the
generator constructs an integer specialization one variable at a time. A
nonzero polynomial of degree at most \(k\) in a variable cannot vanish at all
\(k+1\) integers in `0..k`, so this procedure is deterministic and complete.
The artifact stores the specialization and its nonzero integer determinant.
The verifier reconstructs the weights, level sets, roots, exterior signs, and
integer matrix before replaying that determinant.

The simplest nontrivial certificate is

\[
H=(-1,-1,1,-1,1,1),\qquad z=-1.
\]

It is equivalent to

\[
\lambda_1+\lambda_2+\lambda_4\le2.
\]

Its below-plane set is the single determinant \(124\), its negative-root set
is the single root \(e_4-e_3\), and

\[
E_{43}|123\rangle=|124\rangle.
\]

The determinant matrix is therefore \([X_{123}]\), and the stored
specialization gives determinant one. Reversing this pair fails the trace
condition, which is pinned as a sign-convention regression.

The validity implication does not need a maximal-dimension hypothesis:
Proposition 3.13 of Vergne and Walter proves that surjectivity of the tangent
map gives a valid inequality unconditionally. The trace equality plus a
nonzero square determinant gives an isomorphism, hence surjectivity.

## Affine equalities and canonical GPC

Three admissible hyperplanes pass the Ressayre test in both orientations.
Both halfspaces are valid, so each pair is an exact equality. After a trace
shift and primitive normalization they are

\[
\lambda_1+\lambda_6=1,\qquad
\lambda_2+\lambda_5=1,\qquad
\lambda_3+\lambda_4=1.
\]

Their augmented row RREF is

\[
\begin{pmatrix}
1&0&0&0&0&1&-1\\
0&1&0&0&1&0&-1\\
0&0&1&1&0&0&-1
\end{pmatrix}.
\]

Every remaining inequality is reduced modulo this equality row space. The 15
one-sided Ressayre certificates collapse to only two affine keys:

\[
-\lambda_6\le0
\]

and

\[
\lambda_4-\lambda_5-\lambda_6\le0.
\]

The first is structural nonnegativity. The second is the unique nontrivial
Borland-Dennis GPC. Facetness is checked separately by the exact affine rank
of its tight outer vertices.

## Independent Schubert cross-check

The implementation carries four finite divided-difference certificates. The
first three certify the raw pair rows

\[
\lambda_1+\lambda_6\le1,\quad
\lambda_2+\lambda_5\le1,\quad
\lambda_3+\lambda_4\le1.
\]

They use test spectrum \(a=(1,1,0,0,0,0)\), induced levels
\(2^4,1^{12},0^4\), and a cyclic target permutation of length four. Their
source permutations and reduced divided-difference words are stored in the
artifact, and every exact coefficient is one.

For the nontrivial row, take

\[
a=(1,1,1,0,0,0),\qquad v=(1,2,4,3,5,6).
\]

The induced exterior spectrum has levels \(3^1,2^9,1^9,0^1\). The pullback
polynomial is \(P=x_1+x_2+x_3\), and the exact divided difference is

\[
\partial_3P=1.
\]

This independently certifies \(\lambda_1+\lambda_2+\lambda_4\le2\), which
reduces to the same affine GPC as the Ressayre route.

This coefficient-one result is a rank-6 fact, not a general validity
criterion. The general Klyachko condition is nonvanishing, and no later-rank
implementation may discard a nonzero coefficient merely because it is not
one.

## Completeness theorem for (3,6)

**Theorem.** The generated equalities, structural chamber inequalities, and
the generated nontrivial GPC cut out exactly the pure-state
\(\wedge^3\mathbb C^6\) occupation polytope.

**Proof.** Let \(P\) be the true moment polytope and \(O\) the generated outer
polytope. Every nonstructural row of \(O\) has a replayed Ressayre determinant,
and both sides of every generated equality are valid Ressayre inequalities.
The Weyl and Pauli rows are standard. Therefore \(P\subseteq O\).

The exact plethysm calculation of \(h_m[e_3]\) for \(m\le6\) produces 20
certified attainable spectra, so their convex hull lies in \(P\). Exact
active-set enumeration of \(O\) gives four vertices:

\[
\begin{aligned}
&(1,1,1,0,0,0),\\
&(1,\tfrac12,\tfrac12,\tfrac12,\tfrac12,0),\\
&(\tfrac34,\tfrac34,\tfrac12,\tfrac12,\tfrac14,\tfrac14),\\
&(\tfrac12,\tfrac12,\tfrac12,\tfrac12,\tfrac12,\tfrac12).
\end{aligned}
\]

Every one is present in the plethysm-certified inner set. The Weyl and Pauli
bounds make \(O\) bounded, hence \(O\) is the convex hull of these vertices.
Thus \(O\subseteq P\), and \(P=O\). \(\square\)

This sandwich avoids applying a maximal-dimension Ressayre completeness
theorem directly to the lower-dimensional Borland-Dennis polytope.

## Complexity and scaling boundary

Let \(W=\binom dN\) and \(r=d-1\).

- Brute-force admissible enumeration considers \(\binom Wr\) subsets. Each
  signed-minor normal costs polynomial time in \(d\), bounded here by
  \(O(d^4)\).
- Each oriented trace check costs \(O(Wd+d^2)\).
- A certificate replay with tangent size \(k\) costs \(O(k^3+Wdk)\) integer
  operations after the specialization is supplied.
- Symbolic determinant discovery can be exponential in \(k\). It is a
  discovery step, while the emitted integer specialization is a polynomial
  verification certificate.
- Exact active-set vertex enumeration in affine dimension \(s\), with \(F\)
  distinct outer rows, costs \(O(\binom Fs s^3)\) rational operations.
- The plethysm inner route grows with the partition counts of \(m\) and
  \(Nm\); it is an independent convergence certificate, not a scalable outer
  generator by itself.

The \(\binom Wr\) term is already unacceptable beyond the first ranks. The
next milestone must add Weyl-orbit canonicalization and a completeness-proved
Bruhat or trace pruning layer before attempting rank 8. Empirical edge height,
depth, inheritance, or cyclic-target patterns may prioritize candidates, but
they may not discard candidates until a theorem proves the pruning safe.

## Reproduction

```sh
uv run python scripts/stage1_generator.py
python scripts/verify_stage1_3_6_standalone.py
uv run pytest tests/test_stage1_generator.py
```

The test suite checks exact counts, every determinant replay, mutation
rejection, exterior-action signs, Schubert agreement, the four outer vertices,
the plethysm sandwich, byte-identical artifact regeneration, and only then the
post-generation agreement with the legacy lookup.
