# Pre-registration: held-out Levi-fusion audit

## Status and separation of discovery from test

The systems `(3,7)`, `(3,8)`, and `(4,8)` were inspected first as a pilot.
Their outcomes are discovery data and are not counted as evidence for the
predictions below.

The held-out systems are the primary published systems at ranks 9 and 10:

- `(3,9)`
- `(3,10)`
- `(4,9)`
- `(4,10)`
- `(5,10)`

The measuring script is `scripts/levi_fusion_audit.py`. This document and the
script should be committed together before the held-out artifact is generated.

## Exact object being measured

For a primitive homogeneous GPC row

\[
\tau\cdot\lambda\le 0,
\]

write the distinct values of \(\tau\) as \(h_1<\cdots<h_r\), with
multiplicities \(m_1,\ldots,m_r\). The associated Levi decomposition is

\[
\mathbb C^d=E_1\oplus\cdots\oplus E_r,
\qquad
\dim E_a=m_a.
\]

The zero-grade exterior module is

\[
\left(\bigwedge^N\mathbb C^d\right)_0
=
\bigoplus_{\substack{0\le n_a\le m_a\\
\sum_a n_a=N\\
\sum_a h_an_a=0}}
\bigotimes_a\bigwedge^{n_a}E_a.
\]

The number of direct summands is the **fusion rank**.

Rows are labeled `PADDING`, `FROZEN-CORE`, `PRIMITIVE`, or `STRUCTURAL` by the
independently generated ordered-representation-stability artifact.

## Pilot findings already observed

These are not predictions:

1. The claim "primitive if and only if multi-sector" is false.
   Every one of the four primitive `(3,7)` facets has fusion rank one.

2. The converse claim "inherited implies single-sector" is also false.
   Every nonstructural `(4,8)` row has fusion rank two, while that population
   contains both inherited and primitive rows.

3. Strong compression remains visible:
   - `(3,7)`: 4 rows, 1 exact zero-grade signature.
   - `(3,8)`: 31 rows, 9 signatures.
   - `(4,8)`: 15 rows, 3 signatures.

4. The largest pilot fusion rank is 5, and the largest number of distinct
   \(\tau\)-levels is \(6=2N\), both at `(3,8)`.

The revised target is therefore bounded **Levi-module complexity**, not a
binary single-sector versus multi-sector classification.

## Held-out predictions

### P1. Bounded number of Levi blocks

For every held-out nonstructural row,

\[
r\le 2N.
\]

A violation is any row with more than \(2N\) distinct values of \(\tau\).

### P2. Mechanism compression, confirmatory rather than predictive

At each held-out primary system, the number of distinct unoriented exact
zero-grade signatures is at most half the number of nonstructural rows.

After this pilot was designed, the independently written facet-index
campaign landed a sorted-`tau` shape census on `main`. It already reports
ratios at or below 0.22 at the rank-9 and rank-10 systems. Because the
zero-grade Levi module is a deterministic function of `N` and the sorted
`tau` shape, P2 is now a confirmatory cross-realization check and not new
evidence.

### P3. Primitive rows carry greater fusion complexity at fixed \(N=3\)

At both `(3,9)` and `(3,10)`, the mean fusion rank of primitive rows is strictly
greater than the mean fusion rank of rows labeled `PADDING` or `FROZEN-CORE`.

### P4. Levi signature alone does not decide primitivity

At every held-out primary system containing both primitive and inherited rows,
at least one unoriented zero-grade signature occurs in both classes.

A failure would make the Levi signature more discriminating than the pilot
suggests.

### P5. Slow growth of \(N=3\) fusion rank

For `(3,d)` with \(d=9,10\),

\[
\max f(\tau)\le d-3.
\]

Thus the bounds are 6 and 7.

### P6. Sparse zero-grade module

For every held-out nonstructural row,

\[
\frac{\dim(\bigwedge^N\mathbb C^d)_0}
{\binom dN}
\le 0.55.
\]

## Interpretation rules

- A failure is retained as a finding.
- Passing through rank 10 is evidence for bounded complexity, not a theorem for
  all \(d\).
- Coefficient permutations are not treated as equivalent valid inequalities.
  The signature groups Levi modules only.
- State-support sector occupancy is reported separately and only in shipped
  frames that pass the attainer gate.
- No local Levi or slice calculation is promoted to a global GPC proof without
  a bound Ressayre certificate.
