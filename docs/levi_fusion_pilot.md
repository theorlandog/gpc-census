# Levi-fusion pilot: first exact falsification and compression results

## Scope

This pilot computes the exact zero-grade Levi module for every published facet
of:

- `(3,7)`
- `(3,8)`
- `(4,8)`

The input rows are the repository's published constraints and the labels from
the ordered-representation-stability campaign.

For a row \(\tau\cdot\lambda\le0\), equal values of \(\tau\) define the Levi
blocks. The active particle allocations satisfy

\[
\sum_a n_a=N,
\qquad
\sum_a h_an_a=0.
\]

Each allocation contributes

\[
\bigotimes_a\bigwedge^{n_a}\mathbb C^{m_a}.
\]

## The naive conjecture is false

### Primitive does not imply multi-sector

All four primitive `(3,7)` facets have fusion rank one. They share the exact
zero-grade module

\[
\mathbb C^3\otimes\bigwedge^2\mathbb C^4,
\]

of dimension 18.

The single-sector primitive phenomenon persists at `(3,8)`: row 11 is
primitive and has

\[
\bigwedge^2\mathbb C^6\otimes\mathbb C^2,
\]

of dimension 30.

### Inherited does not imply single-sector

The `(4,8)` system has 4 padding rows, 4 frozen-core rows, 6 primitive rows, and
1 structural row. Every one of its 14 nonstructural facets has fusion rank two.

Therefore fusion rank alone does not distinguish primitive from inherited
facets in either direction.

## The positive result is mechanism compression

| System | Facets | Exact zero-grade signatures | Ratio |
| --- | ---: | ---: | ---: |
| `(3,7)` | 4 | 1 | 25.0% |
| `(3,8)` | 31 | 9 | 29.0% |
| `(4,8)` | 15 | 3 | 20.0% |

At `(3,8)`, the 27 primitive facets collapse to only eight exact Levi-module
signatures.

An independently written facet-index campaign has now extended the same
compression signal across all twelve published systems: facet counts from
1 to 161 collapse to only 1 to 21 sorted-`tau` shapes, with a shape-to-facet
ratio no larger than 0.29. Since the exact zero-grade Levi decomposition is
determined by `N` and the sorted `tau` shape, this supplies a census-wide
upper bound on the number of Levi mechanisms.

The revised hypothesis is therefore:

> The facet list may not stabilize, but the set of irreducible Levi mechanisms
> grows much more slowly than the list of chamber-facing inequalities.

## Complete `(3,8)` mechanism table

| Rows | Label | Fusion rank | Zero-grade dimension | Levi module sectors |
| --- | --- | ---: | ---: | --- |
| 1-4 | padding | 1 | 24 | \(\mathbb C^4\otimes\wedge^2\mathbb C^4\) |
| 5-10 | primitive | 2 | 21 | \(20+1\) |
| 11 | primitive | 1 | 30 | \(\wedge^2\mathbb C^6\otimes\mathbb C^2\) |
| 12-14 | primitive | 2 | 16 | \(4+12\) |
| 15-18 | primitive | 3 | 13 | \(6+1+6\) |
| 19-20 | primitive | 3 | 14 | \(1+12+1\) |
| 21-22 | primitive | 2 | 19 | \(1+18\) |
| 23-27 | primitive | 4 | 10 | \(2+2+4+2\) |
| 28-31 | primitive | 5 | 9 | \(4+1+1+1+2\) |

The sector dimensions are exact products of binomial coefficients.

## Revised theorem target

The first theorem should not say that primitive facets are exactly multi-sector
fusions. The data already refutes that.

A viable statement is:

> Every primitive fermionic facet is controlled by a bounded-complexity
> zero-grade Levi module, possibly a single tensor-product sector and possibly
> a direct sum of several sectors, together with a relative semi-invariant
> certifying the normal directions.

The next held-out test is ranks 9 and 10. Its predictions are frozen in
`docs/prereg_levi_fusion_heldout.md`.
