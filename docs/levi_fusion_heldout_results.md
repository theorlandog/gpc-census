# Held-out Levi-mechanism audit

## Status

**ROW-LEVEL HELD-OUT TEST COMPLETE.**

The audit covers every published facet at

\[
(3,9),\quad(3,10),\quad(4,9),\quad(4,10),\quad(5,10).
\]

The source is the exact primitive-\(\tau\) and inheritance ledger in
`results/data/ordered_representation_stability.json` at Git blob
`323264c12907f8f9a8f636be8db428a05b8ffc15`.

For a row \(\tau\cdot\lambda\le0\), equal entries of \(\tau\) define the Levi
blocks. If their values and multiplicities are \(h_a,m_a\), the zero-grade
exterior module is

\[
\left(\bigwedge^N\mathbb C^d\right)_0
=
\bigoplus_{\substack{0\le n_a\le m_a\\
\sum_a n_a=N\\
\sum_a h_an_a=0}}
\bigotimes_a\bigwedge^{n_a}\mathbb C^{m_a}.
\]

The **fusion rank** is the number of active allocations \((n_a)\).

## Scorecard

| Prediction | Result | Exact outcome |
| --- | --- | --- |
| P1: number of Levi blocks \(r\le2N\) | **FAIL** | Violated at \((3,10)\) and \((4,10)\) |
| P2: at least 2x mechanism compression | **PASS, confirmatory** | 8/52, 16/93, 12/59, 14/124, 20/160 |
| P3: primitive \(N=3\) rows have higher mean fusion rank | **PASS** | Both held-out ranks |
| P4: primitive and inherited rows share mechanisms | **PASS** | Every held-out system |
| P5: \(f_{\max}(3,d)\le d-3\) | **PASS** | \(5\le6\), \(7=7\) |
| P6: zero-grade density at most 55% | **PASS** | Maximum is \(7/15\) |

The interim live update described the two P1 system violations as two failed
hypotheses. The exact scorecard has **one failed hypothesis, P1, at two
systems**. P2 through P6 pass, with P2 counted only as confirmation because an
independent sorted-\(\tau\) census had already established the compression
pattern.

## Per-system results

| System | Facets | Nonstructural mechanisms | Max blocks \(r\) | Max fusion \(f\) | Max zero-grade dimension | Ambient dimension | Density |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| \((3,9)\) | 52 | 8 | 6 | 5 | 31 | 84 | \(31/84\) |
| \((3,10)\) | 93 | 16 | **8** | 7 | 56 | 120 | \(7/15\) |
| \((4,9)\) | 60 | 12 | 7 | 6 | 45 | 126 | \(5/14\) |
| \((4,10)\) | 125 | 14 | **10** | 9 | 50 | 210 | \(5/21\) |
| \((5,10)\) | 161 | 20 | 10 | 9 | 60 | 252 | \(5/21\) |

## P1 falsification: bounded Levi blocks are not the theorem

At \((3,10)\), four primitive facets have eight distinct \(\tau\)-levels,
exceeding \(2N=6\). One exact witness is

\[
\tau=(3,2,-5,-3,-2,-1,-2,-1,0,1).
\]

Its Levi multiplicities are

\[
(-5)^1,\ (-3)^1,\ (-2)^2,\ (-1)^2,\ 0^1,\ 1^1,\ 2^1,\ 3^1,
\]

with fusion rank 7 and zero-grade dimension 12.

The stronger counterexample is at \((4,10)\):

\[
\tau=(0,1,2,-3,3,-6,-4,-5,-2,-1).
\]

Every entry is distinct, so the Levi has ten one-dimensional blocks:

\[
r=10>2N=8.
\]

Nevertheless the active zero-grade module has only nine one-dimensional
sectors:

\[
f=9,\qquad
\dim\left(\bigwedge^4\mathbb C^{10}\right)_0=9.
\]

Thus the number of Levi blocks can be as large as the orbital rank even while
the active zero-grade representation remains tiny.

## P2: mechanism compression survives strongly

The number of chamber-facing facets grows much faster than the number of exact
zero-grade Levi mechanisms:

\[
52\to8,\qquad
93\to16,\qquad
59\to12,\qquad
124\to14,\qquad
160\to20.
\]

The correct compression object is therefore not a bounded number of Levi
blocks. It is the much smaller set of **active zero-grade allocation
mechanisms**, together with chamber placement.

## P3: primitive rows are more fusion-complex at fixed \(N=3\)

At \((3,9)\),

\[
\overline f_{\mathrm{primitive}}
=
\frac{80}{21}
>
\frac{88}{31}
=
\overline f_{\mathrm{inherited}}.
\]

At \((3,10)\),

\[
\overline f_{\mathrm{primitive}}
=
4
>
\frac{42}{13}
=
\overline f_{\mathrm{inherited}}.
\]

So fusion complexity does not classify primitivity, but it does separate the
classes statistically at both genuinely held-out \(N=3\) ranks.

## P4: Levi type does not determine primitivity

The numbers of zero-grade signatures shared by primitive and inherited rows
are

\[
3,\ 2,\ 3,\ 5,\ 7
\]

for

\[
(3,9),\ (3,10),\ (4,9),\ (4,10),\ (5,10).
\]

This confirms that the same Levi module can support both a transported facet
and a genuinely new chamber-facing facet. The missing datum must include
chamber placement and the normal semi-invariant.

## P5: fusion rank still grows slowly

For the held-out \(N=3\) systems,

\[
f_{\max}(3,9)=5,\qquad
f_{\max}(3,10)=7.
\]

The preregistered bound

\[
f_{\max}(3,d)\le d-3
\]

passes and is saturated at \(d=10\). Two points are not a theorem, but this is
now the more plausible bounded-complexity quantity than the number of Levi
blocks.

## P6: zero-grade modules remain sparse

The largest observed fractions of the ambient exterior representation are

\[
\frac{31}{84},\quad
\frac{7}{15},\quad
\frac{5}{14},\quad
\frac{5}{21},\quad
\frac{5}{21}.
\]

Every system passes the 55% bound. More importantly, the rows with the largest
number of Levi blocks have some of the smallest zero-grade modules. The
\((4,10)\) all-distinct witness has only 9 active dimensions out of 210.

## Revised theorem target

The fixed bounded-block conjecture is false. The viable representation-theory
target is now:

> Every primitive fermionic facet is controlled by a sparse, parametrized
> zero-grade Levi mechanism, together with dominant-chamber placement and a
> relative semi-invariant of the normal representation.

A possible complexity measure is not \(r\), the number of Levi blocks, but a
tuple such as

\[
\left(
f,\quad
\dim W_0,\quad
\text{allocation-hypergraph width},\
\text{normal semi-invariant complexity}
\right).
\]

The data supports a picture in which \(r\) may grow with \(d\), while \(f\) and
\(\dim W_0\) grow much more slowly.

This local mechanism program remains viable even though the separate global
Khovanskii/SAGBI campaign finds that the complete stable highest-weight
semigroup does not have a small rank-independent generating set in the tested
range.

## State-sector phase

The optional support-sector audit on shipped exact states was not run in this
connector-only execution. It is not part of the scored row-level
pre-registration. It requires a local checkout containing the complete
`states.jsonl` and `attainer_gate_audit.json`.

The full repository script already contains that path and labels it correctly:
support-sector occupancy is used only when the shipped state passes the
diagonal-1-RDM attainer gate, and it is not promoted to a basis-free fiber
claim.

## Verification

```sh
python scripts/verify_levi_fusion_heldout_standalone.py
pytest -q tests/test_levi_fusion_audit.py
```

Expected output:

```text
held-out Levi audit: PASS
P1: FAIL at (3,10), (4,10)
P2-P6: PASS
```
