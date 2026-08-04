# The four OPEN (3,11) candidates: family evaluation and design search

> **Superseded on the verdicts.** All four candidates are now REFUTED and the
> bracket is closed at 19 TRUE, 31 REFUTED, 0 OPEN. The four claimed level-five
> rows carry exact Ressayre certificates; see `docs/level5_ressayre_3_11.md`.
> This page is retained because its content is independent of that: it closes
> the *design* route, and for candidate 44 it does so by exact counting, which
> is what killed the falsification branch. Read the verdict language below as
> the state at the time of writing, including its tier language: the four
> level-five rows were CLAIMED then and are certified now.

## What this settles, and what it does not

Nothing here settles any of the four. The (3,11) bracket still stands at 19
TRUE, 27 REFUTED, 4 OPEN.

What it does is close the **design route** for all four, and for candidate 44
it does so with an exact certificate rather than a solver verdict. The
settlement previously recorded that cand 23 admits no one-hop-free design; the
same now holds for 26, 34 and 44, so **all four are INTERFERENCE if
attainable**.

## The four candidates

| idx | spectrum | distinct values | denominator |
|---:|---|---:|---:|
| 23 | `(6/7)^2 (1/7)^9` | 2 | 7 |
| 26 | `(7/9)^2 (11/27)^3 (1/27)^6` | 3 | 27 |
| 34 | `(11/17)^4 (1/17)^7` | 2 | 17 |
| 44 | `(1/2)^5 (1/12)^6` | 2 | 12 |

All four have trace 3. Three of the four are two-valued, so their 1-RDM has
the form `b I + (a - b) P` for a projector `P`.

## Family evaluation

Each candidate was evaluated in exact rational arithmetic against all thirteen
stable Grassmann partial families.

| idx | violates | tier | excess |
|---:|---|---|---|
| 23 | AK-RMK-4.2.1, KLY09-EXT-3 | CLAIMED | `1/7` |
| 26 | AK-RMK-4.2.1, KLY09-EXT-1, KLY09-EXT-3 | CLAIMED | `1/27` |
| 34 | AK-RMK-4.2.1, KLY09-EXT-1, KLY09-EXT-2 | CLAIMED | `1/17` |
| 44 | KLY09-EXT-1 | CLAIMED | `1/12` |

All four satisfy every PROVED family, which is what keeps them open.

Two refinements on the archive's summary, which records only that 23, 26 and
34 violate AK-RMK-4.2.1 while 44 violates the extended quadruple. First, each
candidate violates **every** row it violates by exactly `1/denominator`, so
each sits one unit outside the claimed hyperplane. Second, and more usefully,
**candidate 44 is the only one standing against a single claimed row**; the
others have two or three independent claims against them. Attaining 44 would
therefore falsify exactly one published claim rather than a cluster, which
sharpens the existing reason for treating it as the cheapest target.

## Design search

A support is one-hop-free when no two of its triples share exactly two
orbitals. That is precisely the condition that every pair of orbitals lies in
at most one block, so a one-hop-free support **is a partial Steiner triple
system**. For such a support the 1-RDM is diagonal with

\[
\lambda_p=\sum_{S\ni p}w_S,
\]

so a feasible support together with nonnegative weights summing to one is an
explicit attaining state and would certify the candidate TRUE.

The search is a mixed-integer program over all 165 triples of `[11]` with the
1,980 one-hop conflict pairs as exclusion constraints. It reports **infeasible
for all four**.

That verdict comes from CBC in floating point, so for 23, 26 and 34 it is
evidence rather than a certificate.

## Exact proof for candidate 44

Let `H` be the five orbitals at `1/2` and `L` the six at `1/12`, and let
`a, b, c, d` be the total weight on blocks meeting `H` in 3, 2, 1, 0 orbitals.
Counting incidences three ways gives

\[
3a+2b+c=\tfrac52,\qquad a+b+c+d=1,\qquad b+2c+3d=\tfrac12,
\]

and subtracting twice the second from the first leaves

\[
a-c-2d=\tfrac12,\qquad\text{so}\qquad a\ge\tfrac12 .
\]

Under the Steiner condition at most **two** blocks fit inside a five-set,
verified by enumeration: two 3-subsets of a 5-set meet in one or two points,
and a third block always forces a two-point meeting.

*Two inside blocks.* They meet in one orbital, say `x`. Both contain `x`, so
`lambda_x >= a >= 1/2 = lambda_x`, forcing `a = 1/2` and hence `c = d = 0`.
Every block therefore meets `H` in two or three orbitals, so each block with
`|S ∩ H| = 2` carries exactly one light orbital. The two inside blocks consume
six of the ten heavy pairs, leaving four; so at most four blocks can carry a
light orbital, against six light orbitals all needing `lambda = 1/12 > 0`.
Contradiction.

*One inside block.* Its weight is `a = 1/2`, saturating its three orbitals, so
no other block may touch them. The remaining two heavy orbitals admit only the
single block on the one free heavy pair, which carries one light orbital and
leaves the other five at zero. Contradiction.

*No inside block.* `a = 0` contradicts `a >= 1/2`.

Hence no one-hop-free design realises candidate 44.

## Where the gap stands

Closing the design route does not decide attainability: an interference state
whose 1-RDM is not diagonal in any coordinate basis remains possible, and that
is exactly what the design generator cannot reach. The existing contraction
audit leaves positive numerical floors on all four, which points the other way,
toward refutation, but is likewise evidence and not a certificate.

Stage 1 is not a route here. Extrapolating the measured closed-flat growth of
roughly thirty-fold per rank from the (3,7) and (3,8) lattices puts (3,11) at
order `3e10` flats, so the generator that proves rank 7 cannot reach rank 11.

The two live routes are unchanged: prove the claimed level-five inequality, in
which case all four are REFUTED, or attain candidate 44, which falsifies
KLY09-EXT-1. This layer only narrows what an attaining state would have to
look like.

*Resolved after writing: the first route landed. See
`docs/level5_ressayre_3_11.md`.*

## Reproduction

```sh
uv run python scripts/rank11_open_candidate_designs.py
uv run pytest -q tests/test_rank11_open_candidate_designs.py
```
