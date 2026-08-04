# The (3,11) bracket closes: 19 TRUE, 31 REFUTED, 0 OPEN

## Result

The four inequalities that kept candidates 23, 26, 34 and 44 OPEN carry exact
Ressayre certificates. They were CLAIMED, taken from Klyachko 0904.2009 and
CMP Remark 4.2.1 without proof; they are now PROVED in this repository, from
the repository's own generation package.

| family | row | census tier before | after |
|---|---|---|---|
| AK-RMK-4.2.1 | `l1+l2+l4+l7+l11 <= 2` | CLAIMED | PROVED |
| KLY09-EXT-1 | `l2+l3+l4+l5+l11 <= 2` | CLAIMED | PROVED |
| KLY09-EXT-2 | `l1+l3+l4+l6+l11 <= 2` | CLAIMED | PROVED |
| KLY09-EXT-3 | `l1+l2+l5+l6+l11 <= 2` | CLAIMED | PROVED |

Each of the four remaining candidates violates at least one of them, so each
is outside the polytope:

| idx | spectrum | cut by |
|---:|---|---|
| 23 | `(6/7)^2 (1/7)^9` | AK-RMK-4.2.1, KLY09-EXT-3 |
| 26 | `(7/9)^2 (11/27)^3 (1/27)^6` | AK-RMK-4.2.1, KLY09-EXT-1, KLY09-EXT-3 |
| 34 | `(11/17)^4 (1/17)^7` | AK-RMK-4.2.1, KLY09-EXT-1, KLY09-EXT-2 |
| 44 | `(1/2)^5 (1/12)^6` | KLY09-EXT-1 |

The bracket therefore settles at **19 TRUE, 31 REFUTED, 0 OPEN**, which is the
outcome the archive predicted for the case where the series inequalities hold.

Candidate 44 is worth a separate word. It stood against a single claimed row,
which made it the cheapest falsification target in the whole census: attaining
it would have refuted KLY09-EXT-1. That target is now closed from both sides.
The design route was closed by exact counting in
`docs/rank11_open_candidate_designs.md`, and KLY09-EXT-1 is now certified, so
no state of any kind attains it.

## What a certificate proves

The generation package states its convention plainly: a certificate `(H,z)`
proves `H . lambda >= z`. Certification requires two exact conditions. The
trace count, that the number of weights strictly below the affine hyperplane
equals the number of negative roots the hyperplane makes negative; and a
nonzero tangent determinant, computed by fraction-free Bareiss elimination on
exact integers, so a nonzero value is a proof and never a numerical artifact.

Validity is what refutes a candidate. A spectrum violating a valid inequality
is outside the polytope, full stop. Facetness, irredundancy and completeness of
the (3,11) system are all separate questions, and the screen records each of
them as not established. None of them is needed here.

## Independent replay

The repository's validation law is that a claim is not a certificate, and the
rank-7 checkpoint is explicit that the Ressayre tangent determinants sit
outside its independent verifier's trust boundary. Without a separate check the
rank-11 result would rest on the generator checking its own arithmetic.

`scripts/verify_level5_ressayre_3_11_standalone.py` closes that. It imports
nothing from `gpc_census` and uses the standard library only. From the recorded
`(h, z)` and evaluation point it rebuilds the 165 weights of `wedge^3 C^11`, the
on and below level sets, the negative root set and the tangent matrix, and
recomputes the determinant. Two of those steps are deliberately different
derivations, so agreement is evidence and not tautology: the fermionic sign is
obtained by substituting the created orbital into the ordered occupation tuple
and counting inversions, rather than by tracking annihilation and creation
positions the way the generator does; and the determinant is recomputed by
exact `Fraction` Gaussian elimination rather than by the generator's
fraction-free Bareiss routine. The recorded index sets are compared against the
recomputation but never used in it.

Sixty-four checks pass. The verifier also fails when it should. Perturbing the
recorded determinant, hyperplane, on-set, evaluation point or root set was
tried in turn and caught in every case; the determinant and hyperplane
mutations are asserted in the test suite, since a verifier that cannot fail
proves nothing.

## Why the screen is trusted

Calibrated at the two largest ranks where the fixed-`N=3` system is fully
published, by a test that is complete rather than sampled.

Every subset row of size 2 through 7 with right side 1 or 2 was screened, and
every certified row was then maximized by linear programming over the
published system in `src/gpc_census/data/constraints.json`. A maximum at or
below the right side proves the row is valid on the whole polytope.

| rank | rows screened | certified | invalid |
|---:|---:|---:|---:|
| 9 | 984 | 19 | 0 |
| 10 | 1,914 | 31 | 0 |

Fifty certified rows, none of them false. Four of the nineteen and five of the
thirty-one appear verbatim in the published table; the rest are valid but
redundant, which is exactly what a validity screen with no facetness claim
should produce.

## The negative control

An earlier control fed the screen four obviously false rows and watched all
four die. That control was weak, because all four were rejected by the trace
count before the determinant was ever formed. A screen that rejected
everything would have passed it.

The control that matters uses the same calibration runs. Among the rows
screened at ranks 9 and 10, those whose linear-programming maximum exceeds
their right side are demonstrably false. Split them by the route that rejected
them:

| rank | false rows | inadmissible | trace-rejected | reached the determinant | certified |
|---:|---:|---:|---:|---:|---:|
| 9 | 626 | 36 | 587 | 3 | 0 |
| 10 | 1,231 | 0 | 1,228 | 3 | 0 |

Two of those routes are cheap structural shortcuts. `inadmissible` means the
row is not a codimension-one exterior-weight hyperplane at all; `trace_rejected`
means the count of weights below it differs from the count of negative roots.
Neither forms the determinant.

So the honest reading is this. The screen rejected 1,857 false rows and
certified none of them, which is a real result about the screen as a whole. But
the trace count did nearly all of that work: only six false rows reached the
tangent determinant, and all six were rejected there on an exactly zero
determinant. A separate hunt for false rows that clear the trace count took all
145 published facets at ranks 9 and 10 and lowered each right side by one and
by two. Every such row is false, because a facet is tight. All 290 of them were
rejected before the determinant, 242 on the trace count and 48 as
inadmissible, so that hunt added nothing.

The determinant step is therefore the least stressed part of the pipeline in
this calibration, which is worth stating because the rank-11 result depends on
it. What the calibration does establish is that the trace count and
admissibility, taken together, already exclude essentially every false row of
this shape, and that the determinant has never contradicted them.

## The randomized sweep

Every row above has zero-one coefficients. That is the shape of the rank-11
rows, but it is a narrow slice of the row space, and a screen could in
principle be sound on it and wrong elsewhere. A seeded sweep at rank 9 screens
6,000 rows with coefficients drawn from `{-1,0,1,2}` and right sides from
`{0,1,2,3}`, audited the same way.

| rows | certified | invalid | false rows | reached the determinant | certified false |
|---:|---:|---:|---:|---:|---:|
| 6,000 | 5 | 0 | 3,791 | 1 | 0 |

Five certified, all valid. Nearly two thirds of the sweep is demonstrably
false and none of it was certified. This is a sample, not an enumeration, and
it does not make the calibration complete over general rows; what it rules out
is that the soundness at ranks 9 and 10 is an artifact of only ever having
been shown subset rows.

## Rank-11 consistency

The four certified rows were checked against all nineteen spectra the census
has already certified TRUE at (3,11). None violates any of the four. The
Slater point `(1,1,1,0,...)` and the uniform point `(3/11)^11` also satisfy all
four. The rows cut off the four OPEN candidates and nothing else that is known
to be attainable.

## What the closure changed in the repository

`scripts/settle_bracket_3_11.py` now carries LEVEL5-RESSAYRE as a fifth exact
method, applied last so it decides only what the four older methods left OPEN.
`docs/bracket_3_11_settlement.json` regenerates at 19 / 31 / 0, and the four
refutations each name the certified row they violate and by how much. Nothing
else moved: the other 46 verdicts keep their original certificates, which the
test suite asserts.

`scripts/partial_families.py` gains a RESSAYRE tier for these four rows **at
rank 11 only**. The rank-12 rows share their index sets but live in
`wedge^3 C^12`, a different representation, and a rank-11 certificate says
nothing about them, so they stay CLAIMED. Validity does not propagate upward:
restriction sends a valid (3,12) row down to (3,11), not the other way.

## What is not proved

The calibration is complete at ranks 9 and 10. Soundness of the screen at rank
11 follows from the mathematics of the Ressayre criterion, which is
rank-independent, and the calibration is evidence that the implementation
realises that criterion correctly. It is not an enumeration at rank 11, and the
distinction is worth keeping.

Nothing here enumerates the (3,11) system. The bracket settles because valid
inequalities refute the last four candidates, not because the facet list was
computed. Rank 11 remains out of reach for Stage 1: the measured closed-flat
growth of roughly thirtyfold per rank puts (3,11) at order `3e10` flats.

## Reproduction

```sh
uv run python scripts/level5_ressayre_3_11.py            # ~45 min, mostly the sweep
python scripts/verify_level5_ressayre_3_11_standalone.py # standard library only
uv run pytest -q tests/test_level5_ressayre_3_11.py
```
