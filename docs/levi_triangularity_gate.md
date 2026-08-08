# Gate 4: strict linear triangularity is refuted on the known true facets

**Status:** conjecture test, run on every determinant-nonzero Hall-feasible row
at ranks 7, 8 and 9. Strict triangularity fails for the large majority, and the
failing fraction grows with rank.

**Artifact:** `results/data/levi_triangularity_gate.json`, written by
`scripts/levi_triangularity_gate.py`.

## Why this gate and not the other three

The Levi excitation programme proposes four gates. Gate 1, auditing the
excitation-defect bound on the rank-7 to rank-10 populations, cannot falsify
anything: `delta <= N` holds unconditionally, so the bound is vacuous at those
particle numbers, as recorded in `docs/root_budget_boundary.md`. Gates 2 and 3
are implementation rather than test.

Gate 4 is the one that can fail on data that exists. It asks whether the
tangent system of a true facet solves by back substitution, which is what would
remove the generic Groebner calculation from birationality.

## What is measured

For a determinant-nonzero Hall-feasible row the tangent matrix is square and
carries a perfect matching. Orient it by any perfect matching and take strongly
connected components of the induced digraph on rows: those are the irreducible
blocks of the Dulmage-Mendelsohn decomposition. The system is triangular in the
strict sense exactly when every block is `1 x 1`.

The block multiset is matching-independent, which is the Dulmage-Mendelsohn
theorem. Because the whole verdict rests on it, it is re-checked rather than
assumed: 40 rows per rank against four extra random matchings each, with **zero
disagreements**, plus a separate 191-row sweep at `(3,8)` before this gate was
written.

## Result

| System | True-facet rows | Fully triangular | Fraction | Largest block | Largest matrix order |
| --- | ---: | ---: | ---: | ---: | ---: |
| `(3,7)` | 48 | 14 | **29.2%** | 9 | 14 |
| `(3,8)` | 136 | 20 | **14.7%** | 14 | 22 |
| `(3,9)` | 240 | 23 | **9.6%** | 15 | 36 |

Strict triangularity is the exception, not the rule, and the exception is
getting rarer: 29 percent, then 15, then 10.

A rank-9 witness, where a 16 by 16 tangent matrix has a single irreducible
block of size 15 and one trivial block:

```text
tau = (-1, -1, -1, -1, -1, 2, 2, -4, 2)
matrix order 16, blocks [15, 1]
```

Nothing about this row is exotic. It has three distinct levels, it is
Hall-feasible, and its determinant is nonzero, so it is exactly the kind of row
the conjecture is about.

## What this does and does not refute

**Refuted:** the conjecture read strictly, that a true facet becomes linear
triangular in the sense that back substitution alone solves the tangent system.
That is false for 90 percent of rank-9 true facets.

**Not refuted:** the conjecture read as *block* triangular, with irreducible
blocks solved jointly. That remains open. But the premise it would need is
bounded block size, and the largest block grows across the three ranks, 9 then
14 then 15. Three points are not a trend, and the block does not track the
matrix order monotonically (the largest block at rank 9 occurs at order 16, not
at the largest order 36), so this is a caution rather than a second refutation.

**Untouched:** the excitation theorem itself, which is proved, and the sparse
generic-isotropy argument, which is a separate and correct reduction. Gate 4
constrains only the birationality step.

## Consequence for the programme

The proof strategy sketched for the Levi Birationality Theorem is: pick a
minimal positive profile, use the inversion diagrams to find the root variables
reaching it, use Hall equality to make that block square, use nonvanishing to
make it invertible, remove it and induct. The measured block structure says the
peeled object is usually not a single variable but a large irreducible block, so
the induction has to carry joint solves of blocks whose size is not yet known to
be bounded. That is a materially harder theorem than the sketch suggests, and it
is better to know before the birationality stage is written than after.

## Verification

```sh
uv run scripts/levi_triangularity_gate.py
uv run pytest tests/test_levi_triangularity_gate.py
```
