# Gate 4: strict linear triangularity is refuted, on candidates and on facets

**Status:** conjecture test at ranks 7, 8 and 9, on two populations that are
kept apart. Strict triangularity fails on both.

**SCOPE CORRECTION (2026-08-08).** An earlier version of this document called
the wide population "true facets". It is not. It is every Hall-feasible row
with a nonzero Ressayre tangent determinant, which is a superset of the final
system: 48, 136 and 240 rows against published irredundant systems of 4, 31 and
52. Those rows have had no Farkas reduction, no facet witness and no BDR
birationality test. The published rows are now audited separately, and the
refutation is reported on both populations so it does not rest on the wider
scope.

**Artifact:** `results/data/levi_triangularity_gate.json`, written by
`scripts/levi_triangularity_gate.py`.

## Why this gate and not the other three

The Levi excitation programme proposes four gates. Gate 1, auditing the
excitation-defect bound on the rank-7 to rank-10 populations, cannot falsify
anything: `delta <= N` holds unconditionally, so the bound is vacuous at those
particle numbers, as recorded in `docs/root_budget_boundary.md`. Gates 2 and 3
are implementation rather than test.

Gate 4 is the one that can fail on data that exists. It asks whether the
tangent system of a facet solves by back substitution, which is what would
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

## Result on the published irredundant facets

This is the load-bearing table, because these rows are the actual facet system.

| System | Published rows | Fully triangular | Largest block |
| --- | ---: | ---: | ---: |
| `(3,7)` | 4 | **0** | 4 |
| `(3,8)` | 31 | **5** | 7 |
| `(3,9)` | 52 | **5** | 7 |

Every published `(3,7)` facet fails strict triangularity, with block
decompositions `[4]`, `[4]`, `[3,1]`, `[4]`. At `(3,8)` and `(3,9)` only 5 rows
of 31 and of 52 are strictly triangular.

## Result on determinant-nonzero candidates

The wider population, kept separate because it is not the facet system:

| System | Candidate rows | Fully triangular | Fraction | Largest block | Largest matrix order |
| --- | ---: | ---: | ---: | ---: | ---: |
| `(3,7)` | 48 | 14 | 29.2% | 9 | 14 |
| `(3,8)` | 136 | 20 | 14.7% | 14 | 22 |
| `(3,9)` | 240 | 23 | 9.6% | 15 | 36 |

The largest candidate block at rank 9 sits on

```text
tau = (-1, -1, -1, -1, -1, 2, 2, -4, 2)
matrix order 16, blocks [15, 1]
```

**This row is a candidate, not a facet.** Its sorted shape is not among the
eight sorted-`tau` shapes carried by the published `(3,9)` system, so it must
not be quoted as a facet witness. Among published facets the largest block is
7, not 15.

## What this does and does not refute

**Refuted:** the conjecture read strictly, that a facet becomes linear
triangular in the sense that back substitution alone solves the tangent system.
On the published irredundant systems that fails at 4 of 4, 26 of 31 and 47 of
52 rows, so the no-go does not depend on the wider candidate population.

The sharp statement is that Hall feasibility plus a nonzero determinant does
not imply scalar triangularity. A perfect matching says the support admits one;
a nonzero determinant says the system is generically nonsingular. Neither says
it decomposes into scalar equations.

**Not refuted:** the conjecture read as *block* triangular, with irreducible
blocks solved jointly. That remains open. On the published facets the largest
block is 7 at both `(3,8)` and `(3,9)`, which is markedly better behaved than
the candidate population's 14 and 15, so bounded block size is more plausible on
the facet system than the candidate scan suggested.

**Also not refuted, and not tested here:** BDR's own recursive
linear-triangular filter, which analyses the full fiber-equation graph and
eliminates variables iteratively. It is an optional partial validator in that
algorithm, not its birationality engine. Nothing here refutes BDR; it refutes a
hoped-for replacement of its hardest stage.

> **CORRECTION (2026-08-08).** An earlier version of this paragraph asserted
> that BDR's filter is a *strictly stronger* criterion than the DM `1 x 1`
> test. That was never tested, and it is false. The upstream implementation has
> since been fetched at its pinned revision and run on the published rows of
> `(3,7)`, `(4,7)`, `(3,8)` and `(4,8)`: 6 rows are BDR triangular with a
> single irreducible DM block, and 6 rows are fully DM triangular and rejected
> by BDR. The two criteria are incomparable, because they are computed on
> different objects. See `docs/bdr_linear_triangular_gate.md`. The refutation
> below is unaffected; it is a statement about the tangent matrix, which is the
> object the conjecture is about.

**Untouched:** the excitation theorem itself, which is proved, and the sparse
generic-isotropy argument, which is a separate and correct reduction. Gate 4
constrains only the birationality step.

## What the DM decomposition is still good for

The block multiset remains a canonical factorization of the tangent support
into independent irreducible components, which is useful as a preconditioner:
determinant factorization, modular evaluation, Schur complements, caching
repeated component types, and distributing components. What died is the claim
that every component is a scalar, not the decomposition itself.

## Consequence for the programme

The proof strategy sketched for the Levi Birationality Theorem is: pick a
minimal positive profile, use the inversion diagrams to find the root variables
reaching it, use Hall equality to make that block square, use nonvanishing to
make it invertible, remove it and induct. The measured block structure says the
peeled object is usually not a single variable but an irreducible block, so the
induction has to carry joint solves. That is materially harder than the sketch
suggests, and better known before the birationality stage is written than after.

A second conflation is worth naming, because it is the deeper one. A nonzero
tangent determinant is DIFFERENTIAL information: it certifies that the relevant
map is dominant and generically locally finite. BDR Condition A asks for
BIRATIONALITY, generic degree one, which is a global fiber-degree question. A
locally invertible map can still have several generic preimages, and no
statement about tangent support can exclude that. So even a positive
triangularity result would not have discharged Condition A on its own.

## Verification

```sh
uv run scripts/levi_triangularity_gate.py
uv run pytest tests/test_levi_triangularity_gate.py
```
