# Generator ladder gate 1: (3,7) is EXACT

Scores `docs/prereg_generator_gate1_3_7.md`, predictions frozen at commit
`9c55f1b` before the run. All seven PASS.

## Result

`(3,7)` moves from RELATIVE to EXACT on the claim ladder in
`docs/fixed_n3_generator_methodology.md`, by the independent inner/outer
closure route rather than by candidate enumeration.

| prediction | verdict |
|---|---|
| G1 the H-system has exactly 10 vertices | PASS |
| G2 the vertex set equals the reference set exactly | PASS |
| G3 every vertex carries an exact attaining state | PASS |
| G4 every H-row certified on independent replay | PASS |
| G5 the closure is exact, `P = H` | PASS |
| G6 the four rows are load bearing together | PASS |
| G7 each row is individually load bearing | PASS |

The argument: validity of every H-row gives `P` inside `H`; exact
attainability of every enumerated H-vertex gives `H` inside `P` by convexity;
hence `P = H`. Both inclusions are exact, and neither is a sampling argument.

## What was rechecked rather than read

Nothing is taken from the manifest on trust.

The four Ressayre certificates are replayed from their recorded `(h, z)` and
evaluation point, with the 35 weights, both level sets, the negative root set
and the tangent matrix rebuilt from scratch, and the determinant recomputed by
exact `Fraction` elimination rather than the generator's Bareiss routine. All
four reproduce their recorded values exactly:

| row | tangent determinant |
|---|---:|
| `l2+l3+l4+l5 <= 2` | 48,841 |
| second retained row | 1,928,920 |
| third retained row | 1,455,387 |
| fourth retained row | 4,817,604 |

Each attaining state's 1-RDM spectrum is recomputed from that record's own
support and amplitudes. The record's `integer_form` is deliberately not read,
since it is the field under test; using it would make G3 circular. All ten
match, and all ten states are tier `EXACT-CONSTR`.

## The controls, and why they earned their place

Without G6 and G7 a PASS could have come from the ordered slice alone.

The ordered trace-three slice by itself has **13** vertices. The four
Ressayre rows cut it to 10, so they are doing real work.

Individually, dropping each row changes the polytope:

| dropped row | vertices | vertex set changed |
|---:|---:|---|
| 0 | 11 | yes |
| 1 | 12 | yes |
| 2 | 13 | yes |
| 3 | 10 | yes |

Row 3 is the interesting one and the reason the control is stated as a set
comparison. Dropping it leaves the vertex **count** at 10 while changing
**which** ten they are: relaxing a polytope removes the vertices on the
dropped facet and introduces new ones, so a count test can report no change
where the polytope has in fact changed. A first version of this script tested
counts and wrongly scored row 3 as redundant. That was a defect in the
measuring code, corrected before scoring, and it is recorded here because the
corrected criterion is the one the prediction meant.

Two further defects in the measuring code were found and fixed before scoring,
both mine and neither a property of the data: inequality rows were built in
`<=` orientation where `exact_polytope` carries `>= `, and the root operator
for root `(i, j)` is `E_(j,i) = a_j^dag a_i`, which the first version applied
in the wrong direction and which drove every replayed determinant to zero.

## What this does NOT establish

- Nothing about `(3,8)` or higher. Gate 2 is a separate pre-registration.
- Not that candidate ENUMERATION at `(3,7)` is exhaustive. The closure route
  bypasses that question rather than answering it, so the manifest's
  `system_completeness` field is unchanged.
- Not facetness of any row beyond the existing facet witnesses. G7 shows each
  row changes the polytope, which is necessity, not a facet certificate.
- The ten `(3,7)` states are not ten independent observations. Several are
  transported or padded from `(3,6)`, and the artifact records that so this
  cannot later be cited as ten independent attainability measurements.

## Reproduction

```sh
uv run python scripts/generator_gate1_3_7.py
uv run pytest -q tests/test_generator_gate1_3_7.py
```
