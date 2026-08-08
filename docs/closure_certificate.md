# The closure certificate: what the sandwich does not check

**Status:** exact verification layer over the existing `P = O` proof. Nine
systems, 799 vertices, zero mismatches. One implementation subtlety found and
recorded.

**Artifact:** `results/data/closure_certificate.json`, written by
`scripts/closure_certificate.py`.

## What already existed, and what this adds

The two-sided architecture has a completeness half and a verification half, and
the repository already had the completeness half.
`docs/census_inner_sandwich.md` proves `P = O` at all nine census systems by
exactly the sandwich argument: every stored row replays as an exact Ressayre
certificate, giving `P` inside `O`; every vertex of `O` is the spectrum of a
certified state, giving `O` inside `P`. Candidate exhaustiveness never enters.

What was missing is the layer that catches an implementation slip rather than a
mathematical gap: certified dimension, facetness as distinct from validity,
vertex rank, slack rank, and half-filling parity. Those are built here.

This layer does **not** re-prove `P = O` and does not make the facet counts
theorem-level on its own. It makes them hard to get wrong by accident.

## Results

| System | Certified dim | Hull equations | Valid rows | Codim-1 rows | Facets | Vertices | Slack rank |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `(3,6)` | 3 | 3 | 7 | 6 | **4** | 4 | 4 |
| `(3,7)` | 6 | 1 | 11 | 10 | 10 | 10 | 7 |
| `(3,8)` | 7 | 1 | 39 | 39 | 39 | 38 | 8 |
| `(4,8)` | 7 | 1 | 23 | 23 | 23 | 22 | 8 |
| `(3,9)` | 8 | 1 | 61 | 61 | 61 | 58 | 9 |
| `(4,9)` | 8 | 1 | 69 | 69 | 69 | 103 | 9 |
| `(3,10)` | 9 | 1 | 103 | 103 | 103 | 113 | 10 |
| `(4,10)` | 9 | 1 | 135 | 135 | 135 | 159 | 10 |
| `(5,10)` | 9 | 1 | 171 | 171 | 171 | 292 | 10 |

Every column reproduces `docs/polytope_invariants.json` exactly, including the
minimum and maximum facet counts per vertex, which are computed here
independently. All 799 census vertices pass the rank test. Slack rank equals
`m + 1` at every system, and the slack zero pattern reproduces the incidence
relation at every entry.

## The subtlety: a facet is a face, not a row

The first implementation counted rows that pass the codimension-one test. It
reported **6 facets at `(3,6)`**, which is impossible: a 3-dimensional polytope
with 4 vertices is a tetrahedron and has 4.

The cause is exactly the canonicalization requirement that a facet ledger has to
state and that is easy to leave implicit. At `(3,6)` the affine hull carries
three equations rather than one, the Borland-Dennis relations
`lambda_1 + lambda_6 = lambda_2 + lambda_5 = lambda_3 + lambda_4 = 1`. Two rows
that differ by a relation of the hull cut out the same face:

```text
lambda_1 >= lambda_2   and   lambda_5 >= lambda_6   both give tight set {0,2,3}
lambda_2 >= lambda_3   and   lambda_4 >= lambda_5   both give tight set {0,1,3}
```

Six rows, four faces. Canonicalizing by the tight vertex set gives 4 and
restores agreement.

This is worth recording for two reasons. It is the only system where the
collapse happens, so a check that omitted canonicalization would pass eight of
nine systems and fail silently on the one with a degenerate hull. And it is a
concrete instance of why the affine hull must be certified rather than assumed:
had `m = d - 1` been hardcoded, `(3,6)` would have been wrong from the start,
since its trace slice has dimension 5 and the polytope has dimension 3.

Validity is also not facetness, independently of canonicalization: `(3,7)` has
11 valid rows and 10 facets, and `(3,6)` has 7 valid rows and 6 that are
codimension one. Redundant valid rows exist and are removed by the face test.

## Particle-hole parity applies at three systems, not nine

`Theta(lambda)_i = 1 - lambda_{d+1-i}` maps `(N,d)` to `(d-N,d)`. It is an
involution of a single system only at half filling, `2N = d`, which among the
census systems means `(3,6)`, `(4,8)` and `(5,10)`. Elsewhere it relates two
different systems and yields no within-system parity relation.

Where it does apply, the vertex set is closed under it and the parity relation
holds exactly:

| System | Vertices | `Theta`-fixed | Paired orbits | `V = V_fix + 2 V_pair` |
| --- | ---: | ---: | ---: | --- |
| `(3,6)` | 4 | 4 | 0 | PASS |
| `(4,8)` | 22 | 6 | 8 | PASS |
| `(5,10)` | 292 | 8 | 142 | PASS |

At `(3,6)` every vertex is self-dual, which is why the system carries no parity
information at all despite being half filled.

## Polar duality is not built, deliberately

`f_{m-1}(P) = f_0(P^o)` and `f_0(P) = f_{m-1}(P^o)` are equivalent to the two
statements this layer already establishes directly: no facet is redundant (the
codimension-one test) and every listed vertex is a genuine vertex (the rank
test). Constructing the polar and enumerating its vertices would re-derive those
from the same data rather than test them independently, so it would buy a
number that looks like a cross-check but is not one.

## What this does and does not license

- It does **not** establish completeness. That comes from the sandwich, which
  is already proved and is a separate document.
- It does establish that the shipped H-representation and V-representation are
  mutually consistent, exactly, at every census system: right dimension, right
  facets, right vertices, right incidences, right slack rank.
- It says nothing about `(20,40)` or any system outside the census. The
  reverse-search vertex enumeration and the coverage ledger over
  Bulois-Denis-Ressayre parameter space, which are the two pieces a rank-40
  closure would need, are not built here and are not implied by anything here.
- The BDR coverage ledger in particular remains unbuilt and unverifiable in
  this environment, per the INCOMPLETE label in
  `docs/bdr_constructor_comparison.md`.

## Verification

```sh
uv run scripts/closure_certificate.py
uv run pytest tests/test_closure_certificate.py
```
