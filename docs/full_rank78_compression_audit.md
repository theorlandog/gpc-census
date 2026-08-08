# Full rank-7/rank-8 compression audit

**Status:** exact full-population audit.  
**Source commit:** `058e0a34834b2ff3eacfab258d8030264db7d2b3`.

## Results

| Test | `(3,7)` | `(3,8)` | Outcome |
|---|---:|---:|---|
| Trace survivors | 549 | 6,607 | exact population gates |
| Hall zeros | 474 | 6,414 | exact |
| Hall-feasible | 75 | 193 | exact |
| Survivors with 2+ arithmetic holes | 0 | 0 | none |
| One-hole survivors | 0 | 1,019 | interval/one-hole grammar covers all |
| One-hole Hall-feasible | 0 | 7 | exact |
| One-hole determinant-nonzero | 0 | 0 | all seven cancel exactly |
| Nonstructural singleton collisions | 0/547 | 0/6,605 | injective |
| Nonstructural pair collisions | 0/547 | 0/6,605 | injective but unnecessary |
| Maximum optimal fusion width | 3 | 4 | small |
| Maximum balanced fusion width | 3 | 4 | balanced is already optimal |

## Arithmetic-level result

The full survivor populations have hole distributions

```text
(3,7): 0 holes = 549
(3,8): 0 holes = 5,588; 1 hole = 1,019; 2+ holes = 0
```

Only seven one-hole rows survive Hall at rank 8. All have normalized level set
`[0,1,2,4]`, with `3` missing. Exact signed sparse determinant expansion gives
the zero polynomial in all seven cases. Therefore every determinant-nonzero
rank-7/rank-8 candidate has an exact arithmetic-progression level set.

The optional hole must remain in a future grammar because the earlier bounded
rank-9 trial produced exact one-hole nonvanishing examples.

## Incidence tomography

For the positive-side determinant family, the labeled singleton-incidence
vector distinguishes every nonstructural trace survivor:

```text
(3,7): 547 survivors -> 547 singleton signatures
(3,8): 6,605 survivors -> 6,605 singleton signatures
```

Pair incidence is also injective but resolves no additional collision.

## Projection/fusion width

Target-compatible partial states `(particle count, weight sum)` remain tiny:

```text
(3,7): widths 1/2/3 occur 11/73/465 times
(3,8): widths 1/2/3/4 occur 17/222/3,964/2,404 times
```

Balanced trees attain the same width distribution as unrestricted optimal
trees. This is a strong evaluation-compression result, but not a facetness
discriminator: 2,393 rank-8 Hall zeros also have width four.

## Supported constructor architecture

```text
singleton-incidence key
  -> arithmetic interval / optional one-hole mechanism
  -> target-aware fusion-tree dynamic program
  -> Hall obstruction
  -> exact determinant or local semi-invariant
  -> global redundancy reduction
```

The next frontier test is the same streaming audit at rank 9.
