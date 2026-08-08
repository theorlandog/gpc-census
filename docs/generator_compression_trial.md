# Early go/no-go trials for compressed GPC generation

**Status:** post-hoc exploratory trial.

**Source repository commit:** `058e0a34834b2ff3eacfab258d8030264db7d2b3`.

## Executive result

| Idea | Verdict | What the trial established |
| --- | --- | --- |
| Target-aware projection/fusion trees | **GO** | All 83 known system-specific mechanisms have exact optimal block-tree width at most **5**; recursively balanced trees need at most **6** states. |
| Exact arithmetic progressions only | **NO-GO** | One-hole normals pass Hall and three carry exact modular nonzero determinant certificates. |
| Arithmetic progression with bounded holes | **TENTATIVE GO** | In the bounded search, all **135** placements with two or three holes fail Hall. |
| Singleton-incidence tomography | **GO FOR A THEOREM ATTEMPT** | It distinguishes all 83 known mechanisms and every nonstructural placement in the bounded search, with zero cross-shape collisions. |
| Complete generator replacement | **NOT YET** | The bounded normal window and zero-grade calculation do not cover the full candidate universe, the entire tangent determinant problem, or global redundancy. |

## Trial A: exact projection width

For a subset `S` of Levi blocks, retain each partial state `(k,s)` that is achievable in `S` and can be completed through the complement to the global target `(N,0)`. A fusion tree's width is the largest number of retained states at any node.

Across all 83 mechanisms:

```text
optimal tree width distribution:
{"1": 12, "2": 21, "3": 26, "4": 19, "5": 5}

balanced tree width distribution:
{"1": 12, "2": 21, "3": 26, "4": 18, "5": 4, "6": 2}

maximum optimal width: 5
maximum balanced width: 6
maximum optimal one-block chain width: 6
maximum naive sorted-chain width: 7
```

The unpruned full `(k,s)` state set reaches 105; target-aware forward/backward pruning is therefore doing real work.

### Synthetic high-rank stress test

These are constructive upper bounds for one consecutive-level normal family, not claims that the synthetic normals are GPC facets.

| System | Balanced contiguous width | Symmetric-extremes width |
| --- | ---: | ---: |
| (3,10) | 7 | 7 |
| (3,20) | 17 | 17 |
| (3,40) | 37 | 37 |
| (4,10) | 9 | 9 |
| (4,20) | 29 | 27 |
| (4,40) | 69 | 57 |
| (5,10) | 9 | 11 |
| (5,20) | 34 | 43 |
| (5,40) | 94 | 103 |
| (10,20) | 60 | 118 |
| (20,40) | 449 | 1267 |

The synthetic `(20,40)` target needs 449 states under the balanced split rather than materializing all 20-subsets of 40 orbitals. This supports commodity-scale verification of a specified high-rank normal, not generation of every normal.

## Trial B: bounded candidate arithmetic grammar

The search is exhaustive among primitive sorted `N=3` normals with entries in `[-5,5]` at ranks 7, 8, and 9, subject to:

1. zero-level exterior weights of affine rank `d-2`;
2. explicit trace-compatible multiset placements;
3. Hall matching for every trace placement.

| System | Candidate shapes | Trace placements | Hall-feasible | Trace holes | Hall-feasible holes |
| --- | ---: | ---: | ---: | --- | --- |
| (3,7) | 23 | 399 | 69 | 0 holes: 399 | 0 holes: 69 |
| (3,8) | 45 | 3677 | 161 | 0 holes: 3382, 1 holes: 295 | 0 holes: 154, 1 holes: 7 |
| (3,9) | 98 | 61106 | 307 | 0 holes: 40743, 1 holes: 20228, 2 holes: 22, 3 holes: 113 | 0 holes: 276, 1 holes: 31, 2 holes: 0, 3 holes: 0 |

```text
trace placements: 65182
Hall-feasible placements: 537
multi-hole trace placements: 135
multi-hole Hall-feasible placements: 0
one-hole Hall-feasible placements: 38
```

Exact arithmetic progressions alone are therefore too restrictive. The surviving stage-specific hypothesis is:

> Hall-feasible `N=3` candidates have arithmetic-progression levels with at most one missing point.

That remains a bounded-window conjecture.

## Exact one-hole nonvanishing witnesses

The tangent-matrix convention was checked first against all 49 nonstructural pilot facets. All 49 satisfy the trace equality, have perfect matchings, and yield deterministic modular nonzero determinant evaluations.

### Certificate 1

```text
tau: (-5, -3, -2, -1, 0, 1, 2, 3, 4)
prime: 1000003
matrix size: 36
variables: 8
determinant residue: 156335
```

The distinct levels have exactly one missing arithmetic-progression point. A nonzero determinant residue at this integer evaluation proves that the symbolic determinant is nonzero over the integers.
### Certificate 2

```text
tau: (-4, 4, -3, -2, -2, -1, 0, 3, 1)
prime: 1000003
matrix size: 27
variables: 8
determinant residue: 154405
```

The distinct levels have exactly one missing arithmetic-progression point. A nonzero determinant residue at this integer evaluation proves that the symbolic determinant is nonzero over the integers.
### Certificate 3

```text
tau: (0, 2, -2, -1, -1, -1, -1, 0, 0)
prime: 1000003
matrix size: 15
variables: 10
determinant residue: 336353
```

The distinct levels have exactly one missing arithmetic-progression point. A nonzero determinant residue at this integer evaluation proves that the symbolic determinant is nonzero over the integers.

These witnesses prove that a one-hole grammar is needed before global redundancy elimination.

## Trial C: singleton-incidence tomography

For the positive Ressayre side, record only the number of below determinants incident to each orbital.

Results:

- sorted singleton incidence distinguishes all **83/83** known system-specific mechanism shapes;
- labeled singleton incidence distinguishes every nonstructural trace placement in the bounded search;
- the only collisions are structural zero-below families;
- nonstructural cross-shape collisions: **0**.

This is not yet a theorem. It is strong evidence for using the incidence vector as the first canonical hash and for attempting a fixed-weight Chow-parameter theorem.

## Generator consequence

```text
incidence shadow
  -> arithmetic one-hole mechanism
  -> target-aware fusion tree
  -> Hall / modular determinant
  -> global reduction
```

The next decisive test is the same classification over the repository's complete rank-7 and rank-8 trace-survivor sets. The critical falsifiers are:

1. a Hall-feasible normal with two or more arithmetic holes;
2. two nonstructural trace placements with the same singleton incidence vector but different normals;
3. rapidly growing fusion-tree width in the rejected-candidate population.

## Reproduction

```sh
python scripts/generator_compression_trial.py
python scripts/verify_generator_compression_trial.py
```
