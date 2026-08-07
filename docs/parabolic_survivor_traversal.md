# Parabolic survivor traversal

Status: **REFUTED as a speedup, with the mathematics confirmed.** The
reformulation is correct, the incremental update is valid and produces
bit-identical partitions, and it saves 3 to 7 percent of the arithmetic while
running slower in wall time. Two structural facts explain why, and both are
worth keeping.

Certifying artifacts: `scripts/parabolic_survivor_traversal.py`,
`tests/test_parabolic_survivor_traversal.py`.

## The reformulation, which is correct

Let the primitive normal have multiplicities `h = (a_1^{m_1}, ..., a_r^{m_r})`.
The distinct arrangements of `h` are the parabolic quotient
`S_d / (S_{m_1} x ... x S_{m_r})`, the trace condition's inversion number is the
Coxeter length on minimal coset representatives, and the Gaussian multinomial
already used by `trace_survivor_count` is that quotient's Poincare polynomial.

Verified directly rather than assumed, and one half needs a convention stated:

| claim | verdict |
|---|---|
| Young subgroup order is `prod m_i!` | confirmed |
| coset count is `d! / prod m_i!` | confirmed |
| Poincare polynomial equals the Gaussian multinomial | confirmed, convention free |
| minimal coset length equals the arrangement's inversion number | **only with `h` sorted INCREASINGLY** |

With `h` sorted decreasingly, neither the minimum nor the maximum coset length
equals the arrangement's inversion number. The decreasing arrangement is the one
with MAXIMAL inversions, so it is the longest representative, not the shortest.

## The first obstruction: a fixed-length Gray code cannot exist

The proposal was to traverse length-`B` representatives so that consecutive
candidates differ by one adjacent transposition. That is impossible. An adjacent
transposition changes Coxeter length by exactly one, so two distinct length-`B`
representatives are never adjacent.

Measured, not merely argued: across four multisets, **562 adjacent swaps between
arrangements of the same multiset, and zero preserve the inversion number.**

Any traversal must therefore detour through length `B +- 1`, so consecutive
emitted survivors are at least two swaps apart.

## The increment that does exist, and does work

The DFS that already generates the survivors supplies a valid increment without
any Gray code. Consecutive survivors share a prefix, assigning position `k`
determines exactly the `C(k,2)` triples `{i,j,k}` with `i<j<k`, and backtracking
retracts them. Carrying the on/below partition as a stack is then exact.

It is exact: the partitions it produces are identical to a full rebuild on every
survivor at ranks 7 and 8.

## The second obstruction: survivors do not share prefixes

| | (3,7) | (3,8) |
|---|---|---|
| survivors | 549 | 6,607 |
| partitions identical to rebuild | yes | yes |
| triple sums, full rebuild | 19,215 | 369,992 |
| triple sums, incremental | 18,674 | 345,040 |
| **ratio** | **1.03** | **1.07** |
| predicted ratio `d/3` | 2.33 | 2.67 |
| wall time | 0.01 vs 0.01 s | 0.11 vs **0.17 s** |

The prediction assumed prefix sharing. There is almost none: at rank 8 the
incremental route spends **52.2 triple sums per survivor against 56 for a full
rebuild**, so the DFS tree is close to a star and survivors branch apart near
the root. The bookkeeping then costs more than the arithmetic it saves, which is
why the incremental version is slower despite doing less work.

## Why the stage was not worth optimizing anyway

Survivor generation and partition construction together are **0.11 seconds of a
12.4 second rank-8 run**, under one percent. The stage that dominated screening
was the futile evaluation search, and removing that (see
`docs/determinant_matching_audit.md`) took rank 8 from 66 to 12.4 seconds. After
that fix there was no headroom here to recover.

## What survives

The parabolic identification itself, with its convention pinned. It is the right
language for the survivor set, `trace_survivor_count` is already its Poincare
polynomial, and the correctness of the incremental partition is available if a
future stage needs the partition maintained rather than rebuilt. What is refuted
is the claim that traversing this way is faster.

## Reproducing

```sh
python scripts/parabolic_survivor_traversal.py 7
python scripts/parabolic_survivor_traversal.py 8
uv run pytest tests/test_parabolic_survivor_traversal.py
```
