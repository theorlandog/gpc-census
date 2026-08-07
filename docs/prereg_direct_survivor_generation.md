# Pre-registration: generating trace survivors directly

## What is being tested

`docs/screening_orbit_structure.md` closes with a named but unimplemented step:

> The remaining engineering step is to GENERATE the survivors directly instead
> of filtering them out of the full orbit [...] Without that step every row must
> still be visited, which is what sets the reachable rank.

This registers that step. `orbit_canonical.generate_trace_survivors(h, B)`
yields the arrangements of the multiset `h` with exactly `B` inversions, which
is exactly the Ressayre trace test with `B` the orbit-invariant count of weights
below.

## Status of the numbers before this run, stated plainly

This is a REPLICATION target, not a blind prediction, and it is labelled that
way rather than dressed up:

- 549 survivors from 10,682 oriented rows at rank 7, and 6,607 from 332,840 at
  rank 8, are already published in `results/data/screening_orbit_structure.json`.
- The per-orbit counts are already available in closed form from
  `trace_survivor_count`, the Gaussian multinomial, checked exactly on 35 of 35
  rank-7 orbits.
- The rank-7 agreement was observed BEFORE this file was written. The rank-8 run
  was launched before this file was committed.

What is genuinely open is whether the GENERATOR produces the right SET rather
than the right count, and whether it is output linear in practice.

## Predictions

| id | prediction | falsifier |
|---|---|---|
| P-G-1 | generated set equals filtered set on every rank-7 orbit | any orbit where the two sets differ |
| P-G-2 | rank-8 generated total is 6,607 | any other total |
| P-G-3 | every generated arrangement has inversion number exactly `B` | one that does not |
| P-G-4 | search nodes are bounded by a constant times the output | node count scaling with orbit size on some orbit |
| P-G-5 | generation is at least 20x faster than filtering at rank 8 | under 20x |

P-G-4 is the one that carries the argument. Correctness alone would be
worthless here: a generator that enumerated the orbit internally and filtered
would pass P-G-1 to P-G-3 and change nothing about the reachable rank.

## Why output linearity should hold

Prune a branch unless the residual budget lies in `[0, maxinv(rest)]` with
`maxinv = C(m,2) - sum C(m_i,2)`. That interval is exactly the set of attainable
inversion numbers, because the Gaussian multinomial has strictly positive
coefficients across its whole degree range with no internal zeros. So every
surviving branch completes to at least one output and no node is wasted.

If that claim about the coefficients were false, P-G-4 fails and the method
degrades to a smarter filter. It would still be correct.

## What a pass does not license

It does not reduce the number of orbits, which is set by enumeration. It does
not touch the determinant stage. It does not make screening an orbit invariant,
which is refuted and stays refuted. It only removes the row scan.
