# Generating trace survivors instead of filtering for them

Status: **IMPLEMENTED and VERIFIED set against set.** The row scan, which was
the thing setting the reachable rank, is gone. The Ressayre trace survivors are
now produced directly from one representative per orbit, in time proportional to
the number of survivors rather than to the size of the orbit.

Pre-registration: `docs/prereg_direct_survivor_generation.md`. Certifying
artifacts: `orbit_canonical.generate_trace_survivors`,
`tests/test_orbit_canonical.py`.

## The bottleneck this removes

`docs/screening_orbit_structure.md` established the shape of the problem and
then named what was missing:

> It does not reduce the number of CANDIDATES that must be visited: every
> oriented row still has to be touched, which is what sets the reachable rank.

That was the honest statement of the limit. Every row was visited and almost all
were discarded, at a rate that gets worse with rank: survivor fractions 8.8, 5.1
and 2.0 percent at ranks 6, 7 and 8. At rank 9 the oriented rows number about
20 million and at rank 10 about 1.8 billion, essentially all of them destined to
fail one integer comparison.

## Why generation is possible at all

The trace test is `inv(sigma h) == B` where `B`, the count of exterior weights
strictly below the hyperplane, is an `S_d` orbit INVARIANT and `sigma h` ranges
over the arrangements of the multiset `h`. So the survivors of an orbit are
exactly the arrangements of a multiset with a PRESCRIBED inversion number.

That is a classical object. It is counted by the Gaussian multinomial
`[d]_q! / prod [m_i]_q!`, which `trace_survivor_count` already read off in
closed form. Counting it is not the same as producing it, and producing it is
what removes the scan.

## The algorithm, and why it wastes nothing

Build the arrangement left to right. Placing value `v` next costs exactly the
number of remaining elements strictly less than `v`, since each of those will
land to its right and form an inversion. Recurse on the residual budget.

Prune a branch unless the residual budget lies in `[0, maxinv(rest)]` with
`maxinv = C(m,2) - sum C(m_i,2)`.

**That interval is exactly the set of attainable inversion numbers**, because
the Gaussian multinomial has strictly positive coefficients across its whole
degree range, with no internal zeros. So every surviving branch completes to at
least one output: no node is wasted and the cost is proportional to the OUTPUT,
not to the orbit.

This is the property that matters, and it is the one that is easy to fake. A
generator that enumerated the orbit internally and filtered would produce
exactly the same answers and change nothing at all, so output linearity is
measured directly by instrumenting the prune, not asserted.

## Measured

| system | oriented rows | orbits | survivors, filtered | survivors, generated | filter | generate | speedup |
|---|---|---|---|---|---|---|---|
| (3,7) | 10,682 | 35 | 549 | **549** | 1.13 s | **0.0085 s** | **133x** |
| (3,8) | 332,840 | 109 | 6,607 | **6,607** | 53.44 s | **0.0743 s** | **719x** |

Both totals reproduce `results/data/screening_orbit_structure.json` exactly, so
P-G-1 and P-G-2 pass and P-G-5 passes by a factor of 36.

Ten of the 35 rank-7 orbits and **45 of the 109 rank-8 orbits** are killed
outright by `orbit_is_trace_dead` before a single arrangement is built, and
**326,233 of the 332,840 rank-8 rows are never visited at all**.

The speedup RISES with rank, 133x to 719x, for the same reason the survivor
fraction falls: filtering pays for the whole orbit and generation pays for the
survivors, so the gap is the reciprocal of a fraction that is heading toward
zero. That is a change of shape, not a constant factor.

## Verification

The count agreeing is not enough, so three separate things are checked:

1. **Set against set.** On every one of the 35 rank-7 orbits, the generated set
   equals the set obtained by filtering the orbit's actual rows. A generator
   producing the right NUMBER of wrong arrangements would silently replace the
   candidate list, and a count check would not notice.
2. **Against the closed form, off the lattice.** Sixty random multisets, every
   inversion target, checked against `trace_survivor_count` and against a direct
   inversion count of each arrangement produced. The two routes are independent:
   one is a polynomial division, the other a search.
3. **Output linearity, measured.** The prune is instrumented and the node count
   is asserted to be bounded by a constant times the output, on a multiset whose
   orbit has 5,040 members.

## What this does and does not license

**Does.** Screening no longer needs the oriented row list. Given one
representative per orbit, which the orbit enumerator produces, the candidates
that reach the determinant stage are generated directly. Combined with the
determinant work already in place (Hall matching decides 97.1 percent of rank-8
survivors structurally in 0.08 seconds), the per-rank cost is now driven by the
ORBIT count and the SURVIVOR count, both of which are small, rather than by the
hyperplane count, which is not.

**Does not.** It does not reduce the number of orbits: enumeration still sets
that, and the orbit ratio is climbing, not constant. It does not touch the
determinant stage. It does not make screening an `S_d` invariant, which is
refuted in `docs/screening_orbit_structure.md` and stays refuted. And a rank is
not closed by producing candidates: certification and exact reduction still have
to run.

## Reproducing

```sh
uv run pytest tests/test_orbit_canonical.py -k survivor
```
