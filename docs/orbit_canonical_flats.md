# Weyl-orbit canonical enumeration

Status: **PROVED lossless, IMPLEMENTED, and VERIFIED at every level.** The
orbit-canonical enumerator reproduces the entire flat lattice at ranks 6, 7 and
8 without ever materialising it, and reaches **rank 9, which the materialising
enumerator could not finish**. (3,9) has **10,004,154 admissible hyperplanes**,
computed here for the first time, through 231 orbits.

Pre-registration: `docs/prereg_orbit_canonical_flats.md`, committed with the
measuring script at `799d0c2` before the rank-8 run. Certifying artifacts:
`src/gpc_census/generation/orbit_canonical.py`,
`scripts/orbit_canonical_flats.py`,
`results/data/orbit_canonical_flats.json`, `tests/test_orbit_canonical.py`.

## The correction that comes first

**(3,7) candidate exhaustiveness is CLOSED**, and `AGENT_PLAN.md` T2 step 1 was
stale. The methodology doc records that the closed-flat enumerator gives 5,341
hyperplanes in place of 1,623,160 candidate bases, 10,682 oriented rows,
screening resolving every one with none unresolved (10,133 trace failures, 499
exact determinant zeros, 1 structural Pauli row, 49 nonstructural Ressayre
rows, summing exactly to 10,682), and exact reduction returning the four
published facets. T2 step 1 said "candidate exhaustiveness is bypassed, not
proved", which is true of the gate-1 inner/outer closure and false of the flat
enumerator, so reading T2 alone made rank 7 look open. It is not. The plan is
corrected and the next gate is (3,8).

## Why orbit reduction is the one permitted pruning

The set of admissible hyperplanes is `S_d` STABLE: permuting modes permutes the
exterior weights, hence permutes closed flats, hence permutes hyperplanes. That
is verified here rather than assumed, `images_outside_set` being zero at every
rank. So storing one representative per orbit and expanding loses nothing, and
the enumeration is lossless by a theorem rather than by an empirical pattern,
which is what the standing rule requires.

**Stated carefully, because the stronger version is false.** This justifies
orbit reduction for ENUMERATING candidates only. It does NOT say that two
hyperplanes in one orbit are facets together or not at all: the GPC system is a
list of inequalities facing the DOMINANT CHAMBER, the chamber is a fundamental
domain and so is exactly what the Weyl group fails to preserve, and
`docs/screening_orbit_structure.md` measures the consequence, namely that
screening verdicts are not constant on orbits. An earlier version of this
paragraph asserted the stronger claim and was wrong.

## The measured payoff

| rank | hyperplanes | `S_d` orbits | reduction |
|---|---|---|---|
| 6 | 362 | 8 | 45x |
| 7 | 5,341 | 19 | 281x |
| 8 | 166,420 | 56 | 2,972x |
| 9 | **10,004,154** | **231** | **43,308x** |
| 10 | **889,205,792** | **1,337** | **665,076x** |

The raw count explodes and the orbit count barely moves. Orbits grow 8, 19, 56,
231, 1337, ratios 2.4, 2.9, 4.1, 5.8, against raw ratios 14.8, 31.2, 60.1, 88.9.
The orbit ratio is CLIMBING, so the reduction is not a constant-factor law and
must not be extrapolated as one.

## The algorithm

Hold only orbit representatives at each rank. Extend each representative by the
repository's own one-rank closure step, canonicalize every child, keep one per
canonical form. The full lattice is never built.

**Why it loses nothing.** Let `O'` be an orbit of rank `k+1` flats and pick
`F'` in it. `F'` contains some closed rank-`k` flat `F`, which lies in an orbit
with stored representative `R = sigma F`. Then `sigma F'` contains `R` and is
still in `O'`, so extending `R` reaches `O'`. Every orbit at level `k+1` is
therefore reached from some representative at level `k`.

**Canonicalization.** Brute force over `d!` is hopeless past rank 9, so modes
are split by an iterated incidence signature (colour refinement) and only
permutations respecting that ordered partition are tried. The signature is an
invariant of the flat, so refinement can never separate modes that a symmetry
identifies; it only narrows the search. Correctness does not depend on the
refinement being strong, only on it being invariant.

## The factorial tail, found by measuring and then removed

Refinement alone is not enough, and the failure mode was invisible in the orbit
counts because the answers were right; only the cost was wrong. Instrumenting
the canonicalizer to count leaves gave:

| system, level | flats | worst leaves | mean | `d!` |
|---|---|---|---|---|
| (3,8) top | 56 | 5,040 | 317.8 | 40,320 |
| (3,8) peak | 118 | **40,320** | 512.2 | 40,320 |
| (3,9) top | 231 | **362,880** | 3,724.5 | 362,880 |
| (3,9) peak | 494 | **362,880** | 1,750.6 | 362,880 |

The worst case was not a large cell, it was **exactly `d!`, at every rank**: on
the most symmetric flats refinement splits nothing at all and the search
degenerates to full brute force. At `d = 11` that is 39,916,800 permutations for
a single flat, which is a credible cause of the (3,11) run burning hours before
it died.

The fix is **individualization**, the standard nauty move. Refine; if the
partition is discrete the labelling is determined; otherwise choose a target
cell, and for each mode in it give that mode a colour strictly below its
cell-mates, re-refine, and recurse. The minimum over the leaves is the canonical
form.

**Why it is still canonical.** The target cell is chosen by an
isomorphism-invariant rule, the first cell of minimum size greater than one in
refined-colour order, so the entire search tree is an invariant of the flat and
minimising over its leaves is a well defined function of the isomorphism class.
Individualization only narrows the labellings considered, and every permutation
carrying the flat to the minimum respects refinement, so it still appears as a
leaf, which is what keeps the automorphism count exact.

**Measured after.**

| system, level | worst leaves before | after | mean before | after |
|---|---|---|---|---|
| (3,8) top | 5,040 | 5,040 | 317.8 | 309.5 |
| (3,8) peak | 40,320 | **1,440** | 512.2 | **73.5** |

and the (3,8) orbit enumeration falls from **401 seconds to 16.7 seconds**.

**The new cost is at its floor.** Every automorphism of a flat produces a
distinct leaf, so `leaves >= |Aut|` for a canonicalizer of this shape, and
across all 313 flats of the (3,8) orbit lattice the bound is met on 300 with a
worst overshoot of 3x.

## Isolated modes, which is why the floor itself had to be broken

That floor is real but it is not the right floor, and (3,11) proved it. A mode
lying in NO present weight is swapped with any other such mode by a
transposition that leaves the flat pointwise unchanged. Refinement can never
split them, so individualization walks all `k!` of them for nothing, and `|Aut|`
contains that `k!` as a factor.

At low ranks that is most of the modes. Rank 1 of (3,11) is a single weight: 8
of the 11 modes are isolated, `|Aut| = 3! * 8! = 241,920`, and the run produced
no output at all in 43 seconds of CPU because it was doing this 165 times over.

Isolated modes are therefore handled analytically instead of searched: their
cell is never individualized, their positions are assigned in a fixed order, and
the `k!` they contribute is multiplied back into the automorphism count at the
end. The image cannot depend on how they are ordered among themselves, because
they occur in no present weight. Measured on that same rank-1 flat at d = 11:

| | before | after |
|---|---|---|
| leaves | 241,920 | **6** |
| automorphisms reported | 241,920 | 241,920 |
| orbit size | 165 | 165 |

and the orbit size is `C(11,3) = 165`, which is exactly the number of
single-weight flats, so the analytic factor is checked against a count that is
known independently.

**Where the whole lattice ends up:**

| | (3,6) | (3,7) | (3,8) |
|---|---|---|---|
| flats | 31 | 94 | 313 |
| total leaves, refinement only | 964 | 5,180 | 32,774 |
| total leaves, with isolated modes handled | **900** | **4,672** | **28,898** |
| brute force would be | 22,320 | 473,760 | **12,620,160** |

The gain at rank 8 is modest, 437x rather than 385x, because these flats have
few isolated modes. That is exactly the point: the gain is concentrated at low
rank and high `d`, which is the regime the enumerator has to pass THROUGH to
reach a high rank, and it is where the run was actually stalling.

`tests/test_orbit_canonical.py` pins the leaves-to-automorphism ratio rather
than a timing, so the factorial tail cannot return unnoticed, and the orbit-size
sums (362 and 5,341) are what check the analytic `k!` correction.

## Verification, at every level and not just the top

The decisive check is that expanding each orbit by its size, computed from the
orbit-stabilizer theorem as `d! / |Aut|`, rebuilds the materialised lattice
level for level:

| rank | orbits per level | expands to |
|---|---|---|
| 6 | 1, 1, 3, 7, 12, 8 | 1, 20, 190, 735, 1005, 362 |
| 7 | 1, 1, 3, 10, 25, 36, 19 | 1, 35, 595, 4655, 15505, 19019, 5341 |
| 8 | 1, 1, 3, 11, 37, 87, 118, 56 | 1, 56, 1540, 21420, 147630, 467082, 565208, 166420 |
| 9 | 1, 1, 3, 12, 45, 146, 364, 494, 231 | 1, 84, 3486, 78274, 968121, 6383706, 20605599, 27392341, **10004154** |
| 10 | 1, 1, 3, 12, 49, 189, 691, 1840, 2705, 1337 | ..., 1994358765, **889205792** |

Every one of those expansions matches the materialising enumerator exactly at
ranks 6, 7 and 8. At rank 9 the materialising enumerator DIED at level 6, but
it had checkpointed levels 1 to 5 first, and the orbit expansion reproduces
those five stored counts exactly, 84, 3,486, 78,274, 968,121 and 6,383,706,
before continuing past the point of death.

**Cost.** Rank 8 takes 2.8 seconds against 126 seconds materialised. Rank 9
takes 401 seconds and never exceeds 494 orbits at a level, against a
materialised level 7 of 27,392,341 flats that does not fit in memory.

## Memory was the wrong shape too, and that is what killed rank 11

The enumerator held only orbit representatives BETWEEN levels but built every
child of every parent into one dict WITHIN a level, and only then collapsed
them onto canonical forms. So peak memory scaled with the child count, not the
orbit count, and the level being held was about to be discarded. The (3,11) run
died at level 9 holding 6.2 GB, which is that shape and not a compute limit.

`iter_one_rank_children` yields children parent by parent and the enumerator
collapses them as they arrive, so peak memory is now the orbit count of the next
level. Streaming re-canonicalizes a child once per parent that reaches it, so
there is a memo, and the memo is CAPPED: an uncapped one would reintroduce
exactly the per-child memory the change exists to remove, and dropping it costs
time only.

Rerun of (3,9) on the streaming path: the same orbit ladder
`[1, 1, 3, 12, 45, 146, 364, 494, 231]`, the same expansion ending in
**10,004,154**, in **308 seconds at a peak of 112 MB**. The representative kept
per class is the smallest `(mask, basis)` rather than the first seen, so it does
not depend on the order children happen to be produced in.

## Scorecard

| id | prediction | expectation | verdict |
|---|---|---|---|
| P-O-1 | rank-8 flat lattice reproduces the published vector | PASS | **PASS** |
| P-O-2 | the rank-8 `S_d` action is closed | PASS | **PASS**, 0 outside |
| P-O-3 | rank-8 orbit count in [30, 70] | PASS | **PASS**, 56 |
| P-O-4 | rank-8 orbit count in [40, 60] | PASS | **PASS**, 56 |
| P-O-5 | rank-8 reduction exceeds 2,000x | PASS | **PASS**, 2,972x |

Five for five, including the narrow externally proposed band.

## The trap, recorded because the wrong answer is plausible

Hyperplanes are stored with a sign convention, first nonzero entry positive, so
each unoriented hyperplane appears once. Permuting the coordinates of `h`
therefore leaves the stored set unless the image is re-normalized. Without the
renormalization the action sends 201 of 362 images outside the set at rank 6
and the orbit count comes out **15 instead of 8**, and 35 instead of 19 at rank
7. Those are not obviously wrong numbers, and they turned out not to be
meaningless either: 15 and 35 are exactly the counts of ORIENTED hyperplane
orbits, which is what the screening stage actually consumes, since each
hyperplane yields two oriented tau rows. So the un-normalized union-find was
answering a different and equally real question. The script reports
`images_outside_set` and a test asserts it is zero for the unoriented count,
and a second test pins the trap so a convention change cannot silently remove
the guard.

## The (3,8) replication

`166,420` and the level vector `[1, 56, 1540, 21420, 147630, 467082, 565208,
166420]` now have **three independent derivations**: the repo's own
`stage1_reference_series.json` from the gate-2 work, an external numpy
reimplementation over `GF(2^31 - 1)`, and this run. Against `C(56,7) = 2.3e8`
for brute force over candidate bases.

## What this does and does not license

**Does.** Exhaustive admissible-hyperplane enumeration at rank 9, which was out
of reach, with the count certified by an orbit expansion that agrees with the
materialising method wherever the latter got. The enumerator is lossless by a
theorem, not by a pattern.

**Does not.** Enumerating hyperplanes is not a complete GPC system. Each
representative still needs Ressayre screening and exact reduction.

**CORRECTION, 2026-08.** An earlier version of this paragraph said the
screening cost "scales with the ORBIT count only if screening is orbit
invariant, which is true for validity". That parenthetical was WRONG.
`docs/screening_orbit_structure.md` refutes it: at rank 7, 25 of the 35
oriented-tau orbits carry mixed verdicts, several containing both `certified`
and `trace_rejected`. The dominant chamber is a fundamental domain and is
therefore exactly what the Weyl group does not preserve, so a chamber-facing
condition cannot be an orbit invariant. Orbit reduction is lossless for
ENUMERATING candidates and inapplicable to DECIDING them. Nothing here produces
the (3,9) facets.

**Projection, labelled as such.** Orbit ratios 2.4, 2.9, 4.1 do not fix an
asymptotic from four points. If the ratio stays near 4 then rank 10 lands near
900 orbits and rank 11 near 3,600, against raw counts near 1e9 and 1e11. That
is a projection and the rank-10 run is the test of it.

## Reproducing

```sh
python scripts/orbit_canonical_flats.py --ranks 6 7 8
uv run pytest tests/test_orbit_canonical.py
python -c "from gpc_census.generation.orbit_canonical import \
enumerate_orbits_by_augmentation as e; print(e(3, 9)[0])"
```
