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

The moment polytope is Weyl invariant and the ordered slice is a fundamental
domain for the `S_d` action, so two admissible hyperplanes in one orbit are
facets together or not at all. Enumerating representatives DISCARDS NOTHING.
Every other pruning proposed for this search has been an empirical pattern,
which the standing rule forbids until a theorem makes it lossless. This one has
the theorem, which is why it is implemented and the others are not.

## The measured payoff

| rank | hyperplanes | `S_d` orbits | reduction |
|---|---|---|---|
| 6 | 362 | 8 | 45x |
| 7 | 5,341 | 19 | 281x |
| 8 | 166,420 | 56 | 2,972x |
| 9 | **10,004,154** | **231** | **43,308x** |

The raw count explodes and the orbit count barely moves. Orbits grow 8, 19, 56,
231, ratios 2.4, 2.9, 4.1, against raw ratios 14.8, 31.2, 60.1.

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
are first split by an iterated incidence signature (colour refinement) and only
permutations respecting that ordered partition are tried. The signature is an
invariant of the flat, so refinement can never separate modes that a symmetry
identifies; it only narrows the search. Correctness does not depend on the
refinement being strong, only on it being invariant.

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

Every one of those expansions matches the materialising enumerator exactly at
ranks 6, 7 and 8. At rank 9 the materialising enumerator DIED at level 6, but
it had checkpointed levels 1 to 5 first, and the orbit expansion reproduces
those five stored counts exactly, 84, 3,486, 78,274, 968,121 and 6,383,706,
before continuing past the point of death.

**Cost.** Rank 8 takes 2.8 seconds against 126 seconds materialised. Rank 9
takes 401 seconds and never exceeds 494 orbits at a level, against a
materialised level 7 of 27,392,341 flats that does not fit in memory.

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
7. Those are not obviously wrong numbers: they are valid orbit counts of a
valid group action on a set that simply is not the hyperplane set. The script
reports `images_outside_set` and a test asserts it is zero, and a second test
pins the trap itself so a convention change cannot silently remove the guard.

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
representative still needs Ressayre screening and exact reduction, and the
screening cost scales with the ORBIT count only if screening is orbit
invariant, which is true for validity but has not been re-verified in this
repository at rank 9. Nothing here yet produces the (3,9) facets.

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
