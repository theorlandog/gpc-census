# Profile-native candidate generation: the orbit stage without the lattice

**Status:** three proved reductions, one exact cross-algorithm calibration
against five independently known populations, one refutation, and one quantified
gap between a spread bound that is proved and a spread bound that can be run.

**Artifacts:** `src/gpc_census/generation/profiles.py`,
`scripts/profile_native_frontier.py`,
`results/data/generator_frontier_refinement.json`, `tests/test_profiles.py`,
and for the spread-cost table `scripts/profile_spread_cost.py`,
`results/data/profile_spread_cost.json`, `tests/test_profile_spread_cost.py`.

## What was blocking, and what this removes

`docs/orbit_canonical_flats.md` established that the admissible hyperplane set
is `S_d` stable and reached rank 10 by holding only orbit representatives while
walking the closed-flat lattice level by level. It also named its own limit: the
orbit COUNT is small and barely moves, but the lattice the walk passes through
does not stay small. Rank 9 costs 401 seconds, rank 10 needs a checkpointed run
and peaks at 2,705 orbits on an intermediate level, and rank 11 has not been
finished by that route.

The observation this document turns into an algorithm is one the same document
already made and then used only for canonicalization:

> A codimension-one closed flat IS a hyperplane, so it is fully described by its
> exact primitive normal `(h, z)`. Two of them lie in the same `S_d` orbit
> exactly when one `h` is a permutation of the other with matching `z`.

So an orbit is a MULTISET OF INTEGERS. The lattice is a road to an object that
can be written down directly. This module writes it down directly.

## Conventions

Work with the primitive homogeneous normal `tau` rather than with `(h, z)`. On
the trace slice, `tau . lambda <= 0` and `h . x = z` are related by
`tau_i = N h_i - z`, and the enumerator inverts that with
`h = d tau - (sum tau) 1`, `z = -N (sum tau)`, which satisfies `sum(h) = 0` and
`h . omega - z = d (tau . omega)`, so the two descriptions have identical
incidence sets by construction. `incidence_agrees` checks this weight by weight.

An admissible hyperplane is one whose exterior weights have affine rank `d - 1`,
which is the Ressayre admissibility condition already implemented in
`generation.candidate_pipeline`. In `tau` coordinates the weights on it are

```text
Z(tau) = { I in binom([d], N) : sum_{i in I} tau_i = 0 }
```

and admissibility is `rank_Q { chi_I : I in Z(tau) } = d - 1`.

## Theorem A: admissibility is a statement about the value profile

Let `tau` have distinct values `v_1 > ... > v_r` with multiplicities
`m_1, ..., m_r`. Let

```text
P(m) = { c in Z_{>=0}^r : sum_a c_a = N, c_a <= m_a }
S    = { c in P(m) : sum_a c_a v_a = 0 }
```

**Theorem A.** `rank_Q { chi_I : I in Z(tau) } = d - 1` if and only if

1. `rank_Q(S) = r - 1`, and
2. every class `a` with `m_a >= 2` is MOBILE: some `c` in `S` has
   `1 <= c_a <= m_a - 1`.

*Proof.* Let `L` be the span of `{ chi_I : I in Z(tau) }` and `V` the span of the
differences `chi_I - chi_J`. Every `chi_I` has coordinate sum `N`, so
`dim L = dim V + 1` and the claim is `dim V = d - 2`.

Let `pi : R^d -> R^r` sum coordinates within each value class. Then
`dim V = dim(V n ker pi) + dim pi(V)`.

For the quotient part, `pi(chi_I) = c(I)` is exactly the occupancy pattern of
`I`, and `I` lies in `Z(tau)` precisely when `c(I)` lies in `S`, so
`pi(V) = span{ c - c' : c, c' in S }`, whose dimension is the affine dimension of
`S`. Every element of `S` has coordinate sum `N`, so the origin is off the affine
hull and `dim pi(V) = rank_Q(S) - 1`. That is `r - 2` exactly when condition 1
holds, and it can never exceed `r - 2`, because `v` is a nonzero element of the
kernel of `S`.

For the kernel part, write `ker pi n {sum x = 0} = sum_a W_a`, where `W_a` is the
sum-zero subspace of class `a`, of dimension `m_a - 1`, so the total is `d - r`.
Let `Pi = prod_a S_{m_a}` be the class-preserving subgroup. It acts on `Z(tau)`,
hence on `V`, and the `W_a` are pairwise non-isomorphic irreducible
`Pi`-modules, so `V n ker pi` is the sum of the `W_a` it contains and nothing
else. Now fix a class `a` with `m_a >= 2`. If some `c` in `S` has
`1 <= c_a <= m_a - 1`, then a weight `I` realizing `c` meets class `a` without
exhausting it, so swapping one occupied mode of class `a` for an unoccupied one
gives `I'` in `Z(tau)` with `chi_I - chi_{I'} = e_i - e_{i'}`, which has nonzero
component in `W_a`, hence `W_a` is contained in `V`. Conversely if no such `c`
exists, then every `I` in `Z(tau)` meets class `a` in the empty set or in all of
it, so the class-`a` block of every `chi_I` is `0` or all ones, both of which are
`W_a`-orthogonal, and `V` has zero component in `W_a`. Classes with `m_a = 1`
contribute `W_a = 0`. So `dim(V n ker pi) = d - r` exactly when condition 2
holds, and never exceeds it.

Both parts are bounded above by their target, so `dim V = d - 2` if and only if
both are attained. QED

**Condition 2 is load bearing.** At `d = 6`, `tau = (1,1,0,0,0,-2)` has
`S = {(0,3,0), (2,0,1)}` of rank `2 = r - 1`, so the rank half passes. Class 1
has `m_1 = 2` and every pattern gives `c_1` in `{0, 2}`, so it is not mobile, and
the true weight rank is 2, not 5. `tests/test_profiles.py` keeps this row as a
standing counterexample to the rank-only test.

**What Theorem A buys.** The test moves from `R^d` with `binom(d, N)` weights to
`R^r` with at most `binom(r + N - 1, N)` patterns, and it is an `S_d` invariant
by construction, so there is no canonicalization, no refinement, no
individualization and no orbit deduplication. The rank-9 canonicalizer's
841,328 leaves on the top level are simply not needed.

## Lemma L: the search is finite in one integer

**Lemma L.** Suppose `S` is nonempty and `tau` is primitive. Put
`B = gcd_{a,b}(v_a - v_b)` and `l_a = (v_1 - v_a) / B`. Then the `l_a` are
integers, `l_1 = 0`, they strictly increase, their gcd is one, every `c` in `S`
has the same level sum `t = N v_1 / B`, and `v_a` is proportional to
`t - N l_a`.

*Proof.* Integrality and monotonicity are immediate from the definition of `B`.
For `c` in `S`, `0 = sum_a c_a v_a = N v_1 - B sum_a c_a l_a`, so
`sum_a c_a l_a = N v_1 / B`, which does not depend on `c` and is therefore both
constant and an integer. Substituting `v_a = v_1 - B l_a` and multiplying by
`N / B` gives `(N / B) v_a = t - N l_a`. QED

**Corollary (the divisibility lemma, recovered).** `v_a = v_1 - B l_a` with
`gcd(l) = 1` forces `gcd(v) = gcd(v_1, B)`, which is 1 by primitivity, and
`B | N v_1` then gives `B | N`. For `N = 3` the value gap is 1 or 3, which is
exactly what `docs/arithmetic_level_audit.md` derived for facet mechanisms and
what the enumerated candidates show.

So a profile is named by three finite pieces of data plus one unbounded one: a
composition `m` of `d` into `r` parts, an increasing level set `l` with
`l_1 = 0`, an integer target `t` in `[0, N l_r]`, and the level SPREAD `l_r`.

## Monotonicity, which splits the enumeration in two

`S` grows with `m` and rank is monotone, and mobility gets easier as `m_a` grows.
So admissibility is monotone in the multiplicity vector. Patterns cap at `N`, and
mobility needs at most one spare mode beyond a pattern, so admissibility depends
on `m` only through `min(m_a, N + 1)`.

That gives the implementation its shape. Stage one enumerates level profiles and
is independent of `d`, so it is computed once per `(N, r, spread)` and reused
across ranks. Stage two caps by a composition of `d` and re-tests. On the
uncapped set, `rank(S) = r - 1` is a necessary condition, so a level profile that
fails stage one fails at every rank.

## Calibration: five known populations, reproduced exactly

The repository knows five series by routes that share no code with this one. All
five come out exactly, at level spread bound 14.

| rank | unoriented orbits | oriented | self-paired | trace-alive | trace survivors |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 6 | **8** | **15** | **1** | 12 | **64** |
| 7 | **19** | **35** | **3** | **25** | **549** |
| 8 | **56** | **109** | **3** | **64** | **6,607** |
| 9 | **231** | 455 | 7 | **191** | 135,343 |
| 10 | **1337** | 2,667 | 7 | 655 | 3,723,147 |

Bold entries are the independently known values:
`docs/orbit_canonical_flats.md` for the orbit columns,
`docs/levi_triangularity_gate.md` for the trace-alive column, which that document
calls the mechanism population, and `docs/screening_orbit_structure.md` for the
survivors. The remaining entries are new.

Two stronger checks than counting run as well.

- **Set for set.** At `(3,6)`, `(3,7)` and `(4,7)` the canonical `tau` set equals
  the exhaustive closed-flat enumerator's, with no missing and no extra rows.
- **Against the definition.** Every profile at `(3,6)` through `(3,9)` and at
  `(4,7)`, `(4,8)` is re-tested by materializing its weights and computing the
  affine rank. Theorem A and the definition disagree nowhere.
- **End to end.** `enumerate_profile_survivor_taus` feeds the existing survivor
  stage from profile orbits instead of lattice orbits, and the emitted tau SETS
  are equal to `orbit_native.enumerate_orbit_survivor_taus` at `(3,6)`, `(3,7)`
  and `(3,8)`, 64, 549 and 6,607 rows. The replacement is a replacement, not a
  parallel count.

`(4,7)` and `(4,8)` matter because they are a different particle number, so the
reduction is not an `N = 3` accident.

## Cost, against the two routes `docs/cocircuit_chow_decoder.md` measured

That document benchmarked three enumerators reaching the same codimension-one
objects at `(3,8)`, declined the direct cocircuit route, and named the open
question:

> The cheaper experiment, and the one that does not add a dependency, is to ask
> whether the canonical-augmentation enumerator's known cost concentration, 78
> percent of `(3,9)` time in the top level, can be attacked directly.

This is that experiment, and the answer is that the concentration disappears
rather than being attacked: an orbit is a sorted multiset, so the top level
needs no canonicalization at all. Same machine as that table,
`Linux-6.18.5-fc-v20-x86_64`, Python 3.14.0rc2:

| Enumerator at `(3,8)` | Objects it must hold | Seconds |
| --- | ---: | ---: |
| direct cocircuit reverse search | 2,042,570 search nodes, 166,420 hyperplanes | 278.7 |
| materialising closed-flat lattice | 1,369,357 flats | 140.3 |
| orbit-canonical augmentation | 314 orbit representatives | 10.1 |
| **profile enumeration** | **109 oriented orbits** | **2.4** |

The first three rows are quoted from that document. The fourth is a cold run,
no cross-rank cache.

**The seconds column is the least interesting part of this table, and it should
not be read as the finding.** That document's caveat cuts both ways: its
cocircuit backend was a first implementation, and so is this one. The object
column is what does not depend on optimization effort. The profile route holds
the ANSWER and nothing else, where augmentation holds 314 representatives
spread across every level of a lattice it must still traverse.

That is why the gap widens rather than staying a constant factor. At `(3,9)`,
11.2 seconds cold against 401. At `(3,10)`, 51.7 seconds cold against a 520
second augmentation run that was restored from checkpoint rather than run
straight through, so that pair is indicative rather than a clean comparison. At
`(3,11)` augmentation does not finish, which is the rank the cocircuit document
identifies as the first where a new backend would produce anything the
repository does not already have.

Two things this does NOT license. It is not a claim that the cocircuit work is
superseded: that document's decoder and its two theorems are about a different
object, and its cocircuit backend remains an independent low-rank cross-check
of exactly the population enumerated here. And a faster orbit stage does not
close `(3,11)`, because the spread bound below, though proved, is far too large
to enumerate against, and the arrangement layer past it is untouched.

## A target-rank check that does not need the spread bound

Whatever the spread bound turns out to be, rows the repository has already
certified at rank 11 must appear in the enumeration, and they do.

- The four level-five `(3,11)` rows carrying exact Ressayre certificates
  (`results/data/level5_ressayre_3_11.json`) are all present, at level spread 2.
- Of the 13 Altunbulak-Klyachko partial-family rows for `(3,11)`
  (`docs/partial_families_3_11_3_12.json`), 12 are present.

The single absence is a correct rejection rather than a gap. The missing row is
the one tagged WEAK, `lambda_1 + lambda_11 <= 6/5`, whose primitive normal is
`tau = (3,-2,-2,-2,-2,-2,-2,-2,-2,-2,3)`. No triple of its entries sums to zero,
since the available sums are `4`, `-1` and `-6`, so it touches no weight, is not
a supporting hyperplane, and cannot be a Ressayre element. It is a valid
inequality that is not tight anywhere, which is exactly what its WEAK tag says.
Both halves are pinned in `tests/test_profiles.py`.

## Refuted: the bounded-hole arithmetic-level grammar is not a completeness rule

`docs/arithmetic_level_audit.md` found 78 of 83 facet mechanisms to be exact
arithmetic progressions of levels and the other 5 to be progressions with one
hole, and proposed parameterizing a constructor by a progression with at most one
hole. `docs/mechanism_grammar_frontier.md` carried that forward as the surviving
conjecture, with the hole switching on at rank 9.

Scored against the complete CANDIDATE universe rather than against facets, a
bounded hole budget fails, and the failure grows with rank. A budget of `h` means
level spread at most `r - 1 + h`, so `h = 0` is an exact progression.

| rank | h = 0 | h = 1 | h = 2 | h = 3 | complete |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 47 | **56** | 56 | 56 | 56 |
| 9 | 136 | 210 | 225 | **231** | 231 |
| 10 | 398 | 838 | 1,095 | 1,244 | **1,337** |

At rank 10 a one-hole budget misses 499 of 1,337 orbits, a three-hole budget
still misses 93, and a four-hole budget still misses 36. The gap widens with
rank rather than closing. The grammar is a good description of where facets end
up. It is
not a pruning rule, and using it as one would silently delete candidates. The
table is pinned in `tests/test_profiles.py` so the promotion cannot happen by
accident.

This is consistent with, not contrary to, the earlier audits: they measured
facets and mechanisms, this measures candidates, and the candidate population is
the one a coverage theorem has to cover.

## The remaining gap, stated exactly

Lemma L makes the enumeration finite once the level spread is bounded, and
Theorem S below bounds it. The candidate universe at each rank is therefore a
finite object that this enumerator exhausts in principle, and candidate coverage
is decidable rather than open. That discharges the obligation recorded in
`docs/fixed_n3_generator_methodology.md` as a question of principle.

**What is not closed is the distance between the bound that is proved and the
bound that can be run.** Theorem S gives spread at most `2 sqrt(r-1) (N sqrt 2)^(r-2)`,
which is about `2.8e6` at `N = 3, r = 11`, against a realized maximum of 20 and a
run that already costs 1,227 seconds at bound 20. So every rank-11 number below
is still a lower bound. The reason has changed rather than disappeared: the run
stops at spread 20 because 20 is affordable, not because anything proves 20 is
enough.

The protocol used instead is saturation, and it is calibrated where the truth is
known. Raise the bound; if the population stops growing and the realized maximum
spread sits strictly inside the bound, report the plateau.

| rank | realized max spread | plateau equals known orbit count |
| ---: | ---: | --- |
| 6 | 3 | yes, 8 |
| 7 | 4 | yes, 19 |
| 8 | 6 | yes, 56 |
| 9 | 10 | yes, 231 |
| 10 | 14 | yes, 1337 |

At rank 11 the protocol does NOT saturate. The ladder, recorded in the artifact
as `rank_eleven_spread_ladder` and counting oriented orbits, reads

| spread bound | 12 | 14 | 16 | 18 | 20 |
| --- | ---: | ---: | ---: | ---: | ---: |
| oriented orbits | 16,739 | 19,819 | 21,175 | 21,605 | 21,699 |
| unoriented orbits | 8,379 | 9,919 | 10,597 | 10,812 | 10,859 |

and 16 profiles still appear at spread exactly 20, so the realized maximum
reaches the bound rather than sitting inside it. The `(3,11)` candidate universe
is therefore not closed by this route today, and every rank-11 number from it is
a lower bound, conditional on the spread bound in exactly that way.

Two routes out, both well posed:

1. **Sharpen the bound** until it is runnable. Theorem S bounds the spread, but
   five orders of magnitude above the realized maximum. The gap is the work.
2. **Bypass it.** Exact inner and outer closure never needs candidate
   exhaustiveness, and is how ranks 6 through 10 are already closed.

## Theorem S: the pattern set determines the values

The budget argument this document previously proposed is not needed. The
statement follows from Theorem A by a dimension count, in four lines.

**Theorem S.** Let `tau` be admissible with `r >= 2` distinct values `v` and
multiplicities `m`, and let `V = span{ c - c' : c, c' in S }`. Then

```text
V = 1^perp  n  v^perp,   hence   V^perp = span{1, v}.
```

So `S` determines `v` up to the reparametrization `v -> a v + b 1` that Lemma L
already quotients out, and hence determines the primitive level vector `l`
exactly.

*Proof.* Every `c` in `S` has `<1, c> = N`, which is nonzero, so no element of
`S` lies in `V` and `rank_Q(S) = dim V + 1`. Condition (i) of Theorem A is
therefore exactly `dim V = r - 2`. Differences of patterns have coordinate sum
zero, so `V` is contained in `1^perp`; every `c` in `S` has `<c, v> = 0`, so `V`
is contained in `v^perp`. The values are distinct and `r >= 2`, so `v` is not a
multiple of `1` and the two conditions cut out a space of dimension exactly
`r - 2`. A subspace of dimension `r - 2` contained in a space of dimension
`r - 2` is that space. QED

**Corollary (Conjecture S, settled affirmatively).** `S` is a subset of `P(m)`,
which is finite and depends on `m` alone, not on the values. So each
multiplicity vector carries at most `2^|P(m)|` admissible level vectors, and
`sum_a m_a = d` leaves finitely many `m` per rank. The level spread is bounded at
every rank, the enumeration is finite without a spread parameter, and candidate
coverage is decidable.

**Corollary (explicit bound).** Pick `r - 2` independent differences as the rows
of `A` and append the all-ones row to get `A'`. Then `ker A = span{1, v}` by the
theorem, and `ker A' = ker A n 1^perp` is the line through `u = r v - (sum v) 1`,
so Cramer spans it by the vector `g` of signed maximal minors of `A'`. Since
`gcd(u)` divides `r B`, where `B` is the value gap of Lemma L, and `g` is an
integer multiple of the primitive generator, the spread satisfies
`W <= |g_1 - g_r|`. Each pattern row has entries bounded by `N` with absolute
values summing to at most `2N`, hence Euclidean norm at most `N sqrt 2`, and the
appended row has norm `sqrt(r - 1)`. Hadamard then gives

```text
spread  <=  2 sqrt(r - 1) (N sqrt 2)^(r - 2).
```

That is `spread_bound(n, r)` in `profiles.py`. It is true and it is useless as a
search parameter: 12 at `r = 3`, 629,856 at `r = 10` where the realized maximum
is 14, and 2,816,802 at `r = 11` where the realized maximum is at least 20.

**Checked, not assumed.** `level_vector_from_patterns` runs the determination
backwards, recovering `(levels, target)` from the pattern set with the values
never entering. It reproduces the levels of all 3,281 profiles at ranks 6
through 10 and of the 318 at `(4,7)` and `(4,8)`, with no failures. At every one
of those ranks no two profiles sharing a multiplicity vector share a pattern
set, which is the same statement seen from the other side.
`tests/test_profiles.py` pins both.

**What the theorem does not do.** It does not close `(3,11)`. The bound is five
orders of magnitude above what can be enumerated, so the shipped rank-11 numbers
remain lower bounds at the spread actually run. The open problem is no longer
"is the spread bounded" but "is the true bound polynomial in `r`", and the
evidence below is about that question rather than about finiteness.

**Exact, from the complete enumerations.** Maximum spread at the ranks where the
enumeration is closed is 3, 4, 6, 10 and 14 at ranks 6 through 10. That is exact,
not a bound, and it is strictly increasing and far below Theorem S.

**Upper bound on `D`, restricted to `r <= 7` classes.** Writing `D(W)` for the
least ambient rank carrying an admissible profile of level spread `W`, and
minimizing `sum_a min(m_a, N + 1)` over level profiles with at most seven
classes, gives for `W = 1` to `22`:

```text
W     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
D<=7  4  5  6  7  8  8  9  9  9 10 10 10 11 11 11 12 12 12  -  12 12 13
```

with `W = 19` carrying no primitive profile at all in that range. This is an
upper bound on `D`, since a cheaper profile may use more classes, and indeed the
true `D(14) <= 10` while the table reads 11. What it shows is the shape: the cost
of a wide level set climbs with the width, monotonically over the whole measured
range, and it climbs roughly linearly rather than at the exponential rate
Theorem S permits.

That table is regenerated by `scripts/profile_spread_cost.py` into
`results/data/profile_spread_cost.json`, with each row carrying the witness
level set and target that achieves the cost. Two things are checked rather than
asserted. Each stored cost is recomputed from its witness, and for the first
rows every cheaper multiplicity vector is enumerated and shown inadmissible, so
minimality is exhibited and not just reported. And the ceiling is cross-checked
against the exact realized maxima 3, 4, 6, 10, 14 at ranks 6 to 10 from
`generator_frontier_refinement.json`: a profile of spread `W` at rank `d` forces
`D(W) <= d`, so the ceiling must never sit below a realized rank, and it does
not. The two computations share the Theorem A predicate and nothing else.

That gap is now the whole question, and it is a quantitative one. A counterexample
family of unbounded spread at fixed `d` is ruled out by Theorem S, so nothing can
kill the candidate-coverage route outright any more, and the spread-cost table
stops being evidence that a bound exists and becomes evidence about how big it
is. What is still open is whether the true maximum spread at rank `d` is
polynomial, which the exact maxima 3, 4, 6, 10, 14 and the table above both
suggest and neither proves. The budget argument this document previously proposed
as a route to finiteness is still the natural route to a POLYNOMIAL bound: make
`sum_a min(m_a, N + 1) >= f(W)` quantitative and the enumeration becomes
runnable, not merely finite.

## The bottleneck moved, and the measurement says where

With the orbit stage cheap, the next wall is visible in the same run.

| rank | trace-alive orbits | Ressayre trace survivors | survivors per alive orbit |
| ---: | ---: | ---: | ---: |
| 8 | 64 | 6,607 | 103 |
| 9 | 191 | 135,343 | 709 |
| 10 | 655 | 3,723,147 | 5,684 |
| 11 | 2,805 | 108,541,487 | 38,696 |
| 12 | 11,397 | 3,051,122,527 | 267,712 |

Ranks 11 and 12 are conditional on the spread bound, 11 at bound 20 and 12 at
bound 14, so they are lower bounds on the true populations. Rank 13 is reachable
with `--frontier-ranks 11,12,13`; a scratch run put it at 43,719 trace-alive
orbits and 8.6e10 survivors, continuing the trend, but it is not in the shipped
artifact and so is not quoted as a result.

The orbit layer grows by a factor of about four per rank. The ARRANGEMENT layer
grows by a factor of about thirty, because each alive orbit contributes every
arrangement of its multiset with a prescribed inversion number. Screening 108
million rows through a tangent determinant is not a constant-factor problem, and
it is now the binding constraint rather than orbit enumeration.

The honest consequence: `docs/direct_survivor_generation.md` removed the row
SCAN, and it was right that the cost then became the orbit count and the survivor
count. This measurement says the survivor count is the half that matters, and
that the next theorem worth having decides the Ressayre determinant at the orbit
or profile level rather than per arrangement. `docs/screening_orbit_structure.md`
already proves that verdicts are not orbit invariants, so this is a genuinely
open problem and not an engineering task.

## Reproducing

```sh
uv run pytest tests/test_profiles.py tests/test_profile_spread_cost.py
uv run python scripts/profile_native_frontier.py            # calibration only
uv run python scripts/profile_native_frontier.py --frontier # adds ranks 11 and 12
uv run python scripts/profile_spread_cost.py                # the Conjecture S table
```
