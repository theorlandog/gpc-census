# Generator frontier refinement, session report

**Artifact:** `results/data/generator_frontier_refinement.json`, written by
`scripts/profile_native_frontier.py`, guarded by `tests/test_profiles.py`.

**Baseline commit:** `59ba1b5` (tip of `main` at session start), test suite green
at 1021 passed, 22 skipped, 1 xfailed in 1091 seconds on this machine.

This session did not close `(3,11)`. It moved one of the two gates that was
blocking it, measured where the wall actually is, and killed a proposed
compression rule before anything got built on it. What follows says which is
which.

## Headline, in one table

| Claim | Label |
| --- | --- |
| Theorem A: Ressayre admissibility is a condition on the value profile in `R^r`, not on `binom(d,N)` weights in `R^d` | **PROVED**, and checked against the materialized weight rank on every profile at ranks 6 to 9 and at `(4,7)`, `(4,8)` |
| Lemma L: the profile is named by a multiplicity composition, an integer level set and one target, and the value gap divides `N` | **PROVED** |
| Admissibility is monotone in the multiplicity vector and sees `m` only through `min(m_a, N+1)` | **PROVED** |
| The enumerator reproduces all five independently known populations exactly at ranks 6 to 10 | **EXACTLY VERIFIED ON FINITE POPULATION** |
| A bounded-hole arithmetic-level grammar is a candidate-completeness rule | **REFUTED** |
| Theorem S: the pattern set determines the values, so the level spread of an admissible profile is bounded by a function of `d` | **PROVED**, with the explicit bound `2 sqrt(r-1) (N sqrt 2)^(r-2)`, and the recovery checked on every profile at ranks 6 to 10 and at `(4,7)`, `(4,8)` |
| That bound is small enough to enumerate against | **NO**, it is about `2.8e6` at rank 11 against a realized maximum of 20, and closing that gap is the surviving open problem |
| `(3,11)` candidate universe is closed | **NO.** The saturation protocol does not plateau by spread 20 |
| `(3,11)` complete irredundant H-system | **NOT REACHED** |

## Priority 0: reproduction

Pulled `main` at `59ba1b5`. Full suite green as above. `numpy 2.4.1`,
`scipy 1.18.0`, `sympy 1.14.0`, `pulp`, `ortools` per `uv.lock`, Python 3.14.0rc2.

Counts were re-derived rather than inherited from prose. The reference
closed-flat enumerator was re-run at `(3,6)`, `(3,7)` and `(3,8)`: 362, 5,341
and 166,420 hyperplanes in 0.2, 5.6 and 234 seconds, giving 8, 19 and 56
unoriented orbits, matching `docs/orbit_canonical_flats.md`. It was also run at
`(4,7)`, which is not in that table, as the cross-particle-number reference.

## Priority 1 and 7: the orbit stage no longer needs the flat lattice

Full statement, proofs and calibration: `docs/profile_native_generation.md`.

The admissible hyperplane orbits are multisets of integers, and Theorem A tests
one directly in `R^r` where `r` is the number of distinct normal values. Stage
one of the enumeration does not depend on `d` at all. Nothing is canonicalized,
because a sorted multiset is already canonical, so the rank-9 canonicalizer's
841,328 leaves on the top level, which were 78 percent of that run, simply do
not occur.

`docs/cocircuit_chow_decoder.md`, which landed on `main` while this branch was
open, asked exactly this question and left it open:

> The cheaper experiment, and the one that does not add a dependency, is to ask
> whether the canonical-augmentation enumerator's known cost concentration, 78
> percent of `(3,9)` time in the top level, can be attacked directly.

The concentration does not need attacking, it disappears. Extending that
document's `(3,8)` cost table on the same machine: direct cocircuit reverse
search 278.7 s holding 2,042,570 search nodes, materialising flat lattice
140.3 s holding 1,369,357 flats, orbit-canonical augmentation 10.1 s holding 314
representatives, profile enumeration **2.4 s holding 109 oriented orbits**. The
object column is the finding; both cocircuit and profile backends are first
implementations, so neither seconds column is a fair race.

Calibrated against every population the repository already knew:

| rank | unoriented | oriented | self-paired | trace-alive | trace survivors |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 6 | **8** | **15** | **1** | 12 | **64** |
| 7 | **19** | **35** | **3** | **25** | **549** |
| 8 | **56** | **109** | **3** | **64** | **6,607** |
| 9 | **231** | 455 | 7 | **191** | 135,343 |
| 10 | **1337** | 2,667 | 7 | 655 | 3,723,147 |

Bold entries were known independently, from `docs/orbit_canonical_flats.md`,
`docs/levi_triangularity_gate.md` and `docs/screening_orbit_structure.md`. Every
one is reproduced exactly by an algorithm sharing no code with the routes that
produced them. At `(3,6)`, `(3,7)` and `(4,7)` the agreement is set for set, not
just count for count.

`(4,7)` and `(4,8)` are included deliberately: the reduction is not an `N = 3`
artifact.

## Priority 2 and 3: rank-11 status, stated without inflation

**The `(3,11)` candidate universe is not closed.** Lemma L makes the enumeration
finite once the level SPREAD is bounded, and no bound is proved. The saturation
protocol is calibrated at ranks 6 to 10, where raising the bound past the
realized maximum adds nothing and the plateau equals the known orbit count
(realized maxima 3, 4, 6, 10, 14). At rank 11 it does not plateau:

| spread bound | 12 | 14 | 16 | 18 | 20 |
| --- | ---: | ---: | ---: | ---: | ---: |
| oriented orbits | 16,739 | 19,819 | 21,175 | 21,605 | 21,699 |
| unoriented orbits | 8,379 | 9,919 | 10,597 | 10,812 | 10,859 |

Sixteen profiles still appear at spread exactly 20, so the realized maximum
reaches the bound instead of sitting inside it, which is the plateau test
failing. Every rank-11 number below is therefore a LOWER BOUND on the candidate
population, conditional on the spread bound in exactly that way. The ladder is
carried in the artifact as `rank_eleven_spread_ladder`.

Because the candidate stage is open, no `(3,11)` outer H-system was generated,
no vertex enumeration was run, and no attainability search was started. The
self-correcting closure loop of Priority 3 was not exercised. The 19 settled
bracket vertices were not needed and were not used; the bracket settlement is
treated as closed infrastructure per the revised brief.

What this does change for Route B is the cost model, below.

One rank-11 check does hold unconditionally, because it does not depend on the
bound. All four level-five `(3,11)` rows with exact Ressayre certificates appear
in the enumeration, at level spread 2, and 12 of the 13 Altunbulak-Klyachko
partial-family rows appear. The one absence is a correct rejection: the WEAK row
`lambda_1 + lambda_11 <= 6/5` has normal `(3,-2,...,-2,3)`, whose entry triples
sum to 4, -1 or -6 and never to 0, so it touches no weight and is not a
supporting hyperplane at all.

## Priority 10: rank 11 and 12 stage ledger

Measured with the profile enumerator, all conditional on the spread bound
(rank 11 at bound 20, rank 12 at bound 14):

| rank | oriented orbits | unoriented | trace-alive orbits | trace survivors |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 109 | 56 | 64 | 6,607 |
| 9 | 455 | 231 | 191 | 135,343 |
| 10 | 2,667 | 1,337 | 655 | 3,723,147 |
| 11 | 21,699 | 10,859 | 2,805 | 108,541,487 |
| 12 | 135,856 | 67,939 | 11,397 | 3,051,122,527 |

Rank 13 is reachable by the same command with `--frontier-ranks 11,12,13` and
roughly triples the run. It is left out of the shipped artifact deliberately: a
scratch run measured 736,023 oriented orbits, 43,719 trace-alive orbits and
8.6e10 trace survivors there, which is another point on the trend below rather
than a different conclusion, and it is not quoted as a result because nothing
in the shipped artifact backs it. Per-row timings are in the artifact rather
than here, because they are machine dependent.

## The bottleneck has moved, and this is the session's most actionable finding

Survivors per trace-alive orbit: 103, 709, 5,684, 38,696 and 267,712 at ranks 8
through 12.

The ORBIT layer grows by about four per rank. The ARRANGEMENT layer grows by
about thirty, because an alive orbit contributes every arrangement of its
multiset with a prescribed inversion number, and that count is a Gaussian
multinomial coefficient. At rank 11 the pipeline would have to push at least
1.1e8 rows through the Ressayre tangent determinant; at rank 12, at least
3.1e9.

So the binding constraint is no longer candidate enumeration. It is per-row
determinant certification. `docs/direct_survivor_generation.md` removed the row
SCAN and correctly said the cost would then be the orbit count and the survivor
count; this measurement resolves which half matters, and it is the survivors.

`docs/screening_orbit_structure.md` proves that screening verdicts are NOT orbit
invariants, because the dominant chamber is a fundamental domain, so this cannot
be fixed by pushing the determinant up to the orbit. It is an open problem, and
naming it precisely is worth more than another speed run:

> **Next gate.** Find a proof-carrying filter that removes a positive fraction of
> the arrangements of one orbit without evaluating a determinant per arrangement,
> or prove that the number of arrangements reaching a nonzero determinant is
> polynomially smaller than the survivor count.

The Hall structural filter already decides 97.1 percent of rank-8 survivors
structurally, which is the existing partial answer and the natural place to look
for the theorem.

## Refuted: bounded-hole arithmetic levels

`docs/arithmetic_level_audit.md` proposed parameterizing a constructor by an
arithmetic progression of levels with at most one hole, on the evidence that 78
of 83 facet mechanisms are exact progressions and the other 5 have one hole.
`docs/mechanism_grammar_frontier.md` carried it as the surviving conjecture.

Scored against the complete candidate universe, a hole budget `h`, meaning level
spread at most `r - 1 + h`, reaches:

| rank | h=0 | h=1 | h=2 | h=3 | truth |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 47 | **56** | 56 | 56 | 56 |
| 9 | 136 | 210 | 225 | **231** | 231 |
| 10 | 398 | 838 | 1,095 | 1,244 | **1,337** |

One hole misses 499 of 1,337 orbits at rank 10; three holes still miss 93 and
four still miss 36, so the gap widens with rank rather than closing. The
grammar describes where facets land. It is not a candidate pruning rule, and
using it as one deletes candidates silently. The table is pinned as a test.

This does not contradict the earlier audits. They measured facets and mechanisms.
Coverage is a statement about candidates, and that is the population that has to
be covered.

## Priority 4: BDR, not advanced this session

No upstream run was made. The state is as `docs/bdr_constructor_comparison.md`
and `docs/bdr_linear_triangular_gate.md` leave it: the pinned revision
`7ab347d9a7837d68d92bde9fac74606913f395ca` is fetchable, MIT, and has been run
through a throwaway PassageMath environment, reproducing the published
irredundant systems as sets at `(3,7)`, `(4,7)`, `(3,8)`, `(4,8)` modulo the
structural convention. The full-stage ledger, the `PiDominancy` and birationality
stages, and a local exact replay compiler for Condition A all remain unbuilt.
That gate is untouched by this session and is still the principal
representation-theoretic obligation.

Worth recording for whoever picks it up: the profile coordinates of Theorem A
are the natural input for a ramification compiler. A candidate is now named by
`(m, l, t)` rather than by an ambient weight set, which is the compressed
representation Priority 7 asks the ramification machinery to consume.

## What was deliberately not done

- No chemistry, fiber, observable or benchmark work.
- No re-derivation of the rank-11 bracket, which the revised brief marks closed.
- No low-rank performance work on the closed-flat enumerator. It is retained
  unchanged as the completeness oracle that Theorem A is checked against, which
  is now its main job.
- No `(3,11)` H-system, outer vertex enumeration, or attainability search, for
  the reason given above: the candidate stage that would feed them is open.

## Where the next session should start

1. **Make the spread bound runnable.** Conjecture S is settled affirmatively by
   Theorem S in `docs/profile_native_generation.md`: the pattern set determines
   the values, so the spread is bounded at every rank and candidate coverage is
   decidable. The proved bound is `2 sqrt(r-1) (N sqrt 2)^(r-2)`, about `2.8e6`
   at `r = 11` against a realized maximum of 20, so it settles finiteness and
   nothing else. The open problem is now quantitative: is the true maximum spread
   polynomial in `r`? The budget argument that was proposed as a route to
   finiteness is still the natural route here, a lower bound on
   `sum_a min(m_a, N+1)` in terms of the spread. Exact evidence: maximum spread
   is 3, 4, 6, 10, 14 at ranks 6 through 10, strictly increasing, and a
   restricted minimization over profiles with at most seven classes has `D`
   climbing monotonically from 4 to 13 as `W` goes from 1 to 22. A counterexample
   family of unbounded spread at fixed `d` is no longer possible.
2. **Attack the arrangement layer**, per the next gate above. It is the binding
   constraint from rank 11 onward whether or not the spread bound becomes
   runnable.
3. **Route B remains open and is not blocked by either.** Inner and outer closure
   never needs candidate exhaustiveness. A valid outer system for `(3,11)` can be
   assembled from certified rows alone, and the 19 settled bracket vertices are
   inner witnesses for it. The cost model above says nothing against that route;
   it only says the candidate route is not the cheap one.
