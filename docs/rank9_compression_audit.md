# The rank-9 full-population audit: the conjecture's first real test

**Status:** exact full-population audit, ranks 7 to 9. Every determinant
verdict is a proof.

**Artifact:** `results/data/rank9_compression_audit.json`, written by
`scripts/rank9_compression_audit.py`, replayed by
`scripts/verify_rank9_compression_audit.py`.

This is the decisive test named at the end of
`docs/mechanism_grammar_frontier.md`. It does two things: it closes the scope
gap that document declared, and it runs the rank-9 classification the frontier
was waiting on.

## Closing the declared gap

`full_rank78_compression_audit.json` carried population gates that were re-read
from the artifact rather than regenerated, because the generating enumeration
was never shipped with that bundle. Those gates are now rebuilt from the
repository's own machinery, and they agree exactly:

| Column | `(3,7)` | `(3,8)` | Published | Regenerated |
| --- | ---: | ---: | --- | --- |
| trace survivors | 549 | 6,607 | same | same |
| nonstructural survivors | 547 | 6,605 | same | same |
| Hall feasible / zero | 75 / 474 | 193 / 6,414 | same | same |
| hole distribution | `{0: 549}` | `{0: 5588, 1: 1019}` | same | same |
| singleton signatures | 547 | 6,605 | same | same |
| optimal width distribution | `{1:11, 2:73, 3:465}` | `{1:17, 2:222, 3:3964, 4:2404}` | same | same |
| one-hole Hall-feasible rows | 0 | 7 | same `tau` set | same `tau` set |

The orbit counts (35 and 109) and the trace-dead orbit counts (10 and 45) also
reproduce `docs/screening_orbit_structure.md` and
`docs/direct_survivor_generation.md`. Ranks 7 and 8 are therefore a
reproduction test with published answers, run before rank 9 was attempted.

## How rank 9 is reachable at all

Rank 9 has about 20 million oriented rows, which is why the audit stopped at
rank 8. Three existing repository results remove that wall, and this audit is
mostly an assembly of them:

1. the trace test is `inv(sigma h) == B` with `B` an `S_d` invariant, so it
   factors into an orbit invariant and a permutation statistic;
2. `enumerate_orbits_by_augmentation` reaches rank 9, giving 231 unoriented
   top-level orbits, which become 455 oriented sorted-`tau` orbits here;
3. `generate_trace_survivors` produces an orbit's survivors directly, in time
   proportional to the output rather than the orbit.

264 of the 455 rank-9 orbits are killed outright by `orbit_is_trace_dead`
before a single arrangement is built. The whole rank-9 pass takes 25 seconds
once orbit representatives are available.

## The population

| System | Orbits | Trace survivors | Hall feasible | Hall zero | Feasible share |
| --- | ---: | ---: | ---: | ---: | ---: |
| `(3,7)` | 35 | 549 | 75 | 474 | 13.66% |
| `(3,8)` | 109 | 6,607 | 193 | 6,414 | 2.92% |
| `(3,9)` | 455 | 135,343 | 371 | 134,972 | 0.27% |

Survivors grow 12x then 20x while the Hall-feasible count grows 2.6x then 1.9x,
so the feasible share falls by roughly an order of magnitude per rank.

## Result 1: the conjecture faced 5,310 chances to fail and took none

This is the finding, and it is a change in the *kind* of evidence rather than
more of the same.

At ranks 7 and 8 the at-most-one-hole conjecture was close to vacuous on the
full population: there were **no** multi-hole survivors at all, so the Hall
condition was never actually asked the question. Rank 9 is the first rank where
multi-hole rows exist, and they exist in quantity:

```text
(3,7) holes: {0: 549}
(3,8) holes: {0: 5588, 1: 1019}
(3,9) holes: {0: 83885, 1: 46148, 2: 1227, 3: 4083}
```

5,310 rank-9 survivors carry two or three holes. **Hall rejects every one of
them.** The Hall-feasible hole distribution at rank 9 is

```text
{0 holes: 339, 1 hole: 32}
```

so the surviving conjecture is confirmed where it could first have died:

> Every Hall-feasible GPC candidate normal has arithmetic-progression levels
> with at most one missing point.

Note also that the three-hole rows outnumber the two-hole rows, 4,083 to 1,227,
so the multi-hole population is not a thin boundary effect that Hall might
plausibly clip by accident.

## Result 2: the onset is confirmed on the complete population

The frontier document predicted, from a bounded window and from published
facets, that the hole becomes load-bearing at rank 9. The full population
confirms it:

| System | One-hole Hall-feasible | Determinant nonzero |
| --- | ---: | ---: |
| `(3,7)` | 0 | 0 |
| `(3,8)` | 7 | **0** |
| `(3,9)` | 32 | **3** |

Rank 8 has one-hole rows that pass Hall and every one has identically zero
tangent determinant. Rank 9 is the first rank where a one-hole row survives to
a nonzero determinant. The exact-progression grammar is therefore complete
through rank 8 and incomplete from rank 9, now measured on the complete
population rather than inferred.

### The three nonvanishing rows are exactly the bounded search's certificates

```text
tau = (-5, -3, -2, -1, 0, 1, 2, 3, 4)   size 36, 8 variables
tau = (-4,  4, -3, -2, -2, -1, 0, 3, 1)  size 27, 8 variables
tau = ( 0,  2, -2, -1, -1, -1, -1, 0, 0) size 15, 10 variables
```

These are, as a set, the three modular certificates of
`docs/generator_compression_trial.md`, found there by an exhaustive search of a
`[-5,5]` window and here by classifying the entire rank-9 population. The
bounded window missed nothing at rank 9, which is a much stronger statement
about that window than the trial could make on its own.

## Every verdict is a proof

Determinant verdicts are decided in an order chosen so that neither branch
rests on sampling:

1. modular evaluation at random integer points first, where a nonzero residue
   proves the symbolic determinant is nonzero over the integers;
2. exact signed sparse expansion only when every sample vanishes, which is
   exactly the cheap case, because the expansion prunes to an empty state set
   as soon as a level cancels.

Across all 39 one-hole Hall-feasible rows at ranks 7 to 9:

```text
36 proved identically zero by exact expansion
 3 proved nonzero by modular certificate
 0 undecided
```

An earlier pass used an exact-expansion size cutoff and left four rank-9 rows
of size 21 undecided, with modular sampling returning zero every time. Those
four are now proved identically zero: the expansion at size 21 completes
immediately, because a vanishing determinant is the fast case. Nothing in the
audit is left resting on a sample.

## The full determinant split, and an independent campaign that agrees

The audit classifies **every** Hall-feasible row, not only the one-hole ones:

| System | Hall feasible | Determinant nonzero | Identically zero | Unresolved |
| --- | ---: | ---: | ---: | ---: |
| `(3,7)` | 75 | 50 | 25 | 0 |
| `(3,8)` | 193 | 138 | 55 | 0 |
| `(3,9)` | 371 | **242** | **129** | 0 |

Two of the nonzero rows at each rank are structural: an empty positive side
gives a `0x0` matrix whose determinant is the empty product. They are counted
separately so they never inflate the certified count silently.

An independently written rank-9 direct-grammar campaign
(`rank9_direct_grammar_phase.md`, produced outside this repository) reports the
same population from its own orbit-native reconstruction. The two agree on
every quantity both computed:

```text
trace survivors            135,343   135,343
Hall-feasible rows             371       371
Hall-zero rows             134,972   134,972
determinant nonzero            242       242
exact cancellation zeros       129       129
structural among nonzero         2         2
distinct normal mechanisms     191       191
shapes with a Hall placement    52        52
optimal width distribution  {1:24, 2:580, 3:28726, 4:90399, 5:15614}  identical
```

That is two independent implementations agreeing on nine quantities including
a 135,343-row population and a full width distribution, which is a stronger
check on the rank-9 numbers than either campaign alone.

### One coincidence worth not reading as structure

Both "52 grammar shapes with a Hall placement" and "52 retained facets after
exact reduction" appear in that campaign, and they are **different objects that
happen to share a value**. This audit confirms 52 distinct sorted-`tau` shapes
carry a Hall placement at rank 9. The published `(3,9)` system also has 52
facet rows, but those 52 rows carry only **8** distinct sorted-`tau` shapes,
per `levi_fusion_heldout_results.json`. So the coincidence must not be read as
"each Hall shape yields one facet"; the shape-to-facet map is neither injective
nor surjective in the way the matching numbers suggest.

## Result 3: incidence injectivity holds at 135,341 rows

Labeled singleton incidence separates every nonstructural survivor at every
rank:

```text
(3,7):     547 survivors ->     547 signatures, 0 collisions
(3,8):   6,605 survivors ->   6,605 signatures, 0 collisions
(3,9): 135,341 survivors -> 135,341 signatures, 0 collisions
```

That is a 20x larger population than the previous largest test. It remains
measured injectivity on finite sets, not a theorem, and it remains the argument
for using the incidence vector as a canonical first key.

## Result 4: fusion width grows by one per rank

```text
max optimal width: 3 at (3,7), 4 at (3,8), 5 at (3,9)
(3,9) distribution: {1: 24, 2: 580, 3: 28726, 4: 90399, 5: 15614}
```

The caveat from the frontier document survives unchanged and should travel with
any quotation of these numbers: width is an evaluation cost, not a facetness
discriminator. The rank-9 population is 99.7% Hall zeros and its width
distribution is essentially the same as the whole population's, so a small
width does not select facets at rank 9 any more than it did at rank 8.

## What this does not settle

- Rank 9 is `N = 3` only. Nothing here tests `N = 4` or `N = 5` at any rank,
  and the published one-hole mechanisms at `(4,9)`, `(4,10)` and `(5,10)` are
  outside this population.
- Hall feasibility and a nonzero tangent determinant are necessary conditions
  in this pipeline, not facetness. No `(3,9)` facet system is produced here,
  and no local calculation is promoted to a global GPC proof without a bound
  Ressayre certificate.
- The conjecture is now non-vacuously confirmed at rank 9 and remains
  unproved. The next population that could break it is rank 10, where the
  augmentation enumerator reports 1,337 top orbits.

## Verification

```sh
uv run scripts/rank9_compression_audit.py --ranks 7 8 9
uv run scripts/verify_rank9_compression_audit.py
uv run pytest tests/test_rank9_compression_audit.py
```

The rank-9 orbit stage accepts `--checkpoint-dir` so a killed run resumes from
the last completed level.
