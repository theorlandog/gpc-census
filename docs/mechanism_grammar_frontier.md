# The mechanism-grammar frontier: four campaigns reconciled

**Status:** cross-artifact synthesis. Exact on three finite populations, a
proof of exactly one lemma, and no GPC theorem.

**Artifact:** `results/data/mechanism_grammar_synthesis.json`, written by
`scripts/verify_mechanism_grammar_synthesis.py`.

This document reconciles four campaigns that landed separately and whose
headline sentences read as though they disagree. They do not. They are
separated by rank, and once that is made explicit they compose into a single
candidate grammar for GPC facet normals, with one proved constraint, one
refuted conjecture, and one surviving conjecture.

## The inputs and their exact scopes

Four campaigns, plus the arithmetic audit derived from the first two. Scope
discipline matters here more than the numbers, because each measured a
different population.

| Campaign | Population | Exactness |
| --- | --- | --- |
| `levi_fusion_pilot.md` | every published facet of `(3,7)`, `(3,8)`, `(4,8)` | exact, discovery data |
| `levi_fusion_heldout_results.md` | every published facet of `(3,9)`, `(3,10)`, `(4,9)`, `(4,10)`, `(5,10)` | exact, preregistered |
| `generator_compression_trial.md` | primitive `N=3` normals with entries in `[-5,5]` at ranks 7, 8, 9 | exact within a bounded window |
| `full_rank78_compression_audit.md` | the full `(3,7)` and `(3,8)` trace-survivor populations | exact, full population |
| `arithmetic_level_audit.md` | the 83 distinct Levi mechanisms of the first two | exact, post hoc |

The 83-mechanism population is the set of distinct sorted-`tau` shapes: 13
from the pilot, 70 from the held-out audit. Exactly one of the 83, the `(4,8)`
shape, is STRUCTURAL; the other 82 are nonstructural. Every one-hole
mechanism discussed below is nonstructural and held out. The pilot's
compression ratios count all rows including structural ones, while the
held-out P2 ratios are nonstructural-only, which is why the pilot reports 3
signatures at `(4,8)` where the held-out convention reports 2.

## The apparent contradiction

Two sentences, both exact, that look incompatible:

- the bounded trial: **exact arithmetic progressions only is a NO-GO**, because
  one-hole normals pass Hall and three carry modular nonzero determinant
  certificates;
- the full rank-7/rank-8 audit: **every determinant-nonzero rank-7/rank-8
  candidate has an exact arithmetic-progression level set**.

## The resolution: the hole switches on at rank 9

The two statements are about different ranks, and the transition is sharp.

At ranks 7 and 8, on the complete trace-survivor populations, one-hole level
sets exist in quantity but never survive to a nonzero determinant:

```text
(3,7): 549 survivors, 0 with holes
(3,8): 6,607 survivors, 5,588 with 0 holes, 1,019 with 1 hole, 0 with 2+
       of the 1,019, exactly 7 are Hall-feasible
       all 7 have identically zero tangent determinant
```

So through rank 8 the exact-progression grammar is complete, not by absence of
one-hole candidates but because Hall plus an exact determinant expansion kills
every one of them.

From rank 9 the hole becomes load-bearing, and this is established twice over
by independent routes:

1. the bounded search produces three one-hole rank-9 normals with
   deterministic modular nonzero determinant residues, which proves the
   symbolic determinant is nonzero over the integers;
2. five of the 83 published Levi mechanisms are one-hole progressions, all at
   ranks 9 and 10: `(3,10)`, `(4,9)`, `(4,10)`, and `(5,10)` twice.

The second route is the stronger one and was not stated in the source
documents. It does not depend on the bounded window, the modular evaluation,
or any conjecture: these are certified facets of published polytopes whose
level sets are not arithmetic progressions. A grammar that omits the hole
cannot generate the facet system at rank 9.

### A cross-realization agreement at rank 8

Two independently scoped rank-8 searches isolate the *same* one-hole object.
The bounded window `[-5,5]` and the full repository trace-survivor population
both find exactly one one-hole Hall-feasible shape,

```text
tau = (-2, -1, -1, -1, 0, 0, 0, 2)
normalized levels [0, 1, 2, 4], missing 3
7 Hall-feasible placements in both searches
```

The bounded window was not chosen to contain this shape, so the agreement is
evidence that the window is not obviously too narrow at rank 8.

## What is proved: the step divides the particle number

The arithmetic audit reports "the progression step divides `N`" at 83 of 83
mechanisms. That is not evidence for a pattern. It is entailed, and the
distinction matters because a constructor may rely on it.

**Lemma.** Let the distinct levels of a row lie in `a + g*Z` with
`gcd(a, g) = 1`, and suppose some `N`-particle zero-grade allocation exists.
Then `g` divides `N`.

*Proof.* Write the levels as `a + g*j`. A zero-grade allocation has
`sum_j n_j = N` and `sum_j n_j (a + g*j) = 0`, so `N*a + g*S = 0` with
`S = sum_j j*n_j` an integer. Hence `g | N*a`, and `gcd(a, g) = 1` gives
`g | N`. QED

The coprimality hypothesis is automatic for a published row. Normalizing by
`g` makes the gcd of the normalized levels 1, so `gcd` of the distinct levels
equals `gcd(a, g)`; primitivity of `tau` makes that 1. The verifier checks all
three hypotheses (primitive level set, `gcd(a, g) = 1`, at least one
zero-grade allocation) at all 83 mechanisms, so the 83/83 count is a pipeline
consistency check.

Consequence for a constructor: the divisor `g | N` is a *free* parameter
restriction, not a conjecture to be validated. For `N` prime this is severe,
leaving only `g = 1` and `g = N`, which is visible in the measured step
distribution: at `N = 5` the only steps observed are 1 and 5, at `N = 3` only
1 and 3, while `N = 4` also admits 2.

## What is refuted: bounded Levi blocks

The preregistered P1, `r <= 2N` for the number of distinct `tau` levels, is
**FALSE**, at `(3,10)` and `(4,10)`. The sharpest witness is a `(4,10)` row
with every entry distinct:

```text
tau = (0, 1, 2, -3, 3, -6, -4, -5, -2, -1)
r = 10 > 8 = 2N
f = 9,  dim W_0 = 9,  ambient dim = 210
```

The synthesis verifier recomputes this descriptor from `tau` alone rather than
quoting it. The row where the block count is maximal, and equal to the orbital
rank, is also a row where the active zero-grade module is 9 dimensional out of
210. So the block count is the wrong complexity measure and its failure does
not damage the program: `f` and `dim W_0` are the quantities that stay small.
P5, `f_max(3,d) <= d-3`, passes and is saturated at `d = 10`, and P6 keeps the
zero-grade density at or below `7/15`.

## What survives: the one-hole grammar

The surviving conjecture, and the only one worth a proof attempt:

> Every Hall-feasible GPC candidate normal has arithmetic-progression levels
> with at most one missing point.

Its strength is that no two-hole survivor appears in any of three
independently scoped populations:

| Population | Size | One hole | Two or more holes |
| --- | ---: | ---: | ---: |
| published Levi mechanisms | 83 | 5 | **0** |
| bounded `N=3` Hall-feasible placements, ranks 7 to 9 | 537 | 38 | **0** |
| full rank-7/rank-8 trace survivors | 7,156 | 1,019 | **0** |

In the bounded search 135 multi-hole *trace* placements do exist and every one
of them is rejected by Hall, so the mechanism doing the work is the matching
condition rather than the enumeration window.

## The supported constructor architecture

```text
singleton-incidence key
  -> arithmetic progression with step g | N and at most one hole
  -> target-aware fusion-tree dynamic program
  -> Hall obstruction
  -> exact or modular tangent determinant
  -> global redundancy reduction
```

Two of these stages carry caveats that must travel with them.

**Fusion width is an evaluation cost, not a facetness test.** Target-aware
pruning is real: unpruned state sets reach 105 where optimal trees need at
most 5 over the 83 mechanisms, and at most 4 on the full rank-7/rank-8
populations. The synthetic stress test gives a constructive upper bound of 449
balanced-tree states for one consecutive-level `(20,40)` normal family, rather
than all 20-subsets of 40 orbitals; those synthetic normals are not claimed to
be GPC facets, so that number prices verification of a specified normal and
nothing else. And 2,393 *rejected* rank-8 Hall zeros also reach width 4, so
small width does not select facets.

**Incidence tomography is a hash, not a theorem.** The labeled
singleton-incidence vector is injective on all three populations, at 83 of 83
mechanisms, 547 of 547 nonstructural `(3,7)` survivors and 6,605 of 6,605 at
`(3,8)`, with zero nonstructural cross-shape collisions and no additional
resolution from pair incidence. That licenses it as the first canonical key
and motivates a fixed-weight Chow-parameter theorem attempt. It is measured
injectivity on finite sets.

## Evidence boundary

Nothing here is a GPC proof, and the local mechanism picture is not promoted
to a global statement. Specifically:

- the one-hole grammar is a bounded-window and published-facet conjecture, not
  established for any rank above 10 or for any `N` above 5;
- the full-population gates at ranks 7 and 8 (549 and 6,607 survivors, the
  Hall split, the incidence injectivity) were re-read from the audit artifact
  and cross-checked for internal consistency, because the generating
  enumeration was not shipped with that bundle. RESOLVED: they are regenerated
  from the repository's own machinery in `docs/rank9_compression_audit.md` and
  agree exactly, column by column;
- the bounded search is exhaustive only for primitive `N=3` normals with
  entries in `[-5,5]` at ranks 7, 8, 9;
- no local Levi or slice calculation is promoted to a global GPC proof without
  a bound Ressayre certificate, per the validation law.

This program also remains logically independent of the global Khovanskii/SAGBI
result that the stable highest-weight semigroup does not compress across rank.
That campaign refutes a uniform global generating set; this one asks whether
each fixed system's normals are drawn from a small local grammar. The
refutation there does not touch the conjecture here.

## What would falsify the synthesis

In priority order:

1. a Hall-feasible normal with two or more arithmetic holes, anywhere;
2. a published facet whose level set is not contained in an arithmetic
   progression with at most one hole;
3. two nonstructural trace placements with equal labeled singleton incidence
   but different normals;
4. fusion-tree width growing quickly on the rejected-candidate population at
   rank 9, which would remove the evaluation-compression argument;
5. a determinant-nonzero one-hole candidate at rank 7 or rank 8, which would
   move the onset rank below 9 and break the reconciliation above.

## The decisive test has now been run: rank 9

`docs/rank9_compression_audit.md` carries out the audit this document was
waiting on, and it also regenerates the rank-7 and rank-8 gates that the scope
note above flagged as re-read rather than recomputed. Those reproduce exactly,
column by column, so the gap is closed.

The rank-9 outcome sharpens two of the claims above and overturns none.

- **Falsifier 1 was genuinely tested for the first time, and held.** At ranks 7
  and 8 there were no multi-hole survivors at all, so Hall was never asked the
  question and the at-most-one-hole conjecture was nearly vacuous on the full
  population. Rank 9 produces 5,310 multi-hole survivors (1,227 with two holes,
  4,083 with three) and Hall rejects every one. The Hall-feasible hole
  distribution at rank 9 is `{0: 339, 1: 32}`.
- **The rank-9 onset is confirmed on the complete population.** 32 one-hole
  Hall-feasible rows at rank 9, of which exactly 3 have nonzero tangent
  determinant, against 7 one-hole rows at rank 8 of which 0 do. The three
  nonvanishing rows are exactly the bounded search's three certificates, so the
  `[-5,5]` window missed nothing at rank 9.
- Falsifier 3 held at 20x the population: incidence is injective on all 135,341
  nonstructural rank-9 survivors. Falsifier 4 held: maximum optimal width grows
  by one per rank, 3 then 4 then 5, and the width caveat is unchanged because
  the rank-9 population is 99.7 percent Hall zeros with the same width profile.
- Falsifier 5 did not fire: rank 8's one-hole rows are still all
  determinant-zero, so the onset stays at 9.

Every determinant verdict at ranks 7 to 9 is now a proof: 36 rows proved
identically zero by exact expansion, 3 proved nonzero by modular certificate,
none undecided.

The next population that could break the conjecture is rank 10, where the
augmentation enumerator reports 1,337 top orbits. Falsifier 2 remains open at
`N = 4` and `N = 5`, which no `(3,d)` audit reaches.

## Verification

```sh
uv run scripts/arithmetic_level_audit.py --check
uv run scripts/verify_mechanism_grammar_synthesis.py
uv run scripts/verify_full_rank78_compression_audit.py
uv run scripts/verify_generator_compression_trial.py
uv run scripts/verify_levi_fusion_heldout_standalone.py
uv run pytest tests/test_mechanism_grammar_synthesis.py
```

The Levi audit itself regenerates from the repository ledger:

```sh
uv run scripts/levi_fusion_audit.py
```

It reproduces the held-out per-system row counts, signature counts, maximum
fusion ranks and maximum level counts from
`results/data/ordered_representation_stability.json` directly, which closes
the provenance gap left by the standalone verifier's embedded input rows.
