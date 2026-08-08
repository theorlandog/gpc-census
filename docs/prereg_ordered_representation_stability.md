# PRE-REGISTRATION: ordered representation stability for GPC facets (2026-08-08)

Written BEFORE the measuring script was run. Standing rules:
`docs/prereg_template.md`. Certifying artifacts:
`scripts/ordered_representation_stability.py`,
`results/data/ordered_representation_stability.json`,
`tests/test_ordered_representation_stability.py`.

## What is being tested

Whether the chamber-facing GPC rows of the census, organized over the rank
variable `d`, form a finitely generated functor on a category of ordered mode
embeddings, in the sense that every row at large `d` is the image of a row of
bounded rank under a composite of the census transport maps.

The target statement, which this campaign either supports or refutes:

> For fixed `N`, every primitive chamber-facing GPC facet is an ordered
> specialization of one of finitely many bounded-rank templates.

## Objects, stated precisely enough to fail

### Rows

A row of the system at `(N,d)` is an affine inequality `a . lambda <= b` with
`a` in `Z^d` and `b` in `Z`, read on the trace slice `sum(lambda) = N` inside
the dominant chamber `1 >= lambda_1 >= ... >= lambda_d >= 0`. The affine form
carries a trace gauge, `(a,b) -> (a + t*1, b + t*N)`, so the canonical
gauge-free representative is the primitive integral homogeneous vector

```text
tau = prim(N*a - b*1),      row reads  tau . lambda <= 0
```

which is exactly `gpc_census.generation.bdr.tau_from_constraint`. Every count
below is a count of primitive taus, never of affine representatives.

### Width and height

```text
width(tau)   = min k such that tau_{k+1} = tau_{k+2} = ... = tau_d
spread(tau)  = max_i tau_i - min_i tau_i
height(tau)  = max(max_i a_i, b) for the canonical affine form of
               bdr.constraint_from_tau (min coefficient gauged to zero)
```

`width` is gauge-free: `width(tau) <= k` holds if and only if the row has an
affine representative supported on the first `k` coordinates. `width` and
`spread` are the two quantities that decide finite generation, because every
covariant map below fixes the leading entries of `tau` up to a single positive
rescaling.

### The five maps

Write `V(N,d)` for the set of valid chamber-facing rows and `F(N,d)` for the
published irredundant system (`gpc_census.constraints`, ranks 6 to 10, nine
systems). `Delta(N,d)` is the moment polytope, `S(N,d)` its dominant-chamber
ordered slice.

1. **Terminal restriction** `rho: V(N,d+1) -> V(N,d)`, contravariant.
   `tau -> prim(tau_1, ..., tau_d)`. Well defined because the face
   `{lambda_{d+1} = 0}` of `S(N,d+1)` is exactly the image of `S(N,d)` under
   zero padding.

2. **Tight padding extension** `pi: V(N,d) -> V(N,d+1)`, covariant.
   `pi(tau) = prim(q*tau_1, ..., q*tau_d, p)` where `t* = p/q` in lowest terms
   with `q > 0` and

   ```text
   t* = min { (- sum_{i<=d} tau_i v_i) / v_{d+1} : v a vertex of S(N,d+1), v_{d+1} > 0 }
   ```

   This is the strongest valid extension of `tau` by one empty mode. It is the
   canonical covariant candidate; the naive extension `a -> (a, 0)` is gauge
   dependent and is measured separately as a control.

3. **Frozen-core restriction** `phi: V(N+1,d+1) -> V(N,d)`, contravariant.
   On affine forms `(a,b) -> ((a_2,...,a_{d+1}), b - a_1)`. Well defined
   because the face `{lambda_1 = 1}` of `S(N+1,d+1)` is exactly the image of
   `S(N,d)` under the frozen-core lift `lambda -> (1, lambda)`.

4. **Tight frozen-core extension** `psi: V(N,d) -> V(N+1,d+1)`, covariant.
   `(a,b) -> ((c, a), b + c)` with `c` the smallest rational making the row
   valid, namely

   ```text
   c* = max { (sum_{i>=2} a_{i-1} v_i - b) / (1 - v_1) : v a vertex of S(N+1,d+1), v_1 < 1 }
   ```

   then cleared to integers and primitivized.

5. **Particle-hole duality** `delta: V(N,d) -> V(d-N,d)`, the existing
   `gpc_census.constraints` flip `a_i -> -a_{d+1-i}`, `b -> b - sum(a)`.
   Involutive, rank preserving.

6. **Non-terminal restriction** `rho_j`, the control: drop coordinate `j <= d`
   from a row at `(N,d+1)`. This is the order-preserving injection
   `[d] -> [d+1]` whose image omits `j`.

### Systems in scope

The nine published tables are the PRIMARY systems. Each dual `(d-N, d)` that is
not itself published is DERIVED by applying `delta` to rows and
`lambda -> (1 - lambda_d, ..., 1 - lambda_1)` to vertices, which is exact and is
already the census convention. That adds `(4,7)`, `(5,8)`, `(5,9)`, `(6,9)`,
`(6,10)` and `(7,10)`, so the `N=5` series becomes `d = 8, 9, 10` and the
held-out test has a live `N=5` route after all. Derived systems are marked in
the artifact and never counted as independent evidence: a derived row is the
image of a primary row under a rank-preserving bijection.

### Descent classification

Each published row `tau` in `F(N,d)` gets exactly one label, decided in this
fixed order, and every inherited label requires an exact round-trip witness,
never a coefficient resemblance:

```text
STRUCTURAL        tau is a chamber or Pauli row
PADDING           exists sigma in F(N,d-1) with pi(sigma) = tau and rho(tau) = sigma
FROZEN-CORE       exists sigma in F(N-1,d-1) with psi(sigma) = tau and phi(tau) = sigma
PARTICLE-HOLE     exists sigma in F(d-N,d) with delta(sigma) = tau
PRIMITIVE         none of the above
```

The **core** of a row is the result of iterating PADDING and FROZEN-CORE
descent until a PRIMITIVE row is reached; its rank is the **core rank**.
Particle-hole does not change `d` and therefore never lowers core rank; it is
recorded but excluded from the descent recursion.

## Holdout and its provenance

The corpus is the nine published census systems, `1 + 4 + 31 + 52 + 93` rows at
`N=3` for `d = 6..10`, `15 + 60 + 125` at `N=4` for `d = 8..10`, and `161` at
`N=5, d=10`. Total 542 rows.

Under standing rule R2 the inherited rows are **not** independent
observations: a PADDING or FROZEN-CORE row is a function of its donor, so it
can confirm only the invariance of the derivation. Every prediction below
therefore reports `scope`, `independent_scope` (rows whose core is not already
in the corpus at a lower rank, that is, the PRIMITIVE rows) and
`independent_distinct_cores`.

The held-out test P8 uses ranks 6 to 9 as the generator set and rank 10 as the
holdout. The `(5,10)` system has no `N=5` ancestor in the census, so its
antecedent is empty by construction for padding and frozen-core descent from
`N=5`; per R3 an empty antecedent is UNINFORMATIVE, and P8 states the `(5,10)`
sub-prediction separately for that reason.

## Reconnaissance already performed, disclosed before scoring

Rule R4 forbids presenting a prediction that data in hand already settles. The
following was run against the published tables before this file was written,
and is disclosed here so that no prediction below can be read as a discovery
when it is not:

1. The naive extension `a -> (a,0)` on published affine representatives is
   valid at the next rank for 0/1, 4/4, 10/31, 44/52 rows in the `N=3` series
   and 8/15, 19/60 in the `N=4` series. So the naive extension is **not** a
   total map. This motivated defining `pi` by tightness instead.
2. Terminal restriction of every published row at `(N,d+1)` was valid at
   `(N,d)` in all five consecutive pairs, and exactly `|F(N,d)|` of the
   rank-`d+1` rows restricted onto published rank-`d` rows.
3. The affine coefficient ranges of the published tables are visible in the
   source data: maximum coefficient `1, 1, 2, 3, 4` across `(3,6)..(3,10)` and
   maximum `7` at `(5,10)`.

Consequently **P1 and P7 are confirmatory and are scored as CONFIRMATORY, not
as discoveries.** P2, P3, P4, P5, P6, P8, P9 are genuinely open at the time of
writing. Widths have not been computed at all.

## Power caveat, stated before scoring

A PASS on any prediction is a statement about ranks 6 to 10 at `N = 3, 4, 5`
only. Nine systems cannot establish a theorem about all `d`, and this campaign
will not claim one. The asymmetry matters and is the reason the campaign is
worth running: a single row at rank 10 that is provably not in the image of the
transport maps from lower rank **refutes** finite generation in the proposed
category outright, whereas no finite table can prove it.

Rank-10 rows are additionally conditional on AK2008 completeness at rank 10
(`results/data/PROVENANCE.md`). A refutation resting only on rank-10 rows
inherits that condition; a refutation available already at rank 9 does not, and
the artifact must record which of the two it is.

## Predictions

**P1 (terminal restriction is total).** Every published row at `(N,d+1)`
restricts along `rho` to a row valid at `(N,d)`, in all five consecutive pairs.
Expectation: PASS, and CONFIRMATORY, because the face
`{lambda_{d+1} = 0}` of the ordered slice is exactly the lower slice, so this
is a theorem rather than a measurement.

**P2 (non-terminal restriction fails).** There exist `(N,d+1)`, a published row
`tau`, and an index `j <= d` such that dropping coordinate `j` yields a row
that is NOT valid at `(N,d)`. Expectation: PASS, because the ordered slice has
no face isomorphic to `S(N,d)` at a non-terminal coordinate. Scored by the
count of `(row, j)` pairs that fail, per system.

**P3 (tight padding is a section).** `rho(pi(sigma)) = sigma` for every
published row of every system, and `pi(sigma)` is valid at `(N,d+1)` by
construction. Expectation: PASS. A FAIL means the tight extension is not
well posed and the covariant candidate dies immediately.

**P4 (extension does not preserve facetness).** The fraction of published rows
at `(N,d)` whose tight extension is a published row at `(N,d+1)` is strictly
below 1 for at least one consecutive pair. Expectation: PASS. This separates
"valid row" from "facet" and is the reason the brief demands the distinction.

**P5 (primitive counts do not vanish).** In the `N=3` series the PRIMITIVE row
count is nonzero at `d = 8, 9, 10` and non-decreasing over those three ranks.
Expectation: PASS, that is, no stabilization. A FAIL, meaning the primitive
count drops to zero, would be direct evidence FOR finite generation and would
be the more interesting outcome.

**P6 (width is unbounded).** The maximum width over published rows strictly
increases at least twice within the `N=3` series `d = 6..10`, and the maximum
width at `(3,10)` is at least 8. Expectation: PASS. Bounded width would be
strong evidence FOR a finite-generation statement, because a row of width `<= w`
is pulled back from rank `w`.

**P7 (height is unbounded).** `spread` and `height` maxima strictly increase
across the `N=3` series. Expectation: PASS, and CONFIRMATORY per the disclosure
above.

**P8 (held-out rank 10).** Take every published row at ranks 6 to 9 as
generators. Close under `pi`, `psi` and `delta` without adding any generator.
Then:

- P8a: coverage of the 93 published `(3,10)` rows is strictly below 100 percent,
  and the pre-registered numeric window is 30 to 80 percent;
- P8b: coverage of the 125 published `(4,10)` rows is strictly below 100
  percent;
- P8c: coverage of the 161 published `(5,10)` rows is below 30 percent.
  Its antecedent is thin by construction, see the holdout section, so a
  low value here is weak evidence and is labeled as such.

Expectation: PASS on P8a, P8b, P8c. A 100 percent coverage anywhere would be a
major positive result and would immediately promote the target theorem to a
serious conjecture.

**P9 (composability).** For every published row `sigma` at `(N,d)` with
`d + 2 <= 10`, `pi(pi(sigma))` is valid at `(N,d+2)` and
`rho(rho(pi(pi(sigma)))) = sigma`. Expectation: PASS. A FAIL means the
covariant maps do not compose and no functor exists even on the poset.

## Scoring rules

PASS / FAIL / UNINFORMATIVE / ALREADY-FALSIFIED / CONFIRMATORY per prediction,
decided by `results/data/ordered_representation_stability.json` alone. Each
verdict reports `scope`, `independent_scope` and `independent_distinct_cores`.
No prediction is edited after the run. Every count in the write-up is read from
the artifact, never retyped from a run log.

## What would count as evidence against finite generation

Any one of these, taken from the artifact:

1. a published row at rank `d` that is PRIMITIVE, that is, not the tight
   image of any lower-rank published row under `pi`, `psi` or `delta`, at three
   or more consecutive ranks with non-decreasing count;
2. unbounded width, since a covariant map that fixes the leading entries cannot
   raise width beyond the donor's rank;
3. unbounded spread under maps that rescale `tau` by a bounded factor;
4. held-out coverage strictly below 100 percent at rank 10 from a rank-9
   generator set.

## What would count as evidence for finite generation

1. primitive counts falling to zero at some rank and staying there;
2. bounded width across the whole census;
3. held-out coverage of 100 percent at rank 10;
4. `pi` preserving facetness exactly, which would make `F(N,-)` a subfunctor
   rather than merely `V(N,-)`.

## SCORED (2026-08-08, after running; predictions above are unedited)

Verdicts are read from `results/data/ordered_representation_stability.json`
(`scorecard` block), which the measuring script computes so that no verdict is
retyped by hand. Independent scope is the PRIMITIVE count, since padding and
frozen-core rows are functions of their donors.

| prediction | verdict | scope | independent scope | measured |
|---|---|---|---|---|
| P1 terminal restriction total | CONFIRMATORY-PASS | 726 | 293 | 724 nontrivial, all valid; 2 trivial |
| P2 non-terminal restriction fails | PASS | 6258 | 293 | 3048 invalid pairs of 6258 |
| P3 tight extension is a section | PASS | 620 | 159 | 620 round trips exact |
| P4 extension does not preserve facetness | PASS | 620 | 159 | 5 misses |
| P5 primitive counts do not vanish | **FAIL** | 176 | 89 | 27, 21, 41 at (3,8..10) |
| P6 width unbounded | PASS | 181 | 94 | 4, 6, 7, 8, 9 |
| P7 height unbounded | **CONFIRMATORY-FAIL** | 181 | 94 | spread 2, 3, 12, 12, 12 |
| P8a holdout (3,10) | PASS | 93 | 41 | 55.91 percent |
| P8b holdout (4,10) | PASS | 125 | 45 | 64.0 percent |
| P8c holdout (5,10) | **FAIL** | 161 | 56 | 64.6 percent |
| P9 composability | PASS | 86 | 69 | 86 exact |

### The three failures, with mechanism

**P5.** The nonzero clause held at every rank; the non-decreasing clause did
not, since the primitive count at `N=3` goes 27, 21, 41 across ranks 8, 9, 10.
Monotonicity was asserted with no mechanism behind it and it is false. The
substantive conclusion, that primitive rows never stop appearing within the
census, is unaffected, but the dip removes any licence to extrapolate.

**P7.** Flagged CONFIRMATORY in advance because the coefficient ranges were
already visible in the tables. Spread stabilizes at 12 over ranks 8 to 10 and
affine height goes 9, 9, 15, so neither is monotone. Height is not the
obstruction, which closes one route rather than opening one.

**P8c.** Predicted below 30 percent, measured 64.6 percent. This is a
pre-registration bookkeeping error, not a surprise about the mathematics: the
prediction was drafted for a corpus without the derived dual systems, the
"Systems in scope" section was then added to include them, and P8c was not
restated to match. Including the duals opens the `psi` route from `(4,9)` and
the `pi` route from the derived `(5,9)`, which is where the coverage comes
from. Scored as written, and failed.

### Predictions that passed but should not be banked

**P6** passed and is worthless as evidence. The artifact records
`max_width_primitive` and `max_width_inherited` as equal at every system, so
width does not separate inherited rows from primitive ones. `pi` raises width
by construction, so unbounded width is a property of the transport maps rather
than an obstruction to generation. The pre-registered rationale was wrong even
though the prediction came true.

**P4** passed on a technicality. All five misses are the structural chamber
rows or the non-full-dimensional Borland-Dennis system. Once those two families
are excluded by name, facet preservation is exact at 612 of 612, which is the
opposite of what P4 was written to detect.

### Unregistered findings

Two results were not predicted because the questions were not asked in advance,
and they are labeled as unregistered in the artifact rather than folded in with
the scored predictions:

1. no nonstructural row is permutation stable, 0 of 901, which is the
   obstruction that closes the `FI` repair;
2. all 907 published rows are genuine facets, certified by exact affine rank on
   the tight vertices, which is what makes the facet counting in the disproof
   legitimate. No prior in-repo result established irredundancy.
