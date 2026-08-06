# Completeness without candidate exhaustiveness

Status: **P = O PROVED at all nine census systems** -- (3,6), (3,7), (3,8),
(3,9), (3,10), (4,8), (4,9), (4,10) and (5,10) -- with candidate exhaustiveness
unused. Ranks 9 and 10 had no in-repo completeness proof before this, and no
`N > 3` system did.

Certifying artifacts: `scripts/census_inner_sandwich.py`,
`results/data/census_inner_sandwich.json`,
`tests/test_census_inner_sandwich.py`. Machinery: `gpc_census.generation`
(screening), `gpc_census.exact_polytope` (vertex enumeration),
`scripts/verify_states_standalone.py` (attainment).

## The reframe

`docs/fixed_n3_generator_methodology.md` states the claim ladder: a complete
GPC system needs "an exhaustive-enumeration proof **or** an independent exact
inner/outer closure". The (3,6) theorem took the second route, and the second
route does not care where the candidate rows came from.

That is the whole point. Exhaustive candidate enumeration is the expensive and
unfinished obligation. It can be sidestepped at any rank where both halves of a
sandwich are available, and both halves are available at every rank the census
covers.

    P subset O    every row of O is a replayed exact Ressayre certificate.
    O subset P    every vertex of O is the spectrum of a certified state.
    ---------------------------------------------------------------
    P = O         and candidate exhaustiveness never entered.

## Why both halves are cheap

**The outer half is cheap because validity is not exhaustiveness.** Certifying
that a SUPPLIED row is an exact Ressayre inequality is screening, at about
0.14 ms per row. Only proving that no OTHER row exists was ever expensive.
Every stored row certifies, at every system:

| system | rows | certified | rejected | unresolved |
|---|---|---|---|---|
| (3,7) | 4 | 4 | 0 | 0 |
| (3,8) | 31 | 31 | 0 | 0 |
| (3,9) | 52 | 52 | 0 | 0 |
| (3,10) | 93 | 93 | 0 | 0 |
| (4,8) | 15 | 15 | 0 | 0 |
| (4,9) | 60 | 60 | 0 | 0 |
| (4,10) | 125 | 125 | 0 | 0 |
| (5,10) | 161 | 161 | 0 | 0 |

**The inner half is free because the census already paid for it.** There are
799 certified extremal states, each an explicit state carrying an exact
characteristic-polynomial certificate that its 1-RDM spectrum is the claimed
vertex. The plethysm route used at (3,6) does not scale, and it does not have
to: the states are already there.

## Non-circularity, stated precisely

This is the one place the argument could be fooled, so it is worth being exact
about what leans on what.

The census VERTEX LISTS were derived from the constraint tables. Checking
"the vertices of O are census vertices" would therefore be circular, and that
check is NOT the inner half.

The attaining STATES were not derived from the tables. Each is an explicit
state whose spectrum is verified by a characteristic-polynomial identity in
exact arithmetic, referencing no constraint system at all. That is what makes
the inner half independent evidence, and it is why the script re-runs the
state certificates rather than trusting the vertex list.

## Result

Every census system closes.

| system | outer rows | vertices | attained | closed | time |
|---|---|---|---|---|---|
| (3,6) | 1 + 3 equalities | 4 | 4 | **yes** | 0.7 s |
| (3,7) | 4 | 10 | 10 | **yes** | 0.1 s |
| (3,8) | 31 | 38 | 38 | **yes** | 1.8 s |
| (3,9) | 52 | 58 | 58 | **yes** | 3.2 s |
| (3,10) | 93 | 113 | 113 | **yes** | 19 s |
| (4,8) | 15 | 22 | 22 | **yes** | 0.3 s |
| (4,9) | 60 | 103 | 103 | **yes** | 7.3 s |
| (4,10) | 125 | 159 | 159 | **yes** | 26 s |
| (5,10) | 161 | 292 | 292 | **yes** | 177 s |

## The wider systems, and the label that got them there

The census also covers (4,8), (4,9), (4,10) and (5,10). The INNER half ran
there unchanged from the start: every vertex of `O` is attained by a certified
state, 22, 103, 159 and 292 respectively.

The OUTER half did not, at first. The screening pipeline hard coded the
particle number, so at `N = 4` and `N = 5` the validity of the rows rested on
their published source rather than on a replayed certificate. The first version
of this script nonetheless reported those systems as closed, which was an
overstatement of exactly the kind the claim ladder exists to prevent. It was
corrected to `inner_half_only`: inner certified, outer cited, excluded from
`closed_systems`.

That label is what made the fix obvious. Nothing beneath the pipeline was ever
special to `N = 3` -- `exterior_weights`, `negative_roots`, the tangent
determinant, `certify_candidate` and the full-dimension witness are all written
for `wedge^n C^d` -- so the restriction was seven literal `3`s in
`candidate_pipeline.py` plus a `combinations(range(d), 3)` in
`full_dimension.py`. Threading a `particle_number` through both closes all four
systems, and `screen_tau_candidates(..., particle_number=n)` is now the only
call the script makes.

Two things had to be checked before believing it, because plumbing that always
says yes proves nothing:

- **it still rejects.** Bumping one coefficient of each (4,8) row rejects
  15 of 15 and certifies 0.
- **the parameter is load bearing.** Screening the (4,8) rows as if they were
  `N = 3` does not certify them all.

Both are tests, not observations. The `inner_half_only` status stays in the
code and keeps its own test: it is now unused, but it is the labelling that a
future system with a citable outer half must fall back to.

The full-dimension witness generalised for free. The BDR moment-cone theorem
needs a nonempty interior, and the same deterministic vector construction gives
zero compact stabiliser at every wider system: skew and symmetric ranks
(28, 36), (36, 45), (45, 55) and (45, 55) against their required values. The
`n = 3` certificates are byte-identical to before, which is the regression that
matters.

## The predicted wall was not the wall

Exact vertex enumeration of a highly degenerate polytope in affine dimension
`d - 2` with tens of rows was expected to be the obstacle, with naive
active-set costing `C(F,s) * s^3`. It was not: the double-description
enumerator in `gpc_census.exact_polytope`, built for the orbit-closure layer,
does every census system and reproduces the certified vertex counts exactly.
No lrs, no pivoting rewrite, no symmetry reduction was needed.

It did start to bite at the wide end, and the fix was inside the existing
enumerator rather than a new one. The adjacency test dominates `cut`, and the
naive form recomputed a vertex's tight set once per PAIR and an exact rank once
per PAIR. Computing each tight set once per vertex, rejecting pairs whose
common tight rows cannot reach rank `n-1` on an integer count, and memoising
the exact rank on the common tight set are all identities, not heuristics:

| system | enumeration before | after |
|---|---|---|
| (3,9) | 12 s | 1.5 s |
| (3,10) | 247 s | 11 s |
| (4,9) | 57 s | 4.3 s |
| (4,10) | 211 s | 22 s |
| (5,10) | over 40 min, abandoned | 166 s |

`(5,10)` is the honest measure of the win: two separate runs were killed after
26 and 40 minutes inside its enumeration, and it now finishes in under three.
A test pins the identity claim by running the whole `(3,8)` enumeration twice,
once through `cut` and once testing every pair with the unmemoised
`_adjacent`, and comparing the cut sets at every step.

That matters for planning: the three attacks proposed for vertex enumeration
(lrs-style pivoting, symmetry-aware orbit enumeration, containment checking
against the known answer) are all still available, and none of them has had to
be spent yet.

## The self-improving loop, and what it would look like

The script is written so that a vertex which resists is not an error but a
pointer. Each vertex of O either receives a certified attaining state, in which
case it is an inner certificate, or it resists, in which case it indicates a
missing facet. The repository already owns the tooling that turns a resisting
vertex into a row: the contraction attack is a distance oracle (a positive
floor means outside), the SOS pipeline supplies non-attainability evidence, and
the level-5 Ressayre machinery turns a separating direction into a proof. Add
the row, re-enumerate, repeat; termination is a proof rather than a stopping
rule, and every intermediate state is honestly a relative system.

**BOUNDED NULL, recorded as such.** At the systems run here the loop never
fired: zero vertices resisted at all nine. So the loop is implemented and
exercised only on its trivial branch, and it has NOT been demonstrated to
recover a facet. Its first real test is a rank where the inner side is thin,
which is (3,11): 19 TRUE vertices from a bracket rather than a census, so a
resisting vertex there is much more likely to be a genuine missing facet than a
hard attainment problem.

## What this does and does not license

**Does.** P = O at all nine census systems, by an argument that uses no
exhaustiveness claim. Ranks 9 and 10, and every `N > 3` system, now have an
in-repo completeness proof rather than a citation of the published tables plus
a census that was built from them.

**Does not.** It says nothing at a rank with no certified attaining states,
which is exactly the frontier: at a NEW rank there is no census to supply the
inner half, and either the loop above must manufacture it vertex by vertex or
the exhaustive route must be finished. It also does not retire the exhaustive
route as a goal, since a rank whose sandwich fails to close needs it.

## The orbit reduction, if the exhaustive route is ever finished

Should exhaustive enumeration be wanted anyway, `C(W, d-1) / d!` is 5.8e3 at
(3,8) and 1.2e5 at (3,9), so canonicalising flats under the `S_d` action is a
several-order win. It is also the one pruning that is provably safe under the
standing rule that empirical patterns may prioritise candidates but may not
discard them: the moment polytope is Weyl-invariant and the ordered slice is a
fundamental domain, so orbit reduction discards nothing. This is recorded as
available and NOT implemented here.

## Reproducing

```sh
python scripts/census_inner_sandwich.py            # the N=3 ladder, about 25 s
python scripts/census_inner_sandwich.py --wider    # all nine, about 4 min
uv run pytest tests/test_census_inner_sandwich.py
```

The artifact is rewritten after every system, so a run that is interrupted in
the middle of `(5,10)` still ships the eight that closed.
