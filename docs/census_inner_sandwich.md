# Completeness without candidate exhaustiveness

Status: **P = O PROVED at (3,6), (3,7), (3,8), (3,9) and (3,10)**, with
candidate exhaustiveness unused. Ranks 9 and 10 had no in-repo completeness
proof before this.

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
Every stored row certifies:

| system | rows | certified | rejected | unresolved |
|---|---|---|---|---|
| (3,7) | 4 | 4 | 0 | 0 |
| (3,8) | 31 | 31 | 0 | 0 |
| (3,9) | 52 | 52 | 0 | 0 |
| (3,10) | 93 | 93 | 0 | 0 |

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

Every N = 3 census system closes.

| system | outer rows | vertices | attained | closed | time |
|---|---|---|---|---|---|
| (3,6) | 1 + 3 equalities | 4 | 4 | **yes** | 0.6 s |
| (3,7) | 4 | 10 | 10 | **yes** | 0.1 s |
| (3,8) | 31 | 38 | 38 | **yes** | 2.2 s |
| (3,9) | 52 | 58 | 58 | **yes** | 15 s |
| (3,10) | 93 | 113 | 113 | **yes** | 254 s |

## The predicted wall was not the wall

Exact vertex enumeration of a highly degenerate polytope in affine dimension
`d - 2` with tens of rows was expected to be the obstacle, with naive
active-set costing `C(F,s) * s^3`. It was not: the double-description
enumerator in `gpc_census.exact_polytope`, built for the orbit-closure layer,
does it in 0.0 s, 0.6 s, 12 s and 247 s at ranks 7, 8, 9 and 10 and reproduces
the certified vertex counts exactly. No lrs, no pivoting rewrite, no symmetry
reduction was needed at these ranks.

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

**BOUNDED NULL, recorded as such.** At the ranks run here the loop never fired:
zero vertices resisted at all five systems. So the loop is implemented and
exercised only on its trivial branch, and it has NOT been demonstrated to
recover a facet. Its first real test is a rank where the inner side is thin,
which is (3,11): 19 TRUE vertices from a bracket rather than a census, so a
resisting vertex there is much more likely to be a genuine missing facet than a
hard attainment problem.

## What this does and does not license

**Does.** P = O at five systems, by an argument that uses no exhaustiveness
claim. Ranks 9 and 10 in particular now have an in-repo completeness proof
rather than a citation of the published tables plus a census that was built
from them.

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
python scripts/census_inner_sandwich.py            # the N=3 ladder, about 5 min
python scripts/census_inner_sandwich.py --wider    # also N=4 and N=5
uv run pytest tests/test_census_inner_sandwich.py
```
