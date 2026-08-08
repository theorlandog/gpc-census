# The facet-index stability probe

Pre-registration: `docs/prereg_facet_index_stability.md`, committed at `ad55d4c`
before any measurement. Generator: `scripts/facet_index_stability.py`. Artifact:
`results/data/facet_index_stability.json`. Guard test:
`tests/test_facet_index_stability.py`.

Motivated by the measured result in `docs/bdr_constructor_comparison.md`: the
stage with the worst long-run scaling is redundancy elimination, because the
floor no external tool can remove is composition and composition grows with the
retained facet count. So the question worth asking is not "can we generate
candidates faster" but "can an intrinsic criterion replace the reduction".

## Headline

❌ **The AK2008 coefficient-one criterion is REFUTED as a facetness test.** Every
facet has cyclic-Schubert coefficient one, but so do 23 of 45 removed rows at
`(3,7)` and 36 of 106 at `(3,8)`. Coefficient one is necessary and not
sufficient. The repository's standing refusal to use it as a facetness or
irredundancy gate is now vindicated by exact counterexamples rather than by
caution.

🟦 **It survives as a SOUND PREFILTER, and that is worth something measurable.**
Because the forward direction holds, filtering the certified rows down to the
coefficient-one ones cannot discard a facet at these ranks. It cuts the
reduction input by 45 percent at `(3,7)` and 51 percent at `(3,8)`, gives
1.66x and 1.79x on the reduction stage, and returns a byte-identical retained
set. It is not gate-eligible until the forward direction is proved.

❌ **Trailing-zero padding is not a functor on facets.** Padded lower-rank rows
are outright VIOLATED at the next rank in 7 of the 8 published maps where
validity is decidable. The naive stability hypothesis is dead.

🟦 **The normals compress hard, across every published system.** Facet counts
run 1 to 161 while the number of distinct sorted-`tau` shapes runs 1 to 21, a
ratio at or below 0.29 in all twelve systems and at or below 0.22 in every
system with more than four facets. The explosion is in placement, not in shape.

## Part A: the padding and shape census (LANDED, not scored)

Rule R4: these tables were computed in a scratch session before the
pre-registration existed. They are landed with a generator and a guard test,
not reported as discoveries. The artifact records the scratch values beside the
landed ones and all three groups agree exactly.

Only P6, the extension to `N != 3`, was predictive.

### Padding is not a valid map

For each published consecutive pair, every lower facet is padded with a trailing
zero and tested exactly over the rationals against the certified vertex list of
the higher system.

| map | lower facets | violated | valid and tight | in higher facet list |
|---|---|---|---|---|
| (3,6) to (3,7) | 1 | **1** | 0 | 0 |
| (3,7) to (3,8) | 4 | 0 | 4 | 4 |
| (3,8) to (3,9) | 31 | **21** | 10 | 10 |
| (3,9) to (3,10) | 52 | **8** | 44 | 37 |
| (4,7) to (4,8) | 4 | **4** | 0 | 0 |
| (4,8) to (4,9) | 15 | **7** | 8 | 8 |
| (4,9) to (4,10) | 60 | **41** | 19 | 19 |
| (5,9) to (5,10) | 60 | **59** | 1 | 0 |

`(5,8)` to `(5,9)` is reported as validity UNAVAILABLE: `(5,9)` carries a
constraint list but no precomputed vertex file, so the exact test cannot run.
Facet membership is still decided there and is 8 of 31.

Two structural facts fall out.

1. Padding fails badly and generically. Borland-Dennis padded to rank 7 is
   violated. This is expected once stated properly: the rank-`d` polytope is a
   FACE of the rank-`d+1` polytope, and a facet of a face need not extend to a
   valid inequality of the whole.
2. **The `valid but strictly interior` column is zero in every single map.** A
   padded facet never becomes merely redundant. It either stays supporting or
   becomes invalid. That is a sharper statement than the violation counts and it
   holds across all three particle numbers.

### Reconciliation with the ordered-stability campaign

`docs/ordered_representation_stability.md` reports a padding extension that is
total and exact, 612 of 612, which reads as a flat contradiction of the table
above until one notices that **the two maps are different**. That campaign's
`pi_pad` appends the largest last entry keeping the row valid. This one appends
zero in affine coordinates. Neither result is wrong and the pair is sharper
than either alone.

Checked at merge time, importing `pi_pad` rather than reimplementing it so the
two campaigns cannot drift:

| map | in higher facet list | zero-padding agrees with tight |
|---|---|---|
| (3,6) to (3,7) | 0 | 0 |
| (3,7) to (3,8) | 4 | 4 |
| (3,8) to (3,9) | 10 | 10 |
| (3,9) to (3,10) | 37 | 37 |
| (4,7) to (4,8) | 0 | 0 |
| (4,8) to (4,9) | 8 | 8 |
| (4,9) to (4,10) | 19 | 19 |
| (5,9) to (5,10) | 0 | 0 |

🟦 CONFIRMED, 8 of 8 maps, three particle numbers: **a trailing-zero padded
facet lands in the higher facet list exactly when the trailing-zero extension
coincides with the tight extension.** The naive map succeeds precisely on the
rows where it was already tight, and fails everywhere else. So the correct
statement is not that padding fails, it is that the trivial extension is the
wrong one and tightening is the repair.

This also corrects an overreach in an earlier draft of this document, which
concluded that a stability theorem must act on something other than the normal
vector. The ordered-stability branch shows the repair does act on normals, by
tightening the new entry rather than by leaving it zero. What is dead is the
trivial extension, not extension of normals as such.

Note the joint scope limit, which neither branch removes: that branch's
Theorem 1 closes the finite-generation route in this category regardless, so an
exact inheritance law is what the pair delivers, not a stable template theorem.

### The normals compress

| system | facets | max abs tau entry | distinct shapes | shape/facet |
|---|---|---|---|---|
| (3,6) | 1 | 1 | 1 | 1.00 |
| (3,7) | 4 | 2 | 1 | 0.25 |
| (3,8) | 31 | 7 | 9 | 0.29 |
| (3,9) | 52 | 7 | 8 | 0.15 |
| (3,10) | 93 | 7 | 16 | 0.17 |
| (4,7) | 4 | 1 | 1 | 0.25 |
| (4,8) | 15 | 3 | 3 | 0.20 |
| (4,9) | 60 | 9 | 13 | 0.22 |
| (4,10) | 125 | 15 | 15 | 0.12 |
| (5,8) | 31 | 13 | 9 | 0.29 |
| (5,9) | 60 | 23 | 13 | 0.22 |
| (5,10) | 161 | 34 | 21 | 0.13 |

The `N = 3` entry bound of 7 holds at ranks 8, 9 and 10 but does **not**
transfer across particle number: the bound grows to 15 at `(4,10)` and 34 at
`(5,10)`. That was predicted in advance and is exactly what `tau_i = N a_i - b`
implies, so no entry bound uniform in `N` is claimed.

The shape ratio is the quantity that does transfer. P6 predicted at most 0.25
for the four `N != 3` systems and measured 0.20, 0.22, 0.12 and 0.13.

Conditionality: the rank-9 and rank-10 lists are conjecturally complete per
AK2008, so every row of those two tables involving rank 9 or 10 is conditional
on that completeness.

## Part B: the coefficient-one converse (PREDICTIVE)

The production pipeline defers exact cyclic-Schubert evaluation until after
redundancy elimination, precisely so it runs on retained rows rather than on
every certified candidate. The consequence is that **the coefficients of the
rows reduction REMOVES had never been computed at any rank.** That gap is the
converse direction of the criterion, and it is what this measures.

State before the run, checked not assumed: `(3,7)` retained rows were already
known to have coefficient one (`results/data/stage1_ressayre_3_7.json`).
`(3,8)` retained coefficients were stored nowhere. Removed-row coefficients
existed nowhere at any rank.

| | `(3,7)` | `(3,8)` |
|---|---|---|
| certified rows | 49 | 137 |
| facets | 4 | 31 |
| facets with coefficient one | **4 of 4** | **31 of 31** |
| removed rows | 45 | 106 |
| **removed rows with coefficient one** | **23** | **36** |
| removed with coefficient not one | 22 | 70 |
| coefficient zero | 0 | 0 |
| ill-formed cyclic problems | 0 | 0 |
| unresolved by budget | 0 | 0 |
| distinct coefficients on removed rows | 1, 2 | 1, 2, 3, 4, 6, 8, 9, 12, 18 |

Smallest counterexample recorded, a removed row with coefficient one at
`(3,7)`: `tau = [-2, 1, 1, 4, -2, -5, -5]`.

So on certified rows, at both ranks:

```text
facet  =>  coefficient one          holds, 35 of 35
coefficient one  =>  facet          FAILS, 59 counterexamples
```

### Why this is the expected answer, stated after the fact but not a surprise

The local test is one selected cyclic coefficient of one divided-difference
recurrence, not the general Belkale-Kumar or Ressayre deformed product. The
pre-registration predicted the refutation for that reason and it scored PASS.
This refutes THAT criterion. It says nothing about the general deformed-product
criterion, which remains untested here.

### The prefilter, priced (added after scoring)

The forward direction holding is what makes filtering safe: at a rank where
every facet has coefficient one, discarding non-coefficient-one rows cannot
lose a facet.

| system | reduction input | after filter | cut | reduction stage | retained set |
|---|---|---|---|---|---|
| (3,7) | 49 | 27 | 45% | 0.239 s to 0.144 s, 1.66x | identical |
| (3,8) | 137 | 67 | 51% | 0.995 s to 0.555 s, 1.79x | identical |

This is the first measured attack on the stage the BDR comparison identified as
the growing floor. It is not integrated, and must not be, because its soundness
rests on an observation at two ranks and not on a theorem.

## Scored pre-registration

| id | prediction | verdict | scope | independent scope |
|---|---|---|---|---|
| P1 | every retained row at `(3,8)` has coefficient 1 | **PASS** | 31 rows | 1 system |
| P2 | some removed row has coefficient 1, refuting the converse | **PASS** | 151 rows | 2 systems |
| P3 | removed rows at `(3,8)` show both values | **PASS** | 106 rows | 1 system |
| P4 | no certified row has coefficient 0 | **PASS** | 186 rows | 2 systems |
| P5 | every certified row is a well-formed cyclic problem | **PASS** | 186 rows | 2 systems |
| P6 | shape/facet ratio at most 0.25 for the four `N != 3` systems | **PASS** | 4 systems | 4 |

Part A's three scratch tables carry the verdict LANDED, which is not a score,
and the artifact records exact agreement between scratch and landed values for
all of padding violations, entry bounds and shape counts.

All six passed, which is itself worth flagging: a preregistration where nothing
fails is weaker evidence of a sharp hypothesis than one where something does.
P2 was the prediction designed to be decisive, and it decided against the
criterion.

## What this establishes and what it does not

**Establishes.**

- Coefficient one is necessary but not sufficient for facetness on certified
  rows at `(3,7)` and `(3,8)`, with 59 exact counterexamples.
- Trailing-zero padding is not validity preserving on facets, in 7 of 8
  decidable published maps, and a padded facet is never merely redundant.
- Facet normals occupy at most 0.29 of their count in distinct shapes across all
  twelve published systems, and at most 0.22 wherever there are more than four
  facets.
- A sound coefficient-one prefilter cuts reduction input roughly in half at both
  ranks with an identical retained set.

**Does not establish.**

- That every facet has coefficient one at any rank beyond 7 and 8, for any
  particle number beyond 3. Two ranks, one `N`. The prefilter is not
  gate-eligible without that theorem. 🔬 OPEN.
- Anything about the general Belkale-Kumar or Ressayre deformed-product
  criterion. Untested here. 🔬 OPEN.
- Any statement uniform in `d`. The shape ratio is a twelve-system observation,
  not a bound.
- Anything unconditional at ranks 9 and 10, whose lists are conjecturally
  complete per AK2008.

## The next theorem, named

> **Forward coefficient-one.** Every facet of the ordered `wedge^N C^d` moment
> polytope has cyclic-Schubert coefficient one.

That is the exact statement the prefilter needs and the only missing piece
between a measured 1.79x on the growing stage and a gate the repository is
allowed to ship. The converse is now known to be false, so this is the whole
remaining question, and it is a forward-validity statement of exactly the kind
`docs/fixed_n3_generator_methodology.md` requires before promoting anything to a
pruning rule.

The cheapest next experiment is the same census at `(3,9)`, `(4,8)` and `(4,9)`.
It tests the forward direction at a third rank and at a second particle number
at once, and it costs one run of an already-written script with
`COEFFICIENT_SYSTEMS` extended.

## Reproducing

```sh
uv run python scripts/facet_index_stability.py
uv run pytest tests/test_facet_index_stability.py
```
